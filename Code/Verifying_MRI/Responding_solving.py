# Calculating currents in each ring on pyfftw way (one dimensional)
import os
import sys

sys.path.insert(1, os.path.join(sys.path[0], '..'))
import numpy as np
from Impedance_matrix import Mnm
from scipy.sparse.linalg import gmres, LinearOperator, bicgstab, minres, lobpcg, cg
from pyfftw import pyfftw
from tqdm import tqdm

# Function for creating circulant vectors
def Circvec(rings_3d_str, rings_3d_col, data):
    Nz_str, Ny_str, Nx_str = rings_3d_str.shape
    Nz_col, Ny_col, Nx_col = rings_3d_col.shape
    nz, ny, nx = Nz_str + Nz_col - 1, Ny_str + Ny_col - 1, Nx_str + Nx_col - 1
    Z_circvecs = pyfftw.empty_aligned((nz, ny, nx), dtype = 'complex128')
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                x_str_id = (nx - x) * (x >= Nx_col)
                x_col_id = x * (x < Nx_col)

                y_str_id = (ny - y) * (y >= Ny_col)
                y_col_id = y * (y < Ny_col)

                z_str_id = (nz - z) * (z >= Nz_col)
                z_col_id = z * (z < Nz_col)
                
                Z_circvecs[z][y][x] = Mnm(rings_3d_str[z_str_id][y_str_id][x_str_id], rings_3d_col[z_col_id][y_col_id][x_col_id], data)
    return Z_circvecs

def fft_dot(I, ZI, FFT_Z_circvecs, i_vecs, ifft_i_vecs):
    Nz, Ny, Nx = I.shape
    nz, ny, nx = i_vecs.shape

    i_vecs[:Nz, :Ny, :Nx] = I

    pyfftw.FFTW(i_vecs, ifft_i_vecs, axes = (0, 1, 2), direction='FFTW_BACKWARD').execute()
    pyfftw.FFTW(FFT_Z_circvecs * ifft_i_vecs/nz/ny/nx, ZI, axes = (0, 1, 2)).execute()
    
    return ZI[:nz - Nz + 1, :ny - Ny + 1, :nx - Nx + 1].ravel()

def solvesystem(rings_4d, Responded_ring, M_0, Omega, Inductance = {}, phi_0z = 1, tol = 1e5):
    # Unpacking parameters

    orientations = rings_4d.keys()
    Number = np.sum([value.size for value in rings_4d.values()]) + 1

    FFT_M_circvecs = {}
    i_vecs = {}
    ifft_i_vecs = {}
    MI_vecs = {}
    Coil_M = []
    # Preparing empty arrays for pyfftw
    print('Cirvecs forming')
    for pos_str in tqdm(orientations):
        rings_str = rings_4d[pos_str]
        FFT_M_circvecs[pos_str] = {}
        i_vecs[pos_str] = {}
        ifft_i_vecs[pos_str] = {}
        MI_vecs[pos_str] = {}
        for pos_col in orientations:
            rings_col = rings_4d[pos_col]
            M_circvecs = Circvec(rings_str, rings_col, Inductance)

            N_circ = np.array(rings_str.shape) + np.array(rings_col.shape) - 1
            i_vecs[pos_str][pos_col] = np.zeros(N_circ, dtype=complex)
            MI_vecs[pos_str][pos_col] = pyfftw.empty_aligned(N_circ, dtype = 'complex128')

            FFT_M_circvecs[pos_str][pos_col] = pyfftw.empty_aligned(N_circ, dtype = 'complex128')
            ifft_i_vecs[pos_str][pos_col] = pyfftw.empty_aligned(N_circ, dtype = 'complex128')
            pyfftw.FFTW(M_circvecs, FFT_M_circvecs[pos_str][pos_col], axes = (0, 1, 2)).execute()
        rings = np.ravel(rings_str)
        for ring in rings:
            Coil_M.append(Mnm(Responded_ring, ring, Inductance))
    Coil_M = np.array(Coil_M)
        


    print('Circvecs: Done')



    # Caclulating current in each ring
    print('FFT solving')
    CURRENTS = []
    Responded_impedance = []
    I_old = np.ones(Number, dtype = np.complex128)/M_0(Omega[0])
    I_old[-1] = 1/(Responded_ring.R/1j/Omega[0] + Responded_ring.L - 1/(Omega[0] ** 2 * Responded_ring.C))
    Phi_0z = np.zeros(Number)
    Phi_0z[-1] = 1
    for omega in tqdm(Omega):
        def LO(Currents):
            MI = M_0(omega) * Currents
            # Calculate current in the coil and its influence on the rings
            Coil_Current = Currents[-1]
            MI[-1] += Coil_M @ Currents[:-1]
            MI[:-1] += Coil_M * Coil_Current

            I = Currents[:-1]

            # Make start and end indexes for each orientation
            start_str = 0
            end_str = 0
            for pos_str in orientations:
                end_str += rings_4d[pos_str].size

                start_col = 0
                end_col = 0
                for pos_col in orientations:
                    end_col += rings_4d[pos_col].size
                    MI
                    MI[start_str: end_str] += fft_dot(I[start_col:end_col].reshape(rings_4d[pos_col].shape),
                                                      MI_vecs[pos_str][pos_col],
                                                      FFT_M_circvecs[pos_str][pos_col],
                                                      i_vecs[pos_str][pos_col],
                                                      ifft_i_vecs[pos_str][pos_col])
                    start_col += rings_4d[pos_col].size
                start_str += rings_4d[pos_str].size
            
            return MI 
        
        M = LinearOperator(dtype = np.complex128, shape=(Number, Number), matvec=LO)
        I, info = bicgstab(M, Phi_0z, x0 = I_old, tol = tol)

        if info != 0:
            print(f'omega = {omega/2/np.pi/1e6} MHz did not converge')
        
        CURRENTS.append(I*phi_0z)
        Responded_impedance.append(Coil_M * 1j * omega @ I[:-1]/I[-1])
        I_old = I

    print(f'FFT solving: Done, shape = {[(pos, rings_4d[pos].shape) for pos in orientations]}')
    Data = {}

    Data['Omega'] = list(Omega)
    Data['Responded_impedance'] = [list(np.real(Responded_impedance)), list(np.imag(Responded_impedance))]
    Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
    Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]
 
    return Data  
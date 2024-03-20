# Calculating currents in each ring on pyfftw way (one dimensional)

import numpy as np
from Impedance_matrix import Mnm
from scipy.sparse.linalg import gmres, LinearOperator, bicgstab, minres, lobpcg, cg
from pyfftw import pyfftw
from tqdm import tqdm


def Circvec(rings_3d, data):
    Nz, Ny, Nx = rings_3d.shape
    nz, ny, nx = 2*Nz-1, 2*Ny-1, 2*Nx-1
    Z_circvecs = pyfftw.empty_aligned((nz, ny, nx), dtype = 'complex128')
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                x_str_id = (nx - x) * (x >= Nx)
                x_col_id = x * (x < Nx)

                y_str_id = (ny - y) * (y >= Ny)
                y_col_id = y * (y < Ny)

                z_str_id = (nz - z) * (z >= Nz)
                z_col_id = z * (z < Nz)
                
                Z_circvecs[z][y][x] = Mnm(rings_3d[z_str_id][y_str_id][x_str_id], rings_3d[z_col_id][y_col_id][x_col_id], data)
    return Z_circvecs

def fft_dot(I, ZI, FFT_Z_circvecs, i_vecs, ifft_i_vecs):
    nz, ny, nx = ZI.shape[0], ZI.shape[1], ZI.shape[2]
    Nz, Ny, Nx = (nz + 1) // 2, (ny + 1) // 2, (nx + 1) // 2
    i_vecs[:Nz, :Ny, :Nx] = I.reshape((Nz, Ny, Nx))

    pyfftw.FFTW(i_vecs, ifft_i_vecs, axes = (0, 1, 2), direction='FFTW_BACKWARD').execute()
    pyfftw.FFTW(FFT_Z_circvecs * ifft_i_vecs/nz/ny/nx, ZI, axes = (0, 1, 2)).execute()
    
    return ZI[:Nz, :Ny, :Nx].reshape(len(I))

def solvesystem(Params, N, Omega, grad = [0, 0, 0], Inductance = {}):
    # Unpacking parameters

    mu_0 = 4 * np.pi * 10 ** -7
    H_0z = 1

    L, C, R, W, Radius, Dz, Dy, Dx, Rectangle_packing = Params['L'], Params['C'], Params['R'], Params['W'], Params['Radius'], Params['Dz'], Params['Dy'], Params['Dx'], Params['Rectangle_packing']
    CURRENTS = []
    Polarisation = []

    nz = N[f"z"] + N[f"z"] - 1
    ny = N[f"y"] + N[f"y"] - 1
    nx = N[f"x"] + N[f"x"] - 1

    # Creating 3D array of rings
    Rings = Rectangle_packing(N, Radius, Dx, Dy, Dz, W)
    Rings_3d = Rings.reshape((N['z'], N['y'], N['x']))
    # Creating circulant
    M_circvecs = Circvec(Rings_3d, Inductance)
    Number = len(Rings)

    print('circvecs:Done')
    
    # Preparing empty arrays for pyfftw
    i_vecs = np.zeros((nz, ny, nx), dtype=complex)
    MI = pyfftw.empty_aligned((nz, ny, nx), dtype = 'complex128')

    FFT_M_circvecs = pyfftw.empty_aligned((nz, ny, nx), dtype = 'complex128')
    ifft_i_vecs = pyfftw.empty_aligned((nz, ny, nx), dtype = 'complex128')
    pyfftw.FFTW(M_circvecs, FFT_M_circvecs, axes = (0, 1, 2)).execute()
    
    # Caclulating diagonal of Z matrix using gradient
    omega_0 = 1/np.sqrt(L * C)
    print('FFT solving')
    I_old = np.zeros(Number)
    z, y, x = np.meshgrid(np.arange(N['x']), np.arange(N['y']), np.arange(N['z']), indexing = 'ij')
    gradx = x * grad[0] / N['x']
    grady = y * grad[1] / N['y']
    gradz = z * grad[2] / N['z']
    
    omega_0 = 1/np.sqrt(L * C)
    Omega_3d = omega_0 * (1 + gradz + grady + gradx)
    Cgrad = (1/L/Omega_3d ** 2).reshape(Number)
    
    for omega in tqdm(Omega):
        M_0 = R/1j/omega + L - 1/(omega ** 2 * Cgrad)

        def LO(I):
            return fft_dot(I, MI, FFT_M_circvecs, i_vecs, ifft_i_vecs) + M_0 * I

        A = LinearOperator(dtype = np.complex128, shape=(Number, Number), matvec=LO)
        E = H_0z  * np.ones(Number) * np.pi * Radius ** 2 * mu_0
        I, info = bicgstab(A, E, x0 = I_old, tol = E.sum()/Number/1e6)

        Volume = N['x'] * N['y'] * (N['z']-1) * Dz * Dy * Dx
        polarisation = np.pi * Radius ** 2 * I.sum()/Volume/(H_0z*omega)
        
        #print(f'Solving: Done - {round(omega/2/np.pi/10 ** 6, 1)}MHz, grad = {grad}')
        CURRENTS.append(I)
        Polarisation.append(polarisation)
        I_old = I



    Data = {}

    Data['N'] = N
    Data['Omega'] = list(Omega)
    Data['grad'] = grad
    Data['params'] = [Dz, Dy, Dx, W, L, C, R, Radius]
    Data['RealPolarisation'] = list(np.real(Polarisation))
    Data['ImagPolarisation'] = list(np.imag(Polarisation))

    DataCurrents = {}
    DataCurrents['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
    DataCurrents['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]
 
    return Data, DataCurrents   
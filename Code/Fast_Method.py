# Base voltage method
import os
import scipy
import numpy as np
import scipy.linalg
from tqdm import tqdm 
from Impedance_matrix import Mnm
from scipy.sparse.linalg import lgmres, gmres, LinearOperator


os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
# Calculating currents in each ring on pyfftw way (one dimensional)
solvers = {
    'gmres': gmres,
    'lgmres': lgmres
}
# Function for creating circulant vectors
def Circvec(rings_3d_str:list, rings_3d_col:list, data:dict) -> np.ndarray:
    """Function for creating circulant vectors that fully
    describe matrix M and allow to calculate M * I
    using FFT and N * log(N) operations. After circulant multiplication,
    then we cutoff the result to the size of the original matrix.

    Parameters
    ----------
    rings_3d_str : list
        3d array with rings of the same orientation
    rings_3d_col : list
        3d array with rings of the same orientation
    data : dict
        dictionary with modeling parameters

    Returns
    -------
    np.ndarray
        3d array with three times circulant vectors
    """    
    Nz_str, Ny_str, Nx_str = rings_3d_str.shape
    Nz_col, Ny_col, Nx_col = rings_3d_col.shape
    nz, ny, nx = Nz_str + Nz_col - 1, Ny_str + Ny_col - 1, Nx_str + Nx_col - 1
    Z_circvecs = np.zeros((nz, ny, nx), dtype = complex)
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

def fft_dot(I:np.ndarray, ZI:np.ndarray, FFT_Z_circvecs:np.ndarray, i_vecs:np.ndarray, ifft_i_vecs:np.ndarray) -> np.ndarray:
    """Function for calculating dot product of Toeplitz matrix M and vector I
    using FFT and N * log(N) operations.

    Parameters
    ----------
    I : np.ndarray
        vector with currents for multiplication
    ZI : np.ndarray
        vector with result of multiplication
    FFT_Z_circvecs : np.ndarray
        3d array with 3d Fourier transform for three times circulant vectors
    i_vecs : np.ndarray
        3d array with extra space for circulant vectors
    ifft_i_vecs : np.ndarray
        3d array with extra space for inverse Fourier transform

    Returns
    -------
    np.ndarray
        _description_
    """    
    Nz, Ny, Nx = I.shape
    nz, ny, nx = i_vecs.shape

    i_vecs[:Nz, :Ny, :Nx] = I

    ifft_i_vecs = scipy.fft.ifftn(i_vecs, axes = (0, 1, 2))
    ZI = scipy.fft.fftn(FFT_Z_circvecs * ifft_i_vecs, axes = (0, 1, 2))

    return ZI[:nz - Nz + 1, :ny - Ny + 1, :nx - Nx + 1].ravel()

def ifft_dot(I:np.ndarray, ZI:np.ndarray, FFT_Z_circvecs:np.ndarray, i_vecs:np.ndarray, ifft_i_vecs:np.ndarray) -> np.ndarray:
    """Function to calculate dot product of inverse circulant matrix based on 
    Toeplitz matrix M and vector I using FFT and N * log(N) operations.

    Parameters
    ----------
    I : np.ndarray
        vector with currents for multiplication
    ZI : np.ndarray
        vector with result of multiplication
    FFT_Z_circvecs : np.ndarray
        3d array with 3d Fourier transform for three times circulant vectors
    i_vecs : np.ndarray
        3d array with extra space for circulant vectors
    ifft_i_vecs : np.ndarray
        3d array with extra space for inverse Fourier transform

    Returns
    -------
    np.ndarray
        _description_
    """    
    Nz, Ny, Nx = I.shape
    nz, ny, nx = i_vecs.shape

    i_vecs[:Nz, :Ny, :Nx] = I

    ifft_i_vecs = scipy.fft.ifftn(i_vecs, axes = (0, 1, 2))
    ZI = scipy.fft.fftn(ifft_i_vecs/FFT_Z_circvecs, axes = (0, 1, 2))

    return ZI[:nz - Nz + 1, :ny - Ny + 1, :nx - Nx + 1].ravel()

def solvesystem(Params:dict, rings_4d:dict, phi_0z_4d:dict, Inductance:dict = {}, find:str = 'Voltage', tol:float = 1e-5)->dict:
    """Function for solving system of equations using Toeplitz structure for 
    impedance matrix M and iterative solver for linear system of equations.
    Highly recommended to use Voltage method for solving system of equations and 
    you must use it for non-rectangle structures because of infinitly conditioned
    matrix for Currents way.

    Parameters
    ----------
    Params : dict
        Dictionary with modeling parameters
    rings_4d : dict
        Dictionary with rings for each orientation
    phi_0z_4d : dict
        Dictionary with external fluxes for each orientation
    Inductance : dict, optional
        Dictionary with mutual inductance, by default {}
    find : str, optional
        Which matrix equation to solve, by default 'Currents'
    tol : float, optional
        tolerance for residials in iterative methods, by default 1e-5

    Returns
    -------
    dict
        Dictionary with results of solving system of equations
    """    
    # Unpacking parameters
    Params['Solver_type'] = 'Fast'
    solve = solvers[Params['Solver_name']]
    Omegas = Params['Omega']    

    Omega = np.linspace(Omegas[0], Omegas[1], Omegas[2])
    
    rings = sum([rings_4d[orientation] for orientation in rings_4d], [])
    phi_0z = np.array(sum([phi_0z_4d[orientation] for orientation in phi_0z_4d], []))

    orientations = rings_4d.keys()
    for orientation in orientations:
        Nz, Ny, Nx = Params['N'][orientation]['nz'], Params['N'][orientation]['ny'], Params['N'][orientation]['nx']
        rings_4d[orientation] = np.array(rings_4d[orientation]).reshape(Nz, Ny, Nx)
        phi_0z_4d[orientation] = np.array(phi_0z_4d[orientation])
    Number = np.sum([value.size for value in rings_4d.values()])

    FFT_M_circvecs = {}
    i_vecs = {}
    ifft_i_vecs = {}
    MI_vecs = {}

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
            MI_vecs[pos_str][pos_col] = np.zeros(N_circ, dtype=complex)

            FFT_M_circvecs[pos_str][pos_col] = scipy.fft.fftn(M_circvecs, axes = (0, 1, 2))
            ifft_i_vecs[pos_str][pos_col] = np.zeros(N_circ, dtype=complex)

    print('Circvecs: Done')

    # Calculating diagonal of M matrix
    L, C, R = [], [], []
    for ring in rings:
        L.append(ring.L)
        C.append(ring.C)
        R.append(ring.R)
    L, C, R = np.array(L), np.array(C), np.array(R)
    
    M_0 = lambda Omega: (R - 1j * Omega * L + 1j/(Omega * C))/1j/Omega

    # Caclulating current in each ring
    print('FFT solving (voltage)')
    CURRENTS = []
    I_old = np.ones(Number, dtype = np.complex128)/M_0(Omega[0])
    Phi_0z = phi_0z
    P = []

    for omega in tqdm(Omega):
        M_diag = M_0(omega)
        if find == 'Voltage':
            def LO(I):
                MI = M_diag * I
                # Make start and end indexes for each orientation
                start_str = 0
                end_str = 0
                for pos_str in orientations:
                    end_str += rings_4d[pos_str].size

                    start_col = 0
                    end_col = 0
                    for pos_col in orientations:
                        end_col += rings_4d[pos_col].size
                        MI[start_str: end_str] -= fft_dot(I[start_col:end_col].reshape(rings_4d[pos_col].shape),
                                                          MI_vecs[pos_str][pos_col],
                                                          FFT_M_circvecs[pos_str][pos_col],
                                                          i_vecs[pos_str][pos_col],
                                                          ifft_i_vecs[pos_str][pos_col])
                        start_col += rings_4d[pos_col].size
                    start_str += rings_4d[pos_str].size
                return MI/M_diag
            M = LinearOperator(dtype = np.complex128, shape=(Number, Number), matvec=LO)
            I, info = solve(M, Phi_0z/M_diag, x0 = I_old, rtol = tol, atol = 0)
        elif find == 'Currents':
            def LO(I):
                MI = M_diag * I
                # Make start and end indexes for each orientation
                start_str = 0
                end_str = 0
                for pos_str in orientations:
                    end_str += rings_4d[pos_str].size

                    start_col = 0
                    end_col = 0
                    for pos_col in orientations:
                        end_col += rings_4d[pos_col].size
                        MI[start_str: end_str] -= fft_dot(I[start_col:end_col].reshape(rings_4d[pos_col].shape),
                                                          MI_vecs[pos_str][pos_col],
                                                          FFT_M_circvecs[pos_str][pos_col],
                                                          i_vecs[pos_str][pos_col],
                                                          ifft_i_vecs[pos_str][pos_col])
                        start_col += rings_4d[pos_col].size
                    start_str += rings_4d[pos_str].size
                return MI
        
            M = LinearOperator(dtype = np.complex128, shape=(Number, Number), matvec=LO)
            I, info = solve(M, Phi_0z, x0 = I_old, rtol = tol, atol = 0)

        if info != 0:
            print(f'f = {omega/2/np.pi/1e6} MHz did not converge')
        
        CURRENTS.append(I)
        start = 0
        p = []
        for pos in orientations:
            end = start + rings_4d[pos].size
            p.append(np.sum(I[start:end])/(end-start))
            start = end
        P.append(p)

        I_old = I
    
    P = np.array(P) * Params['P_0z']

    print(f'FFT solving: Done, shape = {[(pos, rings_4d[pos].shape) for pos in orientations]}')
    Data = {}

    Data['Params'] = Params
    Data['Omega'] = Omega
    Data['Currents'] = CURRENTS
    Data['Polarization'] = P
    Data['Phi_0z'] = list(phi_0z)

    return Data
# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from Impedance_matrix import Matrix
from tqdm import tqdm





DATA = {}

H_0z = 1

Dimensions = ('z', 'y', 'x')
Orientation = ('z')

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

    dim_old = N['z'] * N['z'] * N['y'] * N['y'] * N['x'] * N['x']
    dim_new = (N['z'], N['z'], N['y'], N['y'], N['x'], N['x'])


    Rings = Rectangle_packing(N, Radius, Dx, Dy, Dz, W)
    Number = len(Rings)

    M = Matrix(Rings, Rings, Data = Inductance)
    print('Matrix: Done')

    print('Straight solving')

    x, y, z = np.meshgrid(np.arange(N['x']), np.arange(N['y']), np.arange(N['z']), indexing = 'ij')
    gradx = x * grad[0] / N['x']
    grady = y * grad[1] / N['y']
    gradz = z * grad[2] / N['z']
    
    omega_0 = 1/np.sqrt(L * C)
    Omega_3d = omega_0 * (1 + gradz + grady + gradx)
    Cgrad = (1/L/Omega_3d ** 2).reshape(Number)
    E = H_0z * np.ones(Number) * np.pi * Radius ** 2 * mu_0

    for omega in tqdm(Omega):

        Zdiag = R/1j/omega + L - 1/( omega ** 2 * Cgrad)

        M_0 = np.diag(Zdiag)
        I = solve(M + M_0, E)
        Volume = N['x'] * N['y'] * (N['z']-1) * Dz * Dy * Dx
        polarisation = np.pi * Radius ** 2 * I.sum()/Volume/(H_0z*omega)
        
        #print(f'Solving: Done - {round(omega/2/np.pi/10 ** 6, 1)}MHz, grad = {grad}')
        CURRENTS.append(I)
        Polarisation.append(polarisation)

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
    

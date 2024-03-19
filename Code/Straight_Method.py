# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
import matplotlib.pyplot as plt
import json





DATA = {}

H_0z = 1

Dimensions = ('z', 'y', 'x')
Orientation = 'z'

def solvesystem(n, Omega, grad = 0.1):
    Data = {}
    CURRENTS = []
    MaxCurrents = []
    Polarisation = []
    MinCurrents = []

    N = {}

    N['x'], N['y'], N['z'] = n, n, n

    nz = N[f"z"] + N[f"z"] - 1
    ny = N[f"y"] + N[f"y"] - 1
    nx = N[f"x"] + N[f"x"] - 1
    dim_old = N['z'] * N['z'] * N['y'] * N['y'] * N['x'] * N['x']
    dim_new = (N['z'], N['z'], N['y'], N['y'], N['x'], N['x'])


    Rings = Rectangle_packing(N["x"], N["y"], N['z'], Radius, a, b, c, w)
    Number = len(Rings)
    Rings_3d = Rings.reshape((N['z'], N['y'], N['x']))
    Z_circvecs = np.zeros((nz, ny, nx), dtype = complex)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                x_str_id = (nx - x) * (x >= N['x'])
                x_col_id = x * (x < N['x'])

                y_str_id = (ny - y) * (y >= N['y'])
                y_col_id = y * (y < N['y'])

                z_str_id = (nz - z) * (z >= N['z'])
                z_col_id = z * (z < N['z'])

                Z_circvecs[z][y][x] = Mnm(Rings_3d[z_str_id][y_str_id][x_str_id], Rings_3d[z_col_id][y_col_id][x_col_id],DATA) * a1/a * mu_0
    print('circvecs:Done')
    i_vecs = np.zeros((nz, ny, nx), dtype=complex)

    # def LO(I):
    #     return dot(I, Z_circvecs=Z_circvecs, i_vecs=i_vecs, Z_0 = Z_0, N = N)


    M = Matrix(Rings, Rings, DATA) * a1/a * mu_0
    print('Matrix: Done')
    Zdiag = np.ones((N['z'], N['y'], N['x']), dtype = complex)
    omega_0 = 1/sqrt(L * C)
    print('Straight solving')
    for omega in Omega:
        for z in range(N['z']):
            for y in range(N['y']):
                for x in range(N['x']):
                    gradx = x * grad[0] / N['x']
                    grady = y * grad[1] / N['y']
                    gradz = z * grad[2] / N['z']
                    omegagrad = omega_0 * (1 + gradz + grady + gradx)
                    Cgrad = 1/L/omegagrad ** 2
                    Zdiag[z][y][x] = R + 1j * omega * L + 1/(1j * omega * Cgrad)
        Z_0 = np.diag(Zdiag.reshape(Number))
        E = H_0z  * np.ones(Number) * 1j * omega * pi * Radius1 ** 2 * mu_0
        I = solve(M*1j*omega + Z_0, E)
        Imax = I.max()
        Imin = I.min()
        Volume = N['x'] * N['y'] * (N['z']-1) * a1 * b1 * c1
        polarisation = pi * Radius1 ** 2 * I.sum()/Volume/H_0z/1j
        
        print(f'Solving: Done - {round(omega/2/pi/10 ** 6, 1)}MHz, grad = {grad}')
        CURRENTS.append(I)
        MaxCurrents.append(Imax)
        MinCurrents.append(Imin)
        Polarisation.append(polarisation)


    Data['N'] = N
    Data['Omega'] = list(Omega)
    Data['grad'] = grad
    Data['params'] = [a, a1, b, b1, c, c1, w, L, C, R]
    Data['Matrix'] = [list(m) for m in M]
    Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
    Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]
    Data['RealMaxCurrents'] = list(np.real(MaxCurrents))
    Data['ImagMaxCurrents'] = list(np.imag(MaxCurrents))
    Data['RealMinCurrents'] = list(np.real(MinCurrents))
    Data['ImagMinCurrents'] = list(np.imag(MinCurrents))
    Data['RealPolarisation'] = list(np.real(Polarisation))
    Data['ImagPolarisation'] = list(np.imag(Polarisation))

    return Data   
    

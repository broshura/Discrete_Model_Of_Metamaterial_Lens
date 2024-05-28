# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
from scipy.linalg import solve
from Impedance_matrix import Matrix
from tqdm import tqdm

def solvesystem(rings, M_0, Omega, phi_0z = 1, Inductance = {}, find = 'Currents'):    
    # Unpacking parameters
    # Solve system in currents terms and return currents in each ring
    if find == 'Currents':
        CURRENTS = []

        print('Matrix forming')
        Number = len(rings)
        M = Matrix(rings, Data = Inductance)
        print('Matrix: Done')

        print('Straight solving')

        # External field
        Phi_0z = np.ones(Number)*phi_0z/np.max(abs(phi_0z))
        for omega in tqdm(Omega):
            # Solve equation (diag(Z_0)/jw - M)I = Phi_0z
            I = solve(np.diag(M_0(omega)) - M, Phi_0z)
            CURRENTS.append(I * np.max(abs(phi_0z)))
        print('Straight solving: Done')
        Data = {}

        Data['Omega'] = list(Omega)
        Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
        Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]
    elif find == 'Voltage':
        sigma_0 = lambda Omega: 1/(M_0(Omega) * 1j * omega)
        Currents = []
        print('Matrix forming')
        Number = len(rings)
        M = Matrix(rings, Data = Inductance)
        print('Matrix: Done')

        print('Straight solving')
        # Solve equation 
        for omega in tqdm(Omega):
            V = solve(np.ones(Number)/(1j * omega) - M_0@np.diag(sigma_0(omega)), Phi_0z)
            CURRENTS.append(V/sigma_0(omega))
        print('Straight solving: Done')
        Data = {}

        Data['Omega'] = list(Omega)
        Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
        Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]

    return Data

def effective_mu(Params, frequency = False):
    mu_0 = 4 * np.pi * 10 ** (-7)
    r = Params['Radius']
    R = Params['R']
    C = Params['C']
    L = Params['L']
    a = Params['Dz']
    b = Params['Dy']
    c = Params['Dx']
    Sigma = Params['Sigma']

    Z = lambda Omega : R - 1j * Omega * L + 1j/(Omega * C) - 1j * Omega * mu_0 * r * Sigma
    Const =  lambda Omega: 1j * Omega * mu_0 * np.pi ** 2 * r ** 4 /(a*b*c)
    if frequency:
        return lambda w: (Z(w) + 2/3 * Const(w))/(Z(w) - 1/3 * Const(w)) 
    return lambda w: (Z(w) + 2/3 * Const(w))/(Z(w) - 1/3 * Const(w))
def spherical_chi(mu):
    return 3 * (mu - 1)/(mu + 2)

def disk_chi(mu):
    return 1 - 1/mu

def needle_chi(mu):
    return mu - 1

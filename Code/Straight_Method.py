# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
from scipy.linalg import solve
from Impedance_matrix import Matrix
from tqdm import tqdm

def solvesystem(rings, M_0, Omega, phi_0z = 1, Inductance = {}):    
    # Unpacking parameters

    CURRENTS = []

    print('Matrix forming')
    Number = len(rings)
    M = Matrix(rings, Data = Inductance)
    print('Matrix: Done')

    print('Straight solving')

    Phi_0z = np.ones(Number)*phi_0z/np.max(abs(phi_0z))
    for omega in tqdm(Omega):
        I = solve(np.diag(M_0(omega)) - M, Phi_0z)
        CURRENTS.append(I * np.max(abs(phi_0z)))
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

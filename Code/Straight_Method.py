# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
from scipy.linalg import solve
from Impedance_matrix import Matrix
from tqdm import tqdm

def solvesystem(Params, rings_4d, phi_0z_4d, Inductance = {}, find = 'Voltage', tol = 0):
    Params['Solver_type'] = 'Straight'
    Omegas = Params['Omega']    
    Omega = np.linspace(Omegas[0], Omegas[1], Omegas[2])
    # Unpacking parameters
    # Solve system in currents terms and return currents in each ring
    rings = sum([rings_4d[orientation] for orientation in rings_4d], [])
    phi_0z = np.array(sum([phi_0z_4d[orientation] for orientation in phi_0z_4d], []))

    Number = len(rings)

    # Diagonal part
    L, C, R = [], [], []
    for ring in rings:
        L.append(ring.L)
        C.append(ring.C)
        R.append(ring.R)
    L, C, R = np.array(L), np.array(C), np.array(R)

    M_0 = lambda Omega: (R - 1j * Omega * L + 1j/(Omega * C))/1j/Omega
    P = []
    CURRENTS = []
    # External field
    Phi_0z = phi_0z/np.max(abs(phi_0z))

    print('Matrix forming')
    M = Matrix(rings, Data = Inductance)
    print('Matrix: Done')
    if find == 'Currents':
        print('Straight solving (Currents)')
        for omega in tqdm(Omega):
            # Solve equation (diag(Z_0)/jw - M)I = Phi_0z
            I = solve(np.diag(M_0(omega)) - M, Phi_0z)
            CURRENTS.append(I * np.max(abs(phi_0z)))
            start = 0
            p = []
            for pos in Params['Orientations']:
                end = start + Params['Numbers'][pos]
                p.append(np.sum(I[start:end])/(end-start))
                start = end
            P.append(p)
        P = np.array(P)*np.max(abs(phi_0z)) * Params['P_0z']
                

        print('Straight solving (Currents): Done')

    elif find == 'Voltage':
        print('Straight solving (Voltage)')
        # Solve equation 
        for omega in tqdm(Omega):
            M_diag = M_0(omega)
            # Solve equation (1/jw - M/M_diag)I = Phi_0z/M_diag
            I = solve(np.eye(Number) - np.diag(1/M_diag)@M, Phi_0z/M_diag)
            CURRENTS.append(I)
            start = 0
            p = []
            for pos in Params['Orientations']:
                end = start + Params['Numbers'][pos]
                p.append(np.sum(I[start:end])/(end-start))
                start = end
            P.append(p)
        P = np.array(P)*np.max(abs(phi_0z)) * Params['P_0z']
                
        print('Straight solving (Voltage): Done')
        
    Data = {}
    Data['Params'] = Params
    Data['Omega'] = Omega
    Data['Currents'] = CURRENTS
    Data['Polarization'] = P
    Data['Phi_0z'] = list(phi_0z)
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

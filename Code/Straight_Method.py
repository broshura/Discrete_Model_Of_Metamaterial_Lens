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
    M = Matrix(rings, rings, Data = Inductance)
    print('Matrix: Done')

    print('Straight solving')

    Phi_0z = np.ones(Number)
    for omega in tqdm(Omega):
        I = solve(M + M_0(omega), Phi_0z)
        CURRENTS.append(I * phi_0z)
    print('Straight solving: Done')
    Data = {}

    Data['Omega'] = list(Omega)
    Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
    Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]
 
    return Data
    

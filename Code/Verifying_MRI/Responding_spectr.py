# Calculate responding impedance for MRI Lenz and coil
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Responding_solving import solvesystem
from Parameters_MRI import Params, Radius, N, Responded_ring, R0, L, C, R, R_coil, C_coil, L_coil
from Geometry import Rectangle_packing
from Impedance_matrix import M_0, Mnm
import numpy as np

# Calculating currents in each ring for anizotropic system on straight way

import numpy as np
from scipy.linalg import solve
from Impedance_matrix import Matrix
from tqdm import tqdm

def straight_solvesystem(rings, M_0, Omega, phi_0z = 1, Inductance = {}):    
    # Unpacking parameters

    CURRENTS = []

    print('Matrix forming')
    Number = len(rings)
    M = Matrix(rings, Data = Inductance)
    print('Matrix: Done')

    print('Straight solving')

    Phi_0z = np.zeros(Number)
    Phi_0z[-1] = 1
    for omega in tqdm(Omega):
        I = solve(M + np.diag(M_0(omega)), Phi_0z/1j/omega)
        CURRENTS.append(I * phi_0z)
    print('Straight solving: Done')
    Data = {}

    Data['Omega'] = list(Omega)
    Data['RealCurrents'] = [list(np.real(i).reshape(Number)) for i in CURRENTS]
    Data['ImagCurrents'] = [list(np.imag(i).reshape(Number)) for i in CURRENTS]

    return Data

import json

Inductance = {}
Omega =  np.linspace(61, 67, 1000) * 2 * np.pi * 1e6
H_0z = 1
mu_0 = 4 * np.pi * 1e-7
# Calculate polarisation for different structures

Dz, Dy, Dx = Params['Dz'], Params['Dy'], Params['Dx']

# Initial position of the first ring for each orientation


Rings_4d ={}
for pos in Params['N']:
    Rings_4d[pos] = Rectangle_packing(Params, r0 = R0[pos], orientation=pos).reshape(
                                                 Params['N'][pos]['nz'],
                                                 Params['N'][pos]['ny'],
                                                 Params['N'][pos]['nx'])

# Rings = np.concatenate([np.ravel(Rings_4d[pos]) for pos in Params['N']])

# Rings = np.append(Rings, Responded_ring)
# print(Mnm(Rings[-1], Rings[-2])/mu_0 /Dz * 2)
# print(Mnm(Rings[-1], Rings[-3], Inductance)/mu_0 /Dz * 2)
# print(Mnm(Rings[0], Rings[1], Inductance)/mu_0 /Dz * 2)
# print(Mnm(Rings[0], Rings[2], Inductance)/mu_0 /Dz * 2)


# Currents = straight_solvesystem(Rings, M_0(Params), Omega, phi_0z=H_0z * mu_0 * np.pi * Radius ** 2, Inductance = Inductance)
Currents = solvesystem(Rings_4d, Responded_ring, M_0(Params), Omega, phi_0z=H_0z * mu_0*np.pi * Radius ** 2, Inductance = Inductance, tol = 1e-5)
Nz, Ny, Nx = [N[pos][f'n{pos}'] for pos in N]
Volume =  Nx * Ny * Nz * Dz * Dy * Dx
P_0z = np.pi * Radius ** 2 /Volume/H_0z
Currents['Shape'] = N
Currents['RealPolarisation'] = [P_0z * sum(i) for i in Currents['RealCurrents']]
Currents['ImagPolarisation'] = [P_0z * sum(i) for i in Currents['ImagCurrents']]
dims = ''.join([pos for pos in N])
# with open(f"/Users/shuramakarenko/LocalDocs/Discrete_Model_Of_Metamaterial_Lens/Code/DATA/Form_Data/izotropic/MRICurrentsS", "w") as f:
#     json.dump(Currents, f)
exclude = ['RealCurrents', 'ImagCurrents']
data = {key: value for key, value in Currents.items() if key not in exclude}
# with open(f"/Users/shuramakarenko/LocalDocs/Discrete_Model_Of_Metamaterial_Lens/Code/DATA/Form_Data/izotropic/MRIDataS", "w") as f:
#     json.dump(data, f)
# Coil_M = np.array([Mnm(ring, Responded_ring, Inductance) for ring in Rings[:-1]])
# Responded_impedance = [Coil_M @ I[:-1] /I[-1] * 1j * omega for omega, I in zip(Omega, np.array(Currents['RealCurrents']) + 1j * np.array(Currents['ImagCurrents']))]
# data['Responded_impedance'] = [list(np.real(Responded_impedance)), list(np.imag(Responded_impedance))]
plt.title('Responding impedance ')
plt.figure(figsize=(12, 5))
plt.subplot(121)
plt.plot(Omega/2/np.pi/1e6, data['Responded_impedance'][0])
plt.subplot(122)
plt.plot(Omega/2/np.pi/1e6, data['Responded_impedance'][1])
plt.savefig('/Users/shuramakarenko/LocalDocs/Discrete_Model_Of_Metamaterial_Lens/Code/Verifying_MRI/Responding_impedance.png')
#plt.show()

# Calculate responding impedance for MRI Lenz and coil
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Fast_Method import solvesystem
from Parameters_MRI import Params, Radius, N
from Geometry import Rectangle_packing
from Impedance_matrix import Mnm, Matrix, M_0
import numpy as np

import json

Inductance = {}
Omega =  np.linspace(61, 67, 1000) * 2 * np.pi * 1e6
H_0z = 1
mu_0 = 4 * np.pi * 1e-7
# Calculate polarisation for different structures

Dz, Dy, Dx = Params['Dz'], Params['Dy'], Params['Dx']
Radius, 
# Initial position of the first ring for each orientation
R0 = {
    'z': {'nz': 0, 'ny': Dy/2, 'nx': Dx/2},
    'y': {'nz': Dz/2, 'ny': Dy, 'nx': Dx/2},
    'x': {'nz': Dz/2, 'ny': Dy/2, 'nx': Dx}
}

Rings_4d ={}
for pos in N:
    Rings_4d[pos] = Rectangle_packing(Params, r0 = R0[pos], orientation=pos).reshape(
                                                 N[pos]['nz'],
                                                 N[pos]['ny'],
                                                 N[pos]['nx'])
Currents = solvesystem(Rings_4d, M_0(Params), Omega, phi_0z=H_0z * mu_0*np.pi * Radius ** 2, Inductance = Inductance, tol = 1e-2)
Nz, Ny, Nx = [N[pos][f'n{pos}'] for pos in N]
Volume =  Nx * Ny * Nz * Dz * Dy * Dx
P_0z = np.pi * Radius ** 2 /Volume/H_0z
Currents['Shape'] = N
Currents['RealPolarisation'] = [P_0z * sum(i) for i in Currents['RealCurrents']]
Currents['ImagPolarisation'] = [P_0z * sum(i) for i in Currents['ImagCurrents']]
dims = ''.join([pos for pos in N])
with open(f"../Data/Form_Data/izotropic/MRICurrents", "w") as f:
    json.dump(Currents, f)
exclude = ['RealCurrents', 'ImagCurrents']
data = {key: value for key, value in Currents.items() if key not in exclude}
with open(f"Data/Form_Data/izotropic/MRIData", "w") as f:
    json.dump(data, f)
import os 
import json 
import numpy as np
import matplotlib.pyplot as plt

from Parameters import *
from Geometry import *
from Straight_Method import solvesystem as straight_solvesystem
from Fast_Method import solvesystem as fast_solvesystem

Solvers = {
    'Straight': straight_solvesystem,
    'Fast': fast_solvesystem
}

def save(filename, Params):
    solver = Solvers[Params['Solver_type']]
    packing = Packings[Params['Packing']]

    rings_4d = packing(Params)
    phi_0z_4d = {
        orientation: list(np.ones(Params['Numbers'][orientation]
                             ) * (orientation == 'z'
                                  ) * mu_0*np.pi * Radius ** 2
                                  ) for orientation in Params['Orientations']
        }
    print('Количество колец:', Params['Numbers'])
    Data = solver(Params, rings_4d, phi_0z_4d)

    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    os.makedirs(f'./{filename}/{name}', exist_ok=True)

    # Saving modeling parameters in readable format
    with open(f'./{filename}/{name}/Params.json', 'w') as f:
        json.dump(Params, f)

    # Save full calculated data in npz format
    calc_data = {
        'Currents': Data['Currents'],
        'Omega': Data['Omega'],
        'Polarization': Data['Polarization'],
        'Phi_0z': Data['Phi_0z'],
    }
    np.savez(f'./{filename}/{name}/Currents.npz', **calc_data)
    
    # Save neccesary data for plotting in npz format
    pol_data = {
        'Polarization': Data['Polarization'],
        'Omega': Data['Omega'],
        'Phi_0z': Data['Phi_0z'],
    }
    np.savez(f'./{filename}/{name}/Polarization.npz', **pol_data)

def open_model(filename, Params, Currents = 'False', Polarization = 'True'):
    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    data = {}
    with open(f'./{filename}/{name}/Params.json', 'r') as f:
        data['Params'] = json.load(f)
    
    if Currents == 'True':
        with np.load(f'./{filename}/{name}/Currents.npz') as f:
            data['Currents'] = f['Currents']
            data['Omega'] = f['Omega']
            data['Polarization'] = f['Polarization']
            data['Phi_0z'] = f['Phi_0z']
    elif Polarization == 'True':
        with np.load(f'./{filename}/{name}/Polarization.npz') as f:
            data['Polarization'] = f['Polarization']
            data['Omega'] = f['Omega']
            data['Phi_0z'] = f['Phi_0z']
    return data


Params['Solver_type'] = 'Fast'
Params['Threads'] = 1

for n in range(1, 101):
    Params['N'], Params['shape'] = to3D(n, 1, 1, 'zyx')
    save('DATA', Params)
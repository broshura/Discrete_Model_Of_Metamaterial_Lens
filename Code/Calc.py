# Data functions
import os
import sys 
import json 
import numpy as np

from Parameters import mu_0, Radius
from Geometry import *
from Straight_Method import solvesystem as straight_solvesystem
from Fast_Method import solvesystem as fast_solvesystem

# Allowed types of solvers
Solvers = {
    'Straight': straight_solvesystem,
    'Fast': fast_solvesystem
}

def save(filename:str, Params:dict)->None:
    """Function to save calculated data in npz format
    in DATA or DATALowTol folder (configurable by filename)

    Parameters
    ----------
    filename : str
        full path to the folder where data will be saved
    Params : dict
        dictionary with modeling parameters
    """    
    solver = Solvers[Params['Solver_type']]
    packing = Packings[Params['Packing']]
    
    rings_4d = packing(Params, Fill=True)
    phi_0z_4d = {
        orientation: list(np.ones(Params['Numbers'][orientation]
                             ) * (orientation == 'z'
                                  ) * mu_0*np.pi * Radius ** 2
                                  ) for orientation in Params['Orientations']
        }
    print('Количество колец:', Params['Numbers'])
    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    print(name)
    Data = solver(Params, rings_4d, phi_0z_4d, tol = 1e-3)
    os.makedirs(f'./{filename}/{name}', exist_ok=True)

    # Saving modeling parameters in readable format
    with open(f'./{filename}/{name}/Params.json', 'w') as f:
        json.dump(Params, f)
    
    # Save full calculated data in npz format
    if np.array(Data['Currents']).nbytes < 1024**3 * 5:  # 5 GB limit
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

def open_model(filename:str, Params:dict, Currents:bool = False, Polarization:bool = True)->dict:
    """Function to open saved data in npz format

    Parameters
    ----------
    filename : str
        full path to the folder where data is saved
    Params : dict
        dictionary with modeling parameters
    Currents : bool, optional
        Get Currents only if it necessary, by default False
    Polarization : bool, optional
        Get Polarization only if it necessary, by default True

    Returns
    -------
    dict
        dictionary with modeling parameters and calculated data
    """    
    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    data = {}
    with open(f'./{filename}/{name}/Params.json', 'r') as f:
        data['Params'] = json.load(f)
    
    if Currents == True:
        with np.load(f'./{filename}/{name}/Currents.npz') as f:
            data['Currents'] = f['Currents']
            data['Omega'] = f['Omega']
            data['Polarization'] = f['Polarization']
            data['Phi_0z'] = f['Phi_0z']
    elif Polarization == True:
        with np.load(f'./{filename}/{name}/Polarization.npz') as f:
            data['Polarization'] = f['Polarization']
            data['Omega'] = f['Omega']
            data['Phi_0z'] = f['Phi_0z']
    return data

# Modeling for exact strucutres.
if __name__ == '__main__':

    from Parameters import Params

    Params['Packing'] = 'Rectangle'
    Params['Solver_type'] = 'Fast'
    Params['Solver_name'] = 'lgmres'
    
    for n in [5]:
        Params['N'], Params['shape'] = to3D(12, 1, 1, 'zyx')
        save('DATALowTol', Params)

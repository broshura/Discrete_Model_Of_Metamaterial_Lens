# Data functions
import os
import sys 
import json 
import numpy as np

from Parameters import *
from Geometry import *
from Straight_Method import solvesystem as straight_solvesystem
from Fast_Method import solvesystem as fast_solvesystem

# Allowed types of solvers
Solvers = {
    'Straight': straight_solvesystem,
    'Fast': fast_solvesystem
}

def z_coord(rings_4d, orientation):
    """Extract only z-coordinates of rings for given orientation"""
    z_coords = []
    rings = rings_4d[orientation]
    # Iterate through 3D array of rings
    for i in range(rings.shape[0]):
        for j in range(rings.shape[1]):
            for k in range(rings.shape[2]):
                z_coords.append(rings[i,j,k].z)
    return np.array(z_coords)

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

    Omega_range = Params['Omega']  # [omega_min, omega_max, num_points]
    Omega = np.linspace(Omega_range[0], Omega_range[1], Omega_range[2])
    

    if Params['Scattering'] == 'Mie_False':

        phi_0z_calc = lambda omega:{
            orientation: list(np.ones(Params['Numbers'][orientation]
                                ) * (orientation == 'z'
                                    ) * mu_0 * np.pi * Radius ** 2 * H_0z
                                    ) for orientation in Params['Orientations']
            }
        
    elif Params['Scattering'] == 'Mie_True':

        phi_0z_calc = lambda omega: {
            orientation: list(np.ones(Params['Numbers'][orientation]) * 
                            (orientation == 'z') * 
                            mu_0 * np.pi * Radius**2 * H_0z * 
                            np.exp(-1j * omega/3e8 * z_coord(rings_4d, orientation)))
            for orientation in Params['Orientations']
        }
    
    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}_{Params["Scattering"]}'
    print(name)
    Data = solver(Omega, Params, rings_4d, phi_0z_calc, tol = 1e-3)
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

def open_model(filename:str, Params:dict, Currents:bool = True, Polarization:bool = True)->dict:
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
    name = f'{Params["Packing"]}_NoGrad_{Params["shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}_{Params["Scattering"]}'
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
        Params['N'], Params['shape'] = to3D(3, 3, 3, 'zyx')
        save('DATALowTol', Params)

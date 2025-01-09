# Data functions
import os
import shutil
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
    print('Number of rings:', Params['Numbers'])
    name = f'{Params["Packing"]}_NoGrad_{Params["Shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    print(name)
    Data = solver(Params, rings_4d, phi_0z_4d, tol = Params['tol'])
    os.makedirs(f'./{filename}/{name}', exist_ok=True)

    # Saving modeling parameters in readable format
    with open(f'./{filename}/{name}/Params.json', 'w') as f:
        json.dump(Params, f)
    
    # Save full calculated data in npz format
    if np.array(Data['Currents']).nbytes < 1024**3 * 5:  # 5 GB limit
        print('Saving size:', np.array(Data['Currents']).nbytes/1024**2, 'MB')
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

    '''
    Configure solver parameters

    Params['Packing']: base structure is parallelepiped 
    with Nz x Ny x Nx cells, then for removing other rings for 
    ellipse and cylinder along chosen axes.
    Possible values: 
    'Rectangle', 'Ellipse', 'Cylinder-z,
    'Cylinder-y', 'Cylinder-x'

    Params['Solver_type']: the way to solve matrix equation
    'Fast' for iterative solver with Params['tol'] absolutre tolerance
    Has linear dependency of total number of rings for memory and 
    NlogN for solving

    'Straight' for straight solving with machine exactness.
    Takes N^2 memory and N^3 for solving equatiion 

    Params['solver_name'] is for different iterative methods, lgmres 
    is highly recommened for best time and stability results.

    Params['N'] and Params['shape'] needed to create numbers for each 
    orientation along each axis and remember shape.
    Params['Orientations'] describes rings with wheech orientation 
    will be in structure (so is there isotropic or anisotropic structure)
    Params['Type'] describes way to create metacell and border conditions 
    (is there border layers or not, for more check Rings_visualize.ipynb)
    '''

    Params['Packing'] = 'Rectangle'
    Params['Solver_type'] = 'Fast'
    Params['Solver_name'] = 'lgmres'
    Params['Tol'] = 1e-5
    Params['Type'] = 'border'
    Params['Orientations'] = 'zyx'

    for n in [7]:
        Params['N'], Params['Shape'] = to3D(n, n, n,
                                            Params['Orientations'],
                                            Params['Type'])
        save(f'DATA_{Params["Type"]}', Params)

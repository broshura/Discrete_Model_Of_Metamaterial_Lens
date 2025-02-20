# Data functions
'''
Please, before edditing this file run this command in shell

git update-index --skip-worktree Code/Calc.py

It's necessary to avoid adding changes, based on 
calculations for different structures

If there is neccessary changes in open of save functions,
you could run this command in shell:

git update-index --no-skip-worktree Code/Calc.py

'''
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
    print("Start Modeling")
    print('Number of rings:', Params['Numbers'])
    name = f'{Params["Packing"]}_NoGrad_{Params["Shape"]}_{Params["Orientations"]}_{Params["Solver_type"]}'
    print(name, '\n')
    Data = solver(Params, rings_4d, phi_0z_4d, tol = Params['Tol'])
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
    size = sum([np.array(calc_data[key]).nbytes for key in calc_data])
    Degrees = ['B', 'kB', 'MB', 'GB', 'TB']
    degree = min(int(np.log(size)/np.log(1024)), 4)
    print('\nSaving currents, size:', round(size/1024**degree, 1), f'{Degrees[degree]}')
    np.savez(f'./{filename}/{name}/Currents.npz', **calc_data)
    
    # Save neccesary data for plotting in npz format
    pol_data = {
        'Polarization': Data['Polarization'],
        'Omega': Data['Omega'],
        'Phi_0z': Data['Phi_0z'],
    }
    # Adding sliced currents dict to pol_data
    if Params['IsSlices']:
        pol_data.update(Data['SlicedCurrents'])
    size = sum([np.array(pol_data[key]).nbytes for key in pol_data])
    Degrees = ['B', 'kB', 'MB', 'GB', 'TB']
    degree = min(int(np.log(size)/np.log(1024)), 4)
    print('Saving polarization, size:', round(size/1024**degree, 1), f'{Degrees[degree]}\n')
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
            data.update({key: f[key] for key in f})
    elif Polarization == True:
        with np.load(f'./{filename}/{name}/Polarization.npz') as f:
            data.update({key: f[key] for key in f})
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

    Params['MemLim'] are necessary to avoid extra-large current files 
    for modeling with more then 100.000 rings. You shold be shure that there 
    is enough place for saving. 
    Possible values:
    string type: 'All_free' - no limit (except possible)
    int type: nummber of bytes for limit, if zero - no Currents will be saved 
    even in process of modeling.
    '''

    Params['Packing'] = 'Rectangle'
    Params['Solver_type'] = 'Fast'
    Params['Solver_name'] = 'lgmres'
    Params['Tol'] = 1e-5
    Params['Type'] = 'border'
    Params['Orientations'] = 'zyx'
    Params['MemLim'] = 1024 ** 2 * 5 # 5 Gb limit
    Params['IsSlices'] = True

    # Example way to use
    for n in [3, 5, 7]:
        Params['N'], Params['Shape'] = to3D(n, n, n,
                                            Params['Orientations'],
                                            Params['Type'])
        Params['Slices'] = {
            'MiddleZZ': {
                'z': {'nz': [Params['N']['z']['nz']//2, Params['N']['z']['nz']//2 +1],
                      'ny': [0, Params['N']['z']['ny']],
                      'nx': [0, Params['N']['z']['nx']]
                }
            },
            'MiddleZY': {
                'z': {
                    'nz': [0, Params['N']['z']['nz']],
                    'ny': [Params['N']['z']['ny']//2, Params['N']['z']['ny']//2 +1],
                    'nx': [0, Params['N']['z']['nx']]
                }
            },
            'BottomZZ': {
                'z': {'nz': [0, 1],
                      'ny': [0, Params['N']['z']['ny']],
                      'nx': [0, Params['N']['z']['nx']]
                }
            },
            'BottomZY': {
                'z': {
                    'nz': [0, Params['N']['z']['nz']],
                    'ny': [0, 1],
                    'nx': [0, Params['N']['z']['nx']]
                }
            }
        }
        
        save(f'DATA_{Params["Type"]}', Params)

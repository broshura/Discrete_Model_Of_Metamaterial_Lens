from Ring_Class import Ring
from typing import List, Dict, Tuple
import numpy as np
eps = np.finfo(float).eps

def to3D(Nz:int, Ny:int, Nx:int, orientations:str = 'z', Type:str = 'border') -> Tuple[Dict[str, Dict[str, int]], str]:
    """Turn cell sizes into 3D sizes depending on the type of the cell (based on orientations)

    Parameters
    ----------
    Nz : int
        Number of cells in z direction
    Ny : int
        Number of cells in y direction
    Nx : int
        Number of cells in x direction
    orientations : str, optional
        Direction of normal vector for ring, by default 'z'
    Type : str, optional
        Could be 'border' if rings on border flat needed in modeling
        else 'open' and create "sharp" structure., by default 'border'

    Returns
    -------
    tuple = (dict, str)
        return dictionary with keys as orientations and
          values as dictionaries with keys 'nz', 'ny', 'nx' 
          and string with shape of the system
    """    
    N = {}
    shape = f'{Nz}x{Ny}x{Nx}'
    if Type == 'border':
        for orientation in orientations:
            N[orientation] = {
                'nz': Nz + (orientation == 'z')*(len(orientations) != 1),
                'ny': Ny + (orientation == 'y')*(len(orientations) != 1),
                'nx': Nx + (orientation == 'x')*(len(orientations) != 1)
            }
    elif Type == 'open':
        for orientation in orientations:
            N[orientation] = {
                'nz': Nz - (orientation == 'z')*(len(orientations) != 1),
                'ny': Ny - (orientation == 'y')*(len(orientations) != 1),
                'nx': Nx - (orientation == 'x')*(len(orientations) != 1)
            }
    return N, shape

def Rectangle_packing(Params:dict, Fill:bool = False) -> Dict[str, List[Ring]]:
    """Returns dict for paralellepiped packing of rings in correct
    order with respect to the orientation of the system (zyx)

    Parameters
    ----------
    Params : dict
        Dictionary with parameters of the system
    Fill : bool, optional
        Only to unify with other packing functions, by default False

    Returns
    -------
    Dict[str, List[Ring]]
        Dictionary with keys as orientations and values as lists of rings
    """     
    # Choosing the type of returned list and stracture of the system
    Params['Packing'] = 'Rectangle'
    orientations = Params['Orientations']
    Rings = {}
    Numbers = {}
    for orientation in orientations:
        rings = []
        N = Params['N'][orientation]
        L, C, R, w = Params['L'], Params['C'], Params['R'], Params['W']
        r, delta_x, delta_y, delta_z = Params['Radius'], Params['Dx'], Params['Dy'], Params['Dz'] 
        
        # Starting point of the first ring for each orientation
        r0 = {
            'nz': delta_z/2 * (1-(orientation == 'z')), 
            'ny': delta_y/2 * (1-(orientation == 'y')),
            'nx': delta_x/2 * (1-(orientation == 'x'))
        }


        # Case with anisotropic system and shifted layers

        shift_x = Params['shift_x']
        shift_y = Params['shift_y']
        shift_z = Params['shift_z']

        nz, ny, nx = N['nz'], N['ny'], N['nx']
        for k in range(nz):
            for j in range(ny): 
                for i in range(nx):
                    # Shift preventing rings from getting out of the borders
                    Shift_x = shift_x * (k * (orientation == 'z') + j * (orientation == 'y')) % delta_x
                    Shift_y = shift_y * (k * (orientation == 'z') + i * (orientation == 'x')) % delta_y
                    Shift_z = shift_z * (j * (orientation == 'y') + i * (orientation == 'x')) % delta_z
                    rings.append(
                        Ring(
                             i * delta_x + Shift_x  + r0['nx'],
                             j * delta_y + Shift_y  + r0['ny'],
                             k * delta_z + Shift_z  + r0['nz'],
                            orientation,
                            r, w, L, C, R)
                    )
        Numbers[orientation] = len(rings)
        Rings[orientation] = rings
    Params['Numbers'] = Numbers
    return Rings

def Ellipse_packing(Params:dict, Fill:bool = False) -> Dict[str, List[Ring]]:
    """Returns dict for ellipsoidal packing of rings in correct, where
    extra rings are removed from the system for straight solver and 
    kept but with infinite resistance (zero conductivity) for iterative solver

    Parameters
    ----------
    Params : dict
        Dictionary with parameters of the system
    Fill : bool, optional
        Remove or fill with infinite-resitance rings, by default False

    Returns
    -------
    Dict[str, List[Ring]]
        Dictionary with keys as orientations and values as lists of rings
    """    
    orientations = Params['Orientations']
    Rings = Rectangle_packing(Params, orientations)
    Params['Packing'] = 'Ellipse'

    dz, dy, dx = Params['Dz'], Params['Dy'], Params['Dx']
    Nz, Ny, Nx = Params['N']['z']['nz'], Params['N']['y']['ny'], Params['N']['x']['nx']
    R_z, R_y, R_x = (Nz-1)*dz/2, (Ny-1)*dy/2, (Nx-1)*dx/2

    for orientation in orientations:
        for Ring in Rings[orientation][:]:
            # Finding middle of the cell position
            r_x = Ring.x 
            r_y = Ring.y 
            r_z = Ring.z 

            distance = (r_x-R_x) ** 2/R_x**2 + (r_y - R_y) ** 2/R_y **2 + (r_z-R_z) ** 2/R_z ** 2
            if distance > 1.00:
                if Fill:
                    Ring.R = 1e200
                    Ring.C = 1e200
                    Ring.L = 1e200
                else:
                    Rings[orientation].remove(Ring)
    Params['Numbers'] = {orientation: len(Rings[orientation]) for orientation in orientations}
    return Rings

def Cylinder_packing(Params:dict, Fill:bool = False, axis:str = 'z') -> Dict[str, List[Ring]]:
    """Returns dict for cylindrical packing of rings in correct, where
    extra rings are removed from the system for straight solver and
    kept but with infinite resistance (zero conductivity) for iterative solver

    Parameters
    ----------
    Params : dict
        Dictionary with parameters of the system
    Fill : bool, optional
        Remove or fill with infinite-resitance rings , by default False
    axis : str, optional
        axis along which structure is symmetric, by default 'z'

    Returns
    -------
    Dict[str, List[Ring]]
        Dictionary with keys as orientations and values as lists of rings
    """    
    Params['Packing'] = f'Cylinder-{axis}'
    orientations = Params['Orientations']
    Rings = Rectangle_packing(Params, orientations)
    
    dz, dy, dx = Params['Dz'], Params['Dy'], Params['Dx']
    Nz, Ny, Nx = Params['N']['z']['nz'], Params['N']['y']['ny'], Params['N']['x']['nx']
    R_z, R_y, R_x = (Nz-1)*dz/2, (Ny-1)*dy/2, (Nx-1)*dx/2

    for orientation in orientations:
        for Ring in Rings[orientation][:]:
            # Finding middle of the cell position
            r_x = Ring.x 
            r_y = Ring.y 
            r_z = Ring.z 

            distance = (r_x-R_x) ** 2/R_x**2 * (
                axis != 'x') + (r_y - R_y) ** 2/R_y **2 * (
                axis != 'y') + (r_z-R_z) ** 2/R_z ** 2 * (
                axis != 'z')
            
            if distance > 1.00:
                if Fill:
                    Ring.R = 1e20
                    Ring.C = 1e20
                    Ring.L = 0
                else:
                    Rings[orientation].remove(Ring)
    return Rings

# Dictionary with all possible packings
Packings = {
    'Rectangle': Rectangle_packing,
    'Ellipse': Ellipse_packing,
    'Cylinder-z': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'z'),
    'Cylinder-y': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'y'),
    'Cylinder-x': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'x')
}


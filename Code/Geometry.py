from Ring_Class import Ring
import numpy as np
eps = np.finfo(float).eps

'''
Turn cell sizes into 3D sizes depending on the type of the cell (based on orientations)
'''
def to3D(Nz, Ny, Nx, orientations = 'z', Type = 'border'):
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

'''
Create a structure of the system: dictionary with keys as orientations and values as lists
of rings in this orientation
'''
def Rectangle_packing(Params, Fill = False):
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

def Ellipse_packing(Params, Fill = False):
    Params['Packing'] = 'Ellipse'
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

            distance = (r_x-R_x) ** 2/R_x**2 + (r_y - R_y) ** 2/R_y **2 + (r_z-R_z) ** 2/R_z ** 2
            if distance > 1.00:
                if Fill:
                    Ring.R = np.inf
                    Ring.C = np.inf
                    Ring.L = 0
                else:
                    Rings[orientation].remove(Ring)
    return Rings

def Cylinder_packing(Params, Fill = False, axis = 'z'):
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
                    Ring.R = np.inf
                    Ring.C = np.inf
                    Ring.L = 0
                else:
                    Rings[orientation].remove(Ring)
    return Rings

Packings = {
    'Rectangle': Rectangle_packing,
    'Ellipse': Ellipse_packing,
    'Cylinder-z': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'z'),
    'Cylinder-y': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'y'),
    'Cylinder-x': lambda Params, Fill = False :Cylinder_packing(Params, Fill, 'x')
}


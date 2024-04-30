from Ring_Class import Ring
import numpy as np
eps = np.finfo(float).eps

def Rectangle_packing(Params, r0 = False, orientation = 'z'):
    N = Params['N'][orientation]
    L, C, R, w = Params['L'], Params['C'], Params['R'], Params['W']
    r, delta_x, delta_y, delta_z = Params['Radius'], Params['Dx'], Params['Dy'], Params['Dz']

    if not r0:
        r0 = {
            'nz': delta_z/2 * (1-(orientation == 'z')), 
            'ny': delta_y/2 * (1-(orientation == 'y')),
            'nx': delta_x/2 * (1-(orientation == 'x'))
              }
    # Case with anisotropic system and shifted layers

    shift_x = Params['shift_x']
    shift_y = Params['shift_y']
    shift_z = Params['shift_z']
    rings = []
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
    return np.array(rings, dtype=Ring)

def Ellipse_packing(Params, Fill = False, r0 = False, orientation = 'z'):
    r_z, r_y, r_x = ((Params['N'][orientation][f'nz']-1)*Params['Dz'])/2, (Params['N'][orientation][f'ny']-1)*Params['Dy']/2, (Params['N'][orientation][f'nx']-1)*Params['Dx']/2

    if not r0:
        r0 = {
            'nz': Params['Dz']/2 * (1-(orientation == 'z')), 
            'ny': Params['Dy']/2 * (1-(orientation == 'y')),
            'nx': Params['Dx']/2 * (1-(orientation == 'x'))
              }

    Rings = Rectangle_packing(Params, r0, orientation).tolist()

    for Ring in Rings[:]:
        distance = (Ring.x-(r_x + r0['nx'])) ** 2/r_x**2 + (Ring.y-(r_y + r0['ny'])) ** 2/r_y **2 + (Ring.z-(r_z+r0['nz'])) ** 2/r_z ** 2
        if distance > 1.0:
            if Fill:
                Ring.R = 1e10 * Ring.R
            else:
                Rings.remove(Ring)
    return Rings

def Cilinder_packing(Params, Fill = False, r0 = False, orientation = 'z', axes = 'z'):
    r_z, r_y, r_x = ((Params['N'][orientation][f'nz']-1)*Params['Dz'])/2, (Params['N'][orientation][f'ny']-1)*Params['Dy']/2, (Params['N'][orientation][f'nx']-1)*Params['Dx']/2

    if not r0:
        r0 = {
            'nz': Params['Dz']/2 * (1-(orientation == 'z')), 
            'ny': Params['Dy']/2 * (1-(orientation == 'y')),
            'nx': Params['Dx']/2 * (1-(orientation == 'x'))
              }

    Rings = Rectangle_packing(Params, r0, orientation).tolist()

    for Ring in Rings[:]:
        distance = (axes !='x')*(Ring.x-(r_x + r0['nx'])) ** 2/r_x**2 + (axes != 'y') * (Ring.y-(r_y + r0['ny'])) ** 2/r_y **2 + (axes != 'z') * (Ring.z-(r_z+r0['nz'])) ** 2/r_z ** 2
        if distance > 1.1:
            if Fill:
                Ring.R = len(Rings) * Ring.R/eps
            else:
                Rings.remove(Ring)
    return Rings



from Ring_Class import Ring
import numpy as np


def Rectangle_packing(Params, r0 = 0, orientation = 'z'):
    N = Params['N'][orientation]
    L, C, R, w = Params['L'], Params['C'], Params['R'], Params['W']
    r, delta_x, delta_y, delta_z = Params['Radius'], Params['Dx'], Params['Dy'], Params['Dz']

    # Case with anisotropic system and shifted layers
    
    shift_x = Params['shift_x']
    shift_y = Params['shift_y']
    shift_z = Params['shift_z']
    rings = []
    nz, ny, nx = N['nz'], N['ny'], N['nx']
    for k in range(nz):
        for j in range(ny): 
            for i in range(nx):
                Shift_x = shift_x * (k * (orientation == 'z') + j * (orientation == 'y'))
                Shift_y = shift_y * (k * (orientation == 'z') + i * (orientation == 'x'))
                Shift_z = shift_z * (j * (orientation == 'y') + i * (orientation == 'x'))
                rings.append(
                    Ring(
                        # Prevent rings from getting out of the borders
                        (i * delta_x + Shift_x) % ((nx) * delta_x) + r0['nx'],
                        (j * delta_y + Shift_y) % ((ny) * delta_y) + r0['ny'],
                        (k * delta_z + Shift_z) % ((nz) * delta_z) + r0['nz'],
                        orientation,
                        r, w, L, C, R)
                )
    return np.array(rings, dtype=Ring)


# Coming soon
def Hexagonal_packing(Params):
    N = Params['N']
    L, C, R, w = Params['L'], Params['C'], Params['R'], Params['W']
    r, delta_x, delta_y, delta_z = Params['Radius'], Params['Dx'], Params['Dy'], Params['Dz']

    rings = []
    nz, ny, nx = N[orientation]['nz'], N[orientation]['ny'], N[orientation]['nx']
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                rings.append(
                    Ring(
                        # Prevent rin gs from getting out of the borders
                        (i * delta_x + shift_x * k + delta_x/2 * j) % ((nx) * delta_x),
                        (np.sqrt(3)/2*j * delta_x + shift_y * k) % ((ny) * delta_y),
                        k * delta_z,
                        orientation,
                        r, w, L, C, R)
                )
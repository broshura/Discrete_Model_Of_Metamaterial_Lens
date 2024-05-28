from Ring_Class import Ring
import numpy as np
eps = np.finfo(float).eps

def to3D(Nz = False, Ny = False, Nx = False, Number = False, orientations = 'z', Type = 'bordered'):
    N = {}
    if Type == 'bordered':
        for orientation in orientations:
            N[orientation] = {
                'nz': Nz + (orientation == 'z'),
                'ny': Ny + (orientation == 'y'),
                'nx': Nx + (orientation == 'x')
            }
    if Type == 'open':
        for orientation in orientations:
            N[orientation] = {
                'nz': Nz - (orientation == 'z'),
                'ny': Ny - (orientation == 'y'),
                'nx': Nx - (orientation == 'x') 
            }
    return N

def Rectangle_packing(Params, orientations = 'z', Type = 'border_straight'):
    # Choosing the type of returned list and stracture of the system
    border_type, list_type = Type.split('_')
    if list_type == 'straight':
        rings = []
    elif list_type == 'fast':
        Rings = {}
    
    for orientation in orientations:
        if list_type == 'fast':
            rings = []
        N = Params['N'][orientation]
        L, C, R, w = Params['L'], Params['C'], Params['R'], Params['W']
        r, delta_x, delta_y, delta_z = Params['Radius'], Params['Dx'], Params['Dy'], Params['Dz']

        if border_type == 'bordered':
            r0 = {
                'nz': delta_z/2 * (1-(orientation == 'z')), 
                'ny': delta_y/2 * (1-(orientation == 'y')),
                'nx': delta_x/2 * (1-(orientation == 'x'))
            }
        elif border_type == 'open':
            r0 = {
                'nz': delta_z * (orientation == 'z') + delta_z/2 * (1-(orientation == 'z')), 
                'ny': delta_y * (orientation == 'y') + delta_y/2 * (1-(orientation == 'y')),
                'nx': delta_x * (orientation == 'x') + delta_x/2 * (1-(orientation == 'x'))
            }
        else: 
            r0 = {
                'nz': 0,
                'ny': 0,
                'nx': 0 
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
        if list_type == 'fast':
            Rings[orientation] = rings
    if list_type == 'straight':
        return rings
    elif list_type == 'fast':
        return Rings

def Ellipse_packing(Params, Fill = False, orientations = 'z', Type = 'border_straight'):
    border_type, list_type = Type.split('_') 
    Rings = Rectangle_packing(Params, orientations, Type = Type)
    
    if list_type == 'straight':
        r_z0 = (max([Ring.z for Ring in Rings]) + min([Ring.z for Ring in Rings]))/2
        r_y0 = (max([Ring.y for Ring in Rings]) + min([Ring.y for Ring in Rings]))/2
        r_x0 = (max([Ring.x for Ring in Rings]) + min([Ring.x for Ring in Rings]))/2

        r_z = (max([Ring.z for Ring in Rings]) - min([Ring.z for Ring in Rings]))/2 + Params['Dz'] * (border_type == 'open') / 2
        r_y = (max([Ring.y for Ring in Rings]) - min([Ring.y for Ring in Rings]))/2 + Params['Dy'] * (border_type == 'open') / 2
        r_x = (max([Ring.x for Ring in Rings]) - min([Ring.x for Ring in Rings]))/2 + Params['Dx'] * (border_type == 'open') / 2

        for Ring in Rings[:]:
            distance = (Ring.x-r_x0) ** 2/r_x**2 + (Ring.y-r_y0) ** 2/r_y **2 + (Ring.z-r_z0) ** 2/r_z ** 2
            if distance > 1.01:
                if Fill:
                    Ring.R = np.inf
                    Ring.C = np.inf
                    Ring.L = 0
                else:
                    Rings.remove(Ring)
    elif list_type == 'fast':
        r_z0 = (max([max([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]))/2
        r_y0 = (max([max([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]))/2
        r_x0 = (max([max([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]))/2

        r_z = (max([max([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dz'] * (border_type == 'open') / 2
        r_y = (max([max([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dy'] * (border_type == 'open') / 2
        r_x = (max([max([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dx'] * (border_type == 'open') / 2

        for orientation in orientations:
            for Ring in Rings[orientation][:]:
                distance = (Ring.x-r_x0) ** 2/r_x**2 + (Ring.y-r_y0) ** 2/r_y **2 + (Ring.z-r_z0) ** 2/r_z ** 2
                if distance > 1.01:
                    if Fill:
                        Ring.R = np.inf
                        Ring.C = np.inf
                        Ring.L = 0
                    else:
                        Rings[orientation].remove(Ring)
    return Rings

def Cylinder_packing(Params, Fill = False, Type='border_straight', orientations = 'z', axes = 'z'):
    border_type, list_type = Type.split('_')
    Rings = Rectangle_packing(Params, orientations = orientation, Type = Type)
    
    if list_type == 'straight':
        r_z0 = (max([Ring.z for Ring in Rings]) + min([Ring.z for Ring in Rings]))/2
        r_y0 = (max([Ring.y for Ring in Rings]) + min([Ring.y for Ring in Rings]))/2
        r_x0 = (max([Ring.x for Ring in Rings]) + min([Ring.x for Ring in Rings]))/2

        r_z = (max([Ring.z for Ring in Rings]) - min([Ring.z for Ring in Rings]))/2 + Params['Dz'] * (border_type == 'open') / 2
        r_y = (max([Ring.y for Ring in Rings]) - min([Ring.y for Ring in Rings]))/2 + Params['Dy'] * (border_type == 'open') / 2
        r_x = (max([Ring.x for Ring in Rings]) - min([Ring.x for Ring in Rings]))/2 + Params['Dx'] * (border_type == 'open') / 2

        for Ring in Rings[:]:
            distance = (Ring.x-r_x0) ** 2/r_x**2 * (axes != 'x') + (Ring.y-r_y0) ** 2/r_y **2 * (axes != 'y') + (Ring.z-r_z0) ** 2/r_z ** 2 * (axes != 'z')
            if distance > 1.01:
                if Fill:
                    Ring.R = np.inf
                    Ring.C = np.inf
                    Ring.L = 0
                else:
                    Rings.remove(Ring)
    
    elif list_type == 'fast':
        r_z0 = (max([max([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]))/2
        r_y0 = (max([max([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]))/2
        r_x0 = (max([max([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]) + min([min([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]))/2

        r_z = (max([max([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.z for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dz'] * (border_type == 'open') / 2
        r_y = (max([max([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.y for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dy'] * (border_type == 'open') / 2
        r_x = (max([max([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]) - min([min([Ring.x for Ring in Rings[orientation]]) for orientation in orientations]))/2 + Params['Dx'] * (border_type == 'open') / 2

        for orientation in orientations:
            for Ring in Rings[orientation][:]:
                distance = (Ring.x-r_x0) ** 2/r_x**2 * (axes != 'x') + (Ring.y-r_y0) ** 2/r_y **2 * (axes != 'y') + (Ring.z-r_z0) ** 2/r_z ** 2 * (axes != 'z')
                if distance > 1.01:
                    if Fill:
                        Ring.R = np.inf
                        Ring.C = np.inf
                        Ring.L = 0
                    else:
                        Rings[orientation].remove(Ring)

    return Rings



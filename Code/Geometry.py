from Ring_Class import Ring

def Rectangle_packing(Params):
    if len(Params['Orientations']) == 1:
        shift_x = Params['shift_x']
        shift_y = Params['shift_y']
        shift_z = Params['shift_z']
    else:
        shift_x = 0
        shift_y = 0
        shift_z = 0 
    
    N = Params['N']

    
    rings = []

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                Shift_x = shift_x * (k * (orientation == 'z') + j * (orientation == 'y'))
                Shift_y = shift_y * (k * (orientation == 'z') + i * (orientation == 'x'))
                Shift_z = shift_z * (j * (orientation == 'y') + i * (orientation == 'x'))
                dx = (orientation != 'x') * delta_x/2
                dy = (orientation != 'y') * delta_y/2
                dz = (orientation != 'z') * delta_z/2
                rings.append(
                    Ring(
                        # Prevent rings from getting out of the borders
                        (i * delta_x + Shift_x) % ((nx) * delta_x) + dx,
                        (j * delta_y + Shift_y) % ((ny) * delta_y) + dy,
                        (k * delta_z + Shift_z) % ((nz) * delta_z) + dz,
                        orientation,
                        r,
                        w)
                )
    return np.array(rings, dtype=Ring)

def Hexagonal_packing(N, r, delta_x, delta_y, delta_z, orientation, shift_x = 0, shift_y = 0):
    nz, ny, nx = N['z'], N['y'], N['x']    
    rings = []

    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                rings.append(
                    Ring(
                        # Prevent rin gs from getting out of the borders
                        (i * delta_x + shift_x * k + delta_x/2 * j) % ((nx) * delta_x),
                        (sqrt(3)/2*j * delta_x + shift_y * k) % ((ny) * delta_y),
                        k * delta_z,
                        orientation,
                        r,
                        w)
                )
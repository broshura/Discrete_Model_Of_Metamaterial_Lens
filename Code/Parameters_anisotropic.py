# Anisotropic cube parameters

from math import pi, log, sqrt
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings
# are normalized to cell sizes thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system

def Rectangle_packing(N, r, delta_x, delta_y, delta_z, w, orientation = "z", shift_x = 0, shift_y = 0, shift_z = 0):
    nz, ny, nx = N['z'], N['y'], N['x']
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

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum
L = 7.14*10**-9                        # Self-inductance
C = 800 * 10 ** -12                    # Capacitance
R = 0.0178                             # Resistance
omega_0 = 1/np.sqrt(L * C)             # Self-frequence                

N = {'x':15, 'y':12, 'z':34}           # Number of cell on each dimension
Delta_z = 1.5 * 10 ** -3               # Distance between layers
Delta_y = 3 * Delta_z                  # y-distance
Delta_x = 3 * Delta_z                  # Length of cell
n_0 = 1 / (Delta_x ** 2 * Delta_z)     # Volume concentration per one ring
delta_z = 1                            # Distance between layers in packing units
delta_y = delta_z * (Delta_y/Delta_z)  # y-distance in packing units
delta_x = delta_z * (Delta_x/Delta_z)  # Length of cell in packing units

shift_x = 0                            # Shifting of next layer along x axes
shift_y = 0                            # Shifting of next layer along y axes

Radius = Delta_x/3                     # Mean radius of rings
radius = delta_x/3                     # Mean radius of rings in packing units

W = 0.7*Radius/3                       # Width of strip
w = 0.7*radius/3                       # Width of strip in packing units

name = "Anisotropic-Nonshifted-cube"   # Name of Data file with this parameter set

Params = {'L': L,                            # Self-inductance
    'C': C,                                  # Capacitance
    'R': R,                                  # Resistance
    'W': 0.7*Radius/3,                       # Width of strip
    'Radius': Delta_x/3,                     # Mean radius of rings
    'Dz': Delta_z,                           # Distance between layers
    'Dy': 3 * Delta_z,                       # Length of cell
    'Dx': 3 * Delta_z,                       # y-distance
    'Rectangle_packing': Rectangle_packing,    

    'N': N,                                  # Number of cell on each dimension
    'dz': 1,                                 # Distance between layers in packing units
    'dy': delta_z * (Delta_y/Delta_z),       # Length of cell in packing units
    'dx': delta_z * (Delta_z/Delta_z),       # y-distance in packing units
    'shift_x': 0,                            # Shifting of next layer along x axes
    'shift_y': 0,                            # Shifting of next layer along y axes
    'radius': delta_x/3,                     # Mean radius of rings in packing units
    'w': 0.7*radius/3,                       # Width of strip in packing units
    'name': "Anisotropic-Nonshifted-cube",    # Name of Data file with this parameter set
    'Orientations': ('z'),
    'Self-frequence': omega_0                
}
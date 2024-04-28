# Anisotropic cube parameters

from numpy import pi, log, sqrt
import numpy as np

# Parameters for system used in modeling

L = 7.14*10**-9                        # Self-inductance
C = np.inf                             # Capacitance
R = 0.0178                             # Resistance
omega_0 = 1/np.sqrt(L * C)             # Self-frequence                

N = {}
N['z'] = {'nx':15, 'ny':12, 'nz':34}   # Number of cell on each dimension
Dz = 1e-3               # Distance between layers
Dy = 4.5e-3                  # y-distance
Dx = 4.5e-3                  # Length of cell
n_0 = 1 / (Dx ** 2 * Dz)     # Volume concentration per one ring

shift_x = 0                            # Shifting of next layer along x axes
shift_y = 0                            # Shifting of next layer along y axes
shift_z = 0                            # Shifting of next layer along z axes

Radius = 1.5e-3                         # Mean radius of rings

W = 0.5e-3                              # Width of strip


name = "Anisotropic-Nonshifted-cube"   # Name of Data file with this parameter set

Params = {'L': L,                            # Self-inductance
    'C': C,                                  # Capacitance
    'R': R,                                  # Resistance
    'W': 0.7*Radius/3,                       # Width of strip
    'Radius': Dx/3,                     # Mean radius of rings
    'Dz': Dz,                           # Distance between layers
    'Dy': 3 * Dz,                       # Length of cell
    'Dx': 3 * Dz,                       # y-distance  
    'N': N,                                  # Number of cell on each dimension
    'shift_x': 0,                            # Shifting of next layer along x axes
    'shift_y': 0,                            # Shifting of next layer along y axes
    'shift_z': 0,                            # Shifting of next layer along z axes
    'Orientations': ('z'),                   # List of orientations for system   
    'Self-frequence': omega_0,               # Self-frequency
}
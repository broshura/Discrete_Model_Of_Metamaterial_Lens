# Anisotropic cube parameters

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

Radius = 1.5e-3                        # Mean radius of rings

W = 0.5e-3                             # Width of strip
d = 18e-6                              # Thickness of strip
Sigma = -0.06                          # Lattice constant 

name = "Anisotropic-Nonshifted-cube"   # Name of Data file with this parameter set
H_0z = 1.0                            # External magnetic field
mu_0 = 4 * np.pi * 10 ** -7            # Permeability of free space
Orientations = 'z'
Omega = np.linspace(0.9 * omega_0, 1.1 * omega_0, 1000) # Frequency range
Params = {
    'L': L,                     # Self-inductance
    'C': C,                     # Capacitance
    'R': R,                     # Resistance
    'W': W,                     # Width of strip
    'Radius': Radius,           # Mean radius of rings
    'Dz': Dz,                   # Length of cell
    'Dy': Dz,
    'Dx': Dz,
    'shift_x': 0,               # Shifting of next layer along x axes
    'shift_y': 0,               # Shifting of next layer along y axes
    'shift_z': 0,               # Shifting of next layer along z axes
    'Orientations': Orientations,
    'Sigma': Sigma,
    'H_0z': H_0z,
    'Omega': [Omega[0], Omega[-1], len(Omega)],
    'mu_0': mu_0,
    'omega_0': omega_0,
    'Solver_type': "Fast",
    'Packing': 'Rectangle',
    'P_0z': np.pi * Radius ** 2 /H_0z/Dz/Dy/Dx,
    'N' : {
        'z':{'nz': 11, 'ny': 10, 'nx': 10},
        'y':{'nz': 10, 'ny': 11, 'nx': 10},
        'x':{'nz': 10, 'ny': 10, 'nx': 11}
    },
    'shape': '10x10x10',
    'Numbers': {
        'z': 1100,
        'y': 1100,
        'x': 1100
    },
    'Threads': 1,
    'Solver_name': 'lgmres'
}

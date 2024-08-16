# This file contains all parameters for modeling and geometry of rings

import numpy as np

# Parameters for system used in modeling

L = 13.459 * 10 ** -9                   # Self-inductance
C = 470 * 10 ** -10                     # Capacitance
R = 0.002                               # Resistance
omega_0 = 1 / np.sqrt(L * C)            # Resonance frequency
Omega =  np.linspace(omega_0*0.9, omega_0*1.1, 1000) # Frequency range
H_0z = 1                                # Amplitude of magnetic field
mu_0 = 4 * np.pi * 10 ** -7                # Magnetic constant

Dz = 15 * 10 ** -3                      # Length of cell
Dy = Dz
Dx = Dz

Radius = 4.935 * 10 ** -3               # Mean radius of rings
W = 0.7 * 0.15 * Dz * 0                 # Width of strip

Sigma = -0.06                           # Lattice constant

Orientations = 'zyx'                    # Orientations of rings

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


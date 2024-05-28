# This file contains all parameters for modeling and geometry of rings
import os
import sys

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Ring_Class import Ring
from numpy import pi
import numpy as np


# Parameters for system used in modeling

L = 13.459 * 10 ** -9                   # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = float(62.46 * 1e6 * 2*np.pi)  # Resounance frequency for effective magnetic permeability

Dz = 15 * 10 ** -3                      # Length of cell
Dy = Dz
Dx = Dz

Radius = 4.935 * 10 ** -3               # Mean radius of rings
W = 0.7 * 0.15 * Dz * 0                 # Width of strip

Sigma = -0.06                           # Lattice constant

Orientations = 'zyx'                    # Orientations of rings
# Adding responding ring to identify resonance frequency

name = "MRI-self"

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
    'Self-frequence': omega_0,  # Resonance frequency for effective magnetic permeability
    'Sigma': Sigma,
}


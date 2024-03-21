# This file contains all parameters for modeling and geometry of rings

from numpy import pi
import numpy as np


# Parameters for system used in modeling

L = 13.46 * 10 ** -9                    # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 1/np.sqrt(L * C)              # Self-frequence

N = {}
N["z"] = {'nz':10, 'ny': 9, 'nx': 9}            # Number of cells on each row for z-oriented rings
N["y"] = {'nz':9, 'ny': 10, 'nx': 9}            # Number of cells on each row for y-oriented rings
N["x"] = {'nz':9, 'ny': 9, 'nx': 10}            # Number of cells on each row for x-oriented rings

Dz = 15 * 10 ** -3                      # Length of cell
Dy = Dz
Dx = Dz
Radius = 4.935 * 10 ** -3               # Mean radius of rings
W = 0.7 * 0.15 * Dz                     # Width of strip

R_coil = 1.5                           # Resistance of responding ring
L_coil = 1.8*10**-7                    # Self-inductance of responding  ring
Radius_coil = 3 * 2.54 * 10 ** -2      # Radius of responding ring

Orientations = ('x', 'y', 'z')
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
    'N': N,                     # List of numbers in each dimension for each orientation
    'shift_x': 0,               # Shifting of next layer along x axes
    'shift_y': 0,               # Shifting of next layer along y axes
    'shift_z': 0,               # Shifting of next layer along z axes
    'Orientations': Orientations,
    'Self-frequence': omega_0  # Frequency of resonance in free space
    # Adding responding ring to identify resonance frequency
}


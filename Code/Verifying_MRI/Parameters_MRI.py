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

N = {}
N["z"] = {'nz':17, 'ny': 2, 'nx': 18}            # Number of cells on each row for z-oriented rings
N["y"] = {'nz':18, 'ny': 3, 'nx': 18}            # Number of cells on each row for y-oriented rings
N["x"] = {'nz':18, 'ny': 2, 'nx': 17}            # Number of cells on each row for x-oriented rings

# N["z"] = {'nz':1, 'ny': 2, 'nx': 2}            # Number of cells on each row for z-oriented rings
# N["y"] = {'nz':2, 'ny': 3, 'nx': 2}            # Number of cells on each row for y-oriented rings
# N["x"] = {'nz':2, 'ny': 2, 'nx': 1}            # Number of cells on each row for x-oriented rings


Dz = 15 * 10 ** -3                      # Length of cell
Dy = Dz
Dx = Dz
Radius = 4.935 * 10 ** -3               # Mean radius of rings
W = 0.7 * 0.15 * Dz * 0                 # Width of strip
Sigma = -0.06                           # Lattice constant
R_coil = 1.5                           # Resistance of responding ring
L_coil = 1.8*10**-7                    # Self-inductance of responding  ring
C_coil = 1e10#np.inf                        # Capacitance of responding ring
Radius_coil = 3 * 2.54 * 10 ** -2/2    # Radius of responding ring

R0 = {                                 # Initial position of the first ring for each orientation
    'z': {'nz': Dz, 'ny': Dy/2, 'nx': Dx/2},
    'y': {'nz': Dz/2, 'ny': 0, 'nx': Dx/2},
    'x': {'nz': Dz/2, 'ny': Dy/2, 'nx': Dx}
}
Responded_x = Dx * N['y']['nx']/2
Responded_y = -Dy
Responded_z = Dz * N['y']['nz']/2


Responded_ring = Ring(Responded_x, Responded_y, Responded_z, 'y', Radius_coil, 0, L_coil, C_coil, R_coil)
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
    'Self-frequence': omega_0,  # Frequency of resonance in free space
    'Number': np.sum([N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx'] for pos in N])  + 1, # Number of rings in the system
    'Responded_x' : Responded_x,
    'Responded_y' : Responded_y,
    'Responded_z' : Responded_z,
    'Responded_pos': 'y',
    'Radius_coil': Radius_coil,
    'L_coil': L_coil,
    'C_coil': C_coil,
    'R_coil': R_coil,
    'W_coil': 0,
    'Sigma': Sigma,
}


# This file contains all parameters for modeling and geometry of rings

from math import pi
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings are normalized to cell size thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system for MRI lenz - parallelogram without borders at x and z surfaces.
def Rectangle_packing(n, responding_ring, r, r_coil, w):
    nx_x, ny_x, nz_x = n["xx"], n['yx'], n['zx']
    nx_y, ny_y, nz_y = n["xy"], n['yy'], n['zy']
    nx_z, ny_z, nz_z = n["xz"], n['yz'], n['zz']
    
    rings = []

    # Putting X-oriented rings
    for k in range(nz_x):
        for j in range(ny_x):
            for i in range(nx_x):
                rings.append(Ring(i * 2, j * 2 + 1, k * 2 + 1, "x", r, w))
    # Putting Y-oriented rings
    for k in range(nz_y):
        for j in range(ny_y):
            for i in range(nx_y):
                rings.append(Ring(i * 2 + 1, j * 2, k * 2 + 1, "y", r, w))
    # Putting Z-oriented rings
    for k in range(nz_z):
        for j in range(ny_z):
            for i in range(nx_z):
                rings.append(Ring(i * 2 + 1, j * 2 + 1, k * 2, "z", r, w))

    if responding_ring:
        rings.append(Ring(0, -2*2, 0, "y", r_coil, 0))
    
    return np.array(rings, dtype=Ring)

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7                # Permeability of vacuum

L = 13.46 * 10 ** -9                    # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 63.28 * 10 ** 6               # Frequency of resonance in free space

N = {}

N["xx"], N["yx"], N["zx"] = 10, 9, 9            # Number of cells on each row for x-oriented rings
N["xy"], N["yy"], N["zy"] = 9, 10, 10           # Number of cells on each row for y-oriented rings
N["xz"], N["yz"], N["zz"] = 9, 9, 10            # Number of cells on each row for z-oriented rings

start, end = {}, {}
start["x"], end["x"] = 0, N['xx'] * N['yx'] * N['zx']
start["y"], end["y"] = end["x"], end['x'] + N['xy'] * N['yy'] * N['zy']
start["z"], end["z"] = end["y"], end['y'] + N['xz'] * N['yz'] * N['zz']

Delta_z = 15 * 10 ** -3                      # Length of cell
delta_z = 2                                  # Length of cell in packing units
Radius = 4.935 * 10 ** -3                    # Mean radius of rings
radius = Radius * delta_z/Delta_z            # Mean radius of rings in packing units
W = 0.7 * 0.15 * Delta_z                     # Width of strip
w = 0.7 * 0.15 * delta_z                     # Width of strip in packing units


R_coil = 1.5                           # Resistance of responding ring
L_coil = 1.8*10**-7                    # Self-inductance of responding  ring
Radius_coil = 3 * 2.54 * 10 ** -2      # Radius of responding ring
radius_coil = Radius_coil * delta_z/Delta_z      # Radius of responding ring in packing units

Orientations = ('x', 'y', 'z')
# Adding responding ring to identify resonance frequency

Responding_Ring = False

# Name of Data file for this parameters
if Responding_Ring:
    name = "MRI"                            
else:
    name = "MRI-self"

Params = {
    'L': L,                     # Self-inductance
    'C': C,                     # Capacitance
    'R': R,                     # Resistance
    'Self-frequence': omega_0,  # Frequency of resonance in free space
    'N': N,                     # List of numbers in each dimension for each orientation
    'start': start,             # List of start and end numbers of different orientations
    'end': end,
    'Dz': Delta_z,                   # Length of cell
    'Dy': Delta_z,
    'Dx': Delta_z,
    'dz': delta_z,                   # Length of cell in packing units
    'dy': delta_z,
    'dx': delta_z,
    'Radius': Radius,           # Mean radius of rings
    'radius': radius,           # Mean radius of rings in packing units
    'W': W,                     # Width of strip
    'w': w,                     # Width of strip in packing units

    'R_coil': R_coil,           # Resistance of responding ring
    'L_coil': L_coil,           # Self-inductance of responding  ring
    'Radius_coil': Radius_coil, # Radius of responding ring
    'radius_coil': radius_coil, # Radius of responding ring in packing units

    'Orientations': Orientations,
    # Adding responding ring to identify resonance frequency

    'Responding_Ring': Responding_Ring,
    'name': name
}


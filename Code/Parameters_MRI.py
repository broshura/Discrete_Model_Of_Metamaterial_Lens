# This file contains all parameters for modeling and geometry of rings

from math import pi
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings are normalized to cell size thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system for MRI lenz - parallelogram without borders at x and z surfaces.
def Rectangle_packing(nx_x, ny_x, nz_x, nx_y, ny_y, nz_y, nx_z, ny_z, nz_z, r, w):
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
    return np.array(rings, dtype=Ring)

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7                # Permeability of vacuum

L = 13.46 * 10 ** -9                    # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 63.28 * 10 ** 6               # Frequency of resonance in free space

N = {}

N["xx"], N["yx"], N["zx"] = 10, 9, 9            # Number of cells on each row for x-oriented rings
N["xy"], N["yy"], N["zy"] = 9, 10, 10            # Number of cells on each row for y-oriented rings
N["xz"], N["yz"], N["zz"] = 9, 9, 10            # Number of cells on each row for y-oriented rings

start, end = {}, {}
start["x"], end["x"] = 0, N['xx'] * N['yx'] * N['zx']
start["y"], end["y"] = end["x"], end['x'] + N['xy'] * N['yy'] * N['zy']
start["z"], end["z"] = end["y"], end['y'] + N['xz'] * N['yz'] * N['zz']

a1 = 15 * 10 ** -3                     # Length of cell
a = 2                                  # Length of cell in packing units
Radius1 = 4.935 * 10 ** -3             # Mean radius of rings
Radius = 0.329 * a                     # Mean radius of rings in packing units
w1 = 0.7 * 0.15 * a1                   # Width of strip
w = 0.7 * 0.15 * a                     # Width of strip in packing units


R_coil = 1.5                           # Resistance of responding ring
L_coil = 1.8*10**-7                    # Self-inductance of responding  ring
Radius_coil1 = 3 * 2.54 * 10 ** -2     # Radius of responding ring
Radius_coil = 2.54 * a                 # Radius of responding ring in packing units

# List of Rings with their coordinates
Rings = Rectangle_packing(N["xx"], N["yx"], N["zx"],
                          N["xy"], N["yy"], N["zy"],
                          N["xz"], N["yz"], N["zz"],
                          Radius, w)
Orientations = ('x', 'y', 'z')
# Adding responding ring to identify resonance frequency

Responding_Ring = False

if Responding_Ring:
    Rings = np.append(Rings, Ring(0, -2*a, 0, "y", Radius_coil, 0))
    name = "MRI"                            # Name of Data file for this parameters
else:
    name = "MRI-self"

Number = len(Rings)                         # Number of Rings

V = np.array([0 for x in range(Number-1)] + [1])   # Voltage on each ring (only at responding)

#print(Rings)



# This file contains all parameters for modeling and geometry of rings

from math import pi
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings are normalized to cell size thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system for MRI lenz - parallelogram without borders at x and z surfaces.
def Rectangle_packing(nx, ny, nz, r):
    rings = []

    # Putting X-oriented rings
    for i in range(1, nx):
        for j in range(ny):
            for k in range(nz):
                rings.append(Ring(i * 2 - nx, j * 2 + 1 - ny, k * 2 + 1 - nz, "x", r))
    # Putting Y-oriented rings
    for i in range(nx):
        for j in range(ny + 1):
            for k in range(nz):
                rings.append(Ring(i * 2 + 1 - nx, j * 2 - ny, k * 2 + 1 - nz, "y", r))
    # Purring Z-oriented rings
    for k in range(1, nz):
        for j in range(ny):
            for i in range(nx):
                rings.append(Ring(i * 2 + 1 - nx, j * 2 + 1 - ny, k * 2 - nz, "z", r))
    return np.array(rings, dtype=Ring)

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum

L = 13.46 * 10 ** -9                     # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 63.28 * 10 ** 6               # Frequency of resonance in free space

Nx, Ny, Nz = 18, 2, 18                 # Number of cell on each row
a1 = 15 * 10 ** -3                     # Length of cell
a = 2                                  # Length of cell in packing units
Radius1 = 4.935 * 10 ** -3             # Mean radius of rings
Radius = 0.329 * a                     # Mean radius of rings in packing units
w1 = 0.7 * 0.15 * a1                   # Width of strip
w = 0.7 * 0.15 * a                     # Width of strip in packing units


R_coil = 1.5                           # Resistance of responding ring
L_coil = 1.8*10**-7                    # Self-inductance of responding  ring
#L_coil = 2.6 * 10 ** -7                # Another parameter for responding ring
Radius_coil1 = 3 * 2.54 * 10 ** -2      # Radius of responding ring
Radius_coil = 2.54 * a                  # Radius of responding ring in packing units

Rings = Rectangle_packing(Nx, Ny, Nz, Radius)   # List of Rings with their coordinates

# Adding responding ring to identify resonance frequency

Rings = np.append(Rings, Ring(0, -4, 0, "y", Radius_coil))

Number = len(Rings)                         # Number of Rings

name = "MRI"                                # Name of Data file for this parameters

V = [0 for x in range(Number-1)] + [1]    # Voltage on each ring (only at responding)



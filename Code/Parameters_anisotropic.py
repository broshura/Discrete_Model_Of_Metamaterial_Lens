# This file contains all parameters for modeling and geometry of rings
# Special for anisotropic structures and copper rings

from math import pi, log, sqrt
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings
# are normalized to cell sizes thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system for MRI lenz

def Rectangle_packing(nx, ny, nz, r, a, b, w,orientation = "z", shift_x = 0, shift_y = 0):
    rings = []

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                rings.append(
                    Ring(
                        # Prevent rings from getting out of the borders
                        (i*a + shift_x * k) % ((nx)*a),
                        (j*a + shift_y * k) % ((ny)*a),
                        k*b,
                        orientation,
                        r,
                        w)
                )
    return np.array(rings, dtype=Ring)

def Hexagonal_packing(nx, ny, nz, r, a, b, orientation, shift_x = 0, shift_y = 0):
    rings = []

    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                rings.append(
                    Ring(
                        # Prevent rin gs from getting out of the borders
                        (i * a + shift_x * k + a/2 * j) % ((nx) * a),
                        (sqrt(3)/2*j * a + shift_y * k) % ((ny) * a),
                        k * b,
                        orientation,
                        r,
                        w)
                )

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum

L = 7.14*10**-9                           # Self-inductance
C = "Infinity"                         # Capacitance
R = 0.0178                         # Resistance

Nx, Ny, Nz = 15, 12, 34                # Number of cell on each dimension
b1 = 1 * 10 ** -3                      # Distance between layers
a1 = 4.5 * b1                         # Length of cell
n_0 = 1 / (a1 ** 2 * b1)               # Volume concentration per one ring
b = 1                                  # Distance between layers in packing units
a = b * (a1/b1)                      # Length of cell in packing units
shift_x = 0                            # Shifting of next layer along x axes
shift_y = 0                            # Shifting of next layer along y axes

Radius1 = a1/3                         # Mean radius of rings
Radius = a/3                           # Mean radius of rings in packing units
w1 = 0.7*Radius/3                          # Width of strip
w = 0.7*Radius/3                       # Width of strip in packing units

# List of Rings with their coordinates and its number

Rings = Rectangle_packing(Nx, Ny, Nz, Radius, a, b, w, "z", shift_x, shift_y)
Number = len(Rings)

name = "Anisotropic-Nonshifted-cube"   # Name of Data file with this parameter set

H_0z = 10
E_w = np.ones(Number) * pi * Radius1 ** 2 * H_0z * mu_0

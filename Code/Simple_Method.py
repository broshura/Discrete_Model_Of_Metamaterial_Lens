import numpy as np
import scipy
from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import linalg
from scipy import special

from Ring_Class import Ring
from Parameters import Number, Thickness, Radius, mu_0, Z_0, j, omega

K = special.ellipkinc       #  Incomplete elliptic integral of the first kind
E = special.ellipeinc       #  Incomplete elliptic integral of the second kind


def L_parallel(dx, dy, dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(pi/2, kappa) - 2 * E(pi/2, kappa))/kappa
        return A * r_2 * (r_2 + db * cos(alpha))/dp
    L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r_1, r_2))
    return L

def L_orthogonal(dx, dy ,dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha))))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(pi/2, kappa) - 2 * E(pi/2, kappa))/kappa
        return A * r_2 * dy * cos(alpha) / dp
    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r_1, r_2))
    return L

def Mnm(First_ring, Second_ring):
    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    if Thickness == 0:
        r = Radius
        if First_ring.pos == Second_ring.pos:
            return mu_0/ (4 * pi) * L_parallel(dx, dy, dz, r, r)
        else:
            return mu_0/ (4 * pi) * L_orthogonal(dx, dy, dz, r, r)
    else:
        R = Radius + Thickness
        r = Radius - Thickness
        if First_ring.pos == Second_ring.pos:
            return mu_0 / (4 * pi) * (L_parallel(dx, dy, dz, r, r) + L_parallel(dx, dy, dz, R, r) + L_parallel(dx, dy, dz, r, R) + L_parallel(dx, dy, dz, R, R))/4
        else:
            return mu_0 / (4 * pi) * (L_orthogonal(dx, dy, dz, r, r) + L_orthogonal(dx, dy, dz, R, r) + L_orthogonal(dx, dy, dz, r, R) + L_orthogonal(dx, dy, dz, R, R))/4

Z = np.eye(Number) * Z_0
Rings = []
for x in range(-1, 2, 2):
    for y in range(-1, 2, 2):
        for z in range(-1, 2, 2):
            Rings.append(Ring(2*x, 2*y, 2*z, "xy"))
Rings = np.array(Rings, dtype=Ring)


for n in range(Number):
    for m in range(Number):
        if n > m:
            R1 = Rings[n]
            R2 = Rings[m]
            Z[n][m] = j * omega * Mnm(R1, R2)
            Z[m][n] = Z[n][m]
V = np.eye(Number)
I = linalg.inv(Z).dot(V)


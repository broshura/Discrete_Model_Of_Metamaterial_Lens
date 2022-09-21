import numpy as np
import scipy
from numpy import sqrt, cos, sin, pi, isnan
from scipy import integrate
from scipy import linalg
from scipy import special

from Ring_Class import Ring
from Parameters import *

K = special.ellipk       #  Сomplete elliptic integral of the first kind
E = special.ellipe       #  Сomplete elliptic integral of the second kind


def L_parallel(dx, dy, dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa) - 2 * E(kappa))/kappa
        return A * r_2 * (r_2 + db * cos(alpha))/dp
    L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r_1, r_2))
    return mu_0 * L / (4 * pi)

def L_orthogonal(dx, dy ,dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa) - 2 * E(kappa))/kappa
        return A * r_2 * dy * cos(alpha) / dp
    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r_1, r_2))
    return mu_0 * L / (4 * pi)

def Mnm(First_ring, Second_ring, Data):
    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    id = (dx, dy, dz, First_ring.pos == Second_ring.pos)
    if id in Data:
        pass
    if w == 0:
        r = Radius
        if First_ring.pos == Second_ring.pos:
            return L_parallel(dx, dy, dz, r, r)
        else:
            return L_orthogonal(dx, dy, dz, r, r)
    else:
        R = Radius + w/2
        r = Radius - w/2
        if First_ring.pos == Second_ring.pos:
            return (L_parallel(dx, dy, dz, r, r) + L_parallel(dx, dy, dz, R, r) + L_parallel(dx, dy, dz, r, R) + L_parallel(dx, dy, dz, R, R))/4
        else:
            return (L_orthogonal(dx, dy, dz, r, r) + L_orthogonal(dx, dy, dz, R, r) + L_orthogonal(dx, dy, dz, r, R) + L_orthogonal(dx, dy, dz, R, R))/4

Z = np.eye(Number) * Z_0
Data = []

for n in range(Number):
    for m in range(Number):
        if n > m:
            R1 = Rings[n]
            R2 = Rings[m]
            Z[n][m] = 1j * omega * Mnm(R1, R2, Data)
            Z[n][m] = round(Z[n][m].real) + 1j * round(Z[n][m].imag)
            Z[m][n] = Z[n][m]




with open("DATA/Result.txt", "w") as res:
    res.write("Calculated eletricity in each ring \n")
    res.write(f"Number of ring: {Number}")
    res.write(f"Parameters of ring:\nResistance: {R} $\Omega$ \nSelf-capacitance: {C} mF \nSelf-Inductance: {L} Hn \n")
    res.write(f"Impedance matrix\n")
    for i in range(Number):
        res.write(" ".join(map(str, Z[i])) + "\n")
    res.write(f"Voltage and electricity on each ring:\n \n")
    I = 1j*linalg.inv(Z).dot(V)
    for i in range(Number):
        ring = Rings[i]
        res.write(f"№ {i + 1} x = {ring.x} y = {ring.y} z = {ring.z} orientation: {ring.pos} U = {V[i]} I = {I[i]} \n")


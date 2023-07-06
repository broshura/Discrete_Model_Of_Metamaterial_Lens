# Calculating geometry matrix M and saving as large data file for optimize next computation
import json

import numpy as np
from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import special
from Parameters_anisotropic import *

K = special.ellipk  # Сomplete elliptic integral of the first kind
E = special.ellipe  # Сomplete elliptic integral of the second kind

# Choose parameter set for computing

Number = len(Rings)


# Computing for parallel-oriented rings
def L_parallel(dx, dy, dz, r1, r2, width = 0):

    # Define function to integrate over first defined parameter

    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa_sq = 4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2)
        kappa = sqrt(kappa_sq)
        A = 1/(2*pi)*sqrt(r_1/dp) * ((2/kappa - kappa) * K(kappa_sq) - 2 * E(kappa_sq)/kappa)
        return A * r_2 * (r_2 + db * cos(alpha))/dp

    #Considering stripe width

    if r1 == r2 and width:
        R = r1 + w / 2
        r = r1 - w / 2
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, r))
        L_2, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
        L_3, err_3 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, R, R))
        L = (L_1 + 2*L_2 + L_3)/4
    elif width:
        id_r1 = r1 == min(r1, r2)
        id_r2 = r2 == min(r1, r2)
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r1 + id_r1 * width/2, r2 + id_r2*width/2))
        L_2, err_2 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r1 - id_r1 * width/2, r2 - id_r2*width/2))
        return (L_1 + L_2)/2
    else:
        L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r1, r2))
    return L

# Computing for orthogonal-oriented rings

def L_orthogonal(dx, dy ,dz, r1, r2, width):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa_sq = 4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2)
        kappa = sqrt(kappa_sq) + 10 ** -7
        A = 1 / (2 * pi) * sqrt(r_1 / (dp + 10 ** -7)) * ((2 / kappa - kappa) * K(kappa_sq) - 2 * E(kappa_sq) / kappa)
        return A * r_2 * dy * cos(alpha) / (dp + 10 ** -7)

    # Considering stripe width

    if r1 == r2 and width:
        R = r1 + w / 2
        r = r1 - w / 2
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, r))
        L_2, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
        L_3, err_3 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, R, R))
        L = (L_1 + 2 * L_2 + L_3) / 4
    elif width:
        id_r1 = r1 == min(r1, r2)
        id_r2 = r2 == min(r1, r2)
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r1 + id_r1 * width/2, r2 + id_r2*width/2))
        L_2, err_2 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r1 - id_r1 * width/2, r2 - id_r2*width/2))
        return (L_1 + L_2)/2
    else:
        L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r1, r2))
    return L

# Computing for any pair

def Mnm(First_ring, Second_ring, Data = {}):
    global cnt_ZY, cnt_YZ, cnt_ZX, cnt_XZ, cnt_XY, cnt_YX

    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    r1 = First_ring.r
    r2 = Second_ring.r

    # To avoid calculating integrals with same params each time
    # there is a dictionary with all parameters and values

    id_1 = f"{dx} {dy} {dz} {r1} {r2} {First_ring.pos}{Second_ring.pos}"
    id_2 = f"{-dx} {-dy} {-dz} {r2} {r1} {Second_ring.pos}{First_ring.pos}"

    if id_1 in Data:
        return Data[id_1]
    elif id_2 in Data:
        return Data[id_2]

    # Consider all types of parallel orientation and symmetry for x-z axes

    if First_ring.pos == Second_ring.pos:
        if First_ring.pos == "z":                           # Z-oriented rings
            l = L_parallel(dx, dy, dz, r1, r2, w)
        elif First_ring.pos == "y":                         # Y-oriented rings
            l = L_parallel(dx, -dz, dy, r1, r2, w)
        else:                                               # X-oriented rings
            l = L_parallel(-dz, dy, dx, r1, r2, w)

    # Consider all types of orthogonal orientation

    else:
        if First_ring.pos == "z":
            if Second_ring.pos == "y":                      # Z-Y oriented pair
                l = L_orthogonal(dx, dy, dz, r1, r2, w)
            else:                                           # Z-X oriented pair
                l = L_orthogonal((dy), (dx), dz, r1, r2, w)
        elif First_ring.pos == "y":
            if Second_ring.pos == "z":                      # Y-Z oriented pair
                l = L_orthogonal(dx, (dz), (dy), r1, r2, w)
            else:                                           # Y-X oriented pair
                l = L_orthogonal(-dz, (dx), (dy), r1, r2,  w)
        elif First_ring.pos == "x":
            if Second_ring.pos == "z":                      # X-Z oriented pair
                l = L_orthogonal(dy, (dz), (dx), r1, r2, w)
            else:                                           # X-Y oriented pair
                l = L_orthogonal((dz), dy, (dx), r1, r2, w)
    Data[id_1], Data[id_2] = l, l
    return l

# Calculating for each pair

M = np.eye(Number) * 0
Data = {}

print(f"Modeling impedance matrix for {name} set of parameters")
print(f"Number of rings: {Number}")


for n in range(Number):
    for m in range(Number):
        if n % (Number//10) == 0 and m % (Number) == 0:
            print(f"Done {n // (Number//10) * 10}%")
        if n != m:
            R1 = Rings[n]
            R2 = Rings[m]
            M[n][m] = Mnm(R1, R2, Data)


# Writing table in Data-file, divide string by \n and elements by " "



print("Saving data...")
with open(f"DATA/Data-{name}.txt", "w") as res:
    for i in range(Number):
        res.write(" ".join(map(str, M[i])) + "\n")

with open(f"DATA/Data-{name}.json", "w") as res:
    res.write(json.dumps(Data))


# Sum of matrix M elements to calculate permeability
with open(f"DATA/SumM-{name}.txt", "w") as res:
    if name.split("-")[0] == "Anisotropic":
        SumM = 0
        rings = Rings.reshape(Nz, Ny, Nx)
        Z, Y, X = Nz//2, Ny // 2, Nx // 2
        print(Z, Y, X)
        for z in [k for k in range(len(rings))]:
            break
            for y in [j for j in range(len(rings[z]))]:
                for x in [X + i for i in range(-len_x, len_x + 1)]:
                    if (Z, Y, X) != (z, y, x):
                        SumM += Mnm(rings[Z][Y][X], rings[z][y][x], Data)

        print(f"Sum of matrix elements:{SumM}")
        res.write(str(SumM))


print("Ended")

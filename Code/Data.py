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
        #print(f"db = {db}")
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        #print(f"dp = {dp}, r_1 = {r_1}, r_2 = {r_2}")
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2))
        #print(f"kappa = {kappa}, K(kappa) = {K(kappa)}, E(kappa) = {E(kappa)}")
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2))/kappa
        #print(f"A = {A}")
        #print("")
        return A * r_2 * (r_2 + db * cos(alpha))/dp
    L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r_1, r_2))
    return L

def L_orthogonal(dx, dy ,dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2))/kappa
        return A * r_2 * dy * cos(alpha) / dp
    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r_1, r_2))
    return L

def Mnm(First_ring, Second_ring, Data = {}):
    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    id = (dx, dy, dz, First_ring.pos, Second_ring.pos)
    if id in Data:
        return Data[id]
    R1 = First_ring.r + w/2
    r1 = First_ring.r - w/2
    R2 = Second_ring.r + w/2
    r2 = Second_ring.r - w/2
    if First_ring.pos == Second_ring.pos:
        if First_ring.pos == "z":
            l = (L_parallel(dx, dy, dz, r1, r2) + L_parallel(dx, dy, dz, R1, r2) + L_parallel(dx, dy, dz, r1, R2) + L_parallel(dx, dy, dz, R1, R2))/4
            Data[id] = l
        elif First_ring.pos == "y":
            l = (L_parallel(dx, -dz, dy, r1, r2) + L_parallel(dx, -dz, dy, R1, r2) + L_parallel(dx, -dz, dy, r1, R2) + L_parallel(dx, -dz, dy, R1, R2)) / 4
            Data[id] = l
        else:
            l = (L_parallel(-dz, dy, dx, r1, r2) + L_parallel(-dz, dy, dx, R1, r2) + L_parallel(-dz, dy, dx, r1, R2) + L_parallel(-dz, dy, dx, R1, R2)) / 4
            Data[id] = l
    else:
        if First_ring.pos == "z":
            if Second_ring.pos == "y":
                l = (L_orthogonal(dx, dy, dz, r1, r2) + L_orthogonal(dx, dy, dz, R1, r2) + L_orthogonal(dx, dy, dz, r1, R2) + L_orthogonal(dx, dy, dz, R1, R2))/4
                Data[id] = l
            else:
                l = (L_orthogonal(-dy, dx, dz, r1, r2) + L_orthogonal(-dy, dx, dz, R1, r2) + L_orthogonal(-dy, dx, dz, r1,R2) + L_orthogonal(-dy, dx, dz, R1, R2)) / 4
                Data[id] = l
        elif First_ring.pos == "y":
            if Second_ring.pos == "z":
                l = (L_orthogonal(-dx, -dy, -dz, r1, r2) + L_orthogonal(-dx, -dy, -dz, R1, r2) + L_orthogonal(-dx, -dy, -dz, r1, R2) + L_orthogonal(-dx, -dy, -dz, R1, R2)) / 4
                Data[id] = l
            else:
                l = (L_orthogonal(dz, dy, dx, r1, r2) + L_orthogonal(dz, dy, dx, R1, r2) + L_orthogonal(dz, dy, dx, r1, R2) + L_orthogonal(dz, dy, dx, R1, R2)) / 4
                #l = L_orthogonal(dz, dx, dy, r1, R2)
                #l = L_orthogonal(dz, dy, dx, r1, R2)
                Data[id] = l
        elif First_ring.pos == "x":
            if Second_ring.pos == "z":
                l = (L_orthogonal(dy, -dx, -dz, r1, r2) + L_orthogonal(dy, -dx, -dz, R1, r2) + L_orthogonal(dy, -dx, -dz, r1, R2) + L_orthogonal(dy, -dx, -dz, R1, R2)) / 4
                Data[id] = l
            else:
                l = (L_orthogonal(dz, -dy, dx, r1, R2) + L_orthogonal(dz, -dy, dx, R1, r2) + L_orthogonal(dz, -dy, dx, r1, R2) + L_orthogonal(dz, -dy, dx, R1, R2)) / 4
                #l = L_orthogonal(dz, -dx, dy, r1, R2)
                #l = L_orthogonal(dz, -dy, dx, r1, R2)
                Data[id] = l
    return l

M = np.eye(Number) * 0
Data = {}
for n in range(Number):
    for m in range(Number):
        if n > m:
            if n % 2000 == 1999 and m%2000 == 1999:
                print(" ")
            R1 = Rings[n]
            R2 = Rings[m]
            M[n][m] = Mnm(R1, R2, Data)
            M[m][n] = M[n][m]

with open("DATA/Data.txt", "w") as res:
    for i in range(Number):
        res.write(" ".join(map(str, M[i])) + "\n")
    pass
# with open("DATA/matrica.txt", "r") as res:
#     RES = res.read()
#     M = [[k for k in x.split(" ")] for x in RES.split("\n")]
# #print(M[1])
# M1 = [[float(M[i][k])*(Radius) ** 1 * 10 ** -7 for k in range(len(M[i]))] for i in range(len(M)-1)]
# #print(M1[1])
# with open("DATA/Data.txt", "r") as res:
#     RES = res.read()
#     M = [[k for k in x.split(" ")] for x in RES.split("\n")]
# #print(M[1])
# M2 = [[float(M[i][k])*(a1/2) ** 1 for k in range(len(M[i]))] for i in range(len(M)-1)]
# #print(M2[1])


# print(Rings[C])
#
# cnt_sign = 0
# cnt_abs = 0
#
# print(" Смотрим на точки где ошиблись в знаке ")
# for i in range(len(M[C])):
#     if M1[C][i] !=0 and M2[C][i]!= 0:
#         id = M2[C][i]/M1[C][i]
#         if id < 0:
#             cnt_sign+=1
#             print(i, id, Rings[C].pos, Rings[i].pos, f"dx = {Rings[i].x - Rings[C].x}, dy = {Rings[i].y - Rings[C].y}, dz = {Rings[i].z - Rings[C].z}")
# print(f": {cnt_sign} колец")
# print(f" Смотрим на точки где ошиблись кратно")
# print("")
# for i in range(len(M[1])):
#     if M1[C][i] !=0 and M2[C][i]!= 0:
#         id = M2[C][i]/M1[C][i]
#         if (id > 10 or id < 9 ) and id > 0:
#             cnt_abs+=1
#             print(i, id, Rings[C].pos, Rings[i].pos,
#                   f"dx = {Rings[i].x - Rings[C].x}, dy = {Rings[i].y - Rings[C].y}, dz = {Rings[i].z - Rings[C].z}")
# print(f"{cnt_abs} колец")
# with open("DATA/Result.txt", "w") as res:
#     res.write("Calculated eletricity in each ring \n")
#     res.write(f"Number of ring: {Number}")
#     res.write(f"Parameters of ring:\nResistance: {R} $\Omega$ \nSelf-capacitance: {C} mF \nSelf-Inductance: {L} Hn \n")
#     res.write(f"Impedance matrix\n")
#     for i in range(Number):
#         res.write(" ".join(map(str, Z[i])) + "\n")
#     res.write(f"Voltage and electricity on each ring:\n \n")
#     I = 1j*linalg.solve(Z, V)
#     for i in range(Number):
#         ring = Rings[i]
#         res.write(f"№ {i + 1} x = {ring.x} y = {ring.y} z = {ring.z} orientation: {ring.pos} U = {V[i]} I = {I[i]} \n")
#

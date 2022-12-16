# This file contains all parameters for modeling and rings

from math import pi
from Ring_Class import Ring
import numpy as np
import matplotlib.pyplot as plt

def Rectangle_packing(nx, ny, nz, r,  type = "closed"):
    rings = []
    if type == "closed":
        for i in range(1, nx):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i * 2 - nx, j * 2 + 1 - ny, k * 2 + 1 - nz, "x", r))
        for i in range(nx):
            for j in range(ny + 1):
                for k in range(nz):
                    rings.append(Ring(i * 2 + 1 - nx, j * 2 - ny, k * 2 + 1 - nz, "y", r))
        for k in range(1, nz):
            for j in range(ny):
                for i in range(nx):
                    rings.append(Ring(i * 2 + 1 - nx, j * 2 + 1 - ny, k * 2 - nz, "z", r))
    elif type == "open":
        for i in range(nx):
            for j in range(ny):
                for k in range(1, nz):
                    rings.append(Ring(i * 2 + 1, j * 2 + 1, k * 2, "xy"))
        for i in range(nx):
            for j in range(1, ny):
                for k in range(nz):
                    rings.append(Ring(i * 2 + 1, j * 2, k * 2 + 1, "xz"))
        for i in range(1,nx):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i * 2, j * 2 + 1, k * 2 + 1, "zy"))
        return np.array(rings, dtype=Ring)
    else:
        return "Unidentified type of packing"
    return np.array(rings, dtype=Ring)
def Sphere_Packing(Nr, type = "closed"):
    if type == "closed":
        pass
    else:
        pass




L = 13.5 * 10 ** -9                     # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 63.28 * 10 ** 6               # Frequency of resonance in free space


Nx, Ny, Nz = 18, 2, 18                    # Number of cell on each row
a1 = 15 * 10 ** -3                       # Length of cell
a = 2
Radius = 0.33*a                    # Mean radius of rings
w = 0.7*0.15 * a                       # Width of strip
mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum

R_coil = 1.5
L_coil = 1.8*10**-7
#L_coil = 2.6*10**-7

Rings = Rectangle_packing(Nx, Ny, Nz, Radius)   # List of Rings with their coordinates
#Rings = [Ring(0, 0, 0, "xy", Radius), Ring(0, 0, a, "xy", Radius)]#, Ring(0, 0, 2*a, "xy", Radius)]
Number = len(Rings)+1                     # Number of Rings
Rings = np.append(Rings, Ring(0, -4, 0, "y", 7.62/1.5/3*5))
V = [0 for x in range(Number-1)] + [1]      # Voltage on each ring


#
# cnt_sign = 0
# cnt_abs = 0

# with open("DATA/matrica.txt", "r") as res:
#     RES = res.read()
#     M = [[k for k in x.split(" ")] for x in RES.split("\n")]
# #print(M[1])
# M1 = [[float(M[i][k]) for k in range(len(M[i]))] for i in range(len(M)-1)]
#print(M1[1])
# with open("DATA/Data.txt", "r") as res:
#     RES = res.read()
#     M = [[k for k in x.split(" ")] for x in RES.split("\n")]
# #print(M[1])
# M2 = [[float(M[i][k])*(a1/2) ** 1/Radius for k in range(len(M[i]))] for i in range(len(M)-1)]
# M1 = M1 + np.eye(len(M1))
# M2 = M2 + np.eye(len(M2)) * 9.6
# print(len(M1), len(M2))
# print(M1[100][101])
# #print(M2[1])
# for i in range(len(M1)):
#     for j in range(len(M1[i])):
#         id = M2[i][j]/M1[i][j]
#         if not (7 < id < 12):
#             cnt_abs += 1
#             print(f"Индексы и расхождение :{i, j, id, cnt_abs}")
#             print(f"Значения которые расходятся: {M2[i][j], M1[i][j]}")
#             print(f'Параметры интегрирования: dx = {Rings[j].x-Rings[i].x}, dy = {Rings[j].y-Rings[i].y}, dz = {Rings[j].z-Rings[i].z}')
#             print(f"Первое кольцо {Rings[i]}")
#             print(f"Второе кольцео: {Rings[j]}")
#             print("")

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
# for i in range(len(M[C])):
#     if M1[C][i] !=0 and M2[C][i]!= 0:
#         id = M2[C][i]/M1[C][i]
#         if (id >10 or id < 9 ) and id > 0:
#             cnt_abs+=1
#             print(i, id, Rings[C].pos, Rings[i].pos,
#                   f"dx = {Rings[i].x - Rings[C].x}, dy = {Rings[i].y - Rings[C].y}, dz = {Rings[i].z - Rings[C].z}")
# print(f"{cnt_abs} колец")
#fig = plt.figure(figsize = (20, 20))

#ax = fig.add_subplot(1, 1, 1, projection = '3d')

#surf = ax.scatter(
#[ring.x for ring in Rings if ring.pos=="xy"], [ring.y for ring in Rings if ring.pos=="xy"], [ring.z for ring in Rings if ring.pos=="xy"],
#color = "b")
#plt.title("Centers of all xy oriented rings", size = 36)
#ax.set_ylim(-18, 18)
#ax.set_xlabel("x", size = 36)
#ax.set_ylabel("y", size = 36)
#ax.set_zlabel("z", size = 36)
#plt.show()
#print(len(Rings))

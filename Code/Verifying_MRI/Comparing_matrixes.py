# Comparing geometric matrix of metamaterial MRI lenz

from Parameters_MRI import delta_z, Dz, radius, Rings
import numpy as np
from numpy import sign

cnt_sign = 0        # Counter of rings with sign errors
cnt_abs = 0         # Counter of rings with large absolute error

# Matrix produced from verified matlab code
# To compare it is necessary to multiply at normalized length
#with open("DATA/matrica.txt", "r") as res:

# Matrix from fast method
with open("DATA/FData-MRI.txt", "r") as res:
    RES = res.read()
    Strings = RES.split("\n")[:len(RES.split("\n"))]
    Nx, Ny, Nz = map(int, Strings[0].split())
    M1 = {}
    for i in range(9):
        M1.update({Strings[2*i + 1]: np.array([float(s) for s in Strings[2*i + 2].split()])})
    for m in M1:
        M1[m].resize(Nx, Ny, Nz)

# Matrix from python - code

with open("DATA/Data-MRI.txt", "r") as res:
    RES = res.read()
    Strings = RES.split("\n")[:len(RES.split("\n")) - 1]
    M2 = [[float(k) for k in x.split(" ")] for x in Strings]

# Ratios between which elements are equal
Ratio_lowest = 0.98
Ratio_highest = 1.02

print("Comparing in absolute ratio \n")
for i in range(len(M2)-1):
    for j in range(len(M2[i])-1):
        # Avoid zero division on diagonal
        if i != j:
            dx, dy, delta_z = Rings[j].x - Rings[i].x, Rings[j].y - Rings[i].y, Rings[j].z - Rings[i].z
            m = Rings[i].pos + Rings[j].pos
            id = M2[i][j]/M1[m][dx//2][dy//2][delta_z//2]
            if not (Ratio_lowest < id < Ratio_highest):
                cnt_abs += 1
                print(f"Indexes and ratio :{i, j, nx, ny, nz, id}")
                print(f"Values: {M2[i][j], M1[m][nx][ny][nz]}")
                print(f'Parameters for integrate: dx = {Rings[j].x-Rings[i].x}, dy = {Rings[j].y-Rings[i].y}, dz = {Rings[j].z-Rings[i].z} orientation = {Rings[j].pos, Rings[i].pos} \n')
print(cnt_abs)

print("Comparing signs of values \n")
for i in range(len(M1)):
    break
    for j in range(len(M1[i])):
        # Avoid zero division on diagonal
        if i != j:
            id = M2[i][j]/M1[i][j]
            if id < 0:
                cnt_sign += 1
                print(f"Indexes and ratio :{i, j, id}")
                print(f"Values: {M2[i][j], M1[i][j]}")
                print(f'Parameters for integrate: dx = {Rings[j].x-Rings[i].x}, dy = {Rings[j].y-Rings[i].y}, dz = {Rings[j].z-Rings[i].z}, orientation = {Rings[j].pos, Rings[i].pos} \n')

print(cnt_sign)

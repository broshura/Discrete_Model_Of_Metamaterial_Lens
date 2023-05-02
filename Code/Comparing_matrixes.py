# Comparing geometric matrix of metamaterial MRI lenz

from Parameters_MRI import a, a1, Radius, Rings

cnt_sign = 0        # Counter of rings with sign errors
cnt_abs = 0         # Counter of rings with large absolute error

# Matrix produced from verified matlab code
# To compare it is necessary to multiply at normalized length
with open("DATA/matrica.txt", "r") as res:
    RES = res.read()
    Strings = RES.split("\n")[:len(RES.split("\n"))-1]
    M1 = [[float(k) * Radius for k in x.split(" ")] for x in Strings]
#M1 = [[float(M[i][k]) for k in range(len(M[i]))] for i in range(len(M)-1)]


# Matrix from python - code

with open("DATA/Data-MRI.txt", "r") as res:
    RES = res.read()
    Strings = RES.split("\n")[:len(RES.split("\n")) - 1]
    M2 = [[float(k)*(a1/a) ** 1 / 1.5 for k in x.split(" ")] for x in Strings]

# Ratios between which elements are equal
Ratio_lowest = 0.98
Ratio_highest = 1.02

print("Comparing in absolute ratio \n")
for i in range(len(M1)):
    for j in range(len(M1[i])):
        # Avoid zero division on diagonal
        if i != j:
            id = M2[i][j]/M1[i][j]
            if not (Ratio_lowest < id < Ratio_highest):
                cnt_abs += 1
                print(f"Indexes and ratio :{i, j, id}")
                print(f"Values: {M2[i][j], M1[i][j]}")
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

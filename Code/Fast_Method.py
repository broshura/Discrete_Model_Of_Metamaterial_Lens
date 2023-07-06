import json
import numpy as np
from scipy.fft import fftn, ifftn
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import LinearOperator
from functools import reduce
from numpy import real, imag
from scipy.linalg import solve
from scipy.linalg import norm

from Parameters_MRI import *

def Z_coil(Omega, Data, Rings, I, Ring_number):
    Z = 0
    Ring = Rings[Ring_number]
    for i in range(len(Rings)):
        ring = Rings[i]
        dx = ring.x - Ring.x
        dy = ring.y - Ring.y
        dz = ring.z - Ring.z
        r1, r2 = ring.r, Ring.r
        pos1, pos2 = ring.pos, Ring.pos
        id = f'{dx} {dy} {dz} {r1} {r2} {pos1}{pos2}'
        Z += Data[id] * I[i]/I[Ring_number]
    return Z * 1j * Omega

def Solve(Z, V):

    I_z = fftn(ifftn( ( V["z"] ) - () ) / fftn(  ) )


# Complete Block-Toeplitz matrix to Block-circulant and remind first row for all orientation


w_max = 75 * 10 ** 6 * 2 * pi       # Maximum frequency in MGz
w_min = 55 * 10 ** 6 * 2 * pi       # Minimum frequency in MGz
n = 100                             # number of dots


print("Reading data..")
with open(f"DATA/Data-{name}.json", "r") as res:
    Data_json = res.read()
    Data = json.loads(Data_json)
    for id in Data:
        Data[id] = complex(Data[id] * (a1 / a) ** 1 * mu_0)


# Adding diagonal element to Data of mutual impedance
Data[f"0 0 0 {Radius} {Radius} xx"] = 0
Data[f"0 0 0 {Radius} {Radius} yy"] = 0
Data[f"0 0 0 {Radius} {Radius} zz"] = 0

# Modeling current in rings for this frequency range

Omega = np.linspace(w_min, w_max, n)                           # Frequency for each solving
Impedance_real = []                                              # Real part of impedance on responding ring
Impedance_imag = []                                              # Imaginary part of impedance on responding ring

print("Converting Data to circulant vector")

m = {}

for pos1 in Orientations:
    for pos2 in Orientations:
        pair = pos1 + pos2

        Max_dx = max(
            [int(data.split()[0]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])
        Max_dy = max(
            [int(data.split()[1]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])
        Max_dz = max(
            [int(data.split()[2]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])

        Min_dx = min(
            [int(data.split()[0]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])
        Min_dy = min(
            [int(data.split()[1]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])
        Min_dz = min(
            [int(data.split()[2]) for data in Data if data.split()[-1:-4:-1] == [str(pair), str(Radius), str(Radius)]])

        dx_0 = Min_dx + N[f"x{pos1}"] * a - a
        dy_0 = Min_dy + N[f"y{pos1}"] * a - a
        dz_0 = Min_dz + N[f"z{pos1}"] * a - a

        DZ = tuple(range(dz_0, Max_dz + a, a)) + tuple(range(Min_dz, dz_0, a))
        DY = tuple(range(dy_0, Max_dy + a, a)) + tuple(range(Min_dy, dy_0, a))
        DX = tuple(range(dx_0, Max_dx + a, a)) + tuple(range(Min_dx, dx_0, a))
        m[pair] = np.array([])

        for dz in DZ:
            for dy in DY:
                for dx in DX:
                    id = f"{dx} {dy} {dz} {Radius} {Radius} {pair}"
                    m[pair] = np.append(m[pair], Data[id])

        nz = N[f"z{pos1}"] + N[f"z{pos2}"] - 1
        ny = N[f"y{pos1}"] + N[f"y{pos2}"] - 1
        nx = N[f"x{pos1}"] + N[f"x{pos2}"] - 1

        m[pair] = m[pair].reshape(nz, ny, nx)






# Reading data file with geometric matrix

print("Optimized modeling")

I0 = np.zeros(Number)

for omega in Omega:
    # Transform mutual M-matrix to impedance matrix
    print(f"Frequency: {round(omega / 2/ pi/10**6, 1)} MGz")
    def func_M(I):
        for pos1 in Orientations:
            m[2*pos1][0][0][0] = (R/1j/omega + L - 1 / omega ** 2 / C)

        i = {}
        for pos1 in Orientations:
            for pos2 in Orientations:
                pair = pos1 + pos2

                nz = N[f"z{pos1}"] + N[f"z{pos2}"] - 1
                ny = N[f"y{pos1}"] + N[f"y{pos2}"] - 1
                nx = N[f"x{pos1}"] + N[f"x{pos2}"] - 1

                I0 = I[start[pos2]: end[pos2]]
                i0 = np.zeros((nz, ny, nx), dtype=complex)
                for j in range(len(I0)):
                    x = j % N[f'x{pos2}']
                    y = (j // N[f'x{pos2}']) % N[f'y{pos2}']
                    z = ((j // N[f'x{pos2}']) // N[f'y{pos2}']) % N[f'z{pos2}']
                    i0[z][y][x] = I0[j]
                i[pair] = i0

        mi = {}
        Mi = {}
        for pos1 in Orientations:
            for pos2 in Orientations:
                pair = pos1 + pos2

                mi[pair] = fftn(fftn(m[pair]) * ifftn(i[pair]))
                Mi[pair] = np.zeros(end[pos1] - start[pos1], dtype=complex)
                for z in range(N[f"z{pos1}"]):
                    for y in range(N[f"y{pos1}"]):
                        for x in range(N[f"x{pos1}"]):
                            Mi[pair][N[f'y{pos1}'] * N[f'x{pos1}'] * z + N[f'x{pos1}'] * y + x] = mi[pair][z][y][x]
            Mi[pos1] = reduce(lambda x, y: x + y, [Mi[f'{pos1}{pos}'] for pos in Orientations])

        MI = reduce(lambda x, y: np.concatenate([x, y]), [Mi[pos] for pos in Orientations])
        return MI

    M_w = LinearOperator((Number, Number), func_M)
    I, num = gmres(M_w, V, x0 = I0, tol = 10 ** -5)

    I0 = I
    Z_0 = Z_coil(omega, Data, Rings, I, -1)

    Impedance_real.append(real(Z_0 + R + L * 1j*omega + 1 / 1j / omega / C))
    Impedance_imag.append(imag(Z_0 + R + L * 1j*omega + 1 / 1j / omega / C))








print("Saving data...")
with open(f"DATA/Responding-{name}-fast.txt", 'w') as res:
    for i in range(len(Omega)):
        res.write(f"{Omega[i]} {Impedance_real[i]} {Impedance_imag[i]}\n")
print("Ended")

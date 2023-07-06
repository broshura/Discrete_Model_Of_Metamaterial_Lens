import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.optimize import curve_fit
from scipy.fft import fftn, ifftn
from functools import reduce
from random import random
from time import time as time
from numpy import sqrt

def N3(N, a, b, c):
    return a * N ** 3 + abs(b) * N ** 2 + abs(c) * N

def N2(N, a, b):
    return a * N ** 2 + b * N

def NlogN(N, a):
    return a * N * np.log(N)

Matrix_time = np.array([], dtype=np.float64)
Delta_Matrix = np.array([], dtype=np.float64())

FFT_time = np.array([], dtype=np.float64)
FFT_time_I = np.array([], dtype=np.float64)
FFT_time_ZI = np.array([], dtype=np.float64)
FFT_time_V = np.array([], dtype=np.float64)
Delta_FFT = np.array([], dtype=np.float64)

Solving_time = np.array([], dtype=np.float64)
Delta_Solving = np.array([], dtype=np.float64)

Numbers = np.array([])


Orientations = ('x', 'y', 'z')

for Nx in range(2, 20):
    for Ny in range(2, 10):
        for Nz in range(2, 10):



            Number = Nx * Ny * Nz * 3

            N = {}

            if (Nx + Ny + Nz) % 4 == 0:
                print(Nx, Ny, Nz)

            N = {}

            N["xx"], N["yx"], N["zx"] = Nx, Ny, Nz  # Number of cells on each row for x-oriented rings
            N["xy"], N["yy"], N["zy"] = Nx, Ny, Nz  # Number of cells on each row for y-oriented rings
            N["xz"], N["yz"], N["zz"] = Nx, Ny, Nz  # Number of cells on each row for y-oriented rings

            start, end = {}, {}
            start["x"], end["x"] = 0, N['xx'] * N['yx'] * N['zx']
            start["y"], end["y"] = end["x"], end['x'] + N['xy'] * N['yy'] * N['zy']
            start["z"], end["z"] = end["y"], end['y'] + N['xz'] * N['yz'] * N['zz']

            Z = np.array([[random() * 10 for i in range(Number)] for j in range(Number)])
            I = np.linspace(Number // 10, Number * 2, Number)
            V = np.linspace(Number // 10, Number, Number)

            matrix_time = 0
            matrix_time2 = 0

            for i in range(100):
                Start = np.float64(time())

                ZI = Z.dot(I)

                End = np.float64(time())

                matrix_time += (End - Start) / 100
                matrix_time2 += (End - Start) ** 2 / 100

            solving_time = 0
            solving_time2 = 0

            for i in range(10):
                Start = np.float64(time())

                Vi = solve(Z, V)

                End = np.float64(time())

                solving_time += (End - Start) / 10
                solving_time2 += (End - Start) ** 2 / 10

            Z = {}

            for pos1 in Orientations:
                for pos2 in Orientations:
                    pair = pos1 + pos2
                    nz = N[f"z{pos1}"] + N[f"z{pos2}"] - 1
                    ny = N[f"y{pos1}"] + N[f"y{pos2}"] - 1
                    nx = N[f"x{pos1}"] + N[f"x{pos2}"] - 1

                    Z[pair] = np.zeros(nx * ny * nz).reshape(nz, ny, nx)

            fft_time = 0
            fft_time2 = 0

            fft_time_i = 0
            fft_time_i2 = 0
            fft_time_zi = 0
            fft_time_zi2 = 0
            fft_time_v = 0
            fft_time_v2 = 0

            for j in range(10):

                Start = np.float64(time())

                Start_i = np.float64(time())

                for pos1 in Orientations:
                    Z[2 * pos1][0][0][0] = 1 + 2 + 3

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

                End_i = np.float64(time())

                Start_zi = np.float64(time())

                zi = {}
                Zi = {}
                for pos1 in Orientations:
                    for pos2 in Orientations:
                        pair = pos1 + pos2

                        zi[pair] = fftn(fftn(Z[pair]) * ifftn(i[pair]))

                End_zi = np.float64(time())

                Start_v = np.float64(time())

                for pos1 in Orientations:
                    for pos2 in Orientations:
                        pair = pos1 + pos2
                        Zi[pair] = np.zeros(end[pos1] - start[pos1], dtype=complex)
                        for z in range(N[f"z{pos1}"]):
                            for y in range(N[f"y{pos1}"]):
                                for x in range(N[f"x{pos1}"]):
                                    Zi[pair][N[f'y{pos1}'] * N[f'x{pos1}'] * z + N[f'x{pos1}'] * y + x] = \
                                    zi[pair][z][y][x]
                    Zi[pos1] = reduce(lambda x, y: x + y, [Zi[f'{pos1}{pos}'] for pos in Orientations])

                MI = reduce(lambda x, y: np.concatenate([x, y]), [Zi[pos] for pos in Orientations])

                End_v = np.float64(time())

                End = np.float64(time())

                fft_time += (End - Start) / 10
                fft_time2 += (End - Start) ** 2 / 10

                fft_time_i += (End_i - Start_i) / 10
                fft_time_i2 += (End_i - Start_i) ** 2 / 10

                fft_time_zi += (End_zi - Start_zi) / 10
                fft_time_zi2 += (End_zi - Start_zi) ** 2 / 10

                fft_time_v += (End_v - Start_v) / 10
                fft_time_v2 += (End_v - Start_v) ** 2 / 10

            Matrix_time = np.append(Matrix_time, matrix_time)
            Delta_Matrix = np.append(Delta_Matrix, sqrt(matrix_time2 - matrix_time ** 2))
            FFT_time = np.append(FFT_time, fft_time)
            Delta_FFT = np.append(Delta_FFT, sqrt(fft_time2 - fft_time ** 2))
            FFT_time_I = np.append(FFT_time_I, fft_time_i)
            Delta_FFT_I = np.append(Delta_FFT, sqrt(fft_time_i2 - fft_time_i ** 2))
            FFT_time_ZI = np.append(FFT_time_ZI, fft_time_zi)
            Delta_FFT_ZI = np.append(Delta_FFT, sqrt(fft_time_zi2 - fft_time_zi ** 2))
            FFT_time_V = np.append(FFT_time_V, fft_time_v)
            Delta_FFT_V = np.append(Delta_FFT, sqrt(fft_time_v2 - fft_time_v ** 2))
            Solving_time = np.append(Solving_time, solving_time)
            Delta_Solving = np.append(Delta_Solving, sqrt(solving_time2 - solving_time ** 2))
            Numbers = np.append(Numbers, Number)

plt.figure(figsize=(12, 8))

a_fft, err_fft = curve_fit(NlogN, Numbers, FFT_time)
a_fft_n2, b_fft_n2 = curve_fit(N2, Numbers, FFT_time)[0]
a_matrix, b_matrix = curve_fit(N2, Numbers, Matrix_time)[0]
a_matrix_n3, b_matrix_n3, c_matrix_n3 = curve_fit(N3, Numbers, Matrix_time)[0]
a_solving_n2, b_solving_n2 = curve_fit(N2, Numbers, Solving_time)[0]
a_solving_n3, b_solving_n3, c_solving_n3 = curve_fit(N3, Numbers, Solving_time)[0]

numbers = np.linspace(2, max(Numbers), 100)

Fit_fft = np.array([NlogN(n, a_fft) for n in numbers])
Fit_fft_n2 = np.array([N2(n, a_fft_n2, b_fft_n2) for n in numbers])
Fit_matrix = np.array([N2(n, a_matrix, b_matrix) for n in numbers])
Fit_matrix_n3 = np.array([N3(n, a_matrix_n3, b_matrix_n3, c_matrix_n3) for n in numbers])
Fit_solving_n2 = np.array([N2(n, a_solving_n2, b_solving_n2) for n in numbers])
Fit_solving_n3 = np.array([N3(n, a_solving_n3, b_solving_n3, c_solving_n3) for n in numbers])


plt.subplot(2, 1, 1)


plt.scatter(Numbers, FFT_time, label = "FFT way", color = 'blue')
plt.scatter(Numbers, Matrix_time, label = "Straight way", color = 'red')
plt.scatter(Numbers, Solving_time, label = "Solving", color = 'yellow')

plt.plot(numbers, Fit_fft, label = f'FFT fit $N \log N$', linestyle = '-', color = 'black')
plt.plot(numbers, Fit_matrix, label = f'Fit Straight $N^2$', linestyle = '-', color = 'pink')
#plt.plot(numbers, Fit_solving_n2, label = 'Fit solving $N^2$', linestyle = '-', color = 'green')
plt.plot(numbers, Fit_solving_n3, label = 'Fit solving $N^3$', linestyle = '-', color = '#00A86B')

plt.legend()
plt.grid(True)
plt.title("Dependence of time on dimension of matrix")
#plt.xlabel('Dimension')
plt.ylabel('Time, s')
plt.xlim(0, max(Numbers)*1.1)
plt.ylim(-Solving_time[-1]/10, Solving_time[-1] * 1.1)


plt.subplot(2, 3, 4)


plt.scatter(Numbers, FFT_time, label = "FFT way", color = 'blue')
plt.scatter(Numbers, Matrix_time, label = "Straight way", color = 'red')
plt.plot(numbers, Fit_fft, label = f'FFT fit $N \log N$', linestyle = '-', color = 'black')
plt.plot(numbers, Fit_matrix, label = f'Fit straight $N^2$', linestyle = '-', color = 'pink')

plt.legend()
plt.grid(True)
plt.title("Comparing multiplying")
plt.xlabel('Dimension')
plt.ylabel('Time, s')
plt.xlim(0, max(Numbers)*1.1)
plt.ylim(-FFT_time[-1]/10, FFT_time[-1] * 1.1)


plt.subplot(2, 3, 5)


plt.scatter(Numbers, FFT_time, label = "FFT way", color = 'blue')
plt.plot(numbers, Fit_fft, label = f'FFT fit $N \log N$', linestyle = '-', color = 'black')
plt.plot(numbers, Fit_fft_n2, label = f'FFT fit $N^2$', linestyle = '-', color = 'yellow')

plt.legend()
plt.grid(True)
plt.title("FFT fitting")
plt.xlabel('Dimension')
plt.ylabel('Time, s')
plt.xlim(0, max(Numbers)*1.1)
plt.ylim(-FFT_time[-1]/10, FFT_time[-1] * 1.1)


plt.subplot(2, 3, 6)


plt.scatter(Numbers, Matrix_time, label = "Straight way", color = 'red')
plt.plot(numbers, Fit_matrix, label = f'Fit straight $N^2$', linestyle = '-', color = 'pink')
plt.plot(numbers, Fit_matrix_n3, label =f'Fit straight $N^3$', linestyle = '-', color = "black")

plt.legend()
plt.grid(True)
plt.title("Straight Multiplying")
plt.xlabel('Dimension')
plt.ylabel('Time, s')
plt.xlim(0, max(Numbers)*1.1)
plt.ylim(-Matrix_time[-1]/10, Matrix_time[-1] * 1.1)


plt.savefig(f"Plots/Time estimation.png")
plt.show()





plt.figure(figsize=(11, 7))

a_i = curve_fit(lambda x, a: a * x, Numbers, FFT_time_I)[0]
a_zi = curve_fit(lambda x, a: a * x, Numbers, FFT_time_ZI)[0]
a_v = curve_fit(lambda x, a: a * x, Numbers, FFT_time_V)[0]

Fit_i = [a_i * n for n in numbers]
Fit_v = [a_v * n for n in numbers]
Fit_zi = [a_zi * n for n in numbers]

plt.scatter(Numbers, FFT_time, label = "Total time", color = 'blue')
plt.scatter(Numbers, FFT_time_I, label = f'Reshape $I$', color = 'red')
plt.scatter(Numbers, FFT_time_ZI, label = f'Multiplying $ZI$', color = 'black')
plt.scatter(Numbers, FFT_time_V, label = f'Reshape $V$', color = 'yellow')
plt.plot(numbers, Fit_i, label = f'Average I', linestyle = '-', color = 'red')
plt.plot(numbers, Fit_zi, label = f'Average ZI', linestyle = '-', color = 'black')
plt.plot(numbers, Fit_v, label = f'Average V', linestyle = '-', color = 'yellow')
#plt.plot(numbers, Fit_fft, label = f'FFT fit $N \log N$', linestyle = '-', color = 'pink')
#plt.plot(numbers, Fit_fft_n2, label = f'FFT fit $N^2$', linestyle = '-', color = '#00A86B')

plt.legend()
plt.grid(True)
plt.title("FFT parts")
plt.xlabel('Dimension')
plt.ylabel('Time, s')
plt.xlim(0, max(Numbers)*1.1)
plt.ylim(-FFT_time[-1]/10, FFT_time[-1] * 1.1)


plt.savefig(f"Plots/FFT parts.png")
plt.show()

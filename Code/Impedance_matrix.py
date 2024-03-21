# Calculating geometry matrix M
import json

import numpy as np
from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import special
from tqdm import tqdm

K = special.ellipk  # Сomplete elliptic integral of the first kind
E = special.ellipe  # Сomplete elliptic integral of the second kind

mu_0 = 4 * pi * 10 ** -7

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
        R = r1 + width / 2
        r = r1 - width / 2
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
        R = r1 + width / 2
        r = r1 - width / 2
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, r))
        L_2, err_2 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
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

    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    r1 = First_ring.r
    r2 = Second_ring.r
    w = First_ring.w

    # To avoid calculating integrals with same params each time
    # there is a dictionary with all parameters and values

    id_1 = f"{dx} {dy} {dz} {r1} {r2} {First_ring.pos}{Second_ring.pos}"
    id_2 = f"{-dx} {-dy} {-dz} {r2} {r1} {Second_ring.pos}{First_ring.pos}"

    if id_1 in Data:
        return Data[id_1]
    elif id_2 in Data:
        return Data[id_2]
    elif dx == 0 and dy == 0 and dz == 0:
        Data[id_1] = 0
        return 0

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
    Data[id_1], Data[id_2] = l * mu_0, l * mu_0
    return l * mu_0

# Calculating mutual inductance for each pair

def Matrix(Rings1, Rings2, Data = {}):
    M = np.zeros((len(Rings1), len(Rings2)))
    for n in tqdm(range(len(Rings1))):
        for m in range(len(Rings2)):
            R1 = Rings1[n]
            R2 = Rings2[m]
            M[n][m] = Mnm(R1, R2, Data)
    return M

# Caclulating diagonal of Z matrix using gradient

def Z_0(Params, grad = [0, 0, 0]):
    L, C, R, N = Params['L'], Params['C'], Params['R'], Params['N']
    Number = np.sum([N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx'] for pos in N])
    Cgrad = np.empty(Number)

    # Making start and end points for each orientation
    start = 0
    end = 0
    for pos, i in zip(N, np.arange(len(N))):
        end += N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx']
        
        z, y, x = np.meshgrid(np.arange(N[pos]['nz']),
                              np.arange(N[pos]['ny']),
                              np.arange(N[pos]['nx']),
                              indexing = 'ij')
        gradx = x * grad[0] / N[pos]['nx']
        grady = y * grad[1] / N[pos]['ny']
        gradz = z * grad[2] / N[pos]['nz']

        Omega_0 = (1 + gradz + grady + gradx)/np.sqrt(L * C).ravel()
        Cgrad[start:end] = 1/L/(Omega_0 ** 2)

        start += N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx']
    return lambda omega: R + 1j * omega * L + 1 / (1j * omega * Cgrad)

# Caclulating diagonal of M matrix using gradient

def M_0(Params, grad = [0, 0, 0]):
    L, C, R, N = Params['L'], Params['C'], Params['R'], Params['N']
    Number = np.sum([N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx'] for pos in N])
    Cgrad = np.empty(Number)

    # Making start and end points for each orientation
    start = 0
    end = 0
    for pos, i in zip(N, np.arange(len(N))):
        end += N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx']
        
        z, y, x = np.meshgrid(np.arange(N[pos]['nz']),
                              np.arange(N[pos]['ny']),
                              np.arange(N[pos]['nx']),
                              indexing = 'ij')
        gradx = x * grad[0] / N[pos]['nx']
        grady = y * grad[1] / N[pos]['ny']
        gradz = z * grad[2] / N[pos]['nz']

        Omega_0 = (1 + gradz + grady + gradx)/np.sqrt(L * C)
        Cgrad[start:end] = 1/L/(Omega_0.ravel() ** 2)

        start += N[pos]['nz'] * N[pos]['ny'] * N[pos]['nx']

    return lambda omega: R / (1j * omega) + L - 1 / (omega ** 2 * Cgrad)

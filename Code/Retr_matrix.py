import json

import numpy as np
from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import special
from tqdm import tqdm
from typing import List
from Ring_Class import Ring

K = special.ellipk  # Сomplete elliptic integral of the first kind
E = special.ellipe  # Сomplete elliptic integral of the second kind

c=3e8
mu_0 = 4 * np.pi * 1e-7

def L_parallel(omega, dx:float, dy:float, dz:float, r1:float, r2:float, width:float = 0):
    k = omega/c
    def A_real(phi_2, dx, dy, dz, r_1, r_2, k):
    #r_1 = R_m, r_2 = R_n

        def dA_x(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2) - r_1 * np.cos(phi_1)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) - r_1 * np.sin(phi_1)) ** 2 + dz ** 2)
            return mu_0 / 4 /np.pi * r_1  * np.cos(phi_1) * np.cos(k*dR) / dR

        def dA_y(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2) - r_1 * np.cos(phi_1)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) - r_1 * np.sin(phi_1)) ** 2 + dz ** 2)
            return  mu_0 / 4 /np.pi * r_1  * np.sin(phi_1) * np.cos(k*dR) / dR

        A_x = integrate.quad(dA_x, 0, 2 * np.pi)[0]
        A_y = integrate.quad(dA_y, 0, 2 * np.pi)[0]

        return np.array([A_x, A_y])

    def A_imag(phi_2, dx, dy, dz, r_1, r_2, k):
    #r_1 = R_m, r_2 = R_n

        def dA_x(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2) - r_1 * np.cos(phi_1)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) - r_1 * np.sin(phi_1)) ** 2 + dz ** 2)
            return mu_0 / 4 /np.pi * r_1  * np.cos(phi_1) * np.sin(k*dR) / dR

        def dA_y(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2) - r_1 * np.cos(phi_1)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) - r_1 * np.sin(phi_1)) ** 2 + dz ** 2)
            return  mu_0 / 4 /np.pi * r_1  * np.sin(phi_1) * np.sin(k*dR) / dR

        A_x = integrate.quad(dA_x, 0, 2 * np.pi)[0]
        A_y = integrate.quad(dA_y, 0, 2 * np.pi)[0]

        return np.array([A_x, A_y])

    def L_my(dx:float, dy:float, dz:float, r1:float, r2:float, k, width:float = 0) -> float:

        def f1(alpha):
            A = A_real(alpha, dx, dy, dz, r1, r2, k)
            return r2 * A[0] * np.cos(alpha) + r2 * A[1] * np.sin(alpha)

        def f2(alpha):
            A = A_imag(alpha, dx, dy, dz, r1, r2, k)
            return r2 * A[0] * np.cos(alpha) + r2 * A[1] * np.sin(alpha)

        imag = integrate.quad(f2, 0.0, 2*pi)[0]
        real = integrate.quad(f1, 0.0, 2*pi)[0]

        return real + 1j * imag
    
    return L_my(dx, dy, dz, r1, r2, k, width = 0)


def L_orthogonal(omega, dx:float, dy:float, dz:float, r1:float, r2:float, width:float = 0):
    k = omega/c
    def A_real(phi_2, dx, dy, dz, r_1, r_2, k):
    #r_1 = R_m, r_2 = R_n

        def dA_x(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) + r_1 * np.sin(phi_1)) ** 2 + (dz+r_2 * np.cos(phi_2)) ** 2)
            return mu_0 / 4 /np.pi * r_1  * np.cos(phi_1) * np.cos(k*dR) / dR

        def dA_y(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) + r_1 * np.sin(phi_1)) ** 2 + (dz+r_2 * np.cos(phi_2)) ** 2)
            return  mu_0 / 4 /np.pi * r_1  * np.sin(phi_1) * np.cos(k*dR) / dR

        A_x = integrate.quad(dA_x, 0, 2 * np.pi)[0]
        A_y = integrate.quad(dA_y, 0, 2 * np.pi)[0]

        return np.array([A_x, A_y])

    def A_imag(phi_2, dx, dy, dz, r_1, r_2, k):
    #r_1 = R_m, r_2 = R_n

        def dA_x(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) + r_1 * np.sin(phi_1)) ** 2 + (dz+r_2 * np.cos(phi_2)) ** 2)
            return mu_0 / 4 /np.pi * r_1  * np.cos(phi_1) * np.sin(k*dR) / dR

        def dA_y(phi_1):
            dR = np.sqrt((dx + r_2 * np.cos(phi_2)) ** 2 + \
            (dy + r_2 * np.sin(phi_2) + r_1 * np.sin(phi_1)) ** 2 + (dz+r_2 * np.cos(phi_2)) ** 2)
            return  mu_0 / 4 /np.pi * r_1  * np.sin(phi_1) * np.sin(k*dR) / dR

        A_x = integrate.quad(dA_x, 0, 2 * np.pi)[0]
        A_y = integrate.quad(dA_y, 0, 2 * np.pi)[0]

        return np.array([A_x, A_y])

    def L_my(dx:float, dy:float, dz:float, r1:float, r2:float, k, width:float = 0) -> float:

        def f1(alpha):
            A = A_real(alpha, dx, dy, dz, r1, r2, k)
            return r2 * A[0] * np.cos(alpha) + r2 * A[1] * np.sin(alpha)

        def f2(alpha):
            A = A_imag(alpha, dx, dy, dz, r1, r2, k)
            return r2 * A[0] * np.cos(alpha) + r2 * A[1] * np.sin(alpha)

        imag = integrate.quad(f2, 0.0, 2*pi)[0]
        real = integrate.quad(f1, 0.0, 2*pi)[0]

        return real + 1j * imag
    
    return L_my(dx, dy, dz, r1, r2, k, width = 0)

def Mnm(omega, First_ring:Ring, Second_ring:Ring, Data:dict = {}) -> float:
    """Calculates mutual inductance between two rings with both parallel 
    and orthogonal orientation. Data is a dictionary with all parameters
    and values to avoid calculating integrals with same params each time.

    Parameters
    ----------
    First_ring : Ring
        first ring
    Second_ring : Ring
        second ring
    Data : dict, optional
        Dictionary with mutual inductance, by default {}

    Returns
    -------
    float
        mutual inductance between two rings
    """    
    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    r1 = First_ring.r
    r2 = Second_ring.r
    w = First_ring.w * 0.7

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

    Data[id_1], Data[id_2] = [l * mu_0] * 2
    return l * mu_0

# Calculating mutual inductance for each pair

def Matrix(omega, rings:List[Ring], Data:dict = {}) -> np.ndarray:
    """Calculates mutual inductance matrix for all rings
    and represents it as a square matrix

    Parameters
    ----------
    rings : List[Ring]
        list of rings
    Data : dict, optional
        Dictionary with mutual inductance, by default {}

    Returns
    -------
    np.ndarray
        mutual inductance matrix with zeros on the diagonal
    """    
    M = np.zeros((len(rings), len(rings)))
    for n in tqdm(range(len(rings))):
        for m in range(n, len(rings)):
            R1 = rings[n]
            R2 = rings[m]
            M[n][m] = Mnm(omega, R1, R2, Data)
            M[m][n] = M[n][m]
    return M
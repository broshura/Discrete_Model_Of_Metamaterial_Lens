import numpy as np
import scipy
from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import linalg
from scipy import special

from Ring_Class import Ring

K = special.ellipkinc       #  Incomplete elliptic integral of the first kind
E = special.ellipe          #  Incomplete elliptic integral of the second kind

def L_parallel(alpha, dx, dy, dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa2 = 4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2)
        A = 1
        return A * r_2 * (r_2 + db * cos(alpha))/dp
    L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r_1, r_2))
    return L
def L_orthogonal(alpha, dx, dy ,dz, r_1, r_2):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa2 = 4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)))
        A = 1
        return A * r_2 * dy * cos(alpha) / dp
    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r_1, r_2))
    return L

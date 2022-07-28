import numpy as np
import scipy
from numpy import sqrt, cos
from scipy import integrate
from scipy import linalg
from scipy import special

from Ring_Class import Ring

K = special.ellipkinc       #  Incomplete elliptic integral of the first kind
E = special.ellipe          #  Incomplete elliptic integral of the second kind

def Deltap_parallel(r, de):
    pass

def Deltap_orthogonal():
    pass

def Kappa2_parallel(r, dp, dz):
    return 4 * r * dp / ( (dp + r) ** 2 + dz ** 2 )

def Kappa2_orthogonal(alpha, r_1, r_2, dp, dz):
    return  4 * r_1 * dp / ( (dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2 )

def A(alpha,):
    pass



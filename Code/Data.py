# Calculating geometry matrix M and saving as large data file for optimize next computation


from numpy import sqrt, cos, sin, pi
from scipy import integrate
from scipy import special

from Parameters import *

K = special.ellipk       #  Сomplete elliptic integral of the first kind
E = special.ellipe       #  Сomplete elliptic integral of the second kind

# Computing for parallel-oriented rings
def L_parallel(dx, dy, dz, r, width = 0):
    # Define function to integrate over first defined parameter
    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2))/kappa
        return A * r_2 * (r_2 + db * cos(alpha))/dp
    #Considering stripe width
    if width:
        R = r + w / 2
        r = r - w / 2
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r, r))
        L_2, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
        L_3, err_3 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, R, R))
        L = (L_1 + 2*L_2 + L_3)/4
    else:
        L, err = integrate.quad(dl, 0, 2 * pi, args= (dx, dy, dz, r, r))
    return L

# Computing for orthogonal-oriented rings
def L_orthogonal(dx, dy ,dz, r, width):
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2))/kappa
        return A * r_2 * dy * cos(alpha) / dp

    # Considering stripe width
    if width:
        R = r + w / 2
        r = r - w / 2
        L_1, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, r))
        L_2, err_1 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
        L_3, err_3 = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, R, R))
        L = (L_1 + 2 * L_2 + L_3) / 4
    else:
        L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, r))
    return L

# Function for responding ring because of other radius

def L_rp(dx, dy, dz, r, R):
    # Define function to integrate over first defined parameter
    def dl(alpha, dx, dy, dz, r_1, r_2):
        db = sqrt(dx ** 2 + dy ** 2)
        dp = sqrt(r_2 ** 2 + db ** 2 + 2 * r_2 * db * cos(alpha))
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + dz ** 2))
        A = sqrt(r_1 / dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2)) / kappa
        return A * r_2 * (r_2 + db * cos(alpha)) / dp
    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
    return L

def L_ro(dx, dy, dz, r, R):
    # Define function to integrate over first defined parameter
    def dl(alpha, dx, dy, dz, r_1, r_2):
        dp = sqrt((dx - r_2 * sin(alpha)) ** 2 + dy ** 2)
        kappa = sqrt(4 * r_1 * dp / ((dp + r_1) ** 2 + (dz - r_2 * cos(alpha)) ** 2))
        A = sqrt(r_1/dp) * ((2 - kappa ** 2) * K(kappa ** 2) - 2 * E(kappa ** 2))/kappa
        return A * r_2 * dy * cos(alpha) / dp

    L, err = integrate.quad(dl, 0, 2 * pi, args=(dx, dy, dz, r, R))
    return L

# Computing for any pair
def Mnm(First_ring, Second_ring, Data = {}):
    dx = Second_ring.x - First_ring.x
    dy = Second_ring.y - First_ring.y
    dz = Second_ring.z - First_ring.z
    r = First_ring.r
    r2 = Second_ring.r

    # To avoid calculating integrals with same params each time
    # there is a dictionary with all parameters and values

    id = (dx, dy, dz, First_ring.pos, Second_ring.pos)
    if id in Data:
        return Data[id]

    #Consider responding ring
    if r2 != r:
        wid = R
        L_p = L_rp
        L_o = L_ro
    else:
        wid = w
        L_p = L_parallel
        L_o = L_orthogonal

    # Consider all types of parallel orientation and symmetry for x-z axes
    if First_ring.pos == Second_ring.pos:
        if First_ring.pos == "z":                           # Z-oriented rings
            l = L_p(dx, dy, dz, r, wid)
            Data[id] = l
        elif First_ring.pos == "y":                         # Y-oriented rings
            l = L_p(dx, -dz, dy, r, wid)
            Data[id] = l
        else:                                               # X-oriented rings
            l = L_p(-dz, dy, dx, r, wid)
            Data[id] = l

    # Consider all types of orthogonal orientation
    else:
        if First_ring.pos == "z":
            if Second_ring.pos == "y":                      # Z-Y oriented pair
                l = L_o(dx, dy, dz, r, wid)
                Data[id] = l
            else:                                           # Z-X oriented pair
                l = L_o(-dy, dx, dz, r, wid)
                Data[id] = l
        elif First_ring.pos == "y":
            if Second_ring.pos == "z":                      # Y-Z oriented pair
                l = L_o(-dx, -dy, -dz, r, wid)
                Data[id] = l
            else:                                           # Y-X oriented pair
                l = L_o(dz, dy, dx, r, wid)
                Data[id] = l
        elif First_ring.pos == "x":
            if Second_ring.pos == "z":                      # X-Z oriented pair
                l = L_o(dy, -dx, -dz, r, wid)
                Data[id] = l
            else:                                           # X-Y oriented pair
                l = L_o(dz, -dy, dx, r, wid)
                Data[id] = l
    return l

# Calculating for each pair

M = np.eye(Number) * 0
Data = {}
for n in range(Number):
    for m in range(Number):
        if n > m:
            R1 = Rings[n]
            R2 = Rings[m]
            M[n][m] = Mnm(R1, R2, Data)
            M[m][n] = M[n][m]

# Writing table in Data-file, divide string by \n and elements by " "
with open("DATA/Data.txt", "w") as res:
    for i in range(Number):
        res.write(" ".join(map(str, M[i])) + "\n")\

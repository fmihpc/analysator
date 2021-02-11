import pytools as pt
import numpy as np
from scipy.special import legendre
from time import perf_counter

P = tuple(map(lambda n: legendre(n), (0, 1, 2, 3, 4)))

# y and z are simply cyclic rotations of this
def interpolate_x(a, x, y, z):
    B_x = 0
    for i in range(5):
        for j in range(4):
            for k in range(4):
                B_x += a[i][j][k] * P[i](x) * P[j](y) * P[k](z)
    return B_x

def interpolate_y(b, x, y, z):
    return interpolate_x(b, y, z, x)

def interpolate_z(c, x, y, z):
    return interpolate_x(c, z, x, y)

# Returns moments as [ijk, order, x, y, z]
def solve_moments_from_B(fg_b):
    x_moments = solve_x_moments(fg_b[:, :, :, 0])
    y_moments = solve_y_moments(fg_b[:, :, :, 1])
    z_moments = solve_z_moments(fg_b[:, :, :, 2])
    return (x_moments, y_moments, z_moments)

# Returns moments for x as [order, x, y, z]
# With order (0, y, z, yy, yz, zz, yyy, yyz, yzz, zzz)
# y and z as cyclic rotations
def solve_x_moments(B_x):
    x_moments = np.zeros((10,) + np.shape(B_x))
    x_moments[0] = B_x
    x_moments[1:3] = np.gradient(B_x)[1:3]
    start = 3
    for i in range(2):
        j = i + 1
        x_moments[start:start+3-j] = np.gradient(x_moments[1+i])[j:3]
        start += (3-j)
    for i in range(3):
        j = i + 1 if i < 2 else 2
        x_moments[start:start+3-j] = np.gradient(x_moments[3+i])[j:3]
        start += (3-j)
    return x_moments

def solve_y_moments(B_y):
    return np.transpose(solve_x_moments(np.transpose(B_y, (1, 2, 0))), (0, 3, 1, 2))

def solve_z_moments(B_z):
    return np.transpose(solve_x_moments(np.transpose(B_z, (2, 0, 1))), (0, 2, 3, 1))

# Solves a, b, c components with a[x, y, z], b[y, z, x] and c[z, x, y] up to given order
# Input B should be the output of solve_moments_from_B
def solve_coefficients(B_moments, xyz, order = 4):
    abc = np.zeros((3, 5, 4, 4))
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    # 4th order
    if (order > 3):
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if k else y + 1, z if j else z + 1)
            print(coords)
            a = abc[i]
            Bx = B_moments[i]

            a[0][3][0] = 1/2 * (Bx[6][xyz] + Bx[6][coords])
            a[0][2][1] = 1/2 * (Bx[7][xyz] + Bx[7][coords])
            a[0][1][2] = 1/2 * (Bx[8][xyz] + Bx[8][coords])
            a[0][0][3] = 1/2 * (Bx[9][xyz] + Bx[9][coords])

            a[1][3][0] = Bx[6][xyz] - Bx[6][coords]
            a[1][2][1] = Bx[7][xyz] - Bx[7][coords]
            a[1][1][2] = Bx[8][xyz] - Bx[8][coords]
            a[1][0][3] = Bx[9][xyz] - Bx[9][coords]

        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            a = abc[i]
            b = abc[j]
            c = abc[k]
            
            # Should be correct, check later
            a[4][0][0] = -1/4 * (b[1][0][3] + c[1][3][0])
            a[3][1][0] = -7/30 * c[1][2][1]
            a[3][0][1] = -7/30 * b[1][1][2]
            a[2][2][0] = -3/20 * c[1][1][2]
            a[2][0][2] = -3/20 * b[1][2][1]

    # 3rd order
    if (order > 2):
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if k else y + 1, z if j else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[0][2][0] = 1/2 * (Bx[3][xyz] + Bx[3][coords]) - 1/6 * a[2][2][0]
            a[0][1][1] = 1/2 * (Bx[4][xyz] + Bx[4][coords])
            a[0][0][2] = 1/2 * (Bx[5][xyz] + Bx[5][coords]) - 1/6 * a[2][0][2]

            a[1][2][0] = Bx[3][xyz] - Bx[3][coords]
            a[1][1][1] = Bx[4][xyz] - Bx[4][coords]
            a[1][0][2] = Bx[5][xyz] - Bx[5][coords]
        
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            a = abc[i]
            b = abc[j]
            c = abc[k]

            # Should be correct, check later
            a[3][0][0] = -1/3 * (b[1][0][2] + c[1][2][0])
            a[2][1][0] = -1/4 * c[1][1][1]
            a[2][0][1] = -1/4 * b[1][1][1]

    # 2nd order
    if (order > 1):
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if k else y + 1, z if j else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[0][1][0] = 1/2 * (Bx[1][xyz] + Bx[1][coords]) - 1/6 * a[2][1][0]
            a[0][0][1] = 1/2 * (Bx[2][xyz] + Bx[2][coords]) - 1/6 * a[2][0][1]

            a[1][1][0] = (Bx[1][xyz] - Bx[1][coords]) - 1/10 * a[3][1][0]
            a[1][0][1] = (Bx[2][xyz] - Bx[2][coords]) - 1/10 * a[3][0][1]
        
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            a = abc[i]
            b = abc[j]
            c = abc[k]

            # Should be correct, check later
            a[2][0][0] = -1/2 * (b[1][0][1] + c[1][1][0]) - 3/35 * a[4][0][0] - 1/20 * (b[3][0][1] + c[3][1][0])

    # 1st order
    if (order > 0):
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if k else y + 1, z if j else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[1][0][0] = (Bx[0][xyz] - Bx[0][coords]) - 1/10 * a[3][0][0]
    
    # 0th order
    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3
        coords = (x if i else x + 1, y if k else y + 1, z if j else z + 1)
        a = abc[i]
        Bx = B_moments[i]

        a[0][0][0] = 1/2 * (Bx[0, x, y, z] + Bx[0][coords]) - 1/6 * a[2][0][0] - 1/70 * a[4][0][0]

    # Check constraint:
    test = abc[0][1][0][0] + abc[1][1][0][0] + abc[2][1][0][0] + 1/10 * (abc[0][3][0][0] + abc[1][3][0][0] + abc[2][3][0][0])
    #if abs(test) > 0:
    #    print("Something went wrong, sum (17) is " + str(test))

    return abc

def center_value(B_moments, xyz, order=4):
    abc = solve_coefficients(B_moments, xyz, order)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    return [interpolate_x(abc[0], 0.5, 0.5, 0.5), interpolate_y(abc[1], 0.5, 0.5, 0.5), interpolate_z(abc[2], 0.5, 0.5, 0.5)]

def center_values(B_moments, coords, order=4):
    return np.array(map(lambda xyz: center_value(B_moments, xyz, order), coords))

def all_center_values(B_moments, order=4):
    return

# Code for testing
f = pt.vlsvfile.VlsvReader('/wrk/users/lkotipal/TEST/bulk.0000057.vlsv')
fg_b = f.read_fsgrid_variable('fg_b')
print(np.shape(fg_b))
fg_b_vol = f.read_fsgrid_variable('fg_b_vol')
print('File read!')
t = perf_counter()
B_moments = solve_moments_from_B(fg_b)
print(f'B_moments solved in {perf_counter() - t} seconds!')
for i in range(5):
    t = perf_counter()
    print(center_value(B_moments, (50, 50, 50), i))
    #abc = solve_coefficients(B_moments, 5, 5, 5, i)
    #print([interpolate_x(abc[0], 0.5, 0.5, 0.5), interpolate_y(abc[1], 0.5, 0.5, 0.5), interpolate_z(abc[2], 0.5, 0.5, 0.5)])
    print(f'Coefficients up to order {i} solved in {perf_counter() - t} seconds!')
print(fg_b_vol[50, 50, 50])
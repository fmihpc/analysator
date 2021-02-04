import pytools as pt
import numpy as np
from scipy.special import legendre

P = tuple(map(lambda n: legendre(n), (0, 1, 2, 3, 4)))

# y and z are simply cyclic rotations of this
def interpolate_x(a, x, y, z):
    B_x = 0
    for i in range(5):
        for j in range(4):
            for k in range(4):
                B_x += a[i][j][k] * P[i](x) * P[j](y) * P[k](z)
    return B_x

# Returns moments as [ijk, order, x, y, z]
def solve_moments_from_B(fg_b):
    # Cyclically rotate each component and rotate coordinates back
    x_moments = np.transpose(solve_x_moments(np.transpose(fg_b[:, :, :, 0], (0, 1, 2))), (0, 1, 2, 3))
    y_moments = np.transpose(solve_x_moments(np.transpose(fg_b[:, :, :, 1], (1, 2, 0))), (0, 3, 1, 2))
    z_moments = np.transpose(solve_x_moments(np.transpose(fg_b[:, :, :, 2], (2, 0, 1))), (0, 2, 3, 1))
    return (x_moments, y_moments, z_moments)

# Returns moments for x as [order, x, y, z]
# With order (0, y, z, yy, yz, zz, yyy, yyz, yzz, zzz)
# y and z as cyclic rotations
def solve_x_moments(B_x):
    x_first = np.gradient(B_x)[1:3]
    x_second = []
    for i in range(2):
        j = i + 1
        x_second += np.gradient(x_first[i])[j:3]
    x_third = []
    for i in range(3):
        if i < 2:
            j = i+1
        else:
            j = 2
        x_third += np.gradient(x_second[i])[j:3]
    return np.array([B_x] + x_first + x_second + x_third)

# Solves a, b, c components with a[x, y, z], b[y, z, x] and c[z, x, y] up to given order
# Input B should be the output of solve_moments_from_B
def solve_coefficients(fg_b, x, y, z, order = 4):
    B_moments = solve_moments_from_B(fg_b)
    abc = np.zeros((3, 5, 4, 4))

    # 4th order
    if (order > 3):
        print('4th order!')
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if j else y + 1, z if k else z + 1)
            print(coords)
            a = abc[i]
            Bx = B_moments[i]

            a[0][3][0] = 1/2 * (Bx[6, x, y, z] + Bx[6][coords])
            a[0][2][1] = 1/2 * (Bx[7, x, y, z] + Bx[7][coords])
            a[0][1][2] = 1/2 * (Bx[8, x, y, z] + Bx[8][coords])
            a[0][0][3] = 1/2 * (Bx[9, x, y, z] + Bx[9][coords])

            a[1][3][0] = Bx[6, x, y, z] - Bx[6][coords]
            a[1][2][1] = Bx[7, x, y, z] - Bx[7][coords]
            a[1][1][2] = Bx[8, x, y, z] - Bx[8][coords]
            a[1][0][3] = Bx[9, x, y, z] - Bx[9][coords]

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
        print('3rd order!')
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if j else y + 1, z if k else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[0][2][0] = 1/2 * (Bx[3, x, y, z] + Bx[3][coords]) - 1/6 * a[2][2][0]
            a[0][1][1] = 1/2 * (Bx[4, x, y, z] + Bx[4][coords])
            a[0][0][2] = 1/2 * (Bx[5, x, y, z] + Bx[5][coords]) - 1/6 * a[2][0][2]

            a[1][2][0] = Bx[3, x, y, z] - Bx[3][coords]
            a[1][1][1] = Bx[4, x, y, z] - Bx[4][coords]
            a[1][0][2] = Bx[5, x, y, z] - Bx[5][coords]
        
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
        print('2nd order!')
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if j else y + 1, z if k else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[0][1][0] = 1/2 * (Bx[1, x, y, z] + Bx[1][coords]) - 1/6 * a[2][1][0]
            a[0][0][1] = 1/2 * (Bx[2, x, y, z] + Bx[2][coords]) - 1/6 * a[2][0][1]

            a[1][1][0] = (Bx[1, x, y, z] - Bx[1][coords]) - 1/10 * a[3][1][0]
            a[1][0][1] = (Bx[2, x, y, z] - Bx[2][coords]) - 1/10 * a[3][0][1]
        
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
        print('1st order!')
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            coords = (x if i else x + 1, y if j else y + 1, z if k else z + 1)
            a = abc[i]
            Bx = B_moments[i]

            a[1][0][0] = (Bx[0, x, y, z] - Bx[0][coords]) - 1/10 * a[3][0][0]
    
    # 0th order
    for i in range(3):
        print('0th order!')
        j = (i+1) % 3
        k = (i+2) % 3
        coords = (x if i else x + 1, y if j else y + 1, z if k else z + 1)
        a = abc[i]
        Bx = B_moments[i]

        a[0][0][0] = 1/2 * (Bx[0, x, y, z] + Bx[0][coords]) - 1/6 * a[2][0][0] - 1/70 * a[4][0][0]

    # Check constraint:
    test = abc[0][1][0][0] + abc[1][1][0][0] + abc[2][1][0][0] + 1/10 * (abc[0][3][0][0] + abc[1][3][0][0] + abc[2][3][0][0])
    if abs(test) > 0:
        print("Something went wrong, sum (17) is " + str(test))

    return abc

# Code for testing
#f = pt.vlsvfile.VlsvReader('TEST/bulk.0000057.vlsv')
#fg_b = f.read_fsgrid_variable('fg_b')
#fg_b_vol = f.read_fsgrid_variable('fg_b_vol')
#for i in range(5):
#    abc = solve_coefficients(fg_b, 5, 5, 5, i)
#    print([interpolate_x(abc[0], 0.5, 0.5, 0.5), interpolate_x(abc[1], 0.5, 0.5, 0.5), interpolate_x(abc[2], 0.5, 0.5, 0.5)])
#print(fg_b_vol[5, 5, 5])
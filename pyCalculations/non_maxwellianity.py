import numpy as np
from scipy.constants import k, m_p
import pytools as pt

# population mass
m = m_p

# names of velocity moments in the bulk file
n_name = 'proton/vg_rho'         
v_name = 'proton/vg_v'
T_name = 'proton/vg_temperature'

def epsilon_M(f,cell):
    # calculates the 'non-maxwellianity' parameter for a distribution function f_i and its corresponding Maxwellian g_M
    # using the definition by Graham et al. (2021):
    # (1/2n) * integral(|f_i - g_M|) d^3v    
    # which is valued between 0 (Maxwellian) and 1 (complete deviation from a Maxwellian)
    #
    # NOTE that this is different to the definition by Greco et al. (2012):
    # (1/n) * sqrt(integral((f_i - g_M)^2) d^3v)    
    #
    # examples comparing these two definitions can be found in Settino et al. (2021)

    # -- PARAMETERS ---
    # f: bulk file
    # cell: ID of the cell containing a VDF

    # try getting the distribution from the given cell
    size = f.get_velocity_mesh_size()
    distribution = np.zeros(4 * 4 * 4 * int(size[0]) * int(size[1]) * int(size[2]))
    try:
        D = f.read_velocity_cells(cell)
        D_keys = list(D.keys())
        D_vals = list(D.values())
        distribution[D_keys] = D_vals
    except:
        print("Could not get VDF from bulk file!")
        return -1
    
    # generate a velocity space 
    vids = np.arange(4 * 4 * 4 * int(size[0]) * int(size[1]) * int(size[2]))
    v = f.get_velocity_cell_coordinates(vids)

    # get velocity moments from bulk file
    n = f.read_variable(n_name,cellids=cell)
    v0 = f.read_variable(v_name,cellids=cell)
    T = f.read_variable(T_name,cellids=cell)

    # generate a maxwellian distribution
    v_sq = np.sum((v - v0) * (v - v0), axis=-1)
    maxwellian = n * (0.5 * m / (np.pi * k * T)) ** 1.5 * np.exp(-0.5 * m * v_sq / (k * T))

    # calculate non-maxwellianity
    extent = f.get_velocity_mesh_extent()
    dvx = (extent[3] - extent[0]) / (4 * size[0])
    dvy = (extent[4] - extent[1]) / (4 * size[1])
    dvz = (extent[5] - extent[2]) / (4 * size[2])  

    epsilon = np.sum(np.abs(distribution - maxwellian))
    epsilon *= dvx * dvy * dvz
    epsilon /= 2.0 * n

    return epsilon


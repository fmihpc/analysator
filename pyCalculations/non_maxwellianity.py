import numpy as np
import sys, os
sys.path.append('/users/alhomark/analysator-spectra')
from scipy.constants import k, m_e
import pytools as pt
import matplotlib.pyplot as plt
import warnings

import random
# population mass
m = 10*m_e

# names of velocity moments in the bulk file
n_name = 'electron/vg_rho'         
v_name = 'electron/vg_v'
v_para_name = 'electron/vg_v_parallel'
v_perp_name = 'electron/vg_v_perpendicular'
T_name = 'electron/vg_temperature'
T_perp_name = 'electron/vg_t_perpendicular'
T_para_name = 'electron/vg_t_parallel'
bulkfile_for_moments = "/home/mjalho/tempo-vakka-1925/bulk/bulk002.0000050.vlsv"

def epsilon_M(f,cell,pop="proton", bulk=None, B=None,
                model="bimaxwellian",
                normorder=1, norm=2,
                dummy=None):

    ''' Calculates the 'non-maxwellianity' parameter for a distribution function f_i and its corresponding Maxwellian g_M.
    
    :param f:           VlsvReader containing VDF data
    :param cell:        CellID for the queried cell

    :kword pop:         Population to calculate the parameter for
    :kword bulk:        Bulk file name to use for moments (if available and needed)
    :kword model:       VDF model to be used. Available models "maxwellian", "bimaxwellian" (default)
    :kword normorder:   Norm used for model-data distance measure (default: 1)
    :kword norm:        Constant norm (default 2, see below)
    :kword dummy:       If not None, generate dummy data for e.g. integration.

    :returns:           scalar, non-Maxwellianity parameter for given model and norm

    The definition is given by Graham et al. (2021):
    (1/2n) * integral(|f_i - g_M|) d^3v    
    
    valued between 0 (bi-Maxwellian) and 1 (complete deviation from a bi-Maxwellian).
    
    NOTE that this is different to the definition by Greco et al. (2012):
    (1/n) * sqrt(integral((f_i - g_M)^2) d^3v)    
    
    examples comparing these two definitions can be found in Settino et al. (2021).



    
    '''
    if dummy is not None:
        warnings.warn("Generating dummy value for non-Maxwellianity!")
        return random.random()
    # try getting the distribution from the given cell
    size = f.get_velocity_mesh_size(pop)
    distribution = np.zeros(4 * 4 * 4 * int(size[0]) * int(size[1]) * int(size[2]))
    try:
        D = f.read_velocity_cells(cell, pop)
        D_keys = list(D.keys())
        D_vals = list(D.values())
        distribution[D_keys] = D_vals
        threshold = 1e-21
        distribution[distribution<threshold] = 0
    except:
        warnings.warn("Could not get VDF from bulk file!")
        return -1
     # calculate non-maxwellianity

    
    extent = f.get_velocity_mesh_extent(pop)
    print(extent)
    dvx = (extent[3] - extent[0]) / (4 * size[0])
    dvy = (extent[4] - extent[1]) / (4 * size[1])
    dvz = (extent[5] - extent[2]) / (4 * size[2])
    dV = dvx*dvy*dvz  
    # generate a velocity space 
    vids = np.arange(4 * 4 * 4 * int(size[0]) * int(size[1]) * int(size[2]))
    v = f.get_velocity_cell_coordinates(vids, pop)
    #print(v)
    #sys.exit()
    vs = f.get_velocity_cell_coordinates(vids[D_keys],pop)

    if B is None:
        try:
            print("No B vector given, trying a restart reducer from " + f.file_name)
            B = f.read_variable("B",cellids=cell)
        except:
            print("No restart B found, trying for vg_b_vol from " + f.file_name)
            try:
                B = f.read_variable("vg_b_vol",cellids=cell)
            except:
                print("No vg_b_vol available either")
        if bulk is not None and B is None:
            try:
                print("Trying to load vg_b_vol from given bulk file "+bulk)
                B = pt.vlsvfile.reader(bulk.read_variable("vg_b_vol",cellids=cell))
            except:
                print("Could not load vg_b_vol from bulk file "+bulk)
        if B is None:
            warnings.warn("No B found, cannot proceed at cellid "+str(cell))
            return -1
        
    if bulk is None:
        n = np.sum(np.array(D_vals))*dV
        v0 = np.average(vs, axis=0,weights=np.array(D_vals)*dV)
        # Generate a parallel-perpedicular coordinate basis and corresponding velocity array
        bhat = B/np.linalg.norm(B)
        vperp2hat = np.cross(bhat,v0)/np.linalg.norm(np.cross(bhat,v0))
        vperp1hat = np.cross(vperp2hat, bhat)/np.linalg.norm(np.cross(vperp2hat,bhat))
        
        R = np.array([bhat, vperp1hat, vperp2hat])
        vb = np.matmul(R,vs.T).T
        vb_mean = np.matmul(R,v0.T).T
        v0_perp = vb_mean[1]
        cov = np.zeros((3,3))
        for i,vv in enumerate(vb):
            ov = np.outer(vv-vb_mean,vv-vb_mean)
            cov = cov + ov*D_vals[i]*dV
        P = cov*m
        T = np.trace(P)/(3*n)
        T_para = P[0,0]/n
        T_perp = (P[1,1]+P[2,2])/(2*n)
        print("Pb", P)
        Pperp = 0.5*(P[1,1]+P[2,2])
        Ag = (P[0,1]**2+P[0,2]**2+P[1,2]**2)/(Pperp**2+2*Pperp*P[0,0])
        print("Ag", Ag)
        
    else:
        # get velocity moments from bulk file
        bulkfile_for_moments_reader=pt.vlsvfile.VlsvReader(bulk)
        print(bulkfile_for_moments_reader.read_parameter("time"))
        n = bulkfile_for_moments_reader.read_variable(n_name,cellids=cell)
        v0 = bulkfile_for_moments_reader.read_variable(v_name,cellids=cell)
        v0_para = bulkfile_for_moments_reader.read_variable(v_para_name,cellids=cell)
        v0_perp = bulkfile_for_moments_reader.read_variable(v_perp_name,cellids=cell)
        
        T = bulkfile_for_moments_reader.read_variable(T_name,cellids=cell)
        T_para = bulkfile_for_moments_reader.read_variable(T_para_name,cellids=cell)
        T_perp = bulkfile_for_moments_reader.read_variable(T_perp_name,cellids=cell)
        # Generate a parallel-perpedicular coordinate basis and corresponding velocity array
        bhat = B/np.linalg.norm(B)
        vperp2hat = np.cross(bhat,v0)/np.linalg.norm(np.cross(bhat,v0))
        vperp1hat = np.cross(vperp2hat, bhat)/np.linalg.norm(np.cross(vperp2hat,bhat))
        
        R = np.array([bhat, vperp1hat, vperp2hat])
        vb_mean = np.matmul(R,v0.T).T
        vb = np.matmul(R,v.T).T
    


    vT_para = np.sqrt(2 * k * T_para / m)
    print("v0_para", v0_para, "v0_perp", v0_perp)
    print("T", T, "T_para", T_para, "T_perp", T_perp, "vT_para", vT_para, "vT_perp", vT_para*T_perp/T_para)

    print("v0",v0,"vb_mean",vb_mean, "v0_para",v0_para, "v0_perp",v0_perp)

    #print("v0", v0, "v_mean", v_mean, "v_mean_manual", v_mean_manual,"v_mean_manual2", v_mean_manual2, "in bcoords", np.matmul(R, v_mean.T).T)


      
    #print(cov, np.sqrt(np.trace(cov)))
    #cov = cov*m
    #print(cov)
    #print(np.trace(cov)/rho_manual2)
    # generate a maxwellian distribution
    if model == "maxwellian":
        v_sq = np.sum((v - v0) * (v - v0), axis=-1)
        model_f = n * (0.5 * m / (np.pi * k * T)) ** 1.5 * np.exp(-0.5 * m * v_sq / (k * T))
    elif model == "bimaxwellian": # Graham et al 2021
        model_f = n * T_para / (np.pi**1.5 * vT_para**3 * T_perp) * np.exp(
        -(vb[:,0] - v0_para)**2/vT_para**2
        -((vb[:,1] - v0_perp)**2 + vb[:,2]**2)/(vT_para**2 * (T_perp/T_para))
        )
    else:
        raise NameError("Unknown VDF model '"+model+"', aborting")

   

    #epsilon = np.sum(np.abs(distribution - model_f))
    epsilon = np.linalg.norm(distribution - model_f, ord=normorder)
    epsilon *= dV / (norm * n)
    print("epsilons", epsilon)
    return epsilon


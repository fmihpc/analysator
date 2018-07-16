import numpy as np
import math
import pytools as pt
from multiprocessing import Pool

def expr_Tanisotropy(pass_maps):
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]
    # here it's assumed to contain Tpara and Tperp, and the function
    # returns the Tpara/Tperp

    # Verify that time averaging wasn't used
    if type(pass_maps[0]) is list:
        print("expr_Tanisotropy expected a single timestep, but got multiple. Exiting.")
        quit()

    kb = 1.38065e-23 # Boltzmann's constant

    # If 1D maps
    if np.ndim(pass_maps[0])==1:
        rho = pass_maps[0][:]
        B = pass_maps[1][:]
        PDiag = pass_maps[2][:]
        POff = pass_maps[3][:]
        
        # Normalising B
        B0 = B[:,0]
        B1 = B[:,1]
        B2 = B[:,2]
        normB = np.sqrt(B0**2 + B1**2 + B2**2)
        B0n = B0/normB
        B1n = B1/normB
        B2n = B2/normB
        
        # Building the temperature tensor from variables
        TTensor00 = PDiag[:,0]/(rho+1.)/kb
        TTensor01 = POff[:,2]/(rho+1.)/kb
        TTensor02 = POff[:,1]/(rho+1.)/kb
        TTensor10 = POff[:,2]/(rho+1.)/kb
        TTensor11 = PDiag[:,1]/(rho+1.)/kb
        TTensor12 = POff[:,0]/(rho+1.)/kb
        TTensor20 = POff[:,1]/(rho+1.)/kb
        TTensor21 = POff[:,0]/(rho+1.)/kb
        TTensor22 = PDiag[:,2]/(rho+1.)/kb
        
        # Calculating the parallel temperature as dot(Bn,TTensor*Bn)
        Tpara = (B0n*(TTensor00*B0n+TTensor01*B1n+TTensor02*B2n)
               + B1n*(TTensor10*B0n+TTensor11*B1n+TTensor12*B2n)
               + B2n*(TTensor20*B0n+TTensor21*B1n+TTensor22*B2n))
        
        # Calculating the perpendicular temperature as 0.5*(trace(TTensor)-T_para)
        Tperp = 0.5*(TTensor00+TTensor11+TTensor22 - Tpara)
        
        # Calculate temperature anisotropy as the para/perp ratio
        Tanisotropy = Tperp/(Tpara+1.)

    elif np.ndim(pass_maps[0])==2:    
        rho = pass_maps[0][:,:]
        B = pass_maps[1][:,:]
        PDiag = pass_maps[2][:,:]
        POff = pass_maps[3][:,:]
        
        # Normalising B
        B0 = B[:,:,0]
        B1 = B[:,:,1]
        B2 = B[:,:,2]
        normB = np.sqrt(B0**2 + B1**2 + B2**2)
        B0n = B0/normB
        B1n = B1/normB
        B2n = B2/normB
        
        # Building the temperature tensor from variables
        TTensor00 = PDiag[:,:,0]/(rho+1.)/kb
        TTensor01 = POff[:,:,2]/(rho+1.)/kb
        TTensor02 = POff[:,:,1]/(rho+1.)/kb
        TTensor10 = POff[:,:,2]/(rho+1.)/kb
        TTensor11 = PDiag[:,:,1]/(rho+1.)/kb
        TTensor12 = POff[:,:,0]/(rho+1.)/kb
        TTensor20 = POff[:,:,1]/(rho+1.)/kb
        TTensor21 = POff[:,:,0]/(rho+1.)/kb
        TTensor22 = PDiag[:,:,2]/(rho+1.)/kb
        
        # Calculating the parallel temperature as dot(Bn,TTensor*Bn)
        Tpara = (B0n*(TTensor00*B0n+TTensor01*B1n+TTensor02*B2n)
               + B1n*(TTensor10*B0n+TTensor11*B1n+TTensor12*B2n)
               + B2n*(TTensor20*B0n+TTensor21*B1n+TTensor22*B2n))
        
        # Calculating the perpendicular temperature as 0.5*(trace(TTensor)-T_para)
        Tperp = 0.5*(TTensor00+TTensor11+TTensor22 - Tpara)
        
        # Calculate temperature anisotropy as the perp/para ratio
        Tanisotropy = Tperp/(Tpara+1.)

    else:
        print("ERROR: passed maps have dimension >2")
        return


    return Tanisotropy




def extract_Taniso_step(step=None):

    filename = filedir_global+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
    print(filename+" is being processed")
    f=pt.vlsvfile.VlsvReader(filename)

    pass_maps=[]
        
    for mapval in pass_vars:
        pass_map = f.read_variable(mapval,cellids=cellid_global)
        if np.ndim(pass_map)==1:
            pass_map = pass_map.reshape([1,1])
        else:
            pass_map = pass_map.reshape([1,1,len(pass_map[0])])
        pass_maps.append(np.ma.asarray(pass_map))

    Tani = expr_Tanisotropy(pass_maps)
    Taniso = Tani.data[0][0]

    return (Taniso)



def extract_Tanisotropy(filedir=None,
                        start=None,
                        stop=None,
                        cellid=None,
                        outputdir=None,
                        numproc=8
                        ):

    global filedir_global, cellid_global
    global pass_vars

    filedir_global = filedir
    cellid_global = cellid

    pass_vars=['rho','B','PTensorDiagonal','PTensorOffDiagonal']




    # Parallel extraction of the temperature anisotropy
    if __name__ == 'temperature_anisotropy':
        pool = Pool(numproc)
        Taniso = pool.map(extract_Taniso_step, range(start,stop+1))
    else:
        print("didn't enter the loop")

    
    print('Taniso = '+str(Taniso))

    outputname = outputdir+'Tanisotropy_'+str(cellid[0])+'_'+str(start)+'_'+str(stop)

    print('Saving array in '+outputname)

    np.save(outputname,Taniso)

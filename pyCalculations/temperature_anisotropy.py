import numpy as np
import math

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

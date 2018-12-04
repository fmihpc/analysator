# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 


import pytools as pt
import sys, os, socket
import numpy as np

# Custom expression function                                                               
def exprMA_cust(exprmaps, requestvariables=False):
    if requestvariables==True:
        return ['va']
    # exprmaps is a dictionary of numpy arrays
    # Each array has 2 dimensions [xsize, ysize]
    # or 3 dimensions [xsize, ysize, components]
    # here the function returns the M_A with a preset bulk velocity

    # Verify that time averaging wasn't used
    if type(pass_maps) is list:
        print("exprMA_cust expected a single timestep, but got multiple. Exiting.")
        quit()

    custombulkspeed=1500000. # m/s
    va = exprmaps['va'][:,:]
    MA = custombulkspeed/va
    return MA


# Second example of a more involved custom expression function
def expr_cav_cust(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['rho', 'B', 'beta']
    # pass_maps is a dictionary of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of dictionaries
    # Each dictionary contains an entry called 'dstep'
    # Which contains the relative time step position, i.e.
    # a value of 0 indicates no time offset.

    # This custom expression returns a map with values of
    # either 0 (solar wind), 0.5 (caviton), or 1.0 (SHFA), calculated against
    # time-averaged background values. This version doesn't do intelligent checks for the
    # format of the incoming data.
    if type(pass_maps) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("expr_cav_cust expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()

    # Multiple time steps were found
    ntimes = len(pass_maps)
    dsteps = [x['dstep'] for x in pass_maps]
    curri = dsteps.index(0)
    thesemaps = pass_maps[curri]
    
    thisrho = np.ma.masked_less_equal(thesemaps['rho'][:,:], 0)
    thisB = np.ma.masked_less_equal(thesemaps['B'],0)
    thisbeta = np.ma.masked_less_equal(thesemaps['beta'],0)
    thisBmag = np.linalg.norm(thisB, axis=-1)
        
    avgrho = np.zeros(np.array(thisrho.shape))
    avgBmag = np.zeros(np.array(thisrho.shape))
    # avgbeta = np.zeros(np.array(thisrho.shape))
    
    for i in range(ntimes):
        if i==curri: # Exclude current frame from background value averaging
            continue
        nowmaps = pass_maps[i]
        print(nowmaps['dstep'])
        avgrho = np.add(avgrho, nowmaps['rho'])
        avgBcalc = np.linalg.norm(nowmaps['B'], axis=-1)
        avgBmag = np.add(avgBmag, avgBcalc)
        # avgbeta = np.add(avgbeta, nowmaps[2])

    avgrho = np.divide(np.ma.masked_less_equal(avgrho,0), np.array([ntimes-1]))
    avgBmag = np.divide(np.ma.masked_less_equal(avgBmag,0), np.array([ntimes-1]))
    #avgbeta = np.divide(np.ma.masked_less_equal(avgbeta,0), np.array([ntimes-1]))

    rhoratio = np.ma.divide(thisrho, avgrho)
    Bmagratio = np.ma.divide(thisBmag, avgBmag)
    #betaratio = np.divide(thisbeta, avgbeta)

    # Calculations using masked arrays proved problematic so a less-than elegant method is used here.
    empty = np.zeros(np.array(thisrho.shape))
    half = empty + 0.5
    one = empty + 1.0
    caviton = np.add(empty, one, where=(rhoratio<0.8))
    caviton = np.add(caviton, one, where=(Bmagratio<0.8))
    shfa = np.add(caviton, one, where=(thisbeta>10))

    combo = np.add(empty, half, where=(caviton>1.5))
    combo2 = np.add(empty, half, where=(shfa>2.5))
    combo3 = combo+combo2

    # Mask out anything that is inside the bow shock
    bowshock = 2.e6
    combo3 = np.ma.masked_where(thisrho>bowshock, combo3)

    return combo3

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

    # Calculate temperature anisotropy as the para/perp ratio
    Tanisotropy = Tpara/(Tperp+1.)

    return Tanisotropy


# Helper function for drawing fsaved cells
def overplotfsaved(ax, XmeshXY,YmeshXY, pass_maps):
    contour_mirror = ax.contour(XmeshXY,YmeshXY,pass_maps[0],[0.5],
                               linewidths=1, colors='black')



# Helper function for drawing on existing panel
def cavitoncontours(ax, XmeshXY,YmeshXY, extmaps, requestvariables=False):
    if requestvariables==True:
        return ['rho', 'B', 'beta']

    # Check if pass_maps has multiple time steps or just one
    if type(pass_maps) is list:
        dsteps = [x['dstep'] for x in pass_maps]
        curri = dsteps.index(0)
        rho = extmaps[curri]['rho']
        beta = extmaps[curri]['beta']
        # take magnitude of B
        B = np.linalg.norm(extmaps[curri]['B'],axis=-1)
    else:
        rho = extmaps['rho']
        beta = extmaps['beta']
        # take magnitude of B
        B = np.linalg.norm(extmaps['B'],axis=-1)

    # Colours to use
    color_cavitons = '#924900'
    color_SHFAs    = '#B66DFF'
    color_BS       = '#FFFF6D'

    # thresholds
    level_bow_shock = 2.e+6
    level_n_caviton = 0.8e+6
    level_B_caviton = 4.e-9
    level_beta_SHFA = 150
    level_beta_SHFA_SW = 10.

    # mask cavitons
    cavitons = np.ma.masked_greater_equal(B,level_B_caviton)
    cavitons.mask[rho > level_n_caviton] = True
    cavitons.fill_value = 0.
    cavitons[cavitons.mask == False] = 1.

    # mask SHFAs
    SHFAs = np.ma.masked_greater_equal(B,level_B_caviton)
    SHFAs.mask[rho > level_n_caviton] = True
    SHFAs.mask[beta < level_beta_SHFA_SW] = True
    SHFAs.fill_value = 0.
    SHFAs[SHFAs.mask == False] = 1.
 
    # draw contours
    contour_shock = ax.contour(XmeshXY,YmeshXY,rho,[level_bow_shock], 
                               linewidths=1.2, colors=color_BS,label='Bow shock')
    contour_cavitons = ax.contour(XmeshXY,YmeshXY,cavitons.filled(),[0.5], linewidths=1.5, colors=color_cavitons)  
    contour_SHFAs = ax.contour(XmeshXY,YmeshXY,SHFAs.filled(),[0.5], linewidths=1.5, colors=color_SHFAs)           



# Helper function for drawing on existing panel
def insetVDF(ax, XmeshXY,YmeshXY, pass_maps):
    # pass_maps is a list of numpy arrays, not used here.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # 'upper right'  : 1,
    # 'upper left'   : 2,
    # 'lower left'   : 3,
    # 'lower right'  : 4,
    # 'right'        : 5,
    # 'center left'  : 6,
    # 'center right' : 7,
    # 'lower center' : 8,
    # 'upper center' : 9,
    # 'center'       : 10

    # Generate inset axes for VDF #1
    VDFcoord = [10.8,0,-5]
    VDFax = inset_axes(ax, width="25%", height="25%", loc=2, bbox_transform=ax.transAxes)
    pt.plot.plot_vdf(filename=fileLocation+bulkname,coordre=VDFcoord,box=[-2.5e6,2.5e6,-2.5e6,2.5e6], fmin=1e-14, fmax=2e-8,slicethick=0,axes=VDFax, unit=6,nocb=1,title='',noxlabels=1,noylabels=1)
    # Add dot at VDF location
    ax.scatter(VDFcoord[0], VDFcoord[2], color='black',marker='o',s=20)
    ax.scatter(VDFcoord[0], VDFcoord[2], color='gray',marker='o',s=2)

    # Generate inset axes for VDF #2
    VDFcoord = [10.8,0,-15]
    VDFax = inset_axes(ax, width="25%", height="25%", loc=6, bbox_transform=ax.transAxes)
    pt.plot.plot_vdf(filename=fileLocation+bulkname,coordre=VDFcoord,box=[-2.5e6,2.5e6,-2.5e6,2.5e6], fmin=1e-14, fmax=2e-8,slicethick=0,axes=VDFax, unit=6,nocb=1,title='',noxlabels=1,noylabels=1)
    # Add dot at VDF location
    ax.scatter(VDFcoord[0], VDFcoord[2], color='black',marker='o',s=20)
    ax.scatter(VDFcoord[0], VDFcoord[2], color='gray',marker='o',s=2)

    # Generate inset axes for VDF #3
    VDFcoord = [10.8,0,-25]
    VDFax = inset_axes(ax, width="25%", height="25%", loc=3, bbox_transform=ax.transAxes)
    pt.plot.plot_vdf(filename=fileLocation+bulkname,coordre=VDFcoord,box=[-2.5e6,2.5e6,-2.5e6,2.5e6], fmin=1e-14, fmax=2e-8,slicethick=0,axes=VDFax, unit=6,nocb=1,title='',noxlabels=1,noylabels=1)
    # Add dot at VDF location
    ax.scatter(VDFcoord[0], VDFcoord[2], color='black',marker='o',s=20)
    ax.scatter(VDFcoord[0], VDFcoord[2], color='gray',marker='o',s=2)
   

fileLocation="/proj/vlasov/2D/BCQ/bulk/"
fluxLocation="/proj/vlasov/2D/BCQ/flux/"
outputLocation=outputdir=os.path.expandvars('$HOME/Plots/')

# Frame extent for this job given as command-line arguments
if len(sys.argv)==3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    # Plot(s) to be made
#    pt.plot.plot_colormap_with_vdf(filename=bulkname, var='temperature', boxre=[-15, -2, -6.5, 6.5],vmin=5e5, vmax=1e8,step=i/2, colormap='nipy_spectral', cellidplot=


#    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[-20, -8, -6, 6], colormap='RdBu_r',cellidplot=[3601701,3601751,3601801,3751851,3601851,345185

    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[-15, -2, -6.5, 6.5], colormap='hot_desaturated',cellidplot=[3601701,3601751,3601801,3751851,36

#    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[-20, -8, -6, 6], colormap='RdBu_r',cellidplot=[3601601,3601651,3601701,3601751,3601801],step=

#    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[-20, -8, -6, 6], colormap='RdBu_r',cellidplot=[3601601,3601651,3601701,3601751,3601801],step=

#    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[-20, -8, -6, 6], colormap='RdBu_r',cellidplot=[3601601,3601651,3601701,3601751,3601801],step=

#    pt.plot.plot_vdf(filename=fileLocation+bulkname, step=j, cellids=3601751,box=[-3e6,3e6,-3e6,3e6],bpara=1,colormap='nipy_spectral', fmin=1e-15, fmax=1e-11, cbulk=1
#    pt.plot.plot_vdf(filename=fileLocation+bulkname, step=j, cellids=3601751,box=[-3e6,3e6,-3e6,3e6],bpara1=1,colormap='nipy_spectral', fmin=1e-15, fmax=1e-11, cbulk=
#    pt.plot.plot_vdf(filename=fileLocation+bulkname, step=j, cellids=3601751,box=[-3e6,3e6,-3e6,3e6],bperp=1,colormap='nipy_spectral', fmin=1e-15, fmax=1e-11, cbulk=1




# BFB stuff
#    pt.plot.plot_colormap_with_vdf(filename=fileLocation+bulkname,boxre=[3,11,0,8], colormap='seismic',cellidplot=[2702151,2852101,2552151],vdfplotlocation='ttl',vdfp


    # Plot a custom nightside reconnection slippage calculation
    pt.plot.plot_helpers.slippageVA=3000000
    pt.plot.plot_colormap(filename=fileLocation+bulkname, run="BCQ",colormap='seismic',step=j,outputdir=outputLocation+'slippage/', wmark=1, fluxdir=fluxLocation, fluxlines=4, boxre=[-40,-10,-15,15], vmin=1e-2, vmax=1e0, pass_vars=['E','B','V'], expression=expr_Slippage)


    # Plot beam number density with inset VDF
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="rhoBeam",run="BCQ",colormap='bwr',step=j,outputdir=outputLocation+'rhoBeamVDF/',wmark=1, boxre=[-10,20,-30,0], external=insetVDF)



   # Useful flags
   # (Find more by entering pt.plot.plot_colormap? in python
   #
   #   nooverwrite = 1:    Set to only perform actions if the target output file does not yet exist                    
   #   boxm=[x0,x1,y0,y1]  Zoom to this box (size given in metres)
   #   boxre=[x0,x1,y0,y1] Zoom to this box (size given in Earth radii)
   #   usesci=0:           Disable scientific notation for colorbar ticks
   #   symlog=0            Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
   #                       Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
   #                       A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
   #   wmark=1             If set to non-zero, will plot a Vlasiator watermark in the top left corner.
   #   draw=1              Set to nonzero to draw image on-screen instead of saving to file (requires x-windowing)



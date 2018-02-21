import pytools as pt
import sys, os, socket
import numpy as np

                     
# Custom expression function                                                               
def exprMA_cust(exprmaps):
    # exprmaps is a list of numpy arrays
    # Each array has 2 dimensions [xsize, ysize]
    # or 3 dimensions [xsize, ysize, components]
    # here it's assumed to contain va, and the function
    # returns the M_A with a preset velocity
    custombulkspeed=750000. # m/s
    va = exprmaps[0][:,:]
    MA = custombulkspeed/va
    return MA

# Helper function for drawing on existing panel
def extcontour(ax, XmeshXY,YmeshXY, extmaps):
    # extmaps consists of [rho, B, beta]

    # take magnitude of B
    extmaps[1]=np.sum(np.asarray(extmaps[1])**2,axis=-1)**(0.5)

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
    cavitons = np.ma.masked_greater_equal(extmaps[1],level_B_caviton)
    cavitons.mask[extmaps[0] > level_n_caviton] = True
    cavitons.fill_value = 0.
    cavitons[cavitons.mask == False] = 1.

    # mask SHFAs
    SHFAs = np.ma.masked_greater_equal(extmaps[1],level_B_caviton)
    SHFAs.mask[extmaps[0] > level_n_caviton] = True
    SHFAs.mask[extmaps[2] < level_beta_SHFA_SW] = True
    SHFAs.fill_value = 0.
    SHFAs[SHFAs.mask == False] = 1.
 
    # draw contours
    contour_shock = ax.contour(XmeshXY,YmeshXY,extmaps[0],[level_bow_shock], 
                               linewidths=1.2, colors=color_BS,label='Bow shock')
    contour_cavitons = ax.contour(XmeshXY,YmeshXY,cavitons.filled(),[0.5], linewidths=1.5, colors=color_cavitons)  
    contour_SHFAs = ax.contour(XmeshXY,YmeshXY,SHFAs.filled(),[0.5], linewidths=1.5, colors=color_SHFAs)           



fileLocation="/proj/vlasov/2D/BCQ/bulk/"
outputLocation=outputdir=os.path.expandvars('$HOME/Plots/')

# Frame extent for this job given as command-line arguments
timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)

for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    # Five different plots

    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="rho",run="BCQ",colormap='viridis_r',step=j,outputdir=outputLocation+'rho/',wmark=1, external=extcontour, extvals=['rho','B','beta'])

    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="B",run="BCQ",colormap='inferno_r',step=j,outputdir=outputLocation+'Bmag/',wmark=1, external=extcontour, extvals=['rho','B','beta'])

    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="rhoBeam",run="BCQ",colormap='bwr',step=j,outputdir=outputLocation+'rhoBeam/',wmark=1, external=extcontour, extvals=['rho','B','beta'])

    pt.plot.plot_colormap(filename=fileLocation+bulkname,run="BCQ",colormap='RdBu',step=j,outputdir=outputLocation+'By/',lin=1, var='B',op='y')

    pt.plot.plot_colormap(filename=fileLocation+bulkname,run="BCQ",colormap='plasma',step=j,outputdir=outputLocation+'cMA/',lin=1,expression=exprMA_cust, exprvals=['va'], vmin=1, vmax=20)

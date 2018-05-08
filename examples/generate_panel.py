
import pytools as pt
import sys, os, socket
import numpy as np
                     
# Custom expression function                                                               
def exprMA_cust(pass_maps):
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]
    # here it's assumed to contain va, and the function
    # returns the M_A with a preset velocity

    # Verify that time averaging wasn't used
    if type(pass_maps[0]) is list:
        print("exprMA_cust expected a single timestep, but got multiple. Exiting.")
        quit()

    custombulkspeed=750000. # m/s
    va = pass_maps[0][:,:]
    MA = custombulkspeed/va
    return MA

# Second example of a more involved custom expression function
def expr_cav_cust(pass_maps):
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of lists, where the top level
    # is 2N+1 timesteps with the middle one the requested time step

    # This custom expression returns a map with values of
    # either 0 (solar wind), 0.5 (caviton), or 1.0 (SHFA), calculated against
    # time-averaged background values. This version doesn't do intelligent checks for the
    # format of the incoming data.
    if type(pass_maps[0]) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("expr_cav_cust expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()

    # Multiple time steps were found
    ntimes = len(pass_maps)
    curri = (ntimes-1)/2
    thesemaps = pass_maps[curri]
    
    thisrho = np.ma.masked_less_equal(thesemaps[0][:,:], 0)
    thisB = np.ma.masked_less_equal(thesemaps[1],0)
    thisbeta = np.ma.masked_less_equal(thesemaps[2],0)
    thisBmag = np.linalg.norm(thisB, axis=-1)
        
    avgrho = np.zeros(np.array(thisrho.shape))
    avgBmag = np.zeros(np.array(thisrho.shape))
    # avgbeta = np.zeros(np.array(thisrho.shape))
    
    for i in range(ntimes):
        if i==curri: # Exclude current frame from background value averaging
            continue
        nowmaps = pass_maps[i]
        avgrho = np.add(avgrho, nowmaps[0])
        avgBcalc = np.linalg.norm(nowmaps[1], axis=-1)
        avgBmag = np.add(avgBmag, avgBcalc)
        # avgbeta = np.add(avgbeta, nowmaps[2])

    avgrho = np.divide(np.ma.masked_less_equal(avgrho,0), np.array([ntimes-1]))
    avgBmag = np.divide(np.ma.masked_less_equal(avgBmag,0), np.array([ntimes-1]))
    #avgbeta = np.divide(np.ma.masked_less_equal(avgbeta,0), np.array([ntimes-1]))

    rhoratio = np.divide(thisrho, avgrho)
    Bmagratio = np.divide(thisBmag, avgBmag)
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


# Helper function for drawing on existing panel
def cavitoncontours(ax, XmeshXY,YmeshXY, pass_maps):
    # pass_maps consists of [rho, B, beta]

    # Check if pass_maps has multiple time steps or just one
    if type(pass_maps[0]) is list:
        # Multiple time steps
        ntimes = len(pass_maps)
        # pick values from "center" timestep
        curri = (ntimes-1)/2
        rho = pass_maps[curri][0]
        B = pass_maps[curri][1]
        beta = pass_maps[curri][2]
    else:
        rho = pass_maps[0]
        B = pass_maps[1]
        beta = pass_maps[2]

    # take magnitude of B
    #B=np.sum(B**2,axis=-1)**(0.5)
    B=np.linalg.norm(B,axis=-1)
        
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
                               linewidths=0.5, colors=color_BS,label='Bow shock')
    contour_cavitons = ax.contour(XmeshXY,YmeshXY,cavitons.filled(),[0.5], linewidths=0.5, colors=color_cavitons)  
    contour_SHFAs = ax.contour(XmeshXY,YmeshXY,SHFAs.filled(),[0.5], linewidths=0.5, colors=color_SHFAs)           


# -------------------------------------------------------------------------------

fileLocation="/proj/vlasov/2D/BCH/bulk/"
fluxLocation="/proj/vlasov/2D/BCH/flux/"
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
    pt.plot.plot_colormap_with_vdf(filename=bulkname, var='temperature', boxre=[-15, -2, -6.5, 6.5],vmin=5e5, vmax=1e8,step=i/2, colormap='nipy_spectral', cellidplot=[3601851,3601801,3601751,3601701,3751901,3451901], outputdir=outputLocation+'mosaic/',fluxdir=fluxLocation,fluxlines=2)



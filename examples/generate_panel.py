import pytools as pt
import sys, os, socket
import numpy as np

                     
# Custom expression function                                                               
def exprMA_cust(exprmaps):
    # exprmaps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]
    # here it's assumed to contain va, and the function
    # returns the M_A with a preset velocity
    custombulkspeed=750000. # m/s
    va = exprmaps[0][:,:]
    MA = custombulkspeed/va
    return MA

# Second example of a more involved custom expression function
def expr_cav_cust(exprmaps):
    # exprmaps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of lists, where the top level
    # is 2N+1 timesteps with the middle one the requested time step

    # This custom expression returns a map with values of
    # either 0 (solar wind), 0.5 (caviton), or 1.0 (SHFA), calculated against
    # time-averaged background values. This version doesn't do intelligent checks for the
    # format of the incoming data.
    ntimes = len(exprmaps)
    curri = (ntimes-1)/2
    thesemaps = exprmaps[curri]
    
    thisrho = np.ma.masked_less_equal(thesemaps[0][:,:], 0)
    thisB = np.ma.masked_less_equal(thesemaps[1],0)
    thisbeta = np.ma.masked_less_equal(thesemaps[2],0)
    thisBmag = np.linalg.norm(thisB, axis=-1)

    avgrho = np.zeros(np.array(thisrho.shape))
    avgBmag = np.zeros(np.array(thisrho.shape))
    #avgbeta = np.zeros(np.array(thisrho.shape))
    
    for i in range(ntimes):
        if i==curri: # Exclude current frame from background value averaging
            continue
        nowmaps = exprmaps[i]
        avgrho = np.add(avgrho, nowmaps[0])
        avgBcalc = np.linalg.norm(nowmaps[1], axis=-1)
        avgBmag = np.add(avgBmag, avgBcalc)
        #avgbeta = np.add(avgbeta, nowmaps[2])
        
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
fluxLocation="/proj/vlasov/2D/BCQ/flux/"
outputLocation=outputdir=os.path.expandvars('$HOME/Plots/')

# Frame extent for this job given as command-line arguments
timetot = range(int(sys.argv[1]), int(sys.argv[2]), 1)

for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    # Five different plots

    # Plot rho with caviton and SHFA contours
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="rho",run="BCQ",colormap='viridis_r',step=j,outputdir=outputLocation+'rho/',wmark=1, external=extcontour, extvals=['rho','B','beta'])

    # Plot B-magnitude with caviton and SHFA contours
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="B",run="BCQ",colormap='inferno_r',step=j,outputdir=outputLocation+'Bmag/',wmark=1, external=extcontour, extvals=['rho','B','beta'])

    # Plot beam number density with magnetic field lines on top
    pt.plot.plot_colormap(filename=fileLocation+bulkname,var="rhoBeam",run="BCQ",colormap='bwr',step=j,outputdir=outputLocation+'rhoBeam/',wmark=1, fluxdir=fluxLocation)

    # Plot a linear plot of magnetic field y-component
    pt.plot.plot_colormap(filename=fileLocation+bulkname,run="BCQ",colormap='RdBu',step=j,outputdir=outputLocation+'By/',lin=1, var='B',op='y')

    # Plot a custom expression of MA using a pre-set bulk flow speed
    pt.plot.plot_colormap(filename=fileLocation+bulkname,run="BCQ",colormap='plasma',step=j,outputdir=outputLocation+'cMA/',lin=1,expression=exprMA_cust, exprvals=['va'], vmin=1, vmax=20)

    # Plot a custom time-averaged caviton and SHFA colourmap with non-timeaveraged contours on top for comparison. Leaves out the
    # Colour bar.and associated title.
   pt.plot.plot_colormap(filename=fileLocation+bulkfile, run="BCQ",colormap='OrRd',step=j,outputdir=outputLocation+'cSHFA/', lin=1,vmin=0,vmax=1, expression=expr_cav_cust, exprvals=['rho','B','beta'], expr_timeavg=20, external=cavitoncontours, extvals=['rho','B','beta'], nocb=1, cbtitle='')

   # Other useful flags
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


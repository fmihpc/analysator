import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os, sys
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.ticker import LogLocator
import matplotlib.ticker as mtick
import colormaps as cmaps
from matplotlib.cbook import get_sample_data
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import ids3d
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import mpl_toolkits.mplot3d.art3d as art3d


# Create the 3d axes and the coordinate axes for the 3d plot
def axes3d(fig, reflevel,  xr, yr, zr, xsize, ysize, zsize):
    # get 3d axes
#    ax = fig.gca(projection='3d')
    ax = fig.add_axes([.1,.1,.8,.8],projection='3d')


    # create 3d coordinate axes lines which cuts the point at xr, yr and zr
    line=art3d.Line3D(*zip((0, yr, zr), (xsize*2**reflevel, yr, zr)),
                      color='black', linewidth=0.8, alpha=0.25, zorder=100)
    ax.add_line(line)
    line=art3d.Line3D(*zip((xr, 0, zr), (xr, ysize*2**reflevel, zr)),
                      color='black', linewidth=0.8, alpha=0.25, zorder=100)
    ax.add_line(line)
    line=art3d.Line3D(*zip((xr, yr, 0), (xr, yr, ysize*2**reflevel)),
                      color='black', linewidth=0.8, alpha=0.25, zorder=100)
    ax.add_line(line)



    # spacing between ticks
    s = int((xsize*2**reflevel)/20)
    # widths of the ticks
    l = int((xsize*2**reflevel)/100)


    ###################
    # plot ticks STARTS    

    # plot the first ticks which position > xr or in the others yr or zr
    i = 1
    while xr+i*s < xsize*2**reflevel:
      line=art3d.Line3D(*zip((xr+i*s, yr-l, zr), (xr+i*s, yr+l, zr)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr+i*s, yr, zr-l), (xr+i*s, yr, zr+l)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1

    # plot the last ticks which position < xr or in the others yr or zr
    i = 1
    while xr-i*s > 0:
      line=art3d.Line3D(*zip((xr-i*s, yr-l, zr), (xr-i*s, yr+l, zr)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr-i*s, yr, zr-l), (xr-i*s, yr, zr+l)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1


    i = 1
    while yr+i*s < ysize*2**reflevel:
      line=art3d.Line3D(*zip((xr-l, yr+i*s, zr), (xr+l, yr+i*s, zr)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr, yr+i*s, zr-l), (xr, yr+i*s, zr+l)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1

    i = 1
    while yr-i*s > 0:
      line=art3d.Line3D(*zip((xr-l, yr-i*s, zr), (xr+l, yr-i*s, zr)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr, yr-i*s, zr-l), (xr, yr-i*s, zr+l)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1


    i = 1
    while zr+i*s < zsize*2**reflevel:
      line=art3d.Line3D(*zip((xr, yr-l, zr+i*s), (xr, yr+l, zr+i*s)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr-l, yr, zr+i*s), (xr+l, yr, zr+i*s)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1

    i = 1
    while zr-i*s > 0:
      line=art3d.Line3D(*zip((xr, yr-l, zr-i*s), (xr, yr+l, zr-i*s)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      line=art3d.Line3D(*zip((xr-l, yr, zr-i*s), (xr+l, yr, zr-i*s)),
                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
      ax.add_line(line)
      i += 1

    # plot ticks ENDS
    #################



    # set the basic 3d coordinate axes off
    ax.set_axis_off()

    # return axes
    return ax




def plot_threeslice(filename=None,
                  vlsvobj=None,
                  filedir=None, step=None, run=None,
                  outputdir=None, outputfile=None,
                  nooverwrite=False,
                  var=None, op=None, operator=None,
                  colormap=None, vmin=None, vmax=None,
                  symmetric=False, lin=None, symlog=None,
                  usesci=True,
                  wmark=False,wmarkb=False,
                  thick=1.0,scale=1.0,
                  expression=None,
                  vscale=1.0,
                  cutpoint=None,cutpointre=None
                  ):

    ''' Plots a 3d plot constructed of three 2d cut throughs.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword run:         run identifier, used for constructing output filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/ or override with PTOUTPUTDIR)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final part will be used as a prefix for the files.
    :kword outputfile:  Singular output file name

    :kword nooverwrite: Set to only perform actions if the target output file does not yet exist                    

    :kword var:         variable to plot, e.g. rho, RhoBackstream, beta, Temperature, MA, Mms, va, vms,
                        E, B, v, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc
    :kword operator:    Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude. 
    :kword op:          duplicate of operator

           
    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, parula, nipy_spectral, RdBu, bwr
    :kword vmin,vmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword symmetric:   Set the absolute value of vmin and vmax to the greater of the two
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2, but this can
                        result in the innermost tick marks overlapping. In this case, using a larger value for 
                        symlog is suggested.
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: True)
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner. If set to a text
                        string, tries to use that as the location, e.g. "NW","NE","SW","SW"
    :kword wmarkb:      As for wmark, but uses an all-black Vlasiator logo.

    :kword scale:       Scale text size (default=1.0)
    :kword thick:       line and axis thickness, default=1.0

    :kword expression:  Optional function which calculates a custom expression to plot. The function
                        receives the same dictionary of numpy arrays as external, as an argument pass_maps,
                        the contents of which are maps of variables. Each is either of size [ysize,xsize]
                        or for multi-dimensional variables (vectors, tensors) it's [ysize,xsize,dim].
                        If the function accepts a second variable, if set to true, it is expected to 
                        return a list of required variables for pass_maps.

    Important note: the dictionaries of arrays passed to external and expression are of shape [ysize,xzize], so
    for some analysis transposing them is necessary. For pre-existing functions to use and to base new functions
    on, see the plot_helpers.py file.

    :kword vscale:      Scale all values with this before plotting. Useful for going from e.g. m^-3 to cm^-3
                        or from tesla to nanotesla. Guesses correct units for colourbar for some known
                        variables.

    :kword cutpoint:    Coordinates of the point through which all three 2D cuts must pass [m]
    :kword cutpointre:  Coordinates of the point through which all three 2D cuts must pass [rE]

    .. code-block:: python

    '''

    # Input file or object
    if filename!=None:
        f=pt.vlsvfile.VlsvReader(filename)
    elif ((filedir!=None) and (step!=None)):
        filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
        f=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        f=vlsvobj
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    if (outputdir is ''):
        outputdir = './'

    if operator==None:
        if op!=None:
            operator=op

    if colormap==None:
        # Default values
        colormap="hot_desaturated"
        if operator=='x' or operator=='y' or operator=='z':
            colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=8*scale # Colour bar ticks and title
    # Small internal colorbar needs increased font size

    # step, used for file name
    if step is not None:
        stepstr = '_'+str(step).rjust(7,'0')
    else:
        if filename:
            stepstr = '_'+filename[-12:-5]
        else:
            stepstr = ''

    # If run name isn't given, just put "plot" in the output file name
    if run is None:
        run='plot'
        if filename:
            # If working within CSC filesystem, make a guess:
            if filename[0:16]=="/proj/vlasov/3D/":
                run = filename[16:19]

    # Verify validity of operator
    operatorstr=''
    operatorfilestr=''
    if operator:
        # .isdigit checks if the operator is an integer (for taking an element from a vector)
        if type(operator) is int:
            operator = str(operator)
        if not operator in 'xyz' and operator is not 'magnitude' and not operator.isdigit():
            print("Unknown operator "+str(operator))
            operator=None
            operatorstr=''
        if operator in 'xyz':
            # For components, always use linear scale, unless symlog is set
            operatorstr='_'+operator
            operatorfilestr='_'+operator
            if symlog is None:
                lin=True
        # index a vector
        if operator.isdigit():
            operator = str(operator)
            operatorstr='_{'+operator+'}'
            operatorfilestr='_'+operator
        # Note: operator magnitude gets operatorstr=''

    # Output file name
    if expression:
        varstr=expression.__name__.replace("/","_")
    else:        
        if not var:
            # If no expression or variable given, defaults to rhom
            var='vg_rhom'
            if f.check_variable("proton/vg_rho"): # multipop v5
                var = 'proton/vg_rho'
            elif f.check_variable("moments"): # restart
                var = 'vg_restart_rhom'
        varstr=var.replace("/","_")

    if not outputfile: # Generate filename
        if not outputdir: # default initial path
            outputdir=pt.plot.defaultoutputdir
        # Sub-directories can still be defined in the "run" variable
        outputfile = outputdir+run+"_threeSlice_"+varstr+operatorfilestr+stepstr+".png"
    else: 
        if outputdir:
            outputfile = outputdir+outputfile

    # Re-check to find actual target sub-directory
    outputprefixind = outputfile.rfind('/')
    if outputprefixind >= 0:            
        outputdir = outputfile[:outputprefixind+1]

    # Ensure output directory exists
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except:
            pass

    if not os.access(outputdir, os.W_OK):
        print(("No write access for directory "+outputdir+"! Exiting."))
        return

    # Check if target file already exists and overwriting is disabled
    if (nooverwrite and os.path.exists(outputfile)):            
        if os.stat(outputfile).st_size > 0: # Also check that file is not empty
            print(("Found existing file "+outputfile+". Skipping."))
            return
        else:
            print(("Found existing file "+outputfile+" of size zero. Re-rendering."))

    Re = 6.371e+6 # Earth radius in m
    # read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")

    if cutpointre is not None:
        cutpoint = np.asarray(cutpointre) * Re


    ############################################
    # Read data and calculate required variables
    ############################################
    if expression==None:
        # Read data from file
        if operator==None:
            operator="pass"
        datamap_info = f.read_variable_info(var, operator=operator)

        cb_title_use = datamap_info.latex
        # if cb_title_use == "": 
        #     cb_title_use = r""+var.replace("_","\_")
        datamap_unit = datamap_info.latexunits

        # If vscale is in use
        if not np.isclose(vscale,1.):
            datamap_unit=r"${\times}$"+pt.plot.cbfmt(vscale,None)
        # Allow specialist units for known vscale and unit combinations
        if datamap_info.units=="s" and np.isclose(vscale,1.e6):
            datamap_unit = r"$\mu$s"
        if datamap_info.units=="s" and np.isclose(vscale,1.e3):
            datamap_unit = "ms"
        if datamap_info.units=="T" and np.isclose(vscale,1.e9):
            datamap_unit = "nT"
        if datamap_info.units=="K" and np.isclose(vscale,1.e-6):
            datamap_unit = "MK"
        if datamap_info.units=="Pa" and np.isclose(vscale,1.e9):
            datamap_unit = "nPa"
        if datamap_info.units=="1/m3" and np.isclose(vscale,1.e-6):
            datamap_unit = r"$\mathrm{cm}^{-3}$"
        if datamap_info.units=="m/s" and np.isclose(vscale,1.e-3):
            datamap_unit = r"$\mathrm{km}\,\mathrm{s}^{-1}$"
        if datamap_info.units=="V/m" and np.isclose(vscale,1.e3):
            datamap_unit = r"$\mathrm{mV}\,\mathrm{m}^{-1}$"
        if datamap_info.units=="eV/cm3" and np.isclose(vscale,1.e-3):
            datamap_unit = r"$\mathrm{keV}\,\mathrm{cm}^{-3}$"

        # Add unit to colorbar title
        if datamap_unit!="":
            cb_title_use = cb_title_use + " ["+datamap_unit+"]"

        # Verify data shape
        datamap = datamap_info.data
        if np.ndim(datamap)==0:
            print("Error, read only single value from vlsv file!",datamap.shape)
            return -1
        if np.ndim(datamap)==2:
            if len(datamap[0,:])!=3:
                print("Error, expected array of 3-element vectors, found array of shape ",datamap.shape)
                return -1
            # 2-dimensional array: take magnitude of three-element vectors
            datamap = np.linalg.norm(datamap, axis=-1)
        if np.ndim(datamap)!=1:
            # Array dimensions not as expected
            print("Error reading variable "+var+"! Found array of shape ",datamap.shape,". Exiting.")
            return -1


        # sort cellid and datamap
        indexids = cellids.argsort()
        cellids = cellids[indexids]
        datamap = datamap[indexids]

        # find the highest refiment level
        reflevel = ids3d.refinement_level(xsize, ysize, zsize, cellids[-1])

        # coordinates of the point where the all three 2d cut throughs cuts
        xr = abs(xmin) + cutpoint[0]
        yr = abs(ymin) + cutpoint[1]
        zr = abs(zmin) + cutpoint[2]


        #############################################################
        # Create the needed data grids for the all three cut throughs
        idlist, indexlist = ids3d.ids3d(cellids, xr, reflevel, xsize, ysize, zsize, xmin=xmin, xmax=xmax)
        datamap1 = datamap[indexlist]
        datamap1 = ids3d.idmesh3d(idlist, datamap1, reflevel, xsize, ysize, zsize, 0, None) # TODO adapt the last argument to also support vectors and tensors

        idlist, indexlist = ids3d.ids3d(cellids, yr, reflevel, xsize, ysize, zsize, ymin=ymin, ymax=ymax)
        datamap2 = datamap[indexlist]
        datamap2 = ids3d.idmesh3d(idlist, datamap2, reflevel, xsize, ysize, zsize, 1, None) # TODO adapt the last argument to also support vectors and tensors

        idlist, indexlist = ids3d.ids3d(cellids, zr, reflevel, xsize, ysize, zsize, zmin=zmin, zmax=zmax)
        datamap3 = datamap[indexlist]
        datamap3 = ids3d.idmesh3d(idlist, datamap3, reflevel, xsize, ysize, zsize, 2, None) # TODO adapt the last argument to also support vectors and tensors
    else:
        # Expression set, use generated or provided colorbar title
        cb_title_use = expression.__name__.replace("_","\_")

    #If automatic range finding is required, find min and max of array
    # Performs range-finding on a masked array to work even if array contains invalid values
    if vmin is not None:
        vminuse=vmin
    else: 
        vminuse=np.ma.amin([np.ma.amin(datamap1),np.ma.amin(datamap2),np.ma.amin(datamap3)])
    if vmax is not None:
        vmaxuse=vmax
    else:
        vmaxuse=np.ma.amax([np.ma.amax(datamap1),np.ma.amax(datamap2),np.ma.amax(datamap3)])

    # If both values are zero, we have an empty array
    if vmaxuse==vminuse==0:
        print("Error, requested array is zero everywhere. Exiting.")
        return 0

    # If vminuse and vmaxuse are extracted from data, different signs, and close to each other, adjust to be symmetric
    # e.g. to plot transverse field components. Always done for symlog.
    if vmin is None and vmax is None:
        if symmetric or np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2) or symlog is not None:
            absval = max(abs(vminuse),abs(vmaxuse))
            vminuse = -absval
            vmaxuse = absval

    # Ensure that lower bound is valid for logarithmic plots
    if (vminuse <= 0) and (lin is None) and (symlog is None):
        # Drop negative and zero values
        vminuse = np.ma.amin([np.ma.amin(np.ma.masked_less_equal(datamap1,0)),np.ma.amin(np.ma.masked_less_equal(datamap2,0)),np.ma.amin(np.ma.masked_less_equal(datamap3,0))])

    # Special case of very small vminuse values
    if ((vmin is None) or (vmax is None)) and (vminuse > 0) and (vminuse < vmaxuse*1.e-5):
        vminuse = vmaxuse*1e-5
        if lin is not None:
            vminuse = 0
    
    # # If symlog scaling is set:
    if symlog is not None:
        if symlog > 0:
            linthresh = symlog 
        else:
            linthresh = max(abs(vminuse),abs(vmaxuse))*1.e-2

    # Lin or log colour scaling, defaults to log
    if lin is None:
        # Special SymLogNorm case
        if symlog is not None:
            if LooseVersion(matplotlib.__version__) < LooseVersion("3.3.0"):
                norm = SymLogNorm(linthresh=linthresh, linscale = 1.0, vmin=vminuse, vmax=vmaxuse, clip=True)
                print("WARNING: colormap SymLogNorm uses base-e but ticks are calculated with base-10.")
                #TODO: copy over matplotlib 3.3.0 implementation of SymLogNorm into pytools/analysator
            else:
                norm = SymLogNorm(base=10, linthresh=linthresh, linscale = 1.0, vmin=vminuse, vmax=vmaxuse, clip=True)
            maxlog=int(np.ceil(np.log10(vmaxuse)))
            minlog=int(np.ceil(np.log10(-vminuse)))
            logthresh=int(np.floor(np.log10(linthresh)))
            logstep=1
            ticks=([-(10**x) for x in range(logthresh, minlog+1, logstep)][::-1]
                    +[0.0]
                    +[(10**x) for x in range(logthresh, maxlog+1, logstep)] )
        else:
            # Logarithmic plot
            norm = LogNorm(vmin=vminuse,vmax=vmaxuse)
            ticks = LogLocator(base=10,subs=range(10)) # where to show labels
    else:
        # Linear
        levels = MaxNLocator(nbins=255).tick_values(vminuse,vmaxuse)
        norm = BoundaryNorm(levels, ncolors=cmapuse.N, clip=True)
        ticks = np.linspace(vminuse,vmaxuse,num=7)


    # Create the scalar mappable to define the face colouring of the surface elements
    scamap = plt.cm.ScalarMappable(cmap=colormap,norm=norm)
    scamap.set_array([])

#    # find the smallest and the largest values of the all three surface grids
#    mini = [np.ma.amin(datamap1[datamap1 > 0]), np.ma.amin(datamap2[datamap2 > 0]), np.ma.amin(datamap3[datamap3 > 0])]
#    maxi = [np.ma.amax(datamap1[datamap1 > 0]), np.ma.amax(datamap2[datamap2 > 0]), np.ma.amax(datamap3[datamap3 > 0])]


#    # the smallest value of the surface grids
#    vmin = np.ma.amin(mini)
#    # the largest value of the surface grids
#    vmax = np.ma.amax(maxi)


#    # do we have logarithmic or linear norm # NOT FINISHED
#    #norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
#    norm = LogNorm(vmin=vmin, vmax=vmax)


    # coordinates of the point where the all three 2d cut throughs cut in terms of cells
    xr = int(round((xr/(xmax - xmin))*xsize*2**reflevel))
    yr = int(round((yr/(ymax - ymin))*ysize*2**reflevel))
    zr = int(round((zr/(zmax - zmin))*zsize*2**reflevel))


    # create a new figure and a 3d axes with a custom 3d coordinate axes 
    fig = plt.figure(figsize=(9,9))
    ax1 = axes3d(fig, reflevel,  xr, yr, zr, xsize, ysize, zsize)


    #########################
    # PLOT THE 3D PLOT STARTS
    #########################

    # create 2d x, y, z arrays for a one surface element of the 3d plot
    y = np.arange(yr + 1)
    z = np.arange(zr + 1)
    Y,Z = np.meshgrid(y,z)
    X = Y*0 + xr

    # plot a surface element of the 3d plot
    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap1[:zr, :yr]),  # facecolors=plt.cm.plasma(norm(datamap1[:zr, :yr])),
                    shade=False, antialiased=False)


    y = np.arange(yr, int(ysize*2**reflevel + 1))
    z = np.arange(zr + 1)
    Y,Z = np.meshgrid(y,z)
    X = Y*0 + xr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap1[:zr, yr:]),  # facecolors=plt.cm.plasma(norm(datamap1[:zr, yr:])),
                    shade=False, antialiased=False)


    y = np.arange(yr + 1)
    z = np.arange(zr, int(zsize*2**reflevel + 1))
    Y,Z = np.meshgrid(y,z)
    X = Y*0 + xr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap1[zr:, :yr]),  # facecolors=plt.cm.plasma(norm(datamap1[zr:, :yr])),
                    shade=False, antialiased=False)


    y = np.arange(yr, int(ysize*2**reflevel + 1))
    z = np.arange(zr, int(zsize*2**reflevel + 1))
    Y,Z = np.meshgrid(y,z)
    X = Y*0 + xr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap1[zr:, yr:]),  # facecolors=plt.cm.plasma(norm(datamap1[zr:, yr:])),
                    shade=False, antialiased=False)


    x = np.arange(xr + 1)
    z = np.arange(zr + 1)
    X,Z = np.meshgrid(x,z)
    Y = X*0 + yr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap2[:zr, :xr]),  # facecolors=plt.cm.plasma(norm(datamap2[:zr, :xr])), 
                    shade=False, antialiased=False)


    x = np.arange(xr, int(xsize*2**reflevel + 1))
    z = np.arange(zr + 1)
    X,Z = np.meshgrid(x,z)
    Y = X*0 + yr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap2[:zr, xr:]),  # facecolors=plt.cm.plasma(norm(datamap2[:zr, xr:])),
                    shade=False, antialiased=False)


    x = np.arange(xr + 1)
    z = np.arange(zr, int(zsize*2**reflevel + 1))
    X,Z = np.meshgrid(x,z)
    Y = X*0 + yr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap2[zr:, :xr]),  # facecolors=plt.cm.plasma(norm(datamap2[zr:, :xr])),
                    shade=False, antialiased=False)


    x = np.arange(xr, int(xsize*2**reflevel + 1))
    z = np.arange(zr, int(zsize*2**reflevel + 1))
    X,Z = np.meshgrid(x,z)
    Y = X*0 + yr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap2[zr:, xr:]),  # facecolors=plt.cm.plasma(norm(datamap2[zr:, xr:])),
                    shade=False, antialiased=False)



    x = np.arange(xr + 1)
    y = np.arange(yr + 1)
    X,Y = np.meshgrid(x,y)
    Z = X*0 + zr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap3[:yr, :xr]),  # facecolors=plt.cm.plasma(norm(datamap3[:yr, :xr])),
                    shade=False, antialiased=False)


    x = np.arange(xr, int(xsize*2**reflevel + 1))
    y = np.arange(yr + 1)
    X,Y = np.meshgrid(x,y)
    Z = X*0 + zr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap3[:yr, xr:]),  # facecolors=plt.cm.plasma(norm(datamap3[:yr, xr:])),
                    shade=False, antialiased=False)


    x = np.arange(xr + 1)
    y = np.arange(yr, int(ysize*2**reflevel + 1))
    X,Y = np.meshgrid(x,y)
    Z = X*0 + zr

    ax1.plot_surface(X, Y, Z, rstride=1, cstride=1,
                    facecolors=scamap.to_rgba(datamap3[yr:, :xr]),  # facecolors=plt.cm.plasma(norm(datamap3[yr:, :xr])),
                    shade=False, antialiased=False)


    x = np.arange(xr, int(xsize*2**reflevel + 1))
    y = np.arange(yr, int(ysize*2**reflevel + 1))
    X,Y = np.meshgrid(x,y)
    Z = X*0 + zr

    fig1 = ax1.plot_surface(X, Y, Z, rstride=1, cstride=1, # The handle (fig1) is used to draw the colorbar further below
                    facecolors=scamap.to_rgba(datamap3[yr:, xr:]),  # facecolors=plt.cm.plasma(norm(datamap3[yr:, xr:])),
                    shade=False, antialiased=False)

    #######################
    # PLOT THE 3D PLOT ENDS
    #######################



    # Split existing axes to make room for colorbar
#    divider = make_axes_locatable(ax1)                         # There seems to be issues with make_axes_locatable with 3d axes
#    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax = fig.add_axes([0.82,0.2,0.03,0.6])                     # TODO find a cleaner way to deal with this
    cbdir="right"; horalign="left"

    # Colourbar title
    if len(cb_title_use)!=0:
        if os.getenv('PTNOLATEX'):
            cb_title_use.replace('\textbf{','')
            cb_title_use.replace('\mathrm{','')
            cb_title_use.replace('}','')
        else:
            cb_title_use = r"\textbf{"+cb_title_use+"}"   

    # Set flag which affects colorbar decimal precision
    if lin is None:
        pt.plot.cb_linear = False
    else:
        pt.plot.cb_linear = True

    # First draw colorbar
    if usesci:
        cb = plt.colorbar(scamap, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmtsci), cax=cax, drawedges=False)
    else:
        cb = plt.colorbar(scamap, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmt), cax=cax, drawedges=False)
    cb.outline.set_linewidth(thick)
    cb.ax.yaxis.set_ticks_position(cbdir)

    cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
    cb_title = cax.set_title(cb_title_use,fontsize=fontsize3,fontweight='bold', horizontalalignment=horalign)
    cb_title.set_position((0.,1.+0.025*scale)) # avoids having colourbar title too low when fontsize is increased

    # Perform intermediate draw if necessary to gain access to ticks
    if (symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2)) or (not lin and symlog is None):
        fig.canvas.draw() # draw to get tick positions

    # Adjust placement of innermost ticks for symlog if it indeed is (quasi)symmetric
    if symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2):
        cbt=cb.ax.yaxis.get_ticklabels()
        (cbtx,cbty) = cbt[len(cbt)//2-1].get_position() # just below zero
        if abs(0.5-cbty)/scale < 0.1:
            cbt[len(cbt)//2-1].set_va("top")
        (cbtx,cbty) = cbt[len(cbt)//2+1].get_position() # just above zero
        if abs(0.5-cbty)/scale < 0.1:
            cbt[len(cbt)//2+1].set_va("bottom")
        if len(cbt)>=7: # If we have at least seven ticks, may want to adjust next ones as well
            (cbtx,cbty) = cbt[len(cbt)//2-2].get_position() # second below zero
            if abs(0.5-cbty)/scale < 0.15:
                cbt[len(cbt)//2-2].set_va("top")
            (cbtx,cbty) = cbt[len(cbt)//2+2].get_position() # second above zero
            if abs(0.5-cbty)/scale < 0.15:
                cbt[len(cbt)//2+2].set_va("bottom")

    # Adjust precision for colorbar ticks
    thesetickvalues = cb.locator()
    if len(thesetickvalues)<2:
        precision_b=1
    else:
        mintickinterval = abs(thesetickvalues[-1]-thesetickvalues[0])
        # find smallest interval
        for ticki in range(len(thesetickvalues)-1):
            mintickinterval = min(mintickinterval,abs(thesetickvalues[ticki+1]-thesetickvalues[ticki]))
        precision_a, precision_b = '{:.1e}'.format(mintickinterval).split('e')
        # e.g. 9.0e-1 means we need precision 1
        # e.g. 1.33e-1 means we need precision 3?
    pt.plot.decimalprecision_cblin = 1
    if int(precision_b)<1: pt.plot.decimalprecision_cblin = str(1+abs(-int(precision_b)))
    cb.update_ticks()

    # if too many subticks in logarithmic colorbar:
    if not lin and symlog is None:
        nlabels = len(cb.ax.yaxis.get_ticklabels()) #/ ratio # TODO ratio
        valids = ['1','2','3','4','5','6','7','8','9']
        if nlabels > 10:
            valids = ['1','2','3','4','5','6','8']
        if nlabels > 19:
            valids = ['1','2','5']
        if nlabels > 28:
            valids = ['1']
        # for label in cb.ax.yaxis.get_ticklabels()[::labelincrement]:
        for labi,label in enumerate(cb.ax.yaxis.get_ticklabels()):
            labeltext = label.get_text().replace('$','').replace('{','').replace('}','').replace('\mbox{\textbf{--}}','').replace('-','').replace('.','').lstrip('0')
            if not labeltext:
                continue
            firstdigit = labeltext[0]
            if not firstdigit in valids: 
                label.set_visible(False)

    # Add Vlasiator watermark
    if (wmark or wmarkb):
        if wmark:
            wm = plt.imread(get_sample_data(watermarkimage))
        else:
            wmark=wmarkb # for checking for placement
            wm = plt.imread(get_sample_data(watermarkimageblack))
        if type(wmark) is str:
            anchor = wmark
        else:
            anchor="NW"
        # Allowed region and anchor used in tandem for desired effect
        if anchor=="NW" or anchor=="W" or anchor=="SW":
            rect = [0.01, 0.01, 0.3, 0.98]
        elif anchor=="NE" or anchor=="E" or anchor=="SE":
            rect = [0.69, 0.01, 0.3, 0.98]
        elif anchor=="N" or anchor=="C" or anchor=="S":
            rect = [0.35, 0.01, 0.3, 0.98]
        newax = fig.add_axes(rect, anchor=anchor, zorder=1)
        newax.imshow(wm)
        newax.axis('off')

    
    # Adjust layout. Uses tight_layout() but in fact this ensures 
    # that long titles and tick labels are still within the plot area.
    plt.tight_layout(pad=0.01)   # TODO check: a warning says tight_layout() might not be compatible with those axes. Seems to work though...
    savefig_pad=0.01
    bbox_inches='tight'

    # Save output to file
    try:
        plt.savefig(outputfile,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)
    except:
        print("Error with attempting to save figure.")
    print(outputfile+"\n")

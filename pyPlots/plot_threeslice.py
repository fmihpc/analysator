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

import time


# Create the 3d axes and the coordinate axes for the 3d plot
def axes3d(fig, reflevel, cutpoint, simext, boxcoords):
    # Create 3d axes
    ax = fig.add_axes([.1,.1,.64,.8],projection='3d')

#    # create 3d coordinate axes lines which cuts the point at xr, yr and zr
    line=art3d.Line3D(*zip((simext[0], cutpoint[1], cutpoint[2]), (simext[1], cutpoint[1], cutpoint[2])),
                      color='black', linewidth=0.5, alpha=1, zorder=20)
    ax.add_line(line)
    line=art3d.Line3D(*zip((cutpoint[0], simext[2], cutpoint[2]), (cutpoint[0], simext[3], cutpoint[2])),
                      color='black', linewidth=0.5, alpha=1, zorder=20)
    ax.add_line(line)
    line=art3d.Line3D(*zip((cutpoint[0], cutpoint[1], simext[4]), (cutpoint[0], cutpoint[1], simext[5])),
                      color='black', linewidth=0.5, alpha=1, zorder=20)
    ax.add_line(line)

    
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    albedo = np.zeros(u.shape)
    albedo = np.arctan(-50.*x)/np.pi+0.5
    levels = MaxNLocator(nbins=255).tick_values(0,1)
    norm = BoundaryNorm(levels, ncolors=255, clip=True)
    scalarmap = plt.cm.ScalarMappable(cmap='Greys',norm=norm)
    scalarmap.set_array([])

    ax.plot_surface(x, y, z, facecolors=scalarmap.to_rgba(albedo),alpha=1,zorder=30)

    ax.set_xlim([boxcoords[0],boxcoords[1]])
    ax.set_ylim([boxcoords[2],boxcoords[3]])
    ax.set_zlim([boxcoords[4],boxcoords[5]])

    ax.axis('equal')


#    # spacing between ticks
#    s = int((xsize*2**reflevel)/20)
#    # widths of the ticks
#    l = int((xsize*2**reflevel)/100)
#
#
#    ###################
#    # plot ticks STARTS    
#
#    # plot the first ticks which position > xr or in the others yr or zr
#    i = 1
#    while xr+i*s < xsize*2**reflevel:
#      line=art3d.Line3D(*zip((xr+i*s, yr-l, zr), (xr+i*s, yr+l, zr)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr+i*s, yr, zr-l), (xr+i*s, yr, zr+l)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#    # plot the last ticks which position < xr or in the others yr or zr
#    i = 1
#    while xr-i*s > 0:
#      line=art3d.Line3D(*zip((xr-i*s, yr-l, zr), (xr-i*s, yr+l, zr)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr-i*s, yr, zr-l), (xr-i*s, yr, zr+l)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#
#    i = 1
#    while yr+i*s < ysize*2**reflevel:
#      line=art3d.Line3D(*zip((xr-l, yr+i*s, zr), (xr+l, yr+i*s, zr)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr, yr+i*s, zr-l), (xr, yr+i*s, zr+l)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#    i = 1
#    while yr-i*s > 0:
#      line=art3d.Line3D(*zip((xr-l, yr-i*s, zr), (xr+l, yr-i*s, zr)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr, yr-i*s, zr-l), (xr, yr-i*s, zr+l)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#
#    i = 1
#    while zr+i*s < zsize*2**reflevel:
#      line=art3d.Line3D(*zip((xr, yr-l, zr+i*s), (xr, yr+l, zr+i*s)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr-l, yr, zr+i*s), (xr+l, yr, zr+i*s)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#    i = 1
#    while zr-i*s > 0:
#      line=art3d.Line3D(*zip((xr, yr-l, zr-i*s), (xr, yr+l, zr-i*s)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      line=art3d.Line3D(*zip((xr-l, yr, zr-i*s), (xr+l, yr, zr-i*s)),
#                        color='black', linewidth=0.7, alpha=0.25, zorder=100+i)
#      ax.add_line(line)
#      i += 1
#
#    # plot ticks ENDS
#    #################



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
                  boxm=None,boxre=None,
                  colormap=None, vmin=None, vmax=None,
                  symmetric=False, absolute=None,
                  lin=None, symlog=None,
                  cbtitle=None, title=None,
                  usesci=True, axisunit=None,
                  pass_full=None,
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

    :kword boxm:        zoom box extents [x0,x1,y0,y1,z0,z1] in metres (default and truncate to: whole simulation box)
    :kword boxre:       zoom box extents [x0,x1,y0,y1,z0,z1] in Earth radii (default and truncate to: whole simulation box)

    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, parula, nipy_spectral, RdBu, bwr
    :kword vmin,vmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword symmetric:   Set the absolute value of vmin and vmax to the greater of the two
    :kword absolute:    Plot the absolute of the evaluated variable
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2, but this can
                        result in the innermost tick marks overlapping. In this case, using a larger value for 
                        symlog is suggested.
    :kword title:       string to use as plot title instead of time.
                        Special case: Set to "msec" to plot time with millisecond accuracy or "musec"
                        for microsecond accuracy. "sec" is integer second accuracy.
    :kword cbtitle:     string to use as colorbar title instead of map name
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: True)
    :kword axisunit:    Plot axes using 10^{axisunit} m (default: Earth radius R_E)

    :kword pass_full:   Set to anything but None in order to pass the full arrays instead of a zoomed-in section

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

    t0 = time.time()

    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    watermarkimageblack=os.path.join(os.path.dirname(__file__), 'logo_black.png')

    # Switch None-keywords to empty lists (this way subsequent calls get correct empty default values
    if boxm is None:
        boxm=[],
    if boxre is None:
        boxre=[]
#    if pass_vars is None:    # Not yet implemented
#        pass_vars=[]

    # Change certain falsy values:
    if not lin and lin is not 0:
        lin = None
    if not symlog and symlog is not 0:
        symlog = None
    if symlog is True:
        symlog = 0
    if (filedir is ''):
        filedir = './'
    if (outputdir is ''):
        outputdir = './'

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

    if not operator:
        if op:
            operator=op

    if not colormap:
        # Default values
        colormap="hot_desaturated"
        if operator and operator in 'xyz':
            colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=8*scale # Colour bar ticks and title

    # Plot title with time
    timeval=f.read_parameter("time")

    # Plot title with time
    if title is None or title=="msec" or title=="musec":        
        if timeval is None:    
            plot_title = ''
        else:
            timeformat='{:4.1f}'
            if title=="sec": timeformat='{:4.0f}'
            if title=="msec": timeformat='{:4.3f}'
            if title=="musec": timeformat='{:4.6f}'
            plot_title = "t="+timeformat.format(timeval)+' s'
    else:
        plot_title = title

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

    # The plot will be saved in a new figure ('draw' and 'axes' keywords not implemented yet)
    if str(matplotlib.get_backend()) is not pt.backend_noninteractive: #'Agg':
        plt.switch_backend(pt.backend_noninteractive)
    plt.switch_backend('TkAgg')

    Re = 6.371e+6 # Earth radius in m
    # read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    xsize = int(xsize)
    ysize = int(ysize)
    zsize = int(zsize)
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")
    sizes = [xsize,ysize,zsize]

    # Read the FSgrid mesh
    [xsizefg, ysizefg, zsizefg] = f.get_fsgrid_mesh_size()
    xsizefg = int(xsizefg)
    ysizefg = int(ysizefg)
    zsizefg = int(zsizefg)
    [xminfg, yminfg, zminfg, xmaxfg, ymaxfg, zmaxfg] = f.get_fsgrid_mesh_extent()
    cellsizefg = (xmaxfg-xminfg)/xsizefg
    pt.plot.plot_helpers.CELLSIZE = cellsizefg
    
    # sort the cellid and the datamap list
    indexids = cellids.argsort()
    cellids = cellids[indexids]

    # find the highest refiment level
    reflevel = ids3d.refinement_level(xsize, ysize, zsize, cellids[-1])
    for i in range(5): # Check if Vlasov grid doesn't reach maximum (fsgrid) refinement
        if xsize*(2**(reflevel + i)) == xsizefg:
            reflevel += i
            break

    # Verify that FSgrid and spatial grid agree
    if ((xmin!=xminfg) or (xmax!=xmaxfg) or
        (ymin!=yminfg) or (ymax!=ymaxfg) or
        (zmin!=zminfg) or (zmax!=zmaxfg) or
        (xsize*(2**reflevel) !=xsizefg) or (ysize*(2**reflevel) !=ysizefg) or (zsize*(2**reflevel) !=zsizefg)):
        print("FSgrid and vlasov grid disagreement!")
        return -1

    # Simulation domain outer boundaries
    simext=[xmin,xmax,ymin,ymax,zmin,zmax]

    # Coordinates of the point through which the three slices pass
    if not cutpoint:
        if cutpointre:
            cutpoint = np.asarray(cutpointre) * Re
        else: # default to [0,0,0]
            print('No cut point coordinates given, defaulting to origin')
            cutpoint = np.asarray([0.,0.,0.])

    ###################################################
    # Find the cellids corresponding to the 3 slices #
    ###################################################

    # {X = x0} slice
    sliceoffset = abs(xmin) + cutpoint[0]
    fgslice_x = int(sliceoffset/cellsizefg)
    idlist_x, indexlist_x = ids3d.ids3d(cellids, sliceoffset, reflevel, xsize, ysize, zsize, xmin=xmin, xmax=xmax)

    # {Y = y0} slice
    sliceoffset = abs(ymin) + cutpoint[1]
    fgslice_y = int(sliceoffset/cellsizefg)
    idlist_y, indexlist_y = ids3d.ids3d(cellids, sliceoffset, reflevel, xsize, ysize, zsize, ymin=ymin, ymax=ymax)

    # {Z = z0} slice
    sliceoffset = abs(zmin) + cutpoint[2]
    fgslice_z = int(sliceoffset/cellsizefg)
    idlist_z, indexlist_z = ids3d.ids3d(cellids, sliceoffset, reflevel, xsize, ysize, zsize, zmin=zmin, zmax=zmax)


    # Select window to draw
    if len(boxm)==6:
        boxcoords=list(boxm)
    elif len(boxre)==6:
        boxcoords=[i*Re for i in boxre]
    else:
        boxcoords=list(simext)

    # If box extents were provided manually, truncate to simulation extents
    # Also subtract one reflevel0-cell in each direction to hide boundary cells
    boxcoords[0] = max(boxcoords[0],simext[0]+cellsize)
    boxcoords[1] = min(boxcoords[1],simext[1]-cellsize)
    boxcoords[2] = max(boxcoords[2],simext[2]+cellsize)
    boxcoords[3] = min(boxcoords[3],simext[3]-cellsize)
    boxcoords[4] = max(boxcoords[4],simext[4]+cellsize)
    boxcoords[5] = min(boxcoords[5],simext[5]-cellsize)

    # Axes and units (default R_E)
    if axisunit is not None: # Use m or km or other
        if np.isclose(axisunit,0):
            axisunitstr = r'm'
        elif np.isclose(axisunit,3):
            axisunitstr = r'km'
        else:
            axisunitstr = r'$10^{'+str(int(axisunit))+'}$ m'
        axisunit = np.power(10,int(axisunit))
    else:
        axisunitstr = r'$\mathrm{R}_{\mathrm{E}}$'
        axisunit = Re

    # Scale data extent and plot box
    simext=[i/axisunit for i in simext]
    boxcoords=[i/axisunit for i in boxcoords]


    #################################################
    # Find rhom map for use in masking out ionosphere
    #################################################
    if f.check_variable("moments"):
        rhomap = f.read_variable("vg_restart_rhom")
    elif f.check_variable("proton/vg_rho"):
        rhomap = f.read_variable("proton/vg_rho")
    elif f.check_variable("proton/vg_rho"):
        rhomap = f.read_variable("vg_rhom")
    else:
        print("error!")
        quit
              
    rhomap = rhomap[indexids]  # sort
    rhomap_x0 = rhomap[indexlist_x] # find required cells (X cut)
    rhomap_y0 = rhomap[indexlist_y] # find required cells (Y cut)
    rhomap_z0 = rhomap[indexlist_z] # find required cells (Z cut)
    # Create the plotting grid
    rhomap_x0 = ids3d.idmesh3d(idlist_x, rhomap_x0, reflevel, xsize, ysize, zsize, 0, None)
    rhomap_y0 = ids3d.idmesh3d(idlist_y, rhomap_y0, reflevel, xsize, ysize, zsize, 1, None)
    rhomap_z0 = ids3d.idmesh3d(idlist_z, rhomap_z0, reflevel, xsize, ysize, zsize, 2, None)


    ############################################
    # Read data and calculate required variables
    ############################################
    if not expression:
        # Read data from file
        if not operator:
            operator="pass"
        datamap_info = f.read_variable_info(var, operator=operator)

        cb_title_use = datamap_info.latex
        datamap_unit = datamap_info.latexunits
        # Check if vscale results in standard unit
        datamap_unit = pt.plot.scaleunits(datamap_info, vscale)

        # Add unit to colorbar title
        if datamap_unit:
            cb_title_use = cb_title_use + "\,["+datamap_unit+"]"

        datamap = datamap_info.data

        # Verify data shape
        if np.ndim(datamap)==0:
            print("Error, read only single value from vlsv file!",datamap.shape)
            return -1

        if var.startswith('fg_'):
            # fsgrid reader returns array in correct shape but needs to be sliced and transposed
            if np.ndim(datamap)==3: # scalar variable
                datamap_x = datamap[fgslice_x,:,:]
                datamap_y = datamap[:,fgslice_y,:]
                datamap_z = datamap[:,:,fgslice_z]
            elif np.ndim(datamap)==4: # vector variable
                datamap_x = datamap[fgslice_x,:,:,:]
                datamap_y = datamap[:,fgslice_y,:,:]
                datamap_z = datamap[:,:,fgslice_z,:]
            elif np.ndim(datamap)==5:  # tensor variable
                datamap_x = datamap[fgslice_x,:,:,:,:]
                datamap_y = datamap[:,fgslice_y,:,:,:]
                datamap_z = datamap[:,:,fgslice_z,:,:]
            else:
                print("Error in reshaping fsgrid datamap!") 
            datamap_x = np.squeeze(datamap_x)
            datamap_x = np.swapaxes(datamap_x, 0,1)
            datamap_y = np.squeeze(datamap_y)
            datamap_y = np.swapaxes(datamap_y, 0,1)
            datamap_z = np.squeeze(datamap_z)
            datamap_z = np.swapaxes(datamap_z, 0,1)

            print('(FSgrid) Size(datamap_y=({:d},{:d})'.format(datamap_y.shape[0],datamap_y.shape[1]))

        else:
            # vlasov grid, AMR
            datamap = datamap[indexids]      # sort
            datamap_x = datamap[indexlist_x] # find required cells (X cut)
            datamap_y = datamap[indexlist_y] # find required cells (Y cut)
            datamap_z = datamap[indexlist_z] # find required cells (Z cut)
            # Create the plotting grid
            if np.ndim(datamap)==1:   # scalar variable
                datamap_x = ids3d.idmesh3d(idlist_x, datamap_x, reflevel, xsize, ysize, zsize, 0, None)
                datamap_y = ids3d.idmesh3d(idlist_y, datamap_y, reflevel, xsize, ysize, zsize, 1, None)
                datamap_z = ids3d.idmesh3d(idlist_z, datamap_z, reflevel, xsize, ysize, zsize, 2, None)
            elif np.ndim(datamap)==2: # vector variable
                datamap_x = ids3d.idmesh3d(idlist_x, datamap_x, reflevel, xsize, ysize, zsize, 0, datamap.shape[1])
                datamap_y = ids3d.idmesh3d(idlist_y, datamap_y, reflevel, xsize, ysize, zsize, 1, datamap.shape[1])
                datamap_z = ids3d.idmesh3d(idlist_z, datamap_z, reflevel, xsize, ysize, zsize, 2, datamap.shape[1])
            elif np.ndim(datamap)==3: # tensor variable
                datamap_x = ids3d.idmesh3d(idlist_x, datamap_x, reflevel, xsize, ysize, zsize, 0, (datamap.shape[1],datamap.shape[2]))
                datamap_y = ids3d.idmesh3d(idlist_y, datamap_y, reflevel, xsize, ysize, zsize, 1, (datamap.shape[1],datamap.shape[2]))
                datamap_z = ids3d.idmesh3d(idlist_z, datamap_z, reflevel, xsize, ysize, zsize, 2, (datamap.shape[1],datamap.shape[2]))
            else:
                print("Dimension error in constructing 2D AMR slice!")
                return -1
            print('(AMR) Size(datamap_y=({:d},{:d})'.format(datamap_y.shape[0],datamap_y.shape[1]))

    else:
        # Expression set, use generated or provided colorbar title
        cb_title_use = expression.__name__ + operatorstr
        print('WARNING: Expressions have not been implemented yet')

    # Now, if map is a vector or tensor, reduce it down
    if np.ndim(datamap_x)==3: # vector
        if datamap_x.shape[2]!=3:
            print("Error, expected array of 3-element vectors, found array of shape ",datamap_x.shape)
            return -1
        # take magnitude of three-element vectors
        datamap_x = np.linalg.norm(datamap_x, axis=-1)
        datamap_y = np.linalg.norm(datamap_y, axis=-1)
        datamap_z = np.linalg.norm(datamap_z, axis=-1)
    if np.ndim(datamap_x)==4: # tensor
        if datamap_x.shape[2]!=3 or datamap_x.shape[3]!=3:
            # This may also catch 3D simulation fsgrid variables
            print("Error, expected array of 3x3 tensors, found array of shape ",datamap_x.shape)
            return -1
        # take trace
        datamap_x = datamap_x[:,:,0,0]+datamap_x[:,:,1,1]+datamap_x[:,:,2,2]
        datamap_y = datamap_y[:,:,0,0]+datamap_y[:,:,1,1]+datamap_y[:,:,2,2]
        datamap_z = datamap_z[:,:,0,0]+datamap_z[:,:,1,1]+datamap_z[:,:,2,2]
    if np.ndim(datamap_x)>=5: # Too many dimensions
        print("Error, too many dimensions in datamap, found array of shape ",datamap_x.shape)
        return -1
    if np.ndim(datamap_x)!=2: # Too many dimensions
        print("Error, too many dimensions in datamap, found array of shape ",datamap_x.shape)
        return -1
        
    # Scale final generated datamap if requested
    datamap_x = datamap_x * vscale
    datamap_y = datamap_y * vscale
    datamap_z = datamap_z * vscale
    
    # Take absolute
    if (absolute):
        datamap_x = abs(datamap_x)
        datamap_y = abs(datamap_y)
        datamap_z = abs(datamap_z)

    # scale the sizes to the heighest refinement level because
    # plotting is done at that level
    sizes[0] = int(sizes[0]*2**reflevel)
    sizes[1] = int(sizes[1]*2**reflevel)
    sizes[2] = int(sizes[2]*2**reflevel)
    finecellsize = cellsize / 2**reflevel

    # Allow title override
    if cbtitle is not None:
        # Here allow underscores for manual math mode
        cb_title_use = cbtitle

    #If automatic range finding is required, find min and max of array
    # Performs range-finding on a masked array to work even if array contains invalid values
    if vmin is not None:
        vminuse=vmin
    else: 
        vminuse=np.ma.amin([np.ma.amin(datamap_x),np.ma.amin(datamap_y),np.ma.amin(datamap_z)])
    if vmax is not None:
        vmaxuse=vmax
    else:
        vmaxuse=np.ma.amax([np.ma.amax(datamap_x),np.ma.amax(datamap_y),np.ma.amax(datamap_z)])

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
        vminuse = np.ma.amin([np.ma.amin(np.ma.masked_less_equal(datamap_x,0)),np.ma.amin(np.ma.masked_less_equal(datamap_y,0)),np.ma.amin(np.ma.masked_less_equal(datamap_z,0))])

    # Special case of very small vminuse values
    if ((vmin is None) or (vmax is None)) and (vminuse > 0) and (vminuse < vmaxuse*1.e-5):
        vminuse = vmaxuse*1e-5
        if lin is not None:
            vminuse = 0
    
    # If symlog scaling is set:
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


    ###############################################################################
    # Making the 12 meshes corresponding to the 12 elementary surfaces to plot #
    ###############################################################################
    # {X = x0 slice}
    [YmeshYmZm,ZmeshYmZm] = scipy.meshgrid(np.linspace(simext[2],cutpoint[1],num=int(round((cutpoint[1]-ymin)/finecellsize))+1),
                               np.linspace(simext[4],cutpoint[2],num=int(round((cutpoint[2]-zmin)/finecellsize))+1))
    XmeshYmZm = np.ones(YmeshYmZm.shape) * cutpoint[0]

    [YmeshYpZm,ZmeshYpZm] = scipy.meshgrid(np.linspace(cutpoint[1],simext[3],num=int(round((ymax-cutpoint[1])/finecellsize))+1),
                               np.linspace(simext[4],cutpoint[2],num=int(round((cutpoint[2]-zmin)/finecellsize))+1))
    XmeshYpZm = np.ones(YmeshYpZm.shape) * cutpoint[0]

    [YmeshYmZp,ZmeshYmZp] = scipy.meshgrid(np.linspace(simext[2],cutpoint[1],num=int(round((cutpoint[1]-ymin)/finecellsize))+1),
                               np.linspace(cutpoint[2],simext[5],num=int(round((zmax-cutpoint[2])/finecellsize))+1))
    XmeshYmZp = np.ones(YmeshYmZp.shape) * cutpoint[0]

    [YmeshYpZp,ZmeshYpZp] = scipy.meshgrid(np.linspace(cutpoint[1],simext[3],num=int(round((ymax-cutpoint[1])/finecellsize))+1),
                               np.linspace(cutpoint[2],simext[5],num=int(round((zmax-cutpoint[2])/finecellsize))+1))
    XmeshYpZp = np.ones(YmeshYpZp.shape) * cutpoint[0]

    # {Y = y0 slice}
    [XmeshXmZm,ZmeshXmZm] = scipy.meshgrid(np.linspace(simext[0],cutpoint[0],num=int(round((cutpoint[0]-xmin)/finecellsize))+1),
                               np.linspace(simext[4],cutpoint[2],num=int(round((cutpoint[2]-zmin)/finecellsize))+1))
    YmeshXmZm = np.ones(XmeshXmZm.shape) * cutpoint[1]

    [XmeshXpZm,ZmeshXpZm] = scipy.meshgrid(np.linspace(cutpoint[0],simext[1],num=int(round((xmax-cutpoint[0])/finecellsize))+1),
                               np.linspace(simext[4],cutpoint[2],num=int(round((cutpoint[2]-zmin)/finecellsize))+1))
    YmeshXpZm = np.ones(XmeshXpZm.shape) * cutpoint[1]

    [XmeshXmZp,ZmeshXmZp] = scipy.meshgrid(np.linspace(simext[0],cutpoint[0],num=int(round((cutpoint[0]-xmin)/finecellsize))+1),
                               np.linspace(cutpoint[2],simext[5],num=int(round((zmax-cutpoint[2])/finecellsize))+1))
    YmeshXmZp = np.ones(XmeshXmZp.shape) * cutpoint[1]

    [XmeshXpZp,ZmeshXpZp] = scipy.meshgrid(np.linspace(cutpoint[0],simext[1],num=int(round((xmax-cutpoint[0])/finecellsize))+1),
                               np.linspace(cutpoint[2],simext[5],num=int(round((zmax-cutpoint[2])/finecellsize))+1))
    YmeshXpZp = np.ones(XmeshXpZp.shape) * cutpoint[1]

    # {Z = z0 slice}
    [XmeshXmYm,YmeshXmYm] = scipy.meshgrid(np.linspace(simext[0],cutpoint[0],num=int(round((cutpoint[0]-xmin)/finecellsize))+1),
                               np.linspace(simext[2],cutpoint[1],num=int(round((cutpoint[1]-ymin)/finecellsize))+1))
    ZmeshXmYm = np.ones(XmeshXmYm.shape) * cutpoint[2]

    [XmeshXpYm,YmeshXpYm] = scipy.meshgrid(np.linspace(cutpoint[0],simext[1],num=int(round((xmax-cutpoint[0])/finecellsize))+1),
                               np.linspace(simext[2],cutpoint[1],num=int(round((cutpoint[1]-ymin)/finecellsize))+1))
    ZmeshXpYm = np.ones(XmeshXpYm.shape) * cutpoint[2]

    [XmeshXmYp,YmeshXmYp] = scipy.meshgrid(np.linspace(simext[0],cutpoint[0],num=int(round((cutpoint[0]-xmin)/finecellsize))+1),
                               np.linspace(cutpoint[1],simext[3],num=int(round((ymax-cutpoint[1])/finecellsize))+1))
    ZmeshXmYp = np.ones(XmeshXmYp.shape) * cutpoint[2]

    [XmeshXpYp,YmeshXpYp] = scipy.meshgrid(np.linspace(cutpoint[0],simext[1],num=int(round((xmax-cutpoint[0])/finecellsize))+1),
                               np.linspace(cutpoint[1],simext[3],num=int(round((ymax-cutpoint[1])/finecellsize))+1))
    ZmeshXpYp = np.ones(XmeshXpYp.shape) * cutpoint[2]

    # Creating lists of meshes to be called in a for loop
    Xmesh_list = [XmeshYmZm,XmeshYpZm,XmeshYmZp,XmeshYpZp, XmeshXmZm,XmeshXpZm,XmeshXmZp,XmeshXpZp, XmeshXmYm,XmeshXpYm,XmeshXmYp,XmeshXpYp]
    Ymesh_list = [YmeshYmZm,YmeshYpZm,YmeshYmZp,YmeshYpZp, YmeshXmZm,YmeshXpZm,YmeshXmZp,YmeshXpZp, YmeshXmYm,YmeshXpYm,YmeshXmYp,YmeshXpYp]
    Zmesh_list = [ZmeshYmZm,ZmeshYpZm,ZmeshYmZp,ZmeshYpZp, ZmeshXmZm,ZmeshXpZm,ZmeshXmZp,ZmeshXpZp, ZmeshXmYm,ZmeshXpYm,ZmeshXmYp,ZmeshXpYp]

    # coordinates of the point where the all three 2d cut throughs cuts [TODO update this part]
    xr = abs(xmin) + cutpoint[0]
    yr = abs(ymin) + cutpoint[1]
    zr = abs(zmin) + cutpoint[2]

    # coordinates of the point where the all three 2d cut throughs cut in terms of cells
    xr = int(round((xr/(xmax - xmin))*xsize*2**reflevel))
    yr = int(round((yr/(ymax - ymin))*ysize*2**reflevel))
    zr = int(round((zr/(zmax - zmin))*zsize*2**reflevel))

    # Creating lists of datamap_i to be called in that same for loop
    datamap_x_list = [datamap_x[:zr,:yr],datamap_x[:zr,yr:],datamap_x[zr:,:yr],datamap_x[zr:,yr:]]
    datamap_y_list = [datamap_y[:zr,:xr],datamap_y[:zr,xr:],datamap_y[zr:,:xr],datamap_y[zr:,xr:]]
    datamap_z_list = [datamap_z[:yr,:xr],datamap_z[:yr,xr:],datamap_z[yr:,:xr],datamap_z[yr:,xr:]]

    # Creating a new figure and a 3d axes with a custom 3d coordinate axes 
    fig = plt.figure(figsize=(6,5),dpi=300)
    ax1 = axes3d(fig, reflevel,  cutpoint/axisunit, simext, boxcoords)

    # Masking and plotting the elementary surfaces one by one (actually three by three)
    for i in range(0,4):
        print('i = '+str(i)+'. Time since start = {:.2f} s'.format(time.time()-t0))
        XmeshYZ = Xmesh_list[i]
        YmeshYZ = Ymesh_list[i]
        ZmeshYZ = Zmesh_list[i]

        XmeshXZ = Xmesh_list[i+4]
        YmeshXZ = Ymesh_list[i+4]
        ZmeshXZ = Zmesh_list[i+4]

        XmeshXY = Xmesh_list[i+8]
        YmeshXY = Ymesh_list[i+8]
        ZmeshXY = Zmesh_list[i+8]

        datamap_x_i = datamap_x_list[i]
        datamap_y_i = datamap_y_list[i]
        datamap_z_i = datamap_z_list[i]

        # The grid generated by meshgrid has all four corners for each cell.
        # We mask using only the centre values.
        # Calculate offsets for cell-centre coordinates
        XmeshXYCentres = XmeshXY[:-1,:-1] + 0.5*(XmeshXY[0,1]-XmeshXY[0,0])
        YmeshXYCentres = YmeshXY[:-1,:-1] + 0.5*(YmeshXY[1,0]-YmeshXY[0,0])

        XmeshXZCentres = XmeshXZ[:-1,:-1] + 0.5*(XmeshXZ[0,1]-XmeshXZ[0,0])
        ZmeshXZCentres = ZmeshXZ[:-1,:-1] + 0.5*(ZmeshXZ[1,0]-ZmeshXZ[0,0])

        YmeshYZCentres = YmeshYZ[:-1,:-1] + 0.5*(YmeshYZ[0,1]-YmeshYZ[0,0])
        ZmeshYZCentres = ZmeshYZ[:-1,:-1] + 0.5*(ZmeshYZ[1,0]-ZmeshYZ[0,0])

        maskgrid_XY = np.ma.array(XmeshXYCentres)
        maskgrid_XZ = np.ma.array(XmeshXZCentres)
        maskgrid_YZ = np.ma.array(YmeshYZCentres)

        if not pass_full:
            # If zoomed-in using a defined box, and not specifically asking to pass all values:
            # Generate mask for only visible section (with small buffer for e.g. gradient calculations)
            maskboundarybuffer = 2.*cellsize/axisunit
            maskgrid_XY = np.ma.masked_where(XmeshXYCentres<(boxcoords[0]-maskboundarybuffer), maskgrid_XY)
            maskgrid_XY = np.ma.masked_where(XmeshXYCentres>(boxcoords[1]+maskboundarybuffer), maskgrid_XY)
            maskgrid_XY = np.ma.masked_where(YmeshXYCentres<(boxcoords[2]-maskboundarybuffer), maskgrid_XY)
            maskgrid_XY = np.ma.masked_where(YmeshXYCentres>(boxcoords[3]+maskboundarybuffer), maskgrid_XY)

            maskgrid_XZ = np.ma.masked_where(XmeshXZCentres<(boxcoords[0]-maskboundarybuffer), maskgrid_XZ)
            maskgrid_XZ = np.ma.masked_where(XmeshXZCentres>(boxcoords[1]+maskboundarybuffer), maskgrid_XZ)
            maskgrid_XZ = np.ma.masked_where(ZmeshXZCentres<(boxcoords[4]-maskboundarybuffer), maskgrid_XZ)
            maskgrid_XZ = np.ma.masked_where(ZmeshXZCentres>(boxcoords[5]+maskboundarybuffer), maskgrid_XZ)

            maskgrid_YZ = np.ma.masked_where(YmeshYZCentres<(boxcoords[2]-maskboundarybuffer), maskgrid_YZ)
            maskgrid_YZ = np.ma.masked_where(YmeshYZCentres>(boxcoords[3]+maskboundarybuffer), maskgrid_YZ)
            maskgrid_YZ = np.ma.masked_where(ZmeshYZCentres<(boxcoords[4]-maskboundarybuffer), maskgrid_YZ)
            maskgrid_YZ = np.ma.masked_where(ZmeshYZCentres>(boxcoords[5]+maskboundarybuffer), maskgrid_YZ)
    
        if np.ma.is_masked(maskgrid_XY):
            # Save lists for masking
            MaskXY_X = np.where(~np.all(maskgrid_XY.mask, axis=1))[0] # [0] takes the first element of a tuple
            MaskXY_Y = np.where(~np.all(maskgrid_XY.mask, axis=0))[0]
            XmeshXYPass = XmeshXY[MaskXY_X[0]:MaskXY_X[-1]+2,:]
            XmeshXYPass = XmeshXYPass[:,MaskXY_Y[0]:MaskXY_Y[-1]+2]
            YmeshXYPass = YmeshXY[MaskXY_X[0]:MaskXY_X[-1]+2,:]
            YmeshXYPass = YmeshXYPass[:,MaskXY_Y[0]:MaskXY_Y[-1]+2]
            ZmeshXYPass = ZmeshXY[MaskXY_X[0]:MaskXY_X[-1]+2,:]
            ZmeshXYPass = ZmeshXYPass[:,MaskXY_Y[0]:MaskXY_Y[-1]+2]
#            XmeshXYCentres = XmeshXYCentres[MaskXY_X[0]:MaskXY_X[-1]+1,:]
#            XmeshXYCentres = XmeshXYCentres[:,MaskXY_Y[0]:MaskXY_Y[-1]+1]
#            YmeshXYCentres = YmeshXYCentres[MaskXY_X[0]:MaskXY_X[-1]+1,:]
#            YmeshXYCentres = YmeshXYCentres[:,MaskXY_Y[0]:MaskXY_Y[-1]+1]
#            ZmeshXYCentres = ZmeshXYCentres[MaskXY_X[0]:MaskXY_X[-1]+1,:]
#            ZmeshXYCentres = ZmeshXYCentres[:,MaskXY_Y[0]:MaskXY_Y[-1]+1]
        else:
            XmeshXYPass = np.ma.array(XmeshXY)
            YmeshXYPass = np.ma.array(YmeshXY)
            ZmeshXYPass = np.ma.array(ZmeshXY)

        if np.ma.is_masked(maskgrid_XZ):
            # Save lists for masking
            MaskXZ_X = np.where(~np.all(maskgrid_XZ.mask, axis=1))[0] # [0] takes the first element of a tuple
            MaskXZ_Z = np.where(~np.all(maskgrid_XZ.mask, axis=0))[0]
            XmeshXZPass = XmeshXZ[MaskXZ_X[0]:MaskXZ_X[-1]+2,:]
            XmeshXZPass = XmeshXZPass[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+2]
            YmeshXZPass = YmeshXZ[MaskXZ_X[0]:MaskXZ_X[-1]+2,:]
            YmeshXZPass = YmeshXZPass[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+2]
            ZmeshXZPass = ZmeshXZ[MaskXZ_X[0]:MaskXZ_X[-1]+2,:]
            ZmeshXZPass = ZmeshXZPass[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+2]
#            XmeshXZCentres = XmeshXZCentres[MaskXZ_X[0]:MaskXZ_X[-1]+1,:]
#            XmeshXZCentres = XmeshXZCentres[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+1]
#            YmeshXZCentres = YmeshXZCentres[MaskXZ_X[0]:MaskXZ_X[-1]+1,:]
#            YmeshXZCentres = YmeshXZCentres[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+1]
#            ZmeshXZCentres = ZmeshXZCentres[MaskXZ_X[0]:MaskXZ_X[-1]+1,:]
#            ZmeshXZCentres = ZmeshXZCentres[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+1]
        else:
            XmeshXZPass = np.ma.array(XmeshXZ)
            YmeshXZPass = np.ma.array(YmeshXZ)
            ZmeshXZPass = np.ma.array(ZmeshXZ)

        if np.ma.is_masked(maskgrid_YZ):
            # Save lists for masking
            MaskYZ_Y = np.where(~np.all(maskgrid_YZ.mask, axis=1))[0] # [0] takes the first element of a tuple
            MaskYZ_Z = np.where(~np.all(maskgrid_YZ.mask, axis=0))[0]
            XmeshYZPass = XmeshYZ[MaskYZ_Y[0]:MaskYZ_Y[-1]+2,:]
            XmeshYZPass = XmeshYZPass[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+2]
            YmeshYZPass = YmeshYZ[MaskYZ_Y[0]:MaskYZ_Y[-1]+2,:]
            YmeshYZPass = YmeshYZPass[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+2]
            ZmeshYZPass = ZmeshYZ[MaskYZ_Y[0]:MaskYZ_Y[-1]+2,:]
            ZmeshYZPass = ZmeshYZPass[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+2]
#            XmeshYZCentres = XmeshYZCentres[MaskYZ_Y[0]:MaskYZ_Y[-1]+1,:]
#            XmeshYZCentres = XmeshYZCentres[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+1]
#            YmeshYZCentres = YmeshYZCentres[MaskYZ_Y[0]:MaskYZ_Y[-1]+1,:]
#            YmeshYZCentres = YmeshYZCentres[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+1]
#            ZmeshYZCentres = ZmeshYZCentres[MaskYZ_Y[0]:MaskYZ_Y[-1]+1,:]
#            ZmeshYZCentres = ZmeshYZCentres[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+1]
        else:
            XmeshYZPass = np.ma.array(XmeshYZ)
            YmeshYZPass = np.ma.array(YmeshYZ)
            ZmeshYZPass = np.ma.array(ZmeshYZ)
    
        # Crop both rhomap and datamap to view region
        if np.ma.is_masked(maskgrid_XY):
            # Strip away columns and rows which are outside the plot region
            rhomap_z = rhomap_z0[MaskXY_X[0]:MaskXY_X[-1]+1,:]
            rhomap_z = rhomap_z[:,MaskXY_Y[0]:MaskXY_Y[-1]+1]
            # Also for the datamap, unless it was already provided by an expression
            if not expression:
                datamap_z_i = datamap_z_i[MaskXY_X[0]:MaskXY_X[-1]+1,:]
                datamap_z_i = datamap_z_i[:,MaskXY_Y[0]:MaskXY_Y[-1]+1]

        if np.ma.is_masked(maskgrid_XZ):
            # Strip away columns and rows which are outside the plot region
            rhomap_y = rhomap_y0[MaskXZ_X[0]:MaskXZ_X[-1]+1,:]
            rhomap_y = rhomap_y[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+1]
            # Also for the datamap, unless it was already provided by an expression
            if not expression:
                datamap_y_i = datamap_y_i[MaskXZ_X[0]:MaskXZ_X[-1]+1,:]
                datamap_y_i = datamap_y_i[:,MaskXZ_Z[0]:MaskXZ_Z[-1]+1]

        if np.ma.is_masked(maskgrid_YZ):
            # Strip away columns and rows which are outside the plot region
            rhomap_x = rhomap_x0[MaskYZ_Y[0]:MaskYZ_Y[-1]+1,:]
            rhomap_x = rhomap_x[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+1]
            # Also for the datamap, unless it was already provided by an expression
            if not expression:
                datamap_x_i = datamap_x_i[MaskYZ_Y[0]:MaskYZ_Y[-1]+1,:]
                datamap_x_i = datamap_x_i[:,MaskYZ_Z[0]:MaskYZ_Z[-1]+1]
    
        # Mask region outside ionosphere. Note that for some boundary layer cells, 
        # a density is calculated, but e.g. pressure is not, and these cells aren't
        # excluded by this method. Also mask away regions where datamap is invalid
        rhomap_z = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap_z), 0)
        rhomap_z = np.ma.masked_where(~np.isfinite(datamap_z_i), rhomap_z)
        XYmask_z = rhomap_z.mask
        if XYmask_z.any():
            if XYmask_z.all():
                # if everything was masked in rhomap, allow plotting
                XYmask_z[:,:] = False
            else:
                # Mask datamap
                datamap_z_i = np.ma.array(datamap_z_i, mask=XYmask_z)

        rhomap_y = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap_y), 0)
        rhomap_y = np.ma.masked_where(~np.isfinite(datamap_y_i), rhomap_y)
        XZmask_y = rhomap_y.mask
        if XZmask_y.any():
            if XZmask_y.all():
                # if everything was masked in rhomap, allow plotting
                XZmask_y[:,:] = False
            else:
                # Mask datamap
                datamap_y_i = np.ma.array(datamap_y_i, mask=XZmask_y)

        rhomap_x = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap_x), 0)
        rhomap_x = np.ma.masked_where(~np.isfinite(datamap_x_i), rhomap_x)
        YZmask_x = rhomap_x.mask
        if YZmask_x.any():
            if YZmask_x.all():
                # if everything was masked in rhomap, allow plotting
                YZmask_x[:,:] = False
            else:
                # Mask datamap
                datamap_x_i = np.ma.array(datamap_x_i, mask=YZmask_x)

        # Building the face colour maps for plot_surface()
        if np.ma.isMaskedArray(datamap_x_i):
                fcolor_x_i = scamap.to_rgba(datamap_x_i.data)
                fcolor_x_i[YZmask_x] = np.array([0,0,0,0])
        else:
                fcolor_x_i = scamap.to_rgba(datamap_x_i)

        if np.ma.isMaskedArray(datamap_y_i):
                fcolor_y_i = scamap.to_rgba(datamap_y_i.data)
                fcolor_y_i[XZmask_y] = np.array([0,0,0,0])
        else:
                fcolor_y_i = scamap.to_rgba(datamap_y_i)

        if np.ma.isMaskedArray(datamap_z_i):
                fcolor_z_i = scamap.to_rgba(datamap_z_i.data)
                fcolor_z_i[XYmask_z] = np.array([0,0,0,0])
        else:
                fcolor_z_i = scamap.to_rgba(datamap_z_i)

        print('Ready to plot a set of surface elements! Time since start = {:.2f} s'.format(time.time()-t0))
        # Plotting the partial {X = x0} cut
        ax1.plot_surface(XmeshYZPass, YmeshYZPass, ZmeshYZPass, rstride=1, cstride=1,
                    facecolors=fcolor_x_i, shade=False, antialiased=False)

        # Plotting the partial {Y = y0} cut
        ax1.plot_surface(XmeshXZPass, YmeshXZPass, ZmeshXZPass, rstride=1, cstride=1,
                    facecolors=fcolor_y_i, shade=False, antialiased=False)

        # Plotting the partial {Z = z0} cut
        ax1.plot_surface(XmeshXYPass, YmeshXYPass, ZmeshXYPass, rstride=1, cstride=1,
                    facecolors=fcolor_z_i, shade=False, antialiased=False)

#******


    print('Adding the colorbar, Time since start = {:.2f} s'.format(time.time()-t0))

    # Split existing axes to make room for colorbar
#    divider = make_axes_locatable(ax1)                         # There seems to be issues with make_axes_locatable with 3d axes
#    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax = fig.add_axes([0.76,0.2,0.03,0.6])                     # TODO find a cleaner way to deal with this
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

    print('Drawing the colorbar, Time since start = {:.2f} s'.format(time.time()-t0))

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
#    if (symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2)) or (not lin and symlog is None):
#        print('Colorbar intermediate draw, Time since start = {:.2f} s'.format(time.time()-t0))
#        fig.canvas.draw() # draw to get tick positions

    # Adjust placement of innermost ticks for symlog if it indeed is (quasi)symmetric
    if symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2):
        print('Tick adjustment, Time since start = {:.2f} s'.format(time.time()-t0))
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

    print('Tick precision, Time since start = {:.2f} s'.format(time.time()-t0))

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

    print('Tick pruning, Time since start = {:.2f} s'.format(time.time()-t0))

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

    print('Watermark test, Time since start = {:.2f} s'.format(time.time()-t0))

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

    print('Calling tight_layout(), Time since start = {:.2f} s'.format(time.time()-t0))
    
    # Adjust layout. Uses tight_layout() but in fact this ensures 
    # that long titles and tick labels are still within the plot area.
#    plt.tight_layout(pad=0.01)   # TODO check: a warning says tight_layout() might not be compatible with those axes. Seems to work though...
#    ax1.axis('equal')
    savefig_pad=0.01
    bbox_inches='tight'

    print('Saving the figure, Time since start = {:.2f} s'.format(time.time()-t0))

    # Save output to file
    try:
        plt.savefig(outputfile,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)
    except:
        print("Error with attempting to save figure.")
    print(outputfile+"\n")

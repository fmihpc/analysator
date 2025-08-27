import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os, sys
import time
import re
import math
import matplotlib.ticker as mtick
import matplotlib.path as mpath
import matplotlib.patches as patches
import matplotlib.projections as projections
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.ticker import LogLocator
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.cbook import get_sample_data
from distutils.version import LooseVersion, StrictVersion

def scaling1d(data, resolution, limits=None):   
    ''' Transfroms data to indices between [0, resolution-1]
    
    :kword data:       1D data to be plotted using bivariate colormaps
    :kword resolution: Number of unique colors per direction
    :kword limits:     Tuple that overrides the colorspace limits from data limits
    
    Returns a list of integer indices.
    '''
    if len(data)==0:
        return np.array([])
    
    if limits is not None:
        vmin, vmax = limits
        data[np.where(data>vmax)] = vmax
        data[np.where(data<vmin)] = vmin
        comparison = data[:, np.newaxis] >= np.linspace(vmin, vmax, resolution, endpoint=False)
        indices = comparison.shape[1] - 1 - np.argmax(comparison[:, ::-1], axis=1)   
        
        return indices   
    
    comparison = data[:, np.newaxis] >= np.linspace(np.amin(data), np.amax(data), resolution, endpoint=False)
    indices = comparison.shape[1] - 1 - np.argmax(comparison[:, ::-1], axis=1)
                                
    return indices
    

def bivariate(resolution_x, resolution_y, a = 0.5, order = [1,0,2], green_multiplier = 1, white_center=True):   
    '''Creates the RGB space for the 2D colormap
    
    :kword resolution:       Number of unique colors per direction
    :kword a:                Controls the slope of the two off-diagonal surfaces, between 0 and 1
    :kword order:            Changes which color occupies which corner of the 2D map
    :kword green_multiplier: Determines the amount of green in the colormap
    :kword white_center:     True to force the colormap to go through white
    
    Returns a (resolution, resolution) array of RGB values.
    '''
    
#    l = resolution
    rgb = np.zeros((resolution_x,resolution_y,3))
    lspace_x = np.linspace(0,1.0,resolution_x)
    lspace_y = np.linspace(0,1.0,resolution_y)

    rgb[:,:,order[0]] += 0.5*lspace_x[:,np.newaxis]+0.5*lspace_y[np.newaxis,:]

    rgb[:,:,order[1]] += a*lspace_x[:,np.newaxis]+(1-a)*lspace_y[np.newaxis,::-1]
    
    rgb[:,:,order[2]] += (1-a)*lspace_x[::-1,np.newaxis]+a*lspace_y[np.newaxis,:]

    rgb[:,:,1] *= green_multiplier

    if white_center:
        lim = 1/np.sqrt(2)
        lspace_second_x = np.linspace(-lim,lim,resolution_x)**2
        lspace_second_y = np.linspace(-lim,lim,resolution_y)**2
        r = np.sqrt(lspace_second_x[:,np.newaxis] + lspace_second_y[np.newaxis,:])
        rgb[:,:,:] += 0.5*(1-r)[:,:,np.newaxis]
        
    rgb[rgb<0] = 0
    rgb[rgb>1] = 1
    return rgb

def plot_2dcolormap(ax, xlimits, zlimits, resolution_x, resolution_y, a = 0.5, order = [1,0,2], green_multiplier = 1, white_center=True, fontsize=12, thick=1):

    rgbs = bivariate(resolution_x, resolution_y, a, order, green_multiplier, white_center)
    rgbmesh, rgbmesh2 = np.meshgrid(np.arange(resolution_x), np.arange(resolution_y))

    rgb_X, rgb_Y = np.meshgrid(np.linspace(xlimits[0], xlimits[1], resolution_x, endpoint=False)+ (xlimits[1]-xlimits[0])/(2*resolution_x) , np.linspace(zlimits[0], zlimits[1], resolution_y, endpoint=False) + (zlimits[1]-zlimits[0])/(2*resolution_y))
    
    rgb_context = rgbs[rgbmesh, rgbmesh2]
    ax.pcolormesh(rgb_X, rgb_Y, rgb_context)
    ax.set_xlabel(r"$X$ [$R_{\rm E}$]", fontsize=fontsize)
    ax.xaxis.set_label_position("top")
    ax.set_ylabel(r"$Z$ [$R_{\rm E}$]", fontsize=fontsize, labelpad=-40)
    ax.set_xlim(xlimits[0], xlimits[1])
    ax.set_ylim(zlimits[0], zlimits[1])
    ax.set_aspect('equal')
    ax.tick_params(labelsize=fontsize, width=thick, length=3*thick, bottom=False, top=True, labelbottom=False, labeltop=True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(thick)

def plot_ionosphere(filename=None,
                  vlsvobj=None,
                  filedir=None, step=None, run=None,
                  outputdir=None, outputfile=None,
                  nooverwrite=False,
                  var=None, op=None, operator=None,
                  colormap=None, vmin=None, vmax=None,
                  contourvar=None,
                  f2=None,
                  contourlevels=np.arange(0,10),
                  contourvmin=None,
                  contourvmax=None,
                  contourcmap=None,
                  contourlinewidths=1,
                  symmetric=False, absolute=None,
                  usesci=True, log=None,
                  lin=None, symlog=None, nocb=False, internalcb=False,
                  minlatitude=60,
                  cbtitle=None, title=None, cbaxes=None,
                  thick=1.0,scale=1.0,vscale=1.0,
                  wmark=False,wmarkb=False,
                  viewdir=1.0,draw=None,
                  axes=None
      ):

    ''' Plots a 2d projection of an ionosphere variable

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
    :kword absolute:    Plot the absolute of the evaluated variable
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: True)
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2, but this can
                        result in the innermost tick marks overlapping. In this case, using a larger value for 
                        symlog is suggested.
    :kword nocb:        Set to suppress drawing of colourbar
    :kword internalcb:  Set to draw colorbar inside plot instead of outside. If set to a text
                        string, tries to use that as the location, e.g. "NW","NE","SW","SW"
    :kword minlatitude: Minimum plot latitude (default=60 degrees)
    :kword title:       string to use as plot title instead of time.
                        Special case: Set to "msec" to plot time with millisecond accuracy or "musec"
                        for microsecond accuracy. "sec" is integer second accuracy.
    :kword cbtitle:     string to use as colorbar title instead of map name
    :kword viewdir:     view direction onto the sphere. Positive value: North pole. Negative values: South pole.
    :kword scale:       Scale text size (default=1.0)
    :kword vscale:      Scale all values with this before plotting. Useful for going from e.g. m^-3 to cm^-3
                        or from tesla to nanotesla. Guesses correct units for colourbar for some known
                        variables.
    :kword axes:        Provide the routine a set of axes to draw within instead of generating a new image.
    :kword cbaxes:      Provide the routine a set of axes for the colourbar.
    :kword thick:       line and axis thickness, default=1.0
    :kword draw:        Set to anything but None or False in order to draw image on-screen instead of saving to file (requires x-windowing)
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner. If set to a text
                        string, tries to use that as the location, e.g. "NW","NE","SW","SW"
    :kword wmarkb:      As for wmark, but uses an all-black Vlasiator logo.

 
    .. code-block:: python

    '''

    IONOSPHERE_RADIUS=6471e3 # R_E + 100 km
    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    watermarkimageblack=os.path.join(os.path.dirname(__file__), 'logo_black.png')


    # Change certain falsy values:
    if not lin and lin != 0:
        lin = None
    if not (log is None):
        lin = None
    if not symlog and symlog != 0:
        symlog = None
    if symlog is True:
        symlog = 0
    if (filedir == ''):
        filedir = './'
    if (outputdir == ''):
        outputdir = './'

    # Input file or object
    if filename!=None:
        f=pt.vlsvfile.VlsvReader(filename)
    elif ((filedir!=None) and (step!=None)):
        filename = glob.glob(filedir+'bulk*'+str(step).rjust(7,'0')+'.vlsv')[0]
        f=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        f=vlsvobj
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    if operator is None:
        if op is not None:
            operator=op

    if not colormap:
        # Default values
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
            plot_title = "t="+timeformat.format(timeval)+' s '
    else:
        plot_title = ''+title

    if viewdir > 0:
        plot_title = "North" + plot_title
    else:
        plot_title = "South" + plot_title

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
    if operator is not None:
        # .isdigit checks if the operator is an integer (for taking an element from a vector)
        if type(operator) is int:
            operator = str(operator)
        if not operator in 'xyz' and operator != 'magnitude' and not operator.isdigit():
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
    if not var:
        # If no expression or variable given, defaults to FACs
        var='ig_fac'
    varstr=var.replace("/","_")

    if not outputfile: # Generate filename
        if not outputdir: # default initial path
            outputdir=pt.plot.defaultoutputdir
        # Sub-directories can still be defined in the "run" variable
        if viewdir > 0:
            outputpole = "_north"
        else:
            outputpole = "_south"
        outputfile = outputdir+run+"_ionosphere_"+varstr+operatorfilestr+outputpole+stepstr+".png"
    else: 
        if outputdir:
            outputfile = outputdir+outputfile

    # Re-check to find actual target sub-directory
    outputprefixind = outputfile.rfind('/')
    if outputprefixind >= 0:            
        outputdir = outputfile[:outputprefixind+1]

    # Ensure output directory exists
    if axes is None and not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except:
            pass

    if axes is None and not os.access(outputdir, os.W_OK):
        print(("No write access for directory "+outputdir+"! Exiting."))
        return

    # Check if target file already exists and overwriting is disabled
    if axes is None and (nooverwrite and os.path.exists(outputfile)):            
        if os.stat(outputfile).st_size > 0: # Also check that file is not empty
            print(("Found existing file "+outputfile+". Skipping."))
            return
        else:
            print(("Found existing file "+outputfile+" of size zero. Re-rendering."))

    # The plot will be saved in a new figure 
    if axes is None and str(matplotlib.get_backend()) != pt.backend_noninteractive: #'Agg':
        plt.switch_backend(pt.backend_noninteractive)

    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    # If axes are provided, leave backend as-is.
    if axes is None and draw is not None:
        if str(matplotlib.get_backend()) != pt.backend_interactive: #'TkAgg': 
            plt.switch_backend(pt.backend_interactive)
    elif axes is None:
        if str(matplotlib.get_backend()) != pt.backend_noninteractive: #'Agg':
            plt.switch_backend(pt.backend_noninteractive)  

    # Read ionosphere mesh node coordinates
    coords = f.get_ionosphere_node_coords()
    # Read ionosphere mesh connectivity
    elements = f.get_ionosphere_element_corners()
    # Read selected variable values
    if operator is None:
       operator="pass"
    datamap_info= f.read_variable_info(var, operator=operator)
    contourmap_info = f.read_variable_info(contourvar)
    if f2 is not None:
        contourmap_info2 = f2.read_variable_info(contourvar)
    # Correction for broken ionosphere units
    datamap_info.latex = re.sub("\\\\text","\\\\mathrm", datamap_info.latex)
    datamap_info.latexunits = re.sub("\\\\mho","\\\\Omega^{-1}", datamap_info.latexunits)
    cb_title_use = datamap_info.latex
    # Check if vscale results in standard unit
    vscale, _, datamap_unit_latex = datamap_info.get_scaled_units(vscale=vscale)
    values = datamap_info.data*vscale
    if np.ndim(values) == 0:
        print("Error, reading variable '" + str(var) + "' from vlsv file!",values.shape)
        return -1
    contourvalues = contourmap_info.data
    #contourvalues2 = contourmap_info2.data

    # Add unit to colorbar title
    if datamap_unit_latex:
       cb_title_use = cb_title_use + "\,["+datamap_unit_latex+"]"

    if viewdir > 0:
        r=np.degrees(np.arccos(coords[:,2]/IONOSPHERE_RADIUS))
    else:
        r=np.degrees(np.arccos(-coords[:,2]/IONOSPHERE_RADIUS))
    theta=np.arctan2(coords[:,1],coords[:,0])-math.pi/2

    # Project nodes and elements into view plane
    mask = []
    if viewdir > 0:
      for e in elements:
         if coords[e[0],2] > 0 and coords[e[1],2] > 0 and coords[e[2],2] > 0:
            mask+=[False]
         else:
            mask+=[True]
    else:
      for e in elements:
         if coords[e[0],2] < 0 and coords[e[1],2] < 0 and coords[e[2],2] < 0:
            mask+=[False]
         else:
            mask+=[True]

    # Build mesh triangulation
    tri = matplotlib.tri.Triangulation(-r*np.sin(theta), r*np.cos(theta), elements, mask)

    # Allow title override
    if cbtitle is not None:
        # Here allow underscores for manual math mode
        cb_title_use = cbtitle

    # If automatic range finding is required, find min and max of array
    if vmin is not None:
        vminuse=vmin
    else: 
        vminuse = np.ma.amin(values)
    if vmax is not None:
        vmaxuse=vmax
    else:
        vmaxuse = np.ma.amax(values)

    # If automatic range finding is required, find min and max of array
    if contourvmin is not None:
        contourvminuse=contourvmin
    else: 
        contourvminuse = np.ma.amin(contourvalues)
    if contourvmax is not None:
        contourvmaxuse=contourvmax
    else:
        contourvmaxuse = np.ma.amax(contourvalues)

    # If both values are zero, we have an empty array
    if vmaxuse==vminuse==0:
        print("Error, requested array is zero everywhere. Exiting.")
        return 0
    # If both values are zero, we have an empty array
    if contourvmaxuse==contourvminuse==0:
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
        vminuse = min(np.ma.amin(values), 0)

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

    # Lin or log colour scaling, defaults to lin
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
            ticks = LogLocator(base=10,subs=list(range(10))) # where to show labels
    else:
        # Linear
        linticks = 7
        if isinstance(lin, int):
            linticks = abs(lin)
            if linticks==1: # old default was to set lin=1 for seven linear ticks
                linticks = 7
                
        levels = MaxNLocator(nbins=255).tick_values(vminuse,vmaxuse)
        norm = BoundaryNorm(levels, ncolors=cmapuse.N, clip=True)
        ticks = np.linspace(vminuse,vmaxuse,num=linticks)

    # Creating a new figure and axes 
    figsize = (6,5)
    if nocb:
        figsize = (5,5)

    if axes is None:
        fig = plt.figure(figsize=figsize,dpi=150)
        ax_cartesian = fig.add_axes([0.1,0.1,0.9,0.9], xlim=(-(90-minlatitude),(90-minlatitude)), ylim=(-(90-minlatitude),(90-minlatitude)), aspect='equal')
        #ax_polar = fig.add_axes([0.1,0.1,0.9,0.9], polar=True, frameon=False, ylim=(0, 90-minlatitude))
        ax_polar = inset_axes(parent_axes=ax_cartesian, width="100%", height="100%", axes_class = projections.get_projection_class('polar'), borderpad=0)
        ax_polar.set_frame_on(False)
        ax_polar.set_aspect('equal')
    else:
        axes.set_xticklabels([])
        axes.set_yticklabels([])
        axes.axis('off')
        ax_cartesian = inset_axes(parent_axes=axes, width="95%", height="95%", borderpad=1, loc='center left')
        ax_cartesian.set_xlim(-(90-minlatitude),(90-minlatitude))
        ax_cartesian.set_ylim(-(90-minlatitude),(90-minlatitude))
        ax_cartesian.set_aspect('equal')
        ax_cartesian.set_frame_on(False)
        ax_polar = inset_axes(parent_axes=ax_cartesian, width="100%", height="100%", axes_class = projections.get_projection_class('polar'), borderpad=0)
        ax_polar.set_frame_on(False)
        ax_polar.set_aspect('equal')
    ax_cartesian.set_xticklabels([])
    ax_cartesian.set_yticklabels([])
    ax_cartesian.axis('off')


    ## Build a circle to map away the regions we're not interested in
    def make_circle(r):
        t = np.arange(0, np.pi * 2.0, 0.01)
        t = t.reshape((len(t), 1))
        x = r * np.cos(t)
        y = r * np.sin(t)
        return np.hstack((x, y))
    path = mpath.Path
    inside_vertices = make_circle(90-minlatitude)
    outside_vertices = make_circle(1000)
    pathCodes = np.ones(len(inside_vertices), dtype=mpath.Path.code_type) * mpath.Path.LINETO
    pathCodes[0] = mpath.Path.MOVETO
    path = mpath.Path(np.concatenate((inside_vertices[::-1],outside_vertices)), np.concatenate((pathCodes,pathCodes)))
    clippingcircle= patches.PathPatch(path, facecolor="#ffffff00", edgecolor='grey', linewidth=thick)

    Re = 6371000

    ### read field line mapping and process it
    fw_feet = np.ma.masked_array(f.read_variable("vg_connection_coordinates_fw"))
    bw_feet = np.ma.masked_array(f.read_variable("vg_connection_coordinates_bw"))
    vg_fluxrope = f.read_variable("vg_fluxrope")
#    vg_curvature = np.ma.masked_array(f.read_variable("vg_curvature", operator="magnitude"))
    cellids = np.ma.masked_array(f.read_variable("CellID"))
    fw_feet[:,0] = np.ma.masked_where(np.linalg.norm(fw_feet, axis=-1) > 6.5e6, fw_feet[:,0])
    bw_feet[:,0] = np.ma.masked_where(np.linalg.norm(bw_feet, axis=-1) > 6.5e6, bw_feet[:,0])
    fw_feet[:,0] = np.ma.masked_where(vg_fluxrope <= 0.01, fw_feet[:,0])
    bw_feet[:,0] = np.ma.masked_where(vg_fluxrope <= 0.01, bw_feet[:,0])
    fw_feet[:,0] = np.ma.masked_where(vg_fluxrope > 7, fw_feet[:,0])
    bw_feet[:,0] = np.ma.masked_where(vg_fluxrope > 7, bw_feet[:,0])
#    fw_feet[:,0] = np.ma.masked_where(vg_curvature < 1.0 / (20.0*Re), fw_feet[:,0])
#    bw_feet[:,0] = np.ma.masked_where(vg_curvature < 1.0 / (20.0*Re), bw_feet[:,0])

    if viewdir > 0: # North
        fw_feet[:,2] = np.ma.masked_where(fw_feet[:,2] < 0, fw_feet[:,2])
        bw_feet[:,2] = np.ma.masked_where(bw_feet[:,2] < 0, bw_feet[:,2])
    else: # South
        fw_feet[:,2] = np.ma.masked_where(fw_feet[:,2] > 0, fw_feet[:,2])
        bw_feet[:,2] = np.ma.masked_where(bw_feet[:,2] > 0, bw_feet[:,2])

    fw_feet_in=np.ma.mask_rows(fw_feet)
    bw_feet_in=np.ma.mask_rows(bw_feet)
    fw_feet_out=np.ma.mask_rows(fw_feet)
    bw_feet_out=np.ma.mask_rows(bw_feet)

    fw_hits_x = np.ma.masked_where(fw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,0])
    bw_hits_x = np.ma.masked_where(bw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,0])
    fw_hits_y = np.ma.masked_where(fw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,1])
    bw_hits_y = np.ma.masked_where(bw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,1])
    fw_hits_z = np.ma.masked_where(fw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,2])
    bw_hits_z = np.ma.masked_where(bw_feet.mask[:,0], f.get_cell_coordinates(cellids)[:,2])

    fw_feet_in[:,0] = np.ma.masked_where(np.abs(fw_hits_y) > 4.5*Re, fw_feet_in[:,0])
    bw_feet_in[:,0] = np.ma.masked_where(np.abs(bw_hits_y) > 4.5*Re, bw_feet_in[:,0])
    fw_feet_out[:,0] = np.ma.masked_where(np.abs(fw_hits_y) <= 4.5*Re, fw_feet_out[:,0])
    bw_feet_out[:,0] = np.ma.masked_where(np.abs(bw_hits_y) <= 4.5*Re, bw_feet_out[:,0])

    fw_feet_in=np.ma.mask_rows(fw_feet_in)
    bw_feet_in=np.ma.mask_rows(bw_feet_in)
    fw_feet_out=np.ma.mask_rows(fw_feet_out)
    bw_feet_out=np.ma.mask_rows(bw_feet_out)

    fw_hits_in_x = np.ma.masked_where(fw_feet_in.mask[:,0], fw_hits_x)
    bw_hits_in_x = np.ma.masked_where(bw_feet_in.mask[:,0], bw_hits_x)
    fw_hits_out_x = np.ma.masked_where(fw_feet_out.mask[:,0], fw_hits_x)
    bw_hits_out_x = np.ma.masked_where(bw_feet_out.mask[:,0], bw_hits_x)

    fw_hits_in_z = np.ma.masked_where(fw_feet_in.mask[:,0], fw_hits_z)
    bw_hits_in_z = np.ma.masked_where(bw_feet_in.mask[:,0], bw_hits_z)
    fw_hits_out_z = np.ma.masked_where(fw_feet_out.mask[:,0], fw_hits_z)
    bw_hits_out_z = np.ma.masked_where(bw_feet_out.mask[:,0], bw_hits_z)

    fw_feet_in=fw_feet_in.compressed().reshape([-1,3])
    bw_feet_in=bw_feet_in.compressed().reshape([-1,3])
    fw_feet_out=fw_feet_out.compressed().reshape([-1,3])
    bw_feet_out=bw_feet_out.compressed().reshape([-1,3])

    fw_hits_in_x = fw_hits_in_x.compressed()
    bw_hits_in_x = bw_hits_in_x.compressed()
    fw_hits_out_x = fw_hits_out_x.compressed()
    bw_hits_out_x = bw_hits_out_x.compressed()

    fw_hits_in_z = fw_hits_in_z.compressed()
    bw_hits_in_z = bw_hits_in_z.compressed()
    fw_hits_out_z = fw_hits_out_z.compressed()
    bw_hits_out_z = bw_hits_out_z.compressed()

    if viewdir > 0:
        try:
            fw_feet_in_r=np.degrees(np.arccos(fw_feet_in[:,2]/IONOSPHERE_RADIUS))
            fw_feet_out_r=np.degrees(np.arccos(fw_feet_out[:,2]/IONOSPHERE_RADIUS))
        except:
            print("no fw feet North")
        try:
            bw_feet_in_r=np.degrees(np.arccos(bw_feet_in[:,2]/IONOSPHERE_RADIUS))
            bw_feet_out_r=np.degrees(np.arccos(bw_feet_out[:,2]/IONOSPHERE_RADIUS))
        except:
            print("no bw feet North")
    else:
        try:
            fw_feet_in_r=np.degrees(np.arccos(-fw_feet_in[:,2]/IONOSPHERE_RADIUS))
            fw_feet_out_r=np.degrees(np.arccos(-fw_feet_out[:,2]/IONOSPHERE_RADIUS))
        except:
            print("no fw feet South")
        try:
            bw_feet_in_r=np.degrees(np.arccos(-bw_feet_in[:,2]/IONOSPHERE_RADIUS))
            bw_feet_out_r=np.degrees(np.arccos(-bw_feet_out[:,2]/IONOSPHERE_RADIUS))
        except:
            print("no bw feet South")
    try:
        fw_feet_in_theta=np.arctan2(fw_feet_in[:,1],fw_feet_in[:,0]) + math.pi
        fw_feet_out_theta=np.arctan2(fw_feet_out[:,1],fw_feet_out[:,0]) + math.pi
    except:
        print("no fw feet")
    try:
        bw_feet_in_theta=np.arctan2(bw_feet_in[:,1],bw_feet_in[:,0]) + math.pi
        bw_feet_out_theta=np.arctan2(bw_feet_out[:,1],bw_feet_out[:,0]) + math.pi
    except:
        print("no bw feet")

    ### THE ACTUAL PLOT HAPPENS HERE ###
    # new tweaks
    resolution_x = 8    
    resolution_y = 7
    a = 0.8
    order = [1,0,2]
    green_multiplier = 1
    white_center = False
    rgb = bivariate(resolution_x, resolution_y, a, order, green_multiplier, white_center)
            
    zlimits = (-10,10) # Z colormap limits
    xlimits = (-30,10)   # X colormap limits

    if cbaxes is not None:
        plot_2dcolormap(cbaxes, xlimits, zlimits, resolution_x, resolution_y, a=a, white_center=white_center, green_multiplier=green_multiplier, order=order, fontsize=fontsize3, thick=thick)
    
    try: 
        scaled_fw_hits_in_z = scaling1d(fw_hits_in_z/Re, resolution_y, limits=zlimits)
        scaled_fw_hits_in_x = scaling1d(fw_hits_in_x/Re, resolution_x, limits=xlimits)
        scaled_fw_hits_out_z = scaling1d(fw_hits_out_z/Re, resolution_y, limits=zlimits)
        scaled_fw_hits_out_x = scaling1d(fw_hits_out_x/Re, resolution_x, limits=xlimits)


        if len(fw_feet_in_theta) > 0 and len(fw_feet_in_r) > 0: # If data exists
            scat = ax_polar.scatter(fw_feet_in_theta, fw_feet_in_r, c=rgb[scaled_fw_hits_in_x, scaled_fw_hits_in_z], s=2000, marker=1, linewidths=6)
        else:
            scat = ax_polar.scatter(fw_feet_in_theta, fw_feet_in_r) # If not, return an empty version of the same object

        if len(fw_feet_out_theta > 0) and len(fw_feet_out_r) > 0: # If data exists
            scat = ax_polar.scatter(fw_feet_out_theta, fw_feet_out_r, c=rgb[scaled_fw_hits_out_x, scaled_fw_hits_out_z], s=2000, marker=2, linewidths=6)
        else:
            scat = ax_polar.scatter(fw_feet_out_theta, fw_feet_out_r) # If not, return an empty version of the same object
    except Exception as e:
        print(e)
        
    try: # The same as above
        scaled_bw_hits_in_z = scaling1d(bw_hits_in_z/Re, resolution_y, limits=zlimits)
        scaled_bw_hits_in_x = scaling1d(bw_hits_in_x/Re, resolution_x, limits=xlimits)
        scaled_bw_hits_out_z = scaling1d(bw_hits_out_z/Re, resolution_y, limits=zlimits)
        scaled_bw_hits_out_x = scaling1d(bw_hits_out_x/Re, resolution_x, limits=xlimits)

        if len(bw_feet_in_theta) > 0 and len(bw_feet_in_r) > 0:
            scat = ax_polar.scatter(bw_feet_in_theta, bw_feet_in_r, c=rgb[scaled_bw_hits_in_x, scaled_bw_hits_in_z], s=2000, marker=1, linewidths=6)
        else: scat = ax_polar.scatter(bw_feet_in_theta, bw_feet_in_r)

        if len(bw_feet_out_theta > 0) and len(bw_feet_out_r) > 0:
            scat = ax_polar.scatter(bw_feet_out_theta, bw_feet_out_r, c=rgb[scaled_bw_hits_out_x, scaled_bw_hits_out_z], s=2000, marker=2, linewidths=6)
        else: scat = ax_polar.scatter(bw_feet_out_theta, bw_feet_out_r)
    except Exception as e:
        print(e)


    #colormap = ax_cartesian.tripcolor(tri, values, cmap=cmapuse, norm=norm)# YPK , shading='gouraud')
    colormap = scat
#    contours = ax_cartesian.tricontour(tri, contourvalues, cmap=contourcmap, levels=contourlevels, vmin=contourvminuse, vmax=contourvmaxuse, linewidths=contourlinewidths)
    contours = ax_cartesian.tricontour(tri, contourvalues, colors=["black"], levels=contourlevels, vmin=contourvminuse, vmax=contourvmaxuse, linewidths=contourlinewidths)
    #contours2 = ax_cartesian.tricontour(tri, contourvalues2, cmap="brg", levels=contourlevels, vmin=contourvminuse, vmax=contourvmaxuse, linewidths=contourlinewidths-2)


    ax_cartesian.add_patch(clippingcircle)

    # Draw polar grid over it
    ax_polar.grid(True, color="black", linestyle='dashdot', linewidth=thick)
    gridlatitudes = np.arange(0., 90.-minlatitude,10.)
    ax_polar.set_rmax(90.-minlatitude);
#    ax_polar.set_rgrids(gridlatitudes, map(lambda x: "", gridlatitudes),angle=225)
    ax_polar.set_rgrids(gridlatitudes, map(lambda x: str(90.-x)+"°", gridlatitudes),angle=30)
    ax_polar.set_thetagrids(np.linspace(0., 360, 13), ["24","2","4","6","8","10","12","14","16","18","20","22","24"])
    ax_polar.set_theta_zero_location('W', offset=0)
    ax_polar.tick_params(labelsize=fontsize2, pad=0.1)

    # Title and plot limits
    if len(plot_title)!=0:# and axes is None:
#        plot_title = pt.plot.mathmode(pt.plot.bfstring(plot_title))
        plot_title = plot_title
        ax_polar.set_title(plot_title,fontsize=fontsize2,fontweight='bold',fontstyle='normal',x=-0.02,y=1, ha="left", va="top")

    # Colourbar title
    if len(cb_title_use)!=0:
        cb_title_use = pt.plot.mathmode(pt.plot.bfstring(cb_title_use))

    # Set flag which affects colorbar decimal precision
    if lin is None:
        pt.plot.cb_linear = False
    else:
        pt.plot.cb_linear = True

    # Creating colorbar axes
    '''
    if not nocb:
        if cbaxes: 
            # Colorbar axes are provided
            cax = cbaxes
            cbdir="right"; horalign="left"
        elif internalcb:
            # Colorbar within plot area
            cbloc=1; cbdir="left"; horalign="right"
            if type(internalcb) is str:
                if internalcb=="NE":
                    cbloc=1; cbdir="left"; horalign="right"
                if internalcb=="NW":
                    cbloc=2; cbdir="right"; horalign="left"
                if internalcb=="SW": 
                    cbloc=3; cbdir="right"; horalign="left"
                if internalcb=="SE": 
                    cbloc=4; cbdir="left";  horalign="right"
            # borderpad default value is 0.5, need to increase it to make room for colorbar title
            cax = inset_axes(ax_polar, width="5%", height="35%", loc=cbloc, borderpad=1.0,
                             bbox_transform=ax_polar.transAxes, bbox_to_anchor=(0.15,0,0.85,0.92))
        else:
            # Split existing axes to make room for colorbar
            if axes is None:
                cax = fig.add_axes([0.9,0.2,0.03,0.6])
            else:
                cax = axes.inset_axes([0.9,0.2,0.03,0.6])
            cbdir="right"; horalign="left"

        # Set flag which affects colorbar decimal precision
        if lin is None:
            pt.plot.cb_linear = False
        else:
            pt.plot.cb_linear = True

        # First draw colorbar
        if usesci:
            cb = plt.colorbar(colormap, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmtsci), cax=cax, drawedges=False)
        else:
            #cb = plt.colorbar(colormap, ticks=ticks, format=mtick.FormatStrFormatter('%4.2f'), cax=cax, drawedges=False)
            cb = plt.colorbar(colormap, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmt), cax=cax, drawedges=False)
        cb.outline.set_linewidth(thick)
        cb.ax.yaxis.set_ticks_position(cbdir)

        if not cbaxes:
            cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
            cb_title = cax.set_title(cb_title_use,fontsize=fontsize3,fontweight='bold', horizontalalignment=horalign)
            cb_title.set_position((0.,1.+0.025*scale)) # avoids having colourbar title too low when fontsize is increased
        else:
            cb.ax.tick_params(labelsize=fontsize)
            cb_title = cax.set_title(cb_title_use,fontsize=fontsize,fontweight='bold', horizontalalignment=horalign)

        # Perform intermediate draw if necessary to gain access to ticks
        if (axes is None) and (symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2)) or (not lin and symlog is None):
            fig.canvas.draw() # draw to get tick positions

        # Adjust placement of innermost ticks for symlog if it indeed is (quasi)symmetric
        if symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2) and len(cb.ax.yaxis.get_ticklabels()) > 2:
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
            nlabels = len(cb.ax.yaxis.get_ticklabels())
            # Force less ticks for internal colorbars
            if internalcb: 
                nlabels = nlabels * 1.5
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
    '''

    # Add Vlasiator watermark
    if (wmark or wmarkb) and (axes is None):
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
    if axes is None:
        plt.tight_layout(pad=0.01)
    savefig_pad=0.01
    bbox_inches='tight'

    # Save output or draw on-screen
    if not draw and axes is None:
        try:
            plt.savefig(outputfile,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)
        except:
            print("Error attempting to save figure: ", sys.exc_info())
        print('...Done!') 
    elif draw is not None and axes is None:
        # Draw on-screen
        plt.draw()
        plt.show()
        print('Draw complete!')

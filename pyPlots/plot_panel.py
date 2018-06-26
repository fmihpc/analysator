import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
import matplotlib.ticker as mtick
import colormaps as cmaps
from matplotlib.cbook import get_sample_data

# TODO: flag for plotting the panel edge-to-edge of the figure.

# Run TeX typesetting through the full TeX engine instead of python's own mathtext. Allows
# for changing fonts, bold math symbols etc, but may cause trouble on some systems.
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.rcParams['text.dvipnghack'] = 'True' # This hack might fix it on some systems

# Register custom colourmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='viridis_r', cmap=matplotlib.colors.ListedColormap(cmaps.viridis.colors[::-1]))
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='plasma_r', cmap=matplotlib.colors.ListedColormap(cmaps.plasma.colors[::-1]))
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='inferno_r', cmap=matplotlib.colors.ListedColormap(cmaps.inferno.colors[::-1]))
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.register_cmap(name='magma_r', cmap=matplotlib.colors.ListedColormap(cmaps.magma.colors[::-1]))
plt.register_cmap(name='parula', cmap=cmaps.parula)
plt.register_cmap(name='parula_r', cmap=matplotlib.colors.ListedColormap(cmaps.parula.colors[::-1]))
# plt.register_cmap(name='cork',cmap=cork_map)
# plt.register_cmap(name='davos_r',cmap=davos_r_map)
plt.register_cmap(name='hot_desaturated', cmap=cmaps.hot_desaturated_colormap)
plt.register_cmap(name='hot_desaturated_r', cmap=cmaps.hot_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step
plt.register_cmap(name='pale_desaturated', cmap=cmaps.pale_desaturated_colormap)
plt.register_cmap(name='pale_desaturated_r', cmap=cmaps.pale_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step

# Different style scientific format for colour bar ticks
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)

def plot_colormap(filename=None,
                  vlsvobj=None,
                  filedir=None, step=None,
                  outputdir=None, nooverwrite=None,
                  var=None, op=None,
                  title=None, cbtitle=None, draw=None, usesci=None,
                  symlog=None,
                  boxm=[],boxre=[],colormap=None,
                  run=None,wmark=None, nocb=None,
                  unit=None, thick=1.0,scale=1.0,
                  noborder=None, noxlabels=None, noylabels=None,
                  vmin=None, vmax=None, lin=None,
                  external=None, expression=None, 
                  #exprvals=None, #extvals=None, #expr_timeavg=None, # These were consolidated
                  pass_vars=None, pass_times=None,
                  fluxfile=None, fluxdir=None,
                  fluxthick=1.0, fluxlines=1
                  ):

    ''' Plots a coloured plot with axes and a colour bar.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
    :kword nooverwrite: Set to only perform actions if the target output file does not yet exist                    
     
    :kword var:         variable to plot, e.g. rho, rhoBeam, beta, temperature, MA, Mms, va, vms,
                        E, B, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc.
    :kword op:          Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude. 
           
    :kword boxm:        zoom box extents [x0,x1,y0,y1] in metres (default and truncate to: whole simulation box)
    :kword boxre:       zoom box extents [x0,x1,y0,y1] in Earth radii (default and truncate to: whole simulation box)
    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, parula, nipy_spectral, RdBu, bwr
    :kword run:         run identifier, used for constructing output filename
    :kword title:       string to use as plot title instead of time
    :kword cbtitle:     string to use as colorbar title instead of map name
    :kword unit:        Plot axes using 10^{unit} m (default: Earth radius R_E)

    :kwird usesci:      Use scientific notation for colorbar ticks? (default: 1)
    :kword vmin,vmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword draw:        Set to nonzero in order to draw image on-screen instead of saving to file (requires x-windowing)

    :kword noborder:    Plot figure edge-to-edge without borders (default off)
    :kword noxlabels:   Suppress x-axis labels and title
    :kword noylabels:   Suppress y-axis labels and title
    :kword scale:       Scale text size (default=1.0)
    :kword thick:       line and axis thickness, default=1.0
   
    :kword external:    Optional function to use for external plotting of e.g. contours. The function
                        receives the following arguments: ax, XmeshXY,YmeshXY, pass_maps
    :kword expression:  Optional function which calculates a custom expression to plot. The function
                        receives the same list of numpy arrays as external, as an argument pass_maps,
                        the contents of which are maps of variables. Each is either of size [ysize,xsize]
                        or for multi-dimensional variables (vectors, tensors) it's [ysize,xsize,dim].
                        Remember to set vmin and vmax manually.

    :kword pass_vars:   Optional list of map names to pass to the external/expression functions 
                        as a list of numpy arrays. Each is either of size [ysize,xsize] or 
                        for multi-dimensional variables (vectors, tensors) it's [ysize,xsize,dim].
    :kword pass_times:  Integer, how many timesteps in each direction should be passed to external/expression
                        functions in pass_vars (e.g. pass_times=1 passes the values of three timesteps). 
                        This causes pass_vars to become a list of timesteps, with each timestep containing
                        a list of numpy arrays as for regular pass_vars.

    :kword fluxfile:    Filename to plot fluxfunction from
    :kword fluxdir:     Directory in which fluxfunction files can be found
    :kword fluxthick:   Scale fluxfunction line thickness
    :kword fluxlines:   Relative density of fluxfunction contours
                            
    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:
    plot_colormap(filename=fileLocation, var="MA", run="BCQ",
                  colormap='nipy_spectral',step=j, outputdir=outputLocation,
                  lin=1, wmark=1, vmin=2.7, vmax=10, 
                  external=cavitoncontours, extvals=['rho','B','beta'])
    # Where cavitoncontours is an external function which receives the arguments
    #  ax, XmeshXY,YmeshXY, extmaps
    # where extmaps is an array of maps for the requested variables.

    # example (simple) use of expressions:
    def exprMA_cust(exprmaps): #where exprmaps contains va, and the function returns the M_A with a preset velocity
        custombulkspeed=750000. # m/s
        va = exprmaps[0][:,:]
        MA = custombulkspeed/va
        return MA
    plot_colormap(filename=fileLocation, vmin=1 vmax=40,
                  expression=exprMA_cust, extvals=['va'],lin=1)

    '''

    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    # watermarkimage=os.path.expandvars('$HOME/appl_taito/analysator/pyPlot/logo_color.png')

    outputprefix = ''
    if outputdir==None:
        outputdir=os.path.expandvars('$HOME/Plots/')
    outputprefixind = outputdir.rfind('/')
    if outputprefixind >= 0:
        outputprefix = outputdir[outputprefixind+1:]
        outputdir = outputdir[:outputprefixind+1]
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    # Input file or object
    if filename!=None:
        f=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        f=vlsvobj
    elif ((filedir!=None) and (step!=None)):
        filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
        f=pt.vlsvfile.VlsvReader(filename)
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    # Flux function files
    if fluxdir!=None:
        if step != None:
            fluxfile = fluxdir+'flux.'+str(step).rjust(7,'0')+'.bin'
            if not os.path.exists(fluxfile):
                fluxfile = fluxdir+'bulk.'+str(step).rjust(7,'0')+'.bin'
        else:            
            fluxfile = fluxdir+'flux.'+filename[-12:-5]+'.bin'
            if not os.path.exists(fluxfile):
                fluxfile = fluxdir+'bulk.'+filename[-12:-5]+'.bin'

    if fluxfile!=None:
        if not os.path.exists(fluxfile):
            print("Error locating flux function file!")
            fluxfile=None
                
    # Scientific notation for colorbar ticks?
    if usesci==None:
        usesci=1
    
    if colormap==None:
        # Default values
        colormap="hot_desaturated"
        if op!=None:
            colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=5*scale # Colour bar ticks

    # Plot title with time
    timeval=None
    timeval=f.read_parameter("time")
    if timeval==None:
        timeval=f.read_parameter("t")
    if timeval==None:
        print "Unknown time format encountered"

    # Plot title with time
    if title==None:        
        if timeval == None:    
            print "Unknown time format encountered"
            plot_title = ''
        else:
            plot_title = "t="+str(np.int(timeval))+' s'
    else:
        plot_title = title

    # step, used for file name
    if step!=None:
        stepstr = '_'+str(step).rjust(7,'0')
    else:
        stepstr = ''

    # If run name isn't given, just put "plot" in the output file name
    if run==None:
        run='plot'

    # Verify validity of operator
    opstr=''
    if op!=None:
        if op!='x' and op!='y' and op!='z':
            print("Unknown operator "+op)
            op=None            
        else:
            # For components, always use linear scale, unless symlog is set
            opstr='_'+op
            if symlog==None:
                lin=1

    # Output file name
    if expression!=None:
        varstr=expression.__name__ 
    else:        
        if var==None:
            # If no expression or variable given, defaults to rho
            var='rho'
        varstr=var
    savefigname = outputdir+outputprefix+run+"_map_"+varstr+opstr+stepstr+".png"

    # Check if target file already exists and overwriting is disabled
    if (nooverwrite!=None and os.path.exists(savefigname)):
        return

    Re = 6.371e+6 # Earth radius in m
    #read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")
    # xsize = f.read_parameter("xcells_ini")
    # ysize = f.read_parameter("ycells_ini")
    # zsize = f.read_parameter("zcells_ini")
    # xmin = f.read_parameter("xmin")
    # xmax = f.read_parameter("xmax")
    # ymin = f.read_parameter("ymin")
    # ymax = f.read_parameter("ymax")
    # zmin = f.read_parameter("zmin")
    # zmax = f.read_parameter("zmax")

    # Check if ecliptic or polar run
    if ysize==1:
        simext=[xmin,xmax,zmin,zmax]
        sizes=[xsize,zsize]
    if zsize==1:
        simext=[xmin,xmax,ymin,ymax]
        sizes=[xsize,ysize]

    # Select window to draw
    if len(boxm)==4:
        boxcoords=boxm
    elif len(boxre)==4:
        boxcoords=[i*Re for i in boxre]
    else:
        boxcoords=simext

    # If box extents were provided manually, truncate to simulation extents
    boxcoords[0] = max(boxcoords[0],simext[0])
    boxcoords[1] = min(boxcoords[1],simext[1])
    boxcoords[2] = max(boxcoords[2],simext[2])
    boxcoords[3] = min(boxcoords[3],simext[3])

    # Axes and units (default R_E)
    if unit!=None: # Use m or km or other
        if unit==0:
            unitstr = r'm'
        if unit==3:
            unitstr = r'km'
        else:
            unitstr = r'$10^{'+str(int(unit))+'}$ m'
        unit = np.power(10,int(unit))
    else:
        unitstr = r'$\mathrm{R}_{\mathrm{E}}$'
        unit = Re
        
    # Scale data extent and plot box
    simext=[i/unit for i in simext]
    boxcoords=[i/unit for i in boxcoords]

    pass_maps=[]

    ##########
    # Read data and calculate required variables
    ##########
    if expression==None:
        if var == 'rho':
            cb_title = r"$n_\mathrm{p} [\mathrm{m}^{-3}]$"
            datamap = f.read_variable("rho")

        elif var == 'rhoBeam':
            cb_title = r"$n_{\mathrm{beam}} [\mathrm{m}^{-3}]$"
            datamap = f.read_variable("RhoBackstream")

        elif var == 'beta':
            cb_title = r"$\beta$"
            datamap = f.read_variable("beta")

        elif var == 'temperature':
            cb_title = r"$T$ [K]"
            datamap = f.read_variable("Temperature")

        elif var == 'MA':
            cb_title = r"$\mathrm{M}_\mathrm{A}$"
            Vmag = f.read_variable("v",operator='magnitude')
            va = f.read_variable("va")
            datamap = Vmag/va

        elif var == 'Mms':
            cb_title = r"$\mathrm{M}_\mathrm{ms}$"
            Vmag = f.read_variable("v",operator='magnitude')
            vms = f.read_variable("vms")
            datamap = Vmag/vms

        elif var == 'va':
            cb_title = r"$v_\mathrm{A}$"
            datamap = f.read_variable("va")

        elif var == 'vms':
            cb_title = r"$v_\mathrm{ms}$"
            datamap = f.read_variable("vms")

        elif var == 'B':
            if op==None:
                cb_title = r"$|B|$ [T]"
                datamap = f.read_variable("B",operator='magnitude')
            else:
                cb_title = r"$B_"+op+"$ [T]"
                datamap = f.read_variable("B",operator=op)
                # datamap = datamap*1e+9 # could be used to ouptut nanotesla instead of tesla

        elif var == 'E':
            if op==None:
                cb_title = r"$|E|$ [V/m]"
                datamap = f.read_variable("E",operator='magnitude')
            else:
                cb_title = r"$E_"+op+"$ [V/m]"
                datamap = f.read_variable("E",operator=op)

        elif var == 'V':
            if op==None:
                cb_title = r"$|V|\,[\mathrm{m}\,\mathrm{s}^{-1}]$"
                datamap = f.read_variable("V",operator='magnitude')
            else:
                cb_title = r"$V_"+op+"\,[\mathrm{m}\,\mathrm{s}^{-1}]$"
                datamap = f.read_variable("V",operator=op)
                # datamap = datamap*1e-3 # Plot this as km/s instead of m/s

        else:
            # Pipe all other vars directly to analysator
            if op==None:
                cb_title = var
                datamap = f.read_variable(var)
                # If value was vector value, take magnitude
                if np.ndim(datamap) != 1:
                    cb_title = r"$|"+var+"|$"
                    datamap = np.linalg.norm(np.asarray(datamap),axis=-1)
            else:
                cb_title = r" "+var+"$_"+op+"$"
                datamap = f.read_variable(var,operator=op)            
            
        if np.ndim(datamap)!=1:
            print("Error reading variable "+var+"! Exiting.")
            return -1

        # Reshape data to an ordered 2D array that can be plotted
        if np.ndim(datamap) != 2:
            datamap = datamap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        

    # If expression or external routine need variables, read them from the file.
    if pass_vars!=None:
        if pass_times==None:
            # Gather the required variable maps for a single time step
            for mapval in pass_vars:
                pass_map = f.read_variable(mapval)
                if np.ndim(pass_map)==1:
                    pass_map = pass_map[cellids.argsort()].reshape([sizes[1],sizes[0]])
                else:
                    pass_map = pass_map[cellids.argsort()].reshape([sizes[1],sizes[0],len(pass_map[0])])
                pass_maps.append(np.ma.asarray(pass_map))
        else:
            # Or gather over a number of time steps
            currstep = int(filename[-12:-5])
            tavg_step_i = -1
            tavg_step = int(pass_times)
            for avgstep in np.arange(currstep-tavg_step, currstep+tavg_step+1,1):
                tavg_step_i = tavg_step_i+1
                filenamestep = filename[:-12]+str(avgstep).rjust(7,'0')+'.vlsv'
                print(filenamestep)
                fstep=pt.vlsvfile.VlsvReader(filenamestep)
                step_cellids = fstep.read_variable("CellID")
                pass_maps.append([])
                for mapval in pass_vars:
                    pass_map = fstep.read_variable(mapval)
                    if np.ndim(pass_map)==1:
                        pass_map = pass_map[step_cellids.argsort()].reshape([sizes[1],sizes[0]])
                    else:
                        pass_map = pass_map[step_cellids.argsort()].reshape([sizes[1],sizes[0],len(pass_map[0])])
                    pass_maps[tavg_step_i].append(np.ma.asarray(pass_map))


    # Optional user-defined expression used for color panel instead of a single pre-existing var
    if expression!=None:
        datamap = expression(pass_maps)
        if np.ndim(datamap)!=2:
            print("Error calling custom expression "+expression+"! Result was not a 2-dimensional array. Exiting.")
            return -1

    # Find region outside ionosphere. Note that for some boundary layer cells, a density is calculated, but
    # e.g. pressure is not, and these cells aren't excluded by this method.
    if f.check_variable("rho"):
        rhomap = f.read_variable("rho")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhomap = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap), 0)
    elif f.check_variable("rhom"):
        rhomap = f.read_variable("rhom")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhomap = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap), 0)
    elif f.check_variable("proton/rho"):
        rhomap = f.read_variable("proton/rho")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhomap = np.ma.masked_less_equal(np.ma.masked_invalid(rhomap), 0)
    else:
        print("Unable to exclude non-zero mass density region from range finder!")
        rhomap = (np.ma.masked_invalid(datamap), 0)

    # If automatic range finding is required, find min and max of array
    # Performs range-finding on a masked array to work even if array contains invalid values
    if vmin!=None:
        vminuse=vmin
    else: 
        vminuse=np.ma.amin(np.ma.masked_where(np.ma.getmask(rhomap),datamap))
    if vmax!=None:
        vmaxuse=vmax
    else:
        vmaxuse=np.ma.amax(np.ma.masked_where(np.ma.getmask(rhomap),datamap) )                   

    # If vminuse and vmaxuse are extracted from data, different signs, and close to each other, adjust to be symmetric
    # e.g. to plot transverse field components
    if vmin==None and vmax==None:
        if (vminuse*vmaxuse < 0) and (abs(abs(vminuse)-abs(vmaxuse))/abs(vminuse) < 0.4 ) and (abs(abs(vminuse)-abs(vmaxuse))/abs(vmaxuse) < 0.4 ):
            absval = max(abs(vminuse),abs(vmaxuse))
            if vminuse < 0:
                vminuse = -absval
                vmaxuse = absval
            else:
                vminuse = absval
                vmaxuse = -absval

    # Check that lower bound is valid for logarithmic plots
    if (vminuse <= 0) and (lin==None) and (symlog==None):
        # Drop negative and zero values
        vminuse = np.ma.amin(np.ma.masked_less_equal(np.ma.masked_where(np.ma.getmask(rhomap),datamap),0))

    # If symlog scaling is set:
    if symlog!=None:
        if symlog>0:
            linthresh = symlog 
        else:
            linthresh = max(abs(vminuse),abs(vmaxuse))*1.e-2

    # Lin or log colour scaling, defaults to log
    if lin==None:
        # Special SymLogNorm case
        if symlog!=None:
            #norm = SymLogNorm(linthresh=linthresh, linscale = 0.3, vmin=vminuse, vmax=vmaxuse, ncolors=cmapuse.N, clip=True)
            norm = SymLogNorm(linthresh=linthresh, linscale = 0.3, vmin=vminuse, vmax=vmaxuse, clip=True)
            maxlog=int(np.ceil(np.log10(vmaxuse)))
            minlog=int(np.ceil(np.log10(-vminuse)))
            logthresh=int(np.floor(np.log10(linthresh)))
            logstep=1
            ticks=([-(10**x) for x in range(logthresh, minlog+1, logstep)][::-1]
                    +[0.0]
                    +[(10**x) for x in range(logthresh, maxlog+1, logstep)] )
        else:
            norm = LogNorm(vmin=vminuse,vmax=vmaxuse)
            ticks = LogLocator(base=10,subs=range(10)) # where to show labels
    else:
        # Linear
        levels = MaxNLocator(nbins=255).tick_values(vminuse,vmaxuse)
        norm = BoundaryNorm(levels, ncolors=cmapuse.N, clip=True)
        ticks = np.linspace(vminuse,vmaxuse,num=7)            

    # Select ploitting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw!=None:
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')  

    # Select image shape to match plotted area, (with some accounting for axes etc)
    boxlenx = boxcoords[1]-boxcoords[0]
    boxleny = boxcoords[3]-boxcoords[2]
    # Round the values so that image sizes won't wobble when there's e.g. a moving box and numerical inaccuracies.
    # This is only done if the box size is suitable for the unit in use.
    if ((boxlenx > 10) and (boxleny > 10)):
        boxlenx = float( 0.05 * int(boxlenx*20*1.024) ) 
        boxleny = float( 0.05 * int(boxleny*20*1.024) ) 
    ratio = np.sqrt(boxleny/boxlenx)
    ratio = boxleny/boxlenx
    # default for square figure is figsize=[4.0,3.15]
    figsize = [4.0,3.15*ratio]
    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Generates the mesh to map the data to.
    # Note, datamap is still of shape [ysize,xsize] (?)
    [XmeshXY,YmeshXY] = scipy.meshgrid(np.linspace(simext[0],simext[1],num=sizes[0]),np.linspace(simext[2],simext[3],num=sizes[1]))
    fig1 = plt.pcolormesh(XmeshXY,YmeshXY,datamap, cmap=colormap,norm=norm)
    ax1 = plt.gca() # get current axes

    # Title and plot limits
    if len(plot_title)!=0:
        ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')

    plt.xlim([boxcoords[0],boxcoords[1]])
    plt.ylim([boxcoords[2],boxcoords[3]])
    ax1.set_aspect('equal')

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(thick)
    ax1.xaxis.set_tick_params(width=thick,length=3)
    ax1.yaxis.set_tick_params(width=thick,length=3)
    #ax1.xaxis.set_tick_params(which='minor',width=3,length=5)
    #ax1.yaxis.set_tick_params(which='minor',width=3,length=5)

    if noxlabels==None:
        plt.xlabel('X ['+unitstr+']',fontsize=fontsize,weight='black')
        plt.xticks(fontsize=fontsize,fontweight='black')
        ax1.xaxis.offsetText.set_fontsize(fontsize)# set axis exponent offset font sizes
    if noylabels==None:
        if ysize==1: #Polar
            plt.ylabel('Z ['+unitstr+']',fontsize=fontsize,weight='black')
        else: #Ecliptic
            plt.ylabel('Y ['+unitstr+']',fontsize=fontsize,weight='black')
        plt.yticks(fontsize=fontsize,fontweight='black')
        ax1.yaxis.offsetText.set_fontsize(fontsize)# set axis exponent offset font sizes

    # Limit ticks, slightly according to ratio
    # ax1.xaxis.set_major_locator(plt.MaxNLocator(int(7/np.sqrt(ratio))))
    # ax1.yaxis.set_major_locator(plt.MaxNLocator(int(7*np.sqrt(ratio))))

    # Optional external additional plotting routine overlayed on color plot
    # Uses the same pass_maps variable as expressions
    if external!=None:
        extresult=external(ax1, XmeshXY,YmeshXY, pass_maps)

    if cbtitle==None:
        if expression!=None:
            cb_title_use = expression.__name__.replace("_","\_") # replaces underscores so math mode subscript mode isn't activated
        else:
            cb_title_use = cb_title
    else:
        cb_title_use = cbtitle

    # Colourbar title
    if len(cb_title_use)!=0:
        cb_title_locy = 1.0 + 0.05#/ratio
        plt.text(1.0, 1.01, cb_title_use, fontsize=fontsize,weight='black', transform=ax1.transAxes, horizontalalignment='center')

    if nocb==None:
        # Witchcraft used to place colourbar
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        # First draw colorbar
        if usesci==0:        
            cb = plt.colorbar(fig1,ticks=ticks,cax=cax, drawedges=False)
        else:
            cb = plt.colorbar(fig1,ticks=ticks,format=mtick.FuncFormatter(fmt),cax=cax, drawedges=False)
        cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
        cb.outline.set_linewidth(thick)

        # if too many subticks:
        if lin==None and usesci!=0 and symlog==None:
            # Note: if usesci==0, only tick labels at powers of 10 are shown anyway.
            # For non-square pictures, adjust tick count
            nlabels = len(cb.ax.yaxis.get_ticklabels()) / ratio
            valids = ['1','2','3','4','5','6','7','8','9']
            if nlabels > 10:
                valids = ['1','2','3','4','5','6','8']
            if nlabels > 19:
                valids = ['1','2','5']
            if nlabels > 28:
                valids = ['1']
            # for label in cb.ax.yaxis.get_ticklabels()[::labelincrement]:
            for label in cb.ax.yaxis.get_ticklabels():
                # labels will be in format $x.0\times10^{y}$
                if not label.get_text()[1] in valids:
                    label.set_visible(False)

    # add flux function contours
    if fluxfile != None:
        # Read binary flux function data from prepared files
        flux_function = np.fromfile(fluxfile,dtype='double').reshape(sizes[1],sizes[0])

        # Find inflow position values
        cid = f.get_cellid( [xmax-2*cellsize, 0,0] )
        ff_v = f.read_variable("v", cellids=cid)
        ff_b = f.read_variable("B", cellids=cid)

        # Account for movement
        bdirsign = -1.0 
        outofplane = [0,1,0] # For ecliptic runs
        if zsize==1:
            outofplane = [0,0,1]  # For polar runs
        if np.inner(np.cross(ff_v,ff_b), outofplane) < 0:
            bdirsign = 1.0
        flux_function = flux_function - timeval * np.linalg.norm(np.cross(ff_v,ff_b)) * bdirsign

        # Mask away ionosphere
        flux_function = np.ma.masked_where(~np.isfinite(rhomap), flux_function)
        flux_function = np.ma.masked_where(rhomap<=0, flux_function)

        # The flux level contours must be fixed instead of scaled based on min/max values in order
        # to properly account for flux freeze-in and advection with plasma
        flux_levels = np.linspace(-10,10,fluxlines*60)
        fluxcont = ax1.contour(XmeshXY,YmeshXY,flux_function,flux_levels,colors='k',linestyles='solid',linewidths=0.5*fluxthick,zorder=2)

    # Add Vlasiator watermark
    if wmark!=None:        
        wm = plt.imread(get_sample_data(watermarkimage))
        newax = fig.add_axes([0.01, 0.90, 0.3, 0.08], anchor='NW', zorder=-1)
        newax.imshow(wm)
        newax.axis('off')

    if noxlabels!=None:
        for label in ax1.xaxis.get_ticklabels():
            label.set_visible(False)
    if noylabels!=None:
        for label in ax1.yaxis.get_ticklabels():
            label.set_visible(False)       

    if noborder==None:
        # adjust layout
        plt.tight_layout()
        savefig_pad=0.1 # The default is 0.1
        bbox_inches=None
    else:
        # adjust layout
        plt.tight_layout(pad=0.01)
        savefig_pad=0.01
        bbox_inches='tight'

    # Save output or draw on-screen
    if draw==None:
        print(savefigname+"\n")
        #plt.savefig(savefigname,dpi=300)
        plt.savefig(savefigname,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)

    else:
        plt.draw()
        plt.show()
    plt.close()
    plt.clf()

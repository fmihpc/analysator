import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
import matplotlib.ticker as mtick
import colormaps as cmaps
from matplotlib.cbook import get_sample_data
import plot_run_defaults

from rotation import rotateVectorToVector

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
# plt.register_cmap(name='cork',cmap=cork_map)
# plt.register_cmap(name='davos_r',cmap=davos_r_map)

# Different style scientific format for colour bar ticks
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)

# find nearest spatial cell with vspace to cid
def getNearestCellWithVspace(vlsvReader,cid):
    cell_candidates = vlsvReader.read(mesh='SpatialGrid',tag='CELLSWITHBLOCKS')
    cell_candidate_coordinates = [vlsvReader.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]
    cell_coordinates = vlsvReader.get_cell_coordinates(cid)
    norms = np.sum((cell_candidate_coordinates - cell_coordinates)**2, axis=-1)**(1./2)
    norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))
    return cell_candidates[i]

# create a 2-dimensional histogram
def doHistogram(f,VX,VY,Vpara,vxBinEdges,vyBinEdges,vthick,wflux=None):
    # Flux weighting?
    if wflux!=None:
        fw = f*np.sqrt( np.sum(np.square([VX,VY,Vpara]),1) ) # use particle flux as weighting in the histogram
    else:
        fw = f # use particle phase-space density as weighting in the histogram

    # Select cells which are within slice area
    indexes = [(abs(Vpara) <= 0.5*vthick) & (VX > min(vxBinEdges)) & (VX < max(vxBinEdges)) & (VY > min(vyBinEdges)) & (VY < max(vyBinEdges)) ]


    # Gather histogram of values
    (nVhist,VXEdges,VYEdges) = np.histogram2d(VX[indexes],VY[indexes],bins=(vxBinEdges,vyBinEdges),weights=fw[indexes],normed=0)
    # Gather histogram of how many cells were summed for the histogram
    (Chist,VXEdges,VYEdges) = np.histogram2d(VX[indexes],VY[indexes],bins=(vxBinEdges,vyBinEdges),normed=0)
    # Correct for summing multiple cells into one histogram output cell
    nonzero = np.where(Chist != 0)

    # nonzerocount = sum([len(row) for row in nonzero])
    # print("indices "+str(nonzerocount))
    # one = np.where(Chist ==1)
    # one = sum( [len(row) for row in one])
    # two = np.where(Chist ==2)
    # two = sum( [len(row) for row in two])
    # three = np.where(Chist ==3)
    # three = sum( [len(row) for row in three])
    # four = np.where(Chist ==4)
    # four = sum( [len(row) for row in four])
    # five = np.where(Chist ==5)
    # five = sum( [len(row) for row in five])
    # six = np.where(Chist ==6)
    # six = sum( [len(row) for row in six])
    # print("vthick "+str(vthick))
    # print("One "+str(one)+" two "+str(two)+" three "+str(three)+" four "+str(four)+" five "+str(five)+" six "+str(six))

    nVhist[nonzero] = np.divide(nVhist[nonzero],Chist[nonzero])


    # Please note that the histogram does not follow the Cartesian convention where x values are on the abscissa
    # and y values on the ordinate axis. Rather, x is histogrammed along the first dimension of the array (vertical),
    # and y along the second dimension of the array (horizontal). This ensures compatibility with histogramdd.
    nVhist = nVhist.transpose()

    # Flux weighting
    if wflux!=None:
        dV = np.abs(vxBinEdges[-1] - vxBinEdges[-2]) # assumes constant bin size
        nVhist = np.divide(nVhist,(dV*4*np.pi)) # normalization

    return (nVhist,VXEdges,VYEdges)
  

# analyze velocity space in a spatial cell (velocity space reducer)
def vSpaceReducer(vlsvReader, cid, slicetype, normvect, VXBins, VYBins, pop="proton", 
                  slicethick=None, wflux=None, cbulk=None, center=None, keepfmin=None):
    # check if velocity space exists in this cell
    if vlsvReader.check_variable('fSaved'): #restart files will not have this value
        if vlsvReader.read_variable('fSaved',cid) != 1.0:
            return (False,0,0,0)

    # Assume velocity cells are cubes
    [vxsize, vysize, vzsize] = vlsvReader.get_velocity_mesh_size(pop=pop)
    # Account for 4x4x4 cells per block
    vxsize = 4*vxsize
    vysize = 4*vysize
    vzsize = 4*vzsize
    [vxmin, vymin, vzmin, vxmax, vymax, vzmax] = vlsvReader.get_velocity_mesh_extent(pop=pop)
    inputcellsize=(vxmax-vxmin)/vxsize
    print("Input velocity grid cell size "+str(inputcellsize))

    velcells = vlsvReader.read_velocity_cells(cid, pop=pop)
    V = vlsvReader.get_velocity_cell_coordinates(velcells.keys(), pop=pop)
    print("Found "+str(len(V))+" v-space cells")

    if cbulk!=None:
        print("Transforming to plasma frame")
        if vlsvReader.check_variable('moments'):
            # This should be a restart file
            moments = np.array(vlsvReader.read_variable('moments',cid))
            if moments==None:
                print("Error reading moments from assumed restart file!")
                exit()
            if len(moments.shape)==2:
                moments = moments[0]
            if moments[0]>0.0:
                bulkv = moments[1:4]/moments[0]
            else:
                bulkv=[0.,0.,0.]
        elif vlsvReader.check_variable('v'):
            # Multipop file with bulk v saved directly
            bulkv = vlsvReader.read_variable('v',cid)
        elif ( vlsvReader.check_variable('rho') and vlsvReader.check_variable('rho_v')):
            # Older regular bulk file
            rhov = vlsvReader.read_variable('rho_v',cid)
            rho = vlsvReader.read_variable('rho',cid)
            bulkv = rhov/rho
        elif vlsvReader.check_variable(pop+'/V'):
            # Multipop file without V saved, use per-population bulk velocity
            bulkv = vlsvReader.read_variable(pop+'/V',cid)
        else:
            print("Error in finding plasma bulk velocity!")
            exit()
        # shift to velocities plasma frame
        V = V - bulkv
    elif center!=None:
        if len(center)==3:
            print("Transforming to frame travelling at speed "+str(center))
            V - V - center
        else:
            print("Error in shape of center vector! Give in form (vx,vy,vz).")

    f = zip(*velcells.items())
    # check that velocity space has cells
    if(len(f) > 0):
        f = np.asarray(zip(*velcells.items())[1])
    else:
        return (False,0,0,0)

    if keepfmin==None:
        # Drop all velocity cells which are below the sparsity threshold. Otherwise the plot will show buffer cells as well.
        fMin = 1e-16 # default
        if vlsvReader.check_variable('MinValue') == True:
            fMin = vlsvReader.read_variable('MinValue',cid)
        ii_f = np.where(f >= fMin)
        print("Dropping velocity cells under fMin value "+str(fMin))
        if len(ii_f) < 1:
            return (False,0,0,0)
        f = f[ii_f]
        V = V[ii_f,:][0,:,:]

    if slicethick==None:
        # Geometric magic to widen the slice to assure that each cell has some velocity grid points inside it.
        # Might still be incorrect, erring on the side of caution.
        # norm_srt = sorted(abs(normvect))
        # if norm_srt[1] > 0:
        #     temp = norm_srt[0]/norm_srt[1]
        #     aratio = (1.+temp)/np.sqrt( 1+temp**2)
        # else:
        #     aratio = 1.
        # gridratio = aratio * norm_srt[2] * (1. + aratio * np.sqrt(norm_srt[0]**2 + norm_srt[1]**2) / norm_srt[2] )
        # if gridratio>1.0:
        #     # Account for numerical inaccuracy
        #     slicethick=inputcellsize*gridratio*1.01
        # else:
        #     slicethick=inputcellsize

        samplebox=np.array([ [0.0,0.0,0.0], [0.0,0.0,1.0], [0.0,1.0,0.0], [0.0,1.0,1.0], [1.0,0.0,0.0], [1.0,0.0,1.0], [1.0,1.0,0.0], [1.0,1.0,1.0] ])
        sbrot = rotateVectorToVector(samplebox,normvect)
        rotminx=np.amin(sbrot[:,0])
        rotmaxx=np.amax(sbrot[:,0])
        rotminy=np.amin(sbrot[:,1])
        rotmaxy=np.amax(sbrot[:,1])
        rotminz=np.amin(sbrot[:,2])
        rotmaxz=np.amax(sbrot[:,2])
        gridratio = np.amax([ rotmaxx-rotminx, rotmaxy-rotminy, rotmaxz-rotminz ])
        slicethick=inputcellsize*gridratio
    else:
        slicethick=inputcellsize*slicethick
    print("Performing slice with a counting thickness of "+str(slicethick))

    if slicetype=="xy":
        VX = V[:,0]
        VY = V[:,1]
        Vpara = V[:,2]
    elif slicetype=="yz":
        VX = V[:,1]
        VY = V[:,2]
        Vpara = V[:,0]
    elif slicetype=="xz":
        VX = V[:,0]
        VY = V[:,2]
        Vpara = V[:,1]
    elif slicetype=="vecperp":
        # Find velocity components in given frame (e.g. B frame: (vx,vy,vz) -> (vperp2,vperp1,vpar))
        N = np.array(normvect)/np.sqrt(normvect[0]**2 + normvect[1]**2 + normvect[2]**2)
        Vrot = rotateVectorToVector(V,N)
        VX = Vrot[:,0]
        VY = Vrot[:,1]
        Vpara = Vrot[:,2]
    elif slicetype=="vecpara":
        N = np.array(normvect)/np.sqrt(normvect[0]**2 + normvect[1]**2 + normvect[2]**2)
        Vrot = rotateVectorToVector(V,N)
        VX = Vrot[:,2]
        VY = Vrot[:,1]
        Vpara = Vrot[:,0]
    else:
        print("Error finding rotation of v-space!")
        return (False,0,0,0)

    # TODO: better rotation so perpendicular component can be defined?

    # create 2-dimensional histogram of velocity components perpendicular to slice-normal vector
    (binsXY,edgesX,edgesY) = doHistogram(f,VX,VY,Vpara,VXBins,VYBins, slicethick, wflux=wflux)
    return (True,binsXY,edgesX,edgesY)





def plot_colormap_with_vdf(filename=None,
                     vlsvobj=None,
                     filedir=None, step=None,
                     outputdir=None,popvdf="proton",
                     var=None, op=None, title=None,
                     draw=None, usesci=None,
                     symlog=None,varnorm=1.,magarrows=False,
                     boxm=[],boxre=[],colormap=None,
                     run=None,notime=None,wmark=None,
                     notre=None, thick=1.0, cbtitle=None,
                     vmin=None, vmax=None, lin=None,
                     fmin=None, fmax=None, cbulk=None,
                     external=None, extvals=None,
                     expression=None, pass_vars=None,
                     cellcoordplot=None,cellidplot=None,
                     vdfplotlocation=None, vdfplotsize=0.1,
                     xlines=False,
                     fluxfile=None, fluxdir=None,
                     fluxthick=1.0, fluxlines=1
                     ):    

    ''' Plots a coloured plot with VDFs in chosen points.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
     
    :kword popvdf:      population (species) for VDF plots, default "proton"
    :kword var:         variable to plot, e.g. rho, rhoBeam, beta, temperature, MA, Mms, va, vms,
                        E, B, V or others. Accepts any variable known by analysator/pytools.
    :kword op:          Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude. 

    :kword boxm:        zoom box extents [x0,x1,y0,y1] in metres (default and truncate to: whole simulation box)
    :kword boxre:       zoom box extents [x0,x1,y0,y1] in Earth radii (default and truncate to: whole simulation box)
    :kword colormap:    colour scale for plot, use e.g. jet, viridis, plasma, inferno, magma, nipy_spectral, RdBu
    :kword run:         run identifier, used for some default vmin,vmax values and for constructing output filename
    :kword notime:      flag to suppress plotting simulation time in title
    :kword title:       string to use as title instead of map name
    :kword cbtitle:     string to use as colour bar title instead of map name
    :kword notre:       flag to use metres (if ==1) or kilometres as axis unit
    :kword thick:       line and axis thickness, default=1.0
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: 1)
    :kword vmin,vmax:   min and max values for colour scale and colour bar of the colormap. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword fmin,fmax:   min and max values for colour scale and colour bar of the VDFs. If no values are given,
                        min and max values for whole plot are used.
    :kword cbulk:       Center plot on position of bulk velocity (for this population)
    :kword lin:         flag for using linear colour scaling instead of log
    :kword symlog:      use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    :kword varnorm:     Normalisation factor for the colormap variable
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword draw:        Draw image on-screen instead of saving to file (requires x-windowing)
   
    :kword external:    Optional function which receives the image axes in order to do further plotting
    :kword extvals:     Optional array of map names to pass to the external function

    :kword expression:  Optional function which calculates a custom expression to plot. Remember to set
                        vmin and vmax manually.
    :kword pass_vars:   Array of map names to pass to the optional expression function (as np.arrays)

    :kword cellcoordplot:   Coordinates of cells to display as circles in the colormap plot, format [x1,y1,z1,...,xn,yn,zn]
    :kword cellidplot:      List of cellIDs to display as circles in the colormap plot
    :kword vdfplotlocation: String specifying the location of the VDF plot with respect to its cell (t:top, b:bottom, l:left, r:right)
    :kword vdfplotsize:     Relative size of the VDF plot axes compared to the main figure

    :kword xlines:      Plots the location of X-lines
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


    # Earth radius [m]
    R_E = 6.371e6

    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    # watermarkimage=os.path.expandvars('$HOME/appl_taito/analysator/pyPlot/logo_color.png')
    # watermarkimage='/homeappl/home/marbat/appl_taito/analysator/logo_color.png'

    if outputdir==None:
        outputdir=os.path.expandvars('$HOME/Plots/')
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
        #colormap="plasma"
        #colormap="viridis_r"
        #colormap="inferno"
        #colormap="seismic"
        colormap="plasma_r"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8 # Most text
    fontsize2=10 # Time title
    fontsize3=5 # Colour bar ticks

    # Plot title with time
    if notime==None:        
        timeval=f.read_parameter("time")
        if timeval == None:
            timeval=f.read_parameter("t")
            if timeval == None:    
                print "Unknown time format encountered"
                plot_title = ''
        if timeval != None:
            plot_title = "t="+"{:.1f}".format(timeval)+' s'
    else:
        plot_title = ''       

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
    savefigname = outputdir+run+"_map_"+varstr+opstr+stepstr+".png"

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
    unitstr = r'$\mathrm{R}_{\mathrm{E}}$'
    unit = Re
    if notre!=None: # Use m or km instead
        if notre==1:
            unit = 1.0
            unitstr = 'm'
        else:
            unit = 1.e3
            unitstr = 'km'

    # Scale data extent and plot box
    simext=[i/unit for i in simext]
    boxcoords=[i/unit for i in boxcoords]

    ##########
    # Read data and calculate required variables
    ##########
    if expression==None:
        if var == 'rho':
            cb_title = r"$n_\mathrm{p} [\mathrm{m}^{-3}]$"
            datamap = f.read_variable("rho")

        elif var == 'rhoBeam':
            cb_title = r"$\rho_{\mathrm{beam}} [\mathrm{m}^{-3}]$"
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
                # datamap = datamap*1e+9 # could be used to output nanotesla instead of tesla

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
                datamap = f.read_variable("v",operator='magnitude')
            else:
                cb_title = r"$V_"+op+"\,[\mathrm{m}\,\mathrm{s}^{-1}]$"
                datamap = f.read_variable("v",operator=op)
                # datamap = datamap*1e-3 # Plot this as km/s instead of m/s

        else:
            # Pipe all other vars directly to analysator
            if op==None:
                cb_title = var
                datamap = f.read_variable(var)
                # If value was vector value, take magnitude
                if np.ndim(datamap) != 1:
                    cb_title = r"$|"+var+"|$"
                    datamap = np.sum(np.asarray(datamap)**2,axis=-1)**(0.5)
            else:
                cb_title = r+""+var+"$_"+op+"$"
                datamap = f.read_variable(var,operator=op)            
            
        if np.ndim(datamap)!=1:
            print("Error reading variable "+var+"! Exiting.")
            return -1

        # Reshape data to an ordered 2D array that can be plotted
        if np.ndim(datamap) != 2:
            datamap = datamap[cellids.argsort()].reshape([sizes[1],sizes[0]])

    else:
    # Optional user-defined expression overrides the var
    # Optional external additional plotting routine
        exprmaps=[]
        if pass_vars==None:
            print("Error, expression must have some variable maps to work on.")
            return
        else:
            # Gather the required variable maps for the expression function
            for mapval in pass_vars:
                exprmap = f.read_variable(mapval)
                if np.ndim(exprmap)==1:
                    exprmap = exprmap[cellids.argsort()].reshape([sizes[1],sizes[0]])
                else:
                    exprmap = exprmap[cellids.argsort()].reshape([sizes[1],sizes[0],len(exprmap[0])])
                exprmaps.append(np.ma.asarray(exprmap))
        datamap = expression(exprmaps)             
        if np.ndim(datamap)!=2:
            print("Error calling custom expression "+expression+"! Result was not a 2-dimensional array. Exiting.")
            return -1

    # Find region outside ionosphere
    if f.check_variable("rho"):
        rhomap = f.read_variable("rho")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhoindex = np.where(rhomap > 1e-30)
    elif f.check_variable("rhom"):
        rhomap = f.read_variable("rhom")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhoindex = np.where(rhomap > 1.e-30)
    elif f.check_variable("proton/rho"):
        rhomap = f.read_variable("proton/rho")
        rhomap = rhomap[cellids.argsort()].reshape([sizes[1],sizes[0]])
        rhoindex = np.where(rhomap > 1e-30)
    else:
        print("Unable to exclude non-zero mass density region from range finder!")
        rhomap = datamap
        rhoindex = np.where(np.isfinite(rhomap))

    datamap = datamap/varnorm

    # If automatic range finding is required, find min and max of array
    if vmin!=None:
        vminuse=vmin
    else: 
        vminuse=np.amin(datamap[rhoindex])        
        # print("Using deduced minimum value "+str(vminuse))
    if vmax!=None:
        vmaxuse=vmax
    else:
        vmaxuse=np.amax(datamap[rhoindex])
        # print("Using deduced maximum value "+str(vmaxuse))

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

    # Check that lower bound is valid
    if (vminuse <= 0) and (lin==None) and (symlog==None):
        # Assume 5 orders of magnitude is enough?
        print("Vmin value invalid for log scale, defaulting to 1.e-5 of maximum value")
        vminuse = vmaxuse*1.e-5

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

    # Select image shape to match plotted area, at least somewhat.
    # default for square figure is figsize=[4.0,3.15]
    ratio = np.sqrt((boxcoords[3]-boxcoords[2])/(boxcoords[1]-boxcoords[0]))
    figsize = [4.0,3.15*ratio]

    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Generates the mesh to map the data to
    [XmeshXY,YmeshXY] = scipy.meshgrid(np.linspace(simext[0],simext[1],num=sizes[0]),np.linspace(simext[2],simext[3],num=sizes[1]))
    fig1 = plt.pcolormesh(XmeshXY,YmeshXY,datamap, cmap=colormap,norm=norm)
    ax1 = plt.gca() # get current axes

    # Optionally plots the locations of the X-lines
    if xlines:
        X_points, O_points = pt.calculations.find_X_O(fluxfile, filename, step)
        X_points_x = np.array([])
        X_points_z = np.array([])
        for ixline in range(0,len(X_points)):
            X_points_x = np.append(X_points_x,X_points[ixline][0]/Re)
            X_points_z = np.append(X_points_z,X_points[ixline][2]/Re)
        ax1.plot(X_points_x,X_points_z,marker="x",markerfacecolor="None",markeredgecolor='r',markeredgewidth=2,linewidth=0)

    # Title and plot limits
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

    # Limit ticks, slightly according to ratio
    ax1.xaxis.set_major_locator(plt.MaxNLocator(int(7/np.sqrt(ratio))))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(int(7*np.sqrt(ratio))))

    plt.xlabel('X ['+unitstr+']',fontsize=fontsize,weight='black')
    if ysize==1: #Polar
        plt.ylabel('Z ['+unitstr+']',fontsize=fontsize,weight='black')
    else: #Ecliptic
        plt.ylabel('Y ['+unitstr+']',fontsize=fontsize,weight='black')
    plt.xticks(fontsize=fontsize,fontweight='black')
    plt.yticks(fontsize=fontsize,fontweight='black')

    # set axis exponent offset font sizes
    ax1.yaxis.offsetText.set_fontsize(fontsize)
    ax1.xaxis.offsetText.set_fontsize(fontsize)

    # Optional external additional plotting routine overlayed on color plot
    if external!=None:
        extmaps=[]
        if extvals!=None:
            for mapval in extvals:
                extmap = f.read_variable(mapval)
                if np.ndim(extmap)==1:
                    extmap = extmap[cellids.argsort()].reshape([sizes[1],sizes[0]])
                else:
                    extmap = extmap[cellids.argsort()].reshape([sizes[1],sizes[0],len(extmap[0])])
                extmaps.append(extmap)
        extresult=external(ax1, XmeshXY,YmeshXY, extmaps)            

    if cbtitle==None:
        if expression!=None:
            cb_title_use = expression.__name__.replace("_","\_") # replaces underscores so math mode subscript mode isn't activated
        else:
            cb_title_use = cb_title
    else:
        cb_title_use = cbtitle

    # Get coordinates if cellIDs were given as input
    if cellidplot!=None:
        cellcoordplot=[]
        for icell in range(0,len(cellidplot)):
            xCid,yCid,zCid = f.get_cell_coordinates(cellidplot[icell])
            cellcoordplot = cellcoordplot + [xCid/R_E,yCid/R_E,zCid/R_E]
        print('Coordinates: '+str(cellcoordplot))







    # Colourbar title
    cb_title_locy = 1.0 + 0.05/ratio
    plt.text(1.0, 1.05, cb_title_use, fontsize=fontsize,weight='black', transform=ax1.transAxes)

    # Witchcraft used to place colourbar
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad="2%")
    # First draw colorbar
    if usesci==0:        
        cb = plt.colorbar(fig1,ticks=ticks,cax=cax)
    else:
        cb = plt.colorbar(fig1,ticks=ticks,format=mtick.FuncFormatter(fmt),cax=cax)
    cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
    cb.outline.set_linewidth(thick)

    # if too many subticks:
    if lin==None and usesci!=0 and symlog==None:
        # Note: if usesci==0, only tick labels at powers of 10 are shown anyway.
        # For non-square pictures, adjust tick count
        nlabels = len(cb.ax.yaxis.get_ticklabels()) / ratio
        if nlabels > 10:
            valids = ['1','2','3','4','5','6','8']
        if nlabels > 19:
            valids = ['1','2','5']
        if nlabels > 28:
            valids = ['1']
        #for label in cb.ax.yaxis.get_ticklabels()[::labelincrement]:
        if nlabels > 10:
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

    # Loop over all cellIDs for which VDF should be plotted
    if cellcoordplot!=None:
        for aux in range(0,len(cellcoordplot)/3):
            x,y,z = cellcoordplot[3*aux],cellcoordplot[3*aux+1],cellcoordplot[3*aux+2]
            print('cellcoord2plot #' + str(aux) + ': x = ' + str(x) + ', y = ' + str(y)  + ', z = ' + str(z))
            ax1.plot(x,z,marker="o",markerfacecolor="None",markeredgecolor='w')
            if magarrows:
                Bvect = f.read_variable("B", cellidplot[aux])
                Bvect = Bvect / 1e-8 #np.sqrt(Bvect[0]**2 + Bvect[1]**2 + Bvect[2]**2)
                dx,dy,dz = Bvect[0],Bvect[1],Bvect[2]
                Byimp = abs(dy)/np.sqrt(dx**2+dy**2+dz**2)
                ax1.arrow(x,z,dx,dz)   #,color=str(Byimp))
                print("arrow plotted"+str(dx)+' '+str(dz))
                print("By relative importance: "+str(Byimp))
                if aux==0:
                    vdfplotlocation = ''
                if dz<0:
                    vdfplotlocation = vdfplotlocation+'t'
                else:
                    vdfplotlocation = vdfplotlocation+'b'
            xratio = (x-boxcoords[0])/(boxcoords[1]-boxcoords[0])
            yratio = (z-boxcoords[2])/(boxcoords[3]-boxcoords[2])
            pos1 = ax1.get_position()
            print('ax1 position: '+str(pos1))

            # Magic numbers so that the x-axis alignment of VDF plots is correct (Note: those values work for equal X and Y ranges)
            # TODO Try to figure out a way to automatically calculate them whichever the X/Y ratio
            magixconstant = 0.05
            magixratio = 0.815

            # Location of VDF plot axes with respect to their cell
            if vdfplotlocation==None:
                vdfaxloc = 'b'
            else:
                vdfaxloc = vdfplotlocation[aux]
            
            if vdfaxloc=='l':
                xax2 = pos1.x0 + magixconstant + pos1.width*magixratio*xratio - vdfplotsize - 0.01
                yax2 = pos1.y0 + yratio*pos1.height - vdfplotsize/2.
            elif vdfaxloc=='r':
                xax2 = pos1.x0 + magixconstant + pos1.width*magixratio*xratio + 0.01
                yax2 = pos1.y0 + yratio*pos1.height - vdfplotsize/2.
            elif vdfaxloc=='t':
                xax2 = pos1.x0 + magixconstant + pos1.width*magixratio*xratio - vdfplotsize/2.
                yax2 = pos1.y0 + yratio*pos1.height + 0.02
            # Default: bottom
            else:
                xax2 = pos1.x0 + magixconstant + pos1.width*magixratio*xratio - vdfplotsize/2.
                yax2 = pos1.y0 + yratio*pos1.height - vdfplotsize - 0.02
 #               print('ax2 x and y: '+str(xax2)+'  '+str(yax2))
            ax2 = fig.add_axes([xax2,yax2,vdfplotsize,vdfplotsize])
            
            
            subplot_vdf(axis=ax2,filename=filename,cellids=cellidplot[aux],box=[-3.e6,3.e6,-3.e6,3.e6],pop=popvdf,colormap='nipy_spectral',notitle=1,noxlabels=1,noylabels=1,bpara=1,fmin=fmin,fmax=fmax,cbulk=cbulk)
 #           ax2.pcolormesh(XmeshXY,YmeshXY,binsXY, cmap='nipy_spectral',norm=norm)




    # adjust layout
#    plt.tight_layout()

    # Save output or draw on-screen
    if draw==None:
        print(savefigname+"\n")
        plt.savefig(savefigname,dpi=300)
    else:
        plt.draw()
        plt.show()
    plt.close()
    plt.clf()





def subplot_vdf(axis=None,
             filename=None,
             vlsvobj=None,
             filedir=None, step=None,
             cellids=None, coordinates=None,
             pop="proton",
             draw=None,unit=None,cellsize=None,
             colormap=None, box=None, cbar=None,
             run=None, notitle=None, wmark=None, thick=1.0,
             fmin=None, fmax=None, slicethick=None,
             xy=None, xz=None, yz=None, normal=None,
             bpara=None, bperp=None,
             cbulk=None, center=None, wflux=None, keepfmin=None,
             legend=None, noborder=None, scale=1.0,
             biglabel=None, biglabloc=None,
             noxlabels=None, noylabels=None,
             coordswap=None
             ):

    ''' Plots a coloured plot with axes and a colour bar.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
     
    :kword cellids:     list of cell IDs to plot for
    :kword coordinates: list of 3-element coordinates to plot for
    :kword pop:         population (species), default "proton"

    :kword colormap:    colour scale for plot, use e.g. jet, viridis, plasma, inferno, magma, nipy_spectral, RdBu
    :kword run:         run identifier, used for constructing output filename
    :kword title:       string to use as title in lieu of map name
    :kword thick:       line and axis thickness, default=1.0
    :kword fmin,fmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot are used.

    :kword box:         extents of plotted velocity grid as [x0,y0,x1,y1]
    :kword unit:        Plot v-axes using 10^unit m/s instead of km/s
   
    :kword xy:          Perform slice in x-y-direction
    :kword xz:          Perform slice in x-z-direction
    :kword yz:          Perform slice in y-z-direction
    :kword normal:      Perform slice in plane perpendicular to given vector
    :kword bpara:       Perform slice in B_para / B_perp2 plane
    :kword bperp:       Perform slice in B_perp1 / B_perp2 plane
                        If no plane is given, default is simulation plane (for 2D simulations)

    :kword cbulk:       Center plot on position of bulk velocity (for this population)
    :kword center:      Center plot on provided velocity vector position
    :kword wflux:       Plot flux instead of distribution function
    :kword slicethick:  Thickness of slice (default 1 cell)
    :kword keepfmin:    Also draw buffer cells with values below fMin

    :kword draw:        Draw image on-screen instead of saving to file (requires x-windowing)
    :kword cbar:        Plot colourbar legend (default off). 
    :kword biglabel:    Plot large label (in top-left corner)
    :kword biglabloc:   Move large label to: 0: NW 1: NE 2: SE 3: SW corner

    :kword noborder:    Plot figure edge-to-edge without borders (default off)
    :kword noxlabels:   Suppress x-axis labels and title
    :kword noylabels:   Suppress y-axis labels and title
    :kword notitle:     flag to suppress plot and colorbar titles
    :kword scale:       Scale text size (default=1.0)

    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:

    '''


    # Input file or object
    if filename!=None:
        vlsvReader=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        vlsvReader=vlsvobj
    elif ((filedir!=None) and (step!=None)):
        filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
        vlsvReader=pt.vlsvfile.VlsvReader(filename)
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    if colormap==None:
        colormap="hot_desaturated"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=5*scale # Colour bar ticks
    fontsize4=16*scale # Big label

    # Plot title with time
    timeval=vlsvReader.read_parameter("time")
    if timeval==None:
        timeval=vlsvReader.read_parameter("t")

    # step, used for file name
    if step!=None:
        stepstr = '_'+str(step).rjust(7,'0')
    else:
        if timeval != None:
            stepstr = '_t'+str(np.int(timeval))
        else:
            stepstr = ''

    # If run name isn't given, just put "plot" in the output file name
    if run==None:
        run='plot'

    # If population isn't defined i.e. defaults to protons, check if 
    # instead should use old version "avgs"
    if pop=="proton":
       if not vlsvReader.check_population(pop):
           if vlsvReader.check_population("avgs"):
               pop="avgs"
               print("Auto-switched to population avgs")
           else:
               print("Unable to detect population "+pop+" in .vlsv file!")
               exit()
    else:
        if not vlsvReader.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            exit()       

    #read in mesh size and cells in ordinary space
    # xsize = vlsvReader.read_parameter("xcells_ini")
    # ysize = vlsvReader.read_parameter("ycells_ini")
    # zsize = vlsvReader.read_parameter("zcells_ini")
    [xsize, ysize, zsize] = vlsvReader.get_spatial_mesh_size()
    # cellids = vlsvReader.read_variable("CellID")
    # vxsize = vlsvReader.read_parameter("vxblocks_ini")*4
    # vysize = vlsvReader.read_parameter("vyblocks_ini")*4
    # vzsize = vlsvReader.read_parameter("vzblocks_ini")*4
    # vxmin = vlsvReader.read_parameter("vxmin")
    # vxmax = vlsvReader.read_parameter("vxmax")
    # vymin = vlsvReader.read_parameter("vymin")
    # vymax = vlsvReader.read_parameter("vymax")
    # vzmin = vlsvReader.read_parameter("vzmin")
    # vzmax = vlsvReader.read_parameter("vzmax")

    # These apparently work with multipop at least for the default mesh of protons.
    # Also works with older vlsv versions.
    [vxsize, vysize, vzsize] = vlsvReader.get_velocity_mesh_size(pop=pop)
    [vxmin, vymin, vzmin, vxmax, vymax, vzmax] = vlsvReader.get_velocity_mesh_extent(pop=pop)
    inputcellsize=(vxmax-vxmin)/vxsize

    # account for 4x4x4 cells per block
    vxsize = 4*vxsize
    vysize = 4*vysize
    vzsize = 4*vzsize

    Re = 6.371e+6 # Earth radius in m
    # unit of velocity
    velUnit = 1e3
    velUnitStr = '[km/s]'
    if unit!=None:
        velUnit = np.power(10,int(unit))
        if unit==1:
            velUnitStr = r'[m/s]'
        else:
            velUnitStr = r'[$10^{'+str(int(unit))+'}$ m/s]'


    if (cellids==None and coordinates==None):
        print("Error: must provide either cell id's or coordinates")
        return -1

    if coordinates!=None:
        if type(coordinates) is not list:
            coordinates = [coordinates]

        # Calculate cell IDs from given coordinates        
        xReq = np.asarray(coordinates).T[0]
        yReq = np.asarray(coordinates).T[1]
        zReq = np.asarray(coordinates).T[2]
        if xReq.shape == yReq.shape == zReq.shape:
            print('Number of points: ' + str(xReq.shape[0]))
        else:
            print('ERROR: bad coordinate variables given')
            exit()
        cidsTemp = []
        for ii in range(xReq.shape[0]):
            cidRequest = (np.int64)(vlsvReader.get_cellid(np.array([xReq[ii],yReq[ii],zReq[ii]])))
            cidNearestVspace = -1
            if cidRequest > 0:
                cidNearestVspace = getNearestCellWithVspace(vlsvReader,cidRequest)
            else:
                print('ERROR: cell not found')
                exit()
            if (cidNearestVspace <= 0):
                print('ERROR: cell with vspace not found')
                exit()
            xCid,yCid,zCid = vlsvReader.get_cell_coordinates(cidRequest)
            xVCid,yVCid,zVCid = vlsvReader.get_cell_coordinates(cidNearestVspace)
            print('Point: ' + str(ii+1) + '/' + str(xReq.shape[0]))
            print('Requested coordinates : ' + str(xReq[ii]/Re) + ', ' + str(yReq[ii]/Re) + ', ' + str(zReq[ii]/Re))
            print('Nearest spatial cell  : ' + str(xCid/Re)    + ', ' + str(yCid/Re)    + ', ' + str(zCid/Re))
            print('Nearest vspace        : ' + str(xVCid/Re)   + ', ' + str(yVCid/Re)   + ', ' + str(zVCid/Re))
            cidsTemp.append(cidNearestVspace)
        cellids = np.unique(cidsTemp)
        print('Unique cells with vspace found: ' + str(len(cidsTemp)))
    else:
        print('Using given cell ids and assuming vspace is stored in them')


    # Loop over all cell ids
    if type(cellids) is not list:
        cellids = [cellids]

    print(cellids)
    for cellid in cellids:
        # Initialise some values
        fminuse=None
        fmaxuse=None

        x,y,z = vlsvReader.get_cell_coordinates(cellid)
        print('cellid ' + str(cellid) + ', x = ' + str(x) + ', y = ' + str(y)  + ', z = ' + str(z))

        # Check slice to perform (and possibly normal vector)
        normvect=None
        if xy==None and xz==None and yz==None and normal==None and bpara==None and bperp==None:
            # Use default slice for this simulation
            # Check if ecliptic or polar run
            if ysize==1: # polar
                slicetype="xz"
                pltxstr=r"$v_x$ "+velUnitStr
                pltystr=r"$v_z$ "+velUnitStr
                normvect=[0,1,0] # used just for cell size normalisation
            elif zsize==1: # ecliptic
                slicetype="xy"
                pltxstr=r"$v_x$ "+velUnitStr
                pltystr=r"$v_y$ "+velUnitStr
                normvect=[0,0,1] # used just for cell size normalisation
            else:
                print("Problem finding default slice direction")
                slicetype="yz"
                pltxstr=r"$v_y$ "+velUnitStr
                pltystr=r"$v_z$ "+velUnitStr
                normvect=[1,0,0] # used just for cell size normalisation
        elif normal!=None:
            if len(normal)==3:
                slicetype="vecperp"
                normvect=normal
                pltxstr=r"$v_1$ "+velUnitStr
                pltystr=r"$v_2$ "+velUnitStr
            else:
                print("Error parsing slice normal vector!")
                exit()
        elif xy!=None:
            slicetype="xy"
            pltxstr=r"$v_x$ "+velUnitStr
            pltystr=r"$v_y$ "+velUnitStr
            normvect=[0,0,1] # used just for cell size normalisation
        elif xz!=None:
            slicetype="xz"
            pltxstr=r"$v_x$ "+velUnitStr
            pltystr=r"$v_z$ "+velUnitStr
            normvect=[0,1,0] # used just for cell size normalisation
        elif yz!=None:
            slicetype="yz"
            pltxstr=r"$v_y$ "+velUnitStr
            pltystr=r"$v_z$ "+velUnitStr
            normvect=[1,0,0] # used just for cell size normalisation
        elif bpara!=None or bperp!=None:
            # Rotate based on B-vector
            if vlsvReader.check_variable("B"):
                Bvect = vlsvReader.read_variable("B", cellid)
            elif (vlsvReader.check_variable("background_B") and vlsvReader.check_variable("perturbed_B")):
                # used e.g. for restart files
                BGB = vlsvReader.read_variable("background_B", cellid)
                PERBB = vlsvReader.read_variable("perturbed_B", cellid)
                Bvect = BGB+PERBB
            else:
                print("Error finding B vector direction!")
                exit()

            if Bvect.shape==(1,3):
                Bvect = Bvect[0]
            normvect = Bvect

            if bperp!=None:
                # slice in b_perp1/b_perp2
                slicetype="vecperp"
                pltxstr=r"$v_{\perp 1}$ "+velUnitStr
                pltystr=r"$v_{\perp 2}$ "+velUnitStr
            else:
                # means bpara!=None, slice in b_parallel/b_perp2 plane
                slicetype="vecpara"
                pltxstr=r"$v_{\parallel}$ "+velUnitStr
                pltystr=r"$v_{\perp}$ "+velUnitStr

        # Extend velocity space and each cell to account for slice directions oblique to axes
        normvect = np.array(normvect)
        normvect = normvect/np.linalg.norm(normvect)

        # Geometric magic to stretch the grid to assure that each cell has some velocity grid points inside it.
        # Might still be incorrect, erring on the side of caution.
        # norm_srt = sorted(abs(normvect))
        # if cellsize==None:
        #     if norm_srt[1] > 0:
        #         temp = norm_srt[0]/norm_srt[1]
        #         aratio = (1.+temp)/np.sqrt( 1+temp**2)
        #     else:
        #         aratio = 1.
        #     gridratio = aratio * norm_srt[2] * (1. + aratio * np.sqrt(norm_srt[0]**2 + norm_srt[1]**2) / norm_srt[2] )
        #     if gridratio>1.0:
        #         gridratio = gridratio*1.01 # account for numerical inaccuracies
        # else:
        #     gridratio = 1.

        if cellsize==None:
            samplebox=np.array([ [0.0,0.0,0.0], [0.0,0.0,1.0], [0.0,1.0,0.0], [0.0,1.0,1.0], [1.0,0.0,0.0], [1.0,0.0,1.0], [1.0,1.0,0.0], [1.0,1.0,1.0] ])
            sbrot = rotateVectorToVector(samplebox,normvect)
            rotminx=np.amin(sbrot[:,0])
            rotmaxx=np.amax(sbrot[:,0])
            rotminy=np.amin(sbrot[:,1])
            rotmaxy=np.amax(sbrot[:,1])
            rotminz=np.amin(sbrot[:,2])
            rotmaxz=np.amax(sbrot[:,2])
            gridratio = np.amax([ rotmaxx-rotminx, rotmaxy-rotminy, rotmaxz-rotminz ])
        else:
            gridratio = cellsize

        # num must be vxsize+1 or vysize+1 in order to do both edges for each cell
        VXBins = np.linspace(vxmin*gridratio,vxmax*gridratio,num=vxsize+1)
        VYBins = np.linspace(vymin*gridratio,vymax*gridratio,num=vysize+1)            
        
        # Read velocity data into histogram
        (checkOk,binsXY,edgesX,edgesY) = vSpaceReducer(vlsvReader,cellid,slicetype,normvect,VXBins, VYBins,pop=pop,
                                                       slicethick=slicethick, wflux=wflux, cbulk=cbulk, 
                                                       center=center,keepfmin=keepfmin)

        # Check that data is ok and not empty
        if checkOk == False:
            print('ERROR: error from velocity space reducer')
            continue

        # Perform swap of coordinate axes, if requested
        if coordswap!=None:
            temp = edgesX
            edgesX = edgesY
            edgesY = temp
            temp = pltxstr
            pltxstr = pltystr
            pltystr = temp
            binsXY = binsXY.T

        # If no other fmin fmax values are given, take min and max of array
        if fmin!=None:
            fminuse=fmin
        else:
            nzindex = np.where(binsXY > 0)
            if np.any(nzindex):
                fminuse=np.amin(binsXY[nzindex])
            else:
                fminuse = 1e-15
        if fmax!=None:
            fmaxuse=fmax
        else:
            nzindex = np.where(binsXY > 0)
            if np.any(nzindex):
                fmaxuse=np.amax(binsXY[nzindex])
            else:
                fmaxuse = 1e-12

        if vlsvReader.check_variable('MinValue') == True:
            fMinFile = vlsvReader.read_variable('MinValue',cellid)
            print("Active f range is "+str(fminuse)+" to "+str(fmaxuse)+" with a vlsv file fMin value of "+str(fMinFile))
        else:
            print("Active f range is "+str(fminuse)+" to "+str(fmaxuse))

        norm = LogNorm(vmin=fminuse,vmax=fmaxuse)
        ticks = LogLocator(base=10,subs=range(10)) # where to show labels

        if box!=None:  # extents of plotted velocity grid as [x0,y0,x1,y1]
            xvalsrange=[box[0],box[1]]
            yvalsrange=[box[2],box[3]]
        else:
            # Find extent of nonzero data
            xindexrange = [vxsize,0]
            yindexrange = [vysize,0]
            for xi in range(len(edgesX)-1):
                for yi in range(len(edgesY)-1):
                    if binsXY[xi,yi] > 0:
                        xindexrange[0] = np.amin([xindexrange[0],xi])
                        xindexrange[1] = np.amax([xindexrange[1],xi])
                        yindexrange[0] = np.amin([yindexrange[0],yi])
                        yindexrange[1] = np.amax([yindexrange[1],yi])

            # leave some buffer
            xindexrange[0] =  np.max([0, 4 * int(np.floor((xindexrange[0]-2.)/4.)) ])
            xindexrange[1] =  np.min([len(edgesX)-1, 4 * int(np.ceil((xindexrange[1]+2.)/4.)) ])
            yindexrange[0] =  np.max([0, 4 * int((np.floor(yindexrange[0]-2.)/4.)) ])
            yindexrange[1] =  np.min([len(edgesY)-1, 4 * int(np.ceil((yindexrange[1]+2.)/4.)) ])

            # If empty VDF: plot whole v-space
            if ((xindexrange==[vxsize,0]) and (yindexrange==[vysize,0])):
                xindexrange = [0,vxsize]
                yindexrange = [0,vysize]

            xvalsrange = [ edgesX[xindexrange[0]] , edgesX[xindexrange[1]] ]
            yvalsrange = [ edgesY[yindexrange[0]] , edgesY[yindexrange[1]] ]

            # TODO make plot area square if it's almost square?

    
        # Plot the slice         
        [XmeshXY,YmeshXY] = scipy.meshgrid(edgesX/velUnit,edgesY/velUnit) # Generates the mesh to map the data to
        axis.pcolormesh(XmeshXY,YmeshXY,binsXY, cmap=colormap,norm=norm)

        plt.xlim([val/velUnit for val in yvalsrange])
        plt.ylim([val/velUnit for val in xvalsrange])
        axis.set_aspect('equal')

        # Grid
        plt.grid(color='grey',linestyle='-',linewidth=0.5)
        plt.xticks(np.arange(xvalsrange[0]/velUnit,xvalsrange[1]/velUnit,1e3))
        plt.yticks(np.arange(yvalsrange[0]/velUnit,yvalsrange[1]/velUnit,1e3))
 #       plt.minorticks_on()

        for axiss in ['top','bottom','left','right']:
            axis.spines[axiss].set_linewidth(0.5)

        axis.xaxis.set_tick_params(width=0.5,length=4)
        axis.yaxis.set_tick_params(width=0.5,length=4)
 #       axis.xaxis.set_tick_params(which='minor',width=0.5*0.8,length=2)
 #       axis.yaxis.set_tick_params(which='minor',width=0.5*0.8,length=2)

        if notitle==None:
            # Title and plot limits
            axis.set_title(plot_title,fontsize=fontsize2,fontweight='bold')

        if noxlabels==None:
            plt.xlabel(pltxstr,fontsize=5,weight='black')
            plt.xticks(fontsize=5,fontweight='black')
            axis.xaxis.offsetText.set_fontsize(5)
        if noylabels==None:
            plt.ylabel(pltystr,fontsize=5,weight='black')
            plt.yticks(fontsize=5,fontweight='black')
            axis.yaxis.offsetText.set_fontsize(5)

        if biglabel!=None:
            if biglabloc==None:
                biglabloc = 0 # default top-left corner
                
            if biglabloc==0:
                BLcoords=[0.02,0.98]
                BLha = "left"
                BLva = "top"
            elif biglabloc==1:
                BLcoords=[0.98,0.98]
                BLha = "right"
                BLva = "top"
            elif biglabloc==2:
                BLcoords=[0.98,0.02]
                BLha = "right"
                BLva = "bottom"
            elif biglabloc==3:
                BLcoords=[0.02,0.02]
                BLha = "left"
                BLva = "bottom"

            plt.text(BLcoords[0],BLcoords[1],biglabel, fontsize=fontsize4,weight='black', transform=axis.transAxes, ha=BLha, va=BLva)
                
        if cbar!=None:
            # Colourbar title
            cb_title_locy = 1.0 + 0.03/ratio
            if notitle==None:
                # in noborder mode, leave out colourbar descriptor
                plt.text(1.0, cb_title_locy, r"$f(v)$", fontsize=fontsize,weight='black', transform=axis.transAxes)

            # Witchcraft used to place colourbar
            divider = make_axes_locatable(axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            # First draw colorbar
            cb = plt.colorbar(fig1,ticks=ticks,cax=cax)
            cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
            cb.outline.set_linewidth(thick)

            # # if too many subticks:
            # # For non-square pictures, adjust tick count
            # nlabels = len(cb.ax.yaxis.get_ticklabels()) / ratio
            # if nlabels > 10:
            #     valids = ['1','2','3','4','5','6','8']
            # if nlabels > 19:
            #     valids = ['1','2','5']
            # if nlabels > 28:
            #     valids = ['1']
            # # for label in cb.ax.yaxis.get_ticklabels()[::labelincrement]:
            # if nlabels > 10:
            #     for label in cb.ax.yaxis.get_ticklabels():
            #         # labels will be in format $x.0\times10^{y}$
            #         if not label.get_text()[1] in valids:
            #             label.set_visible(False)

        if noxlabels!=None:
            for label in axis.xaxis.get_ticklabels():
                label.set_visible(False)
        if noylabels!=None:
            for label in axis.yaxis.get_ticklabels():
                label.set_visible(False) 

        

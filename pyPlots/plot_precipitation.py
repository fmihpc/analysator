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
from multiprocessing import Pool

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

global filedir_global
global cellid_global
global emin_global, emax_global
global alph0_global, pop_global, hemisphere_global

# Different style scientific format for colour bar ticks
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)

def loss_cone_angle(cellcoord=None,cellcoordre=None,deg=False):
    ''' Calculates the value of the loss cone angle at a given location
        :kword cellcoord:    The coordinates (X,Y,Z) of the cell whose loss cone angle value is calculated [in m]
        :kword cellcoordre:  The coordinates (X,Y,Z) of the cell whose loss cone angle value is calculated [in Re]
        :kword deg:          True if user wants the angle in degrees

        :returns:            The value of the loss cone angle

        .. code-blocks:: python
    '''

    # Earth radius [m]
    Re = 6.371e6

    # Convert coordinates to Re if needed
    if cellcoord!=None:
        X = cellcoord[0]/Re
        Y = cellcoord[1]/Re
        Z = cellcoord[2]/Re

    else:
        X = cellcoordre[0]
        Y = cellcoordre[1]
        Z = cellcoordre[2]

    # Calculation of R and L
    R = np.sqrt(X**2+Y**2+Z**2)                   # Radial distance to Earth centre
    if np.sqrt(X**2+Y**2)!=0:                     # Magnetic latitude
        lat_m = np.arctan(Z/np.sqrt(X**2+Y**2))
    else:
        lat_m = np.sign(Z)*np.pi/2.
    L = R/np.cos(lat_m)**2                        # L-shell


    # Analytical formula for loss cone angle (dipole approximation)
    alph0 = np.arcsin(R**-1.5 * (4*L-3*R)**.25/(4*L-3.)**.25)

    # Conversion to degrees if needed
    if deg:
        alph0 = alph0*180./np.pi

    return alph0



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
#    print("Found "+str(len(V))+" v-space cells")

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
#        print("Dropping velocity cells under fMin value "+str(fMin))
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
        # Find velocity components in give nframe (e.g. B frame: (vx,vy,vz) -> (vperp2,vperp1,vpar))
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



def precipitation_spectrum(vlsvReader=None,cid=None,losscone=None,pop=None,emin=None,emax=None,hemisphere=None):
    
    # check if velocity space exists in this cell
    if vlsvReader.check_variable('fSaved'): #restart files will not have this value
        if vlsvReader.read_variable('fSaved',cid) != 1.0:
            return (False,0,0)

    # Assume velocity cells are cubes
    [vxsize, vysize, vzsize] = vlsvReader.get_velocity_mesh_size(pop=pop)
    # Account for 4x4x4 cells per block
    vxsize = 4*vxsize
    vysize = 4*vysize
    vzsize = 4*vzsize
    [vxmin, vymin, vzmin, vxmax, vymax, vzmax] = vlsvReader.get_velocity_mesh_extent(pop=pop)
    inputcellsize=(vxmax-vxmin)/vxsize

    velcells = vlsvReader.read_velocity_cells(cid, pop=pop)
    V = vlsvReader.get_velocity_cell_coordinates(velcells.keys(), pop=pop)
#    print("Found "+str(len(V))+" v-space cells")

    f = zip(*velcells.items())
    # check that velocity space has cells
    if(len(f) > 0):
        f = np.asarray(zip(*velcells.items())[1])
    else:
        return (False,0,0)


    # Drop all velocity cells which are below the sparsity threshold. Otherwise the plot will show buffer cells as well.
    fMin = 1e-16 # default
    if vlsvReader.check_variable('MinValue') == True:
        fMin = vlsvReader.read_variable('MinValue',cid)
    ii_f = np.where(f >= fMin)
    print("Dropping velocity cells under fMin value "+str(fMin))
    if len(ii_f) < 1:
        return (False,0,0)
    f = f[ii_f]
    V = V[ii_f,:][0,:,:]


    # Rotate based on B-vector
    if vlsvReader.check_variable("B"):
        Bvect = vlsvReader.read_variable("B", cid)
    elif (vlsvReader.check_variable("background_B") and vlsvReader.check_variable("perturbed_B")):
        # used e.g. for restart files
        BGB = vlsvReader.read_variable("background_B", cid)
        PERBB = vlsvReader.read_variable("perturbed_B", cid)
        Bvect = BGB+PERBB
    else:
        print("Error finding B vector direction!")
        exit()
    
    if Bvect.shape==(1,3):
        Bvect = Bvect[0]
    normvect = Bvect
    
    normvect = np.array(normvect)
    normvect = normvect/np.linalg.norm(normvect)

    # Performing the rotation
    N = np.array(normvect)/np.sqrt(normvect[0]**2 + normvect[1]**2 + normvect[2]**2)
    Vrot = rotateVectorToVector(V,N)
    VX = Vrot[:,2]
    VY = Vrot[:,1]
    Vpara = Vrot[:,0]

    # Calculate pitch angles
    pitchangle = 0.*VX
    pitchangle[VX==0.] = np.pi/2.
    pitchangle[VX>0.] = np.arctan(np.sqrt(VY[VX>0.]**2+Vpara[VX>0.]**2)/VX[VX>0.])
    pitchangle[VX<0.] = np.pi-np.arctan(np.sqrt(VY[VX<0.]**2+Vpara[VX<0.]**2)/np.absolute(VX[VX<0.]))

    # Constants
    amu = 1.660539e-27 # kg
    if pop=='proton' or pop=='avgs':
        mass = amu
    elif pop=='helium':
        mass = 4.0026*amu
    elif pop=='oxygen':
        mass = 15.999*amu
    
    qe = 1.602177e-19 # C

    # Calculate energy
    normV = np.sqrt(VX**2+VY**2+Vpara**2)
    energy = 0.5*mass*normV**2/qe/1e3 # keV


    # Make the precipitating particle energy flux spectrum
    # -- limitation: low energies are not resolved due to very small pitch angles and finite velocity grid
    energy_bins = np.array([])
    fluxes = np.array([])

    deltaV = 50.e3 # m/s
    for vel in np.arange(vxmin,vxmax,deltaV):

        # Energy in the middle of the bin
        energy_i = 0.5*mass*(vel+deltaV/2.)**2/qe/1e3

        # Calculate corresponding flux
        solid_angle = 2*np.pi*(1.-np.absolute(np.cos(losscone)))
        deltaE = 0.5*mass*((vel+deltaV)**2-vel**2)/qe/1.e3

        # -- Collect precipitating particles within energy range
        #    -- If cell in southern hemisphere, then losscone = pi - losscone
        if hemisphere=='south':
            indexes = [((pitchangle >= np.pi-losscone) * (normV >= vel) * (normV < vel+deltaV))]
        else:
            indexes = [((pitchangle <= losscone) * (normV >= vel) * (normV < vel+deltaV))]

        # -- Calculation of integral value (for now assuming gyrotropy and low dependence on pitch angle inside loss cone)
        integral = (vel+deltaV/2.)**3*np.sum(f[indexes])*solid_angle

        # -- Put epsilon value if nothing in the loss cone [TODO maybe improve this later]
        if integral==0.:
            integral = 1.e-10

        # -- normalisation to angle and energy
        flux_i = integral/solid_angle/deltaE

        # -- propagation to 1 Re (geometric ratio) [TODO check]
        flux_i = flux_i / np.sin(losscone)**2

        # -- conversion to part / cm2 / s / sr / eV
        flux_i = flux_i*1e-7

        # Append to the output arrays
        energy_bins = np.append(energy_bins,energy_i)
        fluxes = np.append(fluxes,flux_i)
        
  #  (fluxes,energy_bins) = np.histogram(a=energy[indexes],bins=np.logspace(np.log10(emin),np.log10(emax),50),weights=f[indexes])


    return (True,energy_bins,fluxes)


def plot_prec_spectrum(filename=None,
                     vlsvobj=None,pop="proton",
                     filedir=None, step=None,
                     outputdir=None,
                     title=None,
                     draw=None, usesci=None,
                     symlog=None,
                     colormap=None,
                     run=None,notime=None,wmark=None,
                     notre=None, thick=1.0, cbtitle=None,
                     vmin=None, vmax=None, lin=None,
                     fmin=None, fmax=None, cbulk=None,
                     emin=None, emax=None,
                     cellcoordplot=None,cellidplot=None
                     ):    

    ''' Plots a precipitating ion flux energy spectrum.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.

    :kword colormap:    colour scale for plot, use e.g. jet, viridis, plasma, inferno, magma, nipy_spectral, RdBu
    :kword run:         run identifier, used for some default vmin,vmax values and for constructing output filename
    :kword notime:      flag to suppress plotting simulation time in title
    :kword title:       string to use as title instead of map name
    :kword cbtitle:     string to use as colour bar title instead of map name
    :kword notre:       flag to use metres (if ==1) or kilometres as axis unit
    :kword thick:       line and axis thickness, default=1.0
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: 1)
    :kword fmin,fmax:   min and max values for colour scale and colour bar of the VDFs. If no values are given,
                        min and max values for whole plot are used.
    :kword cbulk:       Center plot on position of bulk velocity (for this population)
    :kword lin:         flag for using linear colour scaling instead of log
    :kword symlog:      use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword draw:        Draw image on-screen instead of saving to file (requires x-windowing)

    :kword cellcoordplot:   Coordinates of cells to display as circles in the colormap plot, format [x1,y1,z1,...,xn,yn,zn]
    :kword cellidplot:      List of cellIDs to display as circles in the colormap plot
                            
    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:


    '''

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

    # If population isn't defined i.e. defaults to protons, check if 
    # instead should use old version "avgs"
    if pop=="proton":
       if not f.check_population(pop):
           if f.check_population("avgs"):
               pop="avgs"
               print("Auto-switched to population avgs")
           else:
               print("Unable to detect population "+pop+" in .vlsv file!")
               exit()
    else:
        if not f.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            exit()  


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
            plot_title = "t="+str(np.int(timeval))+' s'
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

    # Output file name
    savefigname = outputdir+run+"_prec_spec"+stepstr+"_"+str(cellidplot)+".png"

    Re = 6.371e+6 # Earth radius in m
    #read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")

    # Check if ecliptic or polar run
    if ysize==1:
        simext=[xmin,xmax,zmin,zmax]
        sizes=[xsize,zsize]
    if zsize==1:
        simext=[xmin,xmax,ymin,ymax]
        sizes=[xsize,ysize]

    vminuse=vmin
    vmaxuse=vmax

    # Lin or log colour scaling, defaults to log
    if lin==None:
        # Special SymLogNorm case
        if symlog!=None:
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

    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw!=None:
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')  


    # Select image shape to match plotted area, at least somewhat.
    # default for square figure is figsize=[4.0,3.15]
    figsize = [4.0,3.15]

    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Get coordinates if cellIDs were given as input
    if cellidplot==None:
        print("ERROR: No cell ID given as input")
        return
    else:
        for cellid in cellidplot:
            xCid,yCid,zCid = f.get_cell_coordinates(cellid)

            # Check whether northern or southern hemisphere
            # -- equatorial plane considered northern
            if zCid < 0.:
                hemisphere='south'
            else:
                hemisphere='north'

            # Calculation of loss cone angle value
            alph0 = loss_cone_angle(cellcoord=[xCid,yCid,zCid],deg=False)

            print("alph0 = "+str(alph0))


            # Reduction of the precipitating particle data
            (wentFine,energy,flux) = precipitation_spectrum(vlsvReader=f,cid=cellid,losscone=alph0,pop=pop,emin=emin,emax=emax,hemisphere=hemisphere)

    # Plots the histogram
    if wentFine:
        fig1=plt.loglog(energy,flux,'o-')
        print("flux = "+str(flux)+"\n energy = "+str(energy))
    else:
        print("There was a problem making the histogram")
        return

    ax1 = plt.gca() # get current axes

    # Title and plot limits
    ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')
    plt.xlim([emin,emax])
    plt.ylim([fmin,fmax])
#    ax1.set_aspect('equal')

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(thick)
    ax1.xaxis.set_tick_params(width=thick,length=3)
    ax1.yaxis.set_tick_params(width=thick,length=3)


    # Limit ticks, slightly according to ratio
#    ax1.xaxis.set_major_locator(plt.MaxNLocator(int(7/np.sqrt(ratio))))
#    ax1.yaxis.set_major_locator(plt.MaxNLocator(int(7*np.sqrt(ratio))))

    plt.xlabel('E [keV]',fontsize=fontsize,weight='black')
    plt.ylabel('Flux [protons/cm2/s/sr/eV]',fontsize=fontsize,weight='black')

    plt.xticks(fontsize=fontsize,fontweight='black')
    plt.yticks(fontsize=fontsize,fontweight='black')

    # Grid
    plt.grid(color='grey',linestyle='-')

    # set axis exponent offset font sizes
    ax1.yaxis.offsetText.set_fontsize(fontsize)
    ax1.xaxis.offsetText.set_fontsize(fontsize)
          
    # Add Vlasiator watermark
    if wmark!=None:        
        wm = plt.imread(get_sample_data(watermarkimage))
        newax = fig.add_axes([0.01, 0.90, 0.3, 0.08], anchor='NW', zorder=-1)
        newax.imshow(wm)
        newax.axis('off')

    # adjust layout
    plt.tight_layout()


    # Save output or draw on-screen
    if draw==None:
        print(savefigname+"\n")
        plt.savefig(savefigname,dpi=300)
    else:
        plt.draw()
        plt.show()
    plt.close()
    plt.clf()




def make_keogram_column(step=None):
    filename = filedir_global+'bulk.'+str(step).rjust(7,'0')+'.vlsv'

#    print(filename+" is being processed")
    f=pt.vlsvfile.VlsvReader(filename)

    # Reduction of the precipitating particle data
    (wentFine,energy,flux) = precipitation_spectrum(vlsvReader=f,cid=cellid_global,losscone=alph0_global,pop=pop_global,
                                                    emin=emin_global,emax=emax_global,hemisphere=hemisphere_global)

#    datamap = np.append(datamap,flux)
#    time_keogram = np.append(time_keogram,f.read_parameter("time"))
    out = [f.read_parameter("time"),flux,energy]

    if not wentFine:
        print("There was a problem making the spectrum, filename: "+filename)

    return (out)




def plot_prec_time_spectrum(filedir=None,
                     pop="proton",
                     start=None, stop=None,
                     outputdir=None,
                     title=None,
                     draw=None, usesci=None,
                     symlog=None,
                     colormap=None,
                     run=None,notime=None,wmark=None,
                     notre=None, thick=1.0, cbtitle=None,
                     vmin=None, vmax=None, lin=None,
                     fmin=None, fmax=None, cbulk=None,
                     emin=None, emax=None,
                     cellcoordplot=None,cellidplot=None,
                     numproc=8
                     ):    

    ''' Plots a precipitating ion flux energy spectrum during time interval (colour plot).

    :kword filedir:     Directory where bulk files are located
    :kword start:       Step for starting the plot
    :kword stop:        Step for ending the plot
    :kword outputdir:   Path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.

    :kword colormap:    colour scale for plot, use e.g. jet, viridis, plasma, inferno, magma, nipy_spectral, RdBu
    :kword run:         run identifier, used for some default vmin,vmax values and for constructing output filename
    :kword notime:      flag to suppress plotting simulation time in title
    :kword title:       string to use as title instead of map name
    :kword cbtitle:     string to use as colour bar title instead of map name
    :kword notre:       flag to use metres (if ==1) or kilometres as axis unit
    :kword thick:       line and axis thickness, default=1.0
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: 1)
    :kword fmin,fmax:   min and max values for colour scale and colour bar of the VDFs. If no values are given,
                        min and max values for whole plot are used.
    :kword cbulk:       Center plot on position of bulk velocity (for this population)
    :kword lin:         flag for using linear colour scaling instead of log
    :kword symlog:      use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword draw:        Draw image on-screen instead of saving to file (requires x-windowing)

    :kword cellcoordplot:   Coordinates of cells to display as circles in the colormap plot, format [x1,y1,z1,...,xn,yn,zn]
    :kword cellidplot:      List of cellIDs to display as circles in the colormap plot
                            
    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:


    '''

    global filedir_global
    global cellid_global
    global emin_global, emax_global
    global alph0_global, pop_global, hemisphere_global

    filedir_global=filedir
    emin_global=emin
    emax_global=emax

    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    # watermarkimage=os.path.expandvars('$HOME/appl_taito/analysator/pyPlot/logo_color.png')
    # watermarkimage='/homeappl/home/marbat/appl_taito/analysator/logo_color.png'

    if outputdir==None:
        outputdir=os.path.expandvars('$HOME/Plots/')
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    # Input file list
    if ((filedir!=None) and (start!=None) and (stop!=None)):
        filelist = []
        for step in range(start,stop):
            filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
            filelist.append(filename)
    else:
        print("ERROR: needs a bulk file directory and start/stop steps")
        return

    # Scientific notation for colorbar ticks?
    if usesci==None:
        usesci=1
    
    if colormap==None:
        colormap="YlOrRd"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8 # Most text
    fontsize2=10 # Time title
    fontsize3=5 # Colour bar ticks


    # Stop and start times
    f=pt.vlsvfile.VlsvReader(filelist[-1])
    tstop=None
    tstop=f.read_parameter("time")
    if tstop==None:
        tstop=f.read_parameter("t")
    if tstop==None:
        print "Unknown time format encountered for stop file"
        return

    f=pt.vlsvfile.VlsvReader(filelist[0])
    tstart=None
    tstart=f.read_parameter("time")
    if tstart==None:
        tstart=f.read_parameter("t")
    if tstart==None:
        print "Unknown time format encountered for start file"
        return


    # Plot title with time
    if notime==None:        
        timeval=f.read_parameter("time")
        if timeval == None:
            timeval=f.read_parameter("t")
            if timeval == None:    
                print "Unknown time format encountered"
                plot_title = ''
        if timeval != None:
            plot_title = "t="+str(np.int(timeval))+' s'
    else:
        plot_title = ''       

    # stepstr, used for file name
    if start!=None:
        stepstr = '_t'+str(np.int(tstart))+'s'
    else:
        stepstr = ''

    # If run name isn't given, just put "plot" in the output file name
    if run==None:
        run='plot'

    # Output file name
    savefigname = outputdir+run+"_prec_spec_time"+stepstr+"_"+str(cellidplot)+".png"


    # If population isn't defined i.e. defaults to protons, check if 
    # instead should use old version "avgs"
    if pop=="proton":
       if not f.check_population(pop):
           if f.check_population("avgs"):
               pop="avgs"
               print("Auto-switched to population avgs")
           else:
               print("Unable to detect population "+pop+" in .vlsv file!")
               exit()
    else:
        if not f.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            exit() 
    pop_global = pop 

    Re = 6.371e+6 # Earth radius in m
    #read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")

    # Check if ecliptic or polar run
    if ysize==1:
        simext=[xmin,xmax,zmin,zmax]
        sizes=[xsize,zsize]
    if zsize==1:
        simext=[xmin,xmax,ymin,ymax]
        sizes=[xsize,ysize]

    vminuse=fmin
    vmaxuse=fmax

    # Lin or log colour scaling, defaults to log
    if lin==None:
        # Special SymLogNorm case
        if symlog!=None:
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

    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw!=None:
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')  

    # Select window to draw
    boxcoords = [tstart,tstop,emin,emax]

    # Select image shape to match plotted area, at least somewhat.
    # default for square figure is figsize=[4.0,3.15]
    figsize = [4.0,3.15]
    ratio = 1.

    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Get coordinates if cellIDs were given as input
    if cellidplot==None:
        print("ERROR: No cell ID given as input")
        return
    else:
        for cellid in cellidplot:
            cellid_global = cellid
            xCid,yCid,zCid = f.get_cell_coordinates(cellid)

            # Check whether northern or southern hemisphere
            # -- equatorial plane considered northern
            if zCid < 0.:
                hemisphere='south'
            else:
                hemisphere='north'
            hemisphere_global = hemisphere

            # Calculation of loss cone angle value
            alph0 = loss_cone_angle(cellcoord=[xCid,yCid,zCid],deg=False)
            alph0_global = alph0

            print("alph0 = "+str(alph0))

            # Building the datamap corresponding to the keogram    
            datamap = np.array([])
            time_keogram = np.array([])

            print("Name: "+__name__)

            # Parallel construction of the spectrum keogram
            if __name__ == 'plot_precipitation':
                pool = Pool(numproc)
                return_array = pool.map(make_keogram_column, range(start,stop))
            else:
                print("didn't enter the loop")

            energy = return_array[0][2]

            for j in return_array:
                time_keogram = np.append(time_keogram,j[0])
                datamap = np.append(datamap,j[1])

            # Determine sizes of the keogram array
            sizes=[time_keogram.size,energy.size]


            # Reshape data to an ordered 2D array that can be plotted
            if np.ndim(datamap) != 2:
                datamap = datamap.reshape([sizes[0],sizes[1]])

            datamap = np.transpose(datamap)

    print("Now making the plot")


    # Generates the mesh to map the data to.
    # Note, datamap is still of shape [ysize,xsize] (?)
    [XmeshXY,YmeshXY] = scipy.meshgrid(time_keogram,energy)

    fig1 = plt.pcolormesh(XmeshXY,YmeshXY,datamap, cmap=colormap,norm=norm)
    ax1 = plt.gca() # get current axes

    # Save the keogram to file in case further processing is needed
    datamap.dump(outputdir+'datamap')
    np.save(outputdir+'time_keogram',time_keogram)
    np.save(outputdir+'energy_scale',energy)


    ax1.set_yscale("log")

    # Title and plot limits
#    ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')
    plt.xlim([time_keogram[0],time_keogram[-1]])
    plt.ylim([emin,emax])

    plt.xlabel('Time [s]',fontsize=fontsize,weight='black')
    plt.ylabel('E [keV]',fontsize=fontsize,weight='black')

    plt.xticks(fontsize=fontsize,fontweight='black')
    plt.yticks(fontsize=fontsize,fontweight='black')

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(thick)
    ax1.xaxis.set_tick_params(width=thick,length=3)
    ax1.yaxis.set_tick_params(width=thick,length=3)


    # Limit ticks, slightly according to ratio
#    ax1.xaxis.set_major_locator(plt.MaxNLocator(int(7/np.sqrt(ratio))))
#    ax1.yaxis.set_major_locator(plt.MaxNLocator(int(7*np.sqrt(ratio))))

    # Colourbar title
    cbtitle = 'Flux\n [ion/cm$^2$/s/sr/eV]'
    plt.text(1.0, 1.04, cbtitle, fontsize=fontsize3+1,weight='black', transform=ax1.transAxes)

    # Witchcraft used to place colourbar
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
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

    # Grid
#    plt.grid(color='grey',linestyle='-')

    # set axis exponent offset font sizes
    ax1.yaxis.offsetText.set_fontsize(fontsize)
    ax1.xaxis.offsetText.set_fontsize(fontsize)
          
    # Add Vlasiator watermark
    if wmark!=None:        
        wm = plt.imread(get_sample_data(watermarkimage))
        newax = fig.add_axes([0.01, 0.90, 0.3, 0.08], anchor='NW', zorder=-1)
        newax.imshow(wm)
        newax.axis('off')

    # adjust layout
    plt.tight_layout(pad=1.5)


    # Save output or draw on-screen
    if draw==None:
        print(savefigname+"\n")
        plt.savefig(savefigname,dpi=300)
    else:
        plt.draw()
        plt.show()
    plt.close()
    plt.clf()

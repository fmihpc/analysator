import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
import matplotlib.ticker as mtick
import colormaps as cmaps
from matplotlib.cbook import get_sample_data
from multiprocessing import Pool

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

global filedir_global, filetype_global 
global cellid_global
global emin_global, emax_global, enum_global
global pop_global 


# Different style scientific format for colour bar ticks
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)


def energy_spectrum(vlsvReader, cid, pop, emin, emax, enum=10, fluxout=False):
    ''' Calculates the energy spectrum of a single file at a single cellid
        
        param: vlsvReader         vlsvReader handle
        param: cid                Cellid number to be used (int)
        param: pop                Name of population (string)
        param: emin               Minimum energy for spectrum [keV] (float)
        param: emax               Maximum energy for spectrum [keV] (float)
        param: enum               Number of energy bins for spectrum (int)
        param: fluxout            If true returns fluxes instead of PSD (bool)
        returns:                  Energy bin centers and PSD/flux for each bin (arrays)
    '''

    # check if velocity space exists in this cell
    if vlsvReader.check_variable('fSaved'): #restart files will not have this value
        if vlsvReader.read_variable('fSaved',cid) != 1.0:
            print('Velocity space not found for this cellID!')
            return (False,0,0)

    # Get velocity data
    velcells = vlsvReader.read_velocity_cells(cid, pop=pop)
    V = vlsvReader.get_velocity_cell_coordinates(velcells.keys(), pop=pop)

    # check that velocity space has cells
    if(len(velcells) > 0):
        f = np.asarray(velcells.values())
    else:
        print('Velocity space cells empty!')
        return (False,0,0)

    # Drop all velocity cells which are below the sparsity threshold. Otherwise the plot will show buffer cells as well.
    if pop=='electron':
        fMin = 1e-21 # default
    else:
        fMin = 1e-16 # default
    if vlsvReader.check_variable('MinValue') == True:
        fMin = vlsvReader.read_variable('MinValue',cid)
    ii_f = np.where(f >= fMin)
    if len(ii_f[0]) < 1:
        print('No velocity cells found above threshold: '+str(fMin))
        return (False,0,0)
    f_sparse = f[ii_f]
    V_sparse = V[ii_f,:][0,:,:]

    # Constants
    amu = 1.660539e-27 # kg
    if pop=='proton' or pop=='avgs':
        mass = amu
    elif pop=='helium':
        mass = 4.0026*amu
    elif pop=='oxygen':
        mass = 15.999*amu
    elif pop=='electron':
        mass = 9.10938e-31 # kg
    else:
        print('Population not known! mass needs to be assigned in function energy_spectrum')
        return(False, 0,0)

    qe = 1.602177e-19 # C

    # Calculate energy
    VX = V_sparse[:,0]
    VY = V_sparse[:,1]
    VZ = V_sparse[:,2]
    normV = np.sqrt(VX**2+VY**2+VZ**2)
    cell_energy = 0.5*mass*normV**2/qe/1.e3 # keV (Non-relativistic)

    # Calculate particle energy flux/PSD spectrum
    try:
        energy_bin_edges = np.logspace(np.log2(emin), np.log2(emax), num=enum+1, base=2.)
    except RuntimeWarning:
        print('emin and emax have to be positive numbers!')
    except ValueError:
        print('emin and emax have to be positive numbers!')
    dataout = np.ones(len(energy_bin_edges)-1)*1.E-30
    energyout = np.zeros_like(dataout)

    for c,el in enumerate(energy_bin_edges[:-1]):
        # Energy in the middle of the bin [keV]
        energy_i = (energy_bin_edges[c] + energy_bin_edges[c+1]) / 2. 
        # Velocity in the middle of the bin (m/s)
        vel_i = np.sqrt(2.*qe*1.e3*energy_i/mass)

        shellmask = (cell_energy > el) & (cell_energy <= energy_bin_edges[c+1])
        
        if all(shellmask == False):
            print('WARNING: No cells found between bin edges ' + str(energy_bin_edges[c]) + ' and ' + str(energy_bin_edges[c+1]) )
            print('Output for this channel will be set to 1.E-30! Consider adjusting the energy settings.')
        else:
            if fluxout:
                # Flux (what is measured by spacecraft) in particles/(cm2 s sr eV)
                dataout[c] = vel_i**2/mass*np.mean(f_sparse[shellmask])*1.e-4*qe
            else:
                # PSD average
                dataout[c] = np.mean(f_sparse[shellmask])

        energyout[c] = energy_i

    if all(dataout==1.E-30):
        print('WARNING: No cells found between the energy limits. Consider increasing them!')
        print('Max, min energy = ' + str(max(cell_energy)) + ', ' + str(min(cell_energy))) 

    return (True, energyout, dataout )


def make_timemap(step):
    ''' Auxiliary function to be used for parallelisation of the time vs. energy spectrum plot

        param: step       Time step of file to be processed 
        returns:          Time of simulation [s], energy of bins [keV] and PSD/flux for each bin

    '''
    # Getting file handle
    filename = filedir_global+filetype_global+'.'+str(step).rjust(7,'0')+'.vlsv' 
    f = pt.vlsvfile.VlsvReader(filename)
    print(filename + " is being processed...")

    # Getting energy spectrum data
    (success, energy, particledata) = energy_spectrum(f, cellid_global, pop_global, emin_global, emax_global, enum=enum_global)

    time = f.read_parameter("time")
    if time is None:      # in BCH, at some point "t" was changed to "time"
        time = f.read_parameter("t")

    out = [time, particledata, energy]

    if success == False:
        sys.exit("There was a problem making the spectrum, filename: "+filename)

    return (out)

def get_energy_spectrum(filedir, filetype, pop, start, stop, cid, emin, emax, enum=16, fluxout=False, numproc=8):    

    ''' Outputs data arrays of time/energy spectrum during time interval.

    :param filedir:         Directory where files are located
    :param filetype:        Type of file to be used [example: 'bulk']
    :param pop:             Name of population
    :param start:           Step for starting the data 
    :param stop:            Step for ending the data
    :param cid:             cellID number
    :param emin,emax,enum:  min and max values (edges) for energy levels [keV] and number of energy bins 
    :kword fluxout:         If True, outputs fluxes instead of PSD (Default: False)
    :kword numproc:         Number of processes for parallelisation (default: 8)

    :returns:               Tuple (time, energy, datamap) containing the time data [1-D array], the energy channels [1-D array]
                            and the data output [2-D array] as PSD or flux.

    
    # Example usage:
    import pytools as pt
    data = pt.plot.get_energy_spectrum('./', 'bulk', 'proton', 0, 100, 10000, 1, 100, enum=12)

    '''

    global filedir_global, filetype_global, cellid_global, pop_global
    global emin_global, emax_global, enum_global

    # check if should use old version "avgs"
    filename = filedir+filetype+'.'+str(start).rjust(7,'0')+'.vlsv' 
    f=pt.vlsvfile.VlsvReader(filename)
    if pop=="proton":
       if not f.check_population(pop):
           if f.check_population("avgs"):
               pop="avgs"
               print("Auto-switched to population avgs")
           else:
               print("Unable to detect population "+pop+" in .vlsv file!")
               sys.exit()
    else:
        if not f.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            sys.exit() 
    
    # TODO do not use global variables, check that variables are valid
    filedir_global = filedir
    filetype_global = filetype
    pop_global = pop
    emin_global = emin
    emax_global = emax
    enum_global = enum
    cellid_global = cid
    
    datamap = np.array([])
    time_ar = np.array([])

    # Parallel construction of the spectrum
    if __name__ == 'plot_energy':
        pool = Pool(numproc)
        return_array = pool.map(make_timemap, range(start,stop+1))
    else:
        print("Problem with parallelization")

    # Creating datamap
    for j in return_array:
        time_ar = np.append(time_ar,j[0])
        datamap = np.append(datamap,j[1])
    energy_ar = return_array[0][2]
    
    #Serial construction of spectrum (for debugging)
    #for step in range(start,stop+1):
    #    func_return = make_timemap(step)
    #    time_ar = np.append(time_ar,func_return[0])
    #    datamap = np.append(datamap,func_return[1])
    #energy_ar = func_return[2]

    # Reshape data to an ordered 2D array that can be plotted
    sizes=[time_ar.size,energy_ar.size]
    if np.ndim(datamap) != 2:
        datamap = datamap.reshape([sizes[0],sizes[1]])

    datamap = np.transpose(datamap)

    return (time_ar, energy_ar, datamap)


def plot_energy_spectrum(filedir=None, filetype='bulk',
                     pop='proton',
                     start=1, stop=20,
                     outputdir=None,
                     emin=None, emax=None, enum=None,
                     colormap=None,
                     title=None,
                     draw=None, usesci=1,
                     run=None, wmark=None,
                     notre=None, thick=1.0, cbtitle=None,
                     lin=None,
                     fmin=None, fmax=None,
                     cellcoordplot=None, cellidplot=None,
                     numproc=8
                     ):    

    ''' Plots a time/energy spectrum during time interval (colour plot).

    :kword filedir:         Directory where files are located
    :kword filetype:        Type of file to be used [default = 'bulk']
    :kword start:           File number (step) for starting the plot
    :kword stop:            File number (step) for ending the plot
    :kword outputdir:       Path to directory where output files are created (default: $HOME/Plots/)
                            If directory does not exist, it will be created. If the string does not end in a
                            forward slash, the final parti will be used as a perfix for the files.

    :kword emin,emax,enum:  min and max values (edges) for energy levels [keV] and number of energy bins 
    :kword colormap:        colour scale for plot, use e.g. jet, viridis, plasma, inferno, magma, nipy_spectral, RdBu
    :kword title:           string to use as title instead of map name
    :kword draw:            Draw image on-screen instead of saving to file (requires x-windowing)
    :kword usesci:          Use scientific notation for colorbar ticks? (default: 1)
    :kword run:             run identifier, used for some default vmin,vmax values and for constructing output filename
    :kword wmark:           If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword cbtitle:         string to use as colour bar title instead of map name
    :kword notre:           flag to use metres (if ==1) or kilometres as axis unit
    :kword thick:           line and axis thickness, default=1.0
    :kword lin:             flag for using linear colour scaling instead of log
    :kword fmin,fmax:       min and max values for colour scale and colour bar of the PSD. 

    :kword cellcoordplot:   Coordinates of cell to be plottedd (3-D array)
    :kword cellidplot:      cellID be plotted (list)
    :kword numproc:         Number of processes for parallelisation (default: 8)

    :returns:               Outputs an image to a file or to the screen.

    
    # Example usage:
    import pytools as pt
    pt.plot.plot_energy_spectrum(filedir='./', start=0, stop=100, emin=1, emax=50, enum=16, cellidplot=[300100])

    '''

    global filedir_global, filetype_global, cellid_global, pop_global
    global emin_global, emax_global, enum_global

    # TODO do not use global variables, check that variables are valid
    filedir_global=filedir
    filetype_global=filetype
    emin_global=emin
    emax_global=emax
    enum_global=enum

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
        for step in range(start,stop+1):
            filename = filedir+filetype+'.'+str(step).rjust(7,'0')+'.vlsv' 
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

    # TODO check if files have VDFs
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

    # Plot title
    plot_title = ''       

    # stepstr, used for file name
    stepstr = '_t'+str(np.int(tstart))+'s'

    # If run name isn't given, just put "plot" in the output file name
    if run==None:
        run='plot'

    # Output file name
    savefigname = outputdir+run+"_ener_time_spec"+stepstr+"_"+str(cellidplot[0])+".png"


    # If population isn't defined i.e. defaults to protons, check if 
    # instead should use old version "avgs"
    if pop=="proton":
       if not f.check_population(pop):
           if f.check_population("avgs"):
               pop="avgs"
               print("Auto-switched to population avgs")
           else:
               print("Unable to detect population "+pop+" in .vlsv file!")
               sys.exit()
    else:
        if not f.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            sys.exit() 
    pop_global = pop 

    #if fmin!=None and fmax!=None:
    # TODO how to handle fmin and fmax
    # Lin or log colour scaling, defaults to log
    if lin==None:
        norm = LogNorm(vmin=fmin,vmax=fmax)
        ticks = LogLocator(base=10,subs=range(10)) # where to show labels
    else:
        # Linear
        levels = MaxNLocator(nbins=255).tick_values(fmin,fmax)
        norm = BoundaryNorm(levels, ncolors=cmapuse.N, clip=True)
        ticks = np.linspace(fmin,fmax,num=7)            
   
    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw!=None:
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')  

    # Checking if CellID exists
    if cellidplot==None and cellcoordplot==None :
        print("ERROR: No cell ID or coordinate given as input")
        return

    #TODO find cellid from coordinates
    cellid_global = cellidplot
    datamap = np.array([])
    time_ar = np.array([])

    # Parallel construction of the spectrum
    if __name__ == 'plot_energy':
        pool = Pool(numproc)
        return_array = pool.map(make_timemap, range(start,stop+1))
    else:
        print("Problem with parallelization")

    # Creating datamap
    for j in return_array:
        time_ar = np.append(time_ar,j[0])
        datamap = np.append(datamap,j[1])
    energy = return_array[0][2]
    
    #Serial construction of spectrum (for debugging)
    #for step in range(start,stop+1):
    #    func_return = make_timemap(step)
    #    time_ar = np.append(time_ar,func_return[0])
    #    datamap = np.append(datamap,func_return[1])
    #energy = func_return[2]

    # Reshape data to an ordered 2D array that can be plotted
    sizes=[time_ar.size,energy.size]
    if np.ndim(datamap) != 2:
        datamap = datamap.reshape([sizes[0],sizes[1]])

    datamap = np.transpose(datamap)

    # Select figure size
    ratio=1. #TODO calculate this
    figsize = [4.0,3.15]

    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Generates the mesh to map the data to.
    # Note, datamap is still of shape [ysize,xsize] (?)
    [XmeshXY,YmeshXY] = scipy.meshgrid(time_ar,energy)

    fig1 = plt.pcolormesh(XmeshXY,YmeshXY,datamap, cmap=colormap,norm=norm)
    ax1 = plt.gca()

    ax1.set_yscale("log")

    # Title and plot limits
    #ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')
    plt.xlim([time_ar[0],time_ar[-1]])
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
    ax1.xaxis.set_major_locator(plt.MaxNLocator(int(7/np.sqrt(ratio))))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(int(7*np.sqrt(ratio))))

    # Colourbar title
    cbtitle = 'Flux\n [1/cm$^2$/s/sr/eV]'
    cbtitle = 'PSD [s$^3$/m$^6$]' # TODO: make this flexible
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
    if lin==None and usesci!=0:
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
    #plt.grid(color='grey',linestyle='-')

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

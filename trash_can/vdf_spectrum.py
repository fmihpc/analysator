import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os, sys

# Run TeX typesetting through the full TeX engine instead of python's own mathtext. Allows
# for changing fonts, bold math symbols etc, but may cause trouble on some systems.
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.rcParams['text.dvipnghack'] = 'True' # This hack might fix it on some systems

# find nearest spatial cell with vspace to cid
def getNearestCellWithVspace(vlsvReader,cid):
    cell_candidates = vlsvReader.read(mesh='SpatialGrid',tag='CELLSWITHBLOCKS')
    if len(cell_candidates)==0:
        print("Error: No velocity distributions found!")
        sys.exit()
    cell_candidate_coordinates = [vlsvReader.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]
    cell_coordinates = vlsvReader.get_cell_coordinates(cid)
    norms = np.sum((cell_candidate_coordinates - cell_coordinates)**2, axis=-1)**(1./2)
    norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))
    return cell_candidates[i]

# analyze velocity space in a spatial cell (velocity space reducer)
def vSpaceReduceToSpectrum(vlsvReader, cid, Bins, vmin, vmax, pop="proton", 
                  wflux=None, cbulk=None, center=None, keepfmin=None, energy=None):
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
                sys.exit()
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
            sys.exit()
        # shift to velocities plasma frame
        V = V - bulkv
    elif center!=None:
        if len(center)==3:
            print("Transforming to frame travelling at speed "+str(center))
            V = V - center
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

    # Correct to cell-averages
    V = V + 0.5*inputcellsize

    # Use enegy for binning instead of velocity?
    if energy!=None:
        mp = 1.672622e-27
        JtoeV = .16021773e-18
        V = 0.5 * mp * np.square(V) / JtoeV

    # Now take magnitude of velocity
    V = np.linalg.norm(V, axis=-1)

    # check extents
    if vmin==None:
        vmin = np.amin(V)
    if vmax==None:
        vmax = np.amax(V)

    # Flux weighting?
    if wflux!=None:
        fw = f*V # use particle flux as weighting in the histogram
    else:
        fw = f # use particle phase-space density as weighting in the histogram

    # Gather histogram of values
    (Vhist,VEdges) = np.histogram(V,bins=Bins, range=(vmin,vmax),weights=fw,normed=0)

    # Flux weighting
    if wflux!=None:
        dV = np.abs(VEdges[-1] - VEdges[-2]) # assumes constant bin size
        nVhist = np.divide(Vhist,(dV*4*np.pi)) # normalization

    return (True,Vhist,VEdges)
  


def plot_spectrum(filename=None,
             vlsvobj=None,
             filedir=None, step=None,
             cellids=None, pop="proton",
             coordinates=None, coordre=None, 
             outputdir=None,
             draw=None,unit=None,title=None,
             run=None, wmark=None, thick=1.0,
             fmin=None, fmax=None, energy=None,xlog=None,
             vmin=None, vmax=None, nBins=50,
             cbulk=None, center=None, wflux=None, keepfmin=None,
             noborder=None, scale=1.0,
             noxlabels=None, noylabels=None,
             ):

    ''' Plots a coloured plot with axes and a colour bar.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
     
    :kword cellids:     list of cell IDs to plot VDF for
    :kword coordinates: list of 3-element spatial coordinates to plot VDF for (given in metres)
    :kword coordre: list of 3-element spatial coordinates to plot VDF for (given in Earth radii)
    :kword pop:         Population to plot, default proton

    :kword run:         run identifier, used for constructing output filename
    :kword title:       string to use as plot title instead of time
    :kword fmin,fmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot are used.

    :kword unit:        Plot x-axis unit (in m/s, default: 1000 i.e. km/s)
    :kword energy:      Plot x-axis in energy instead of velocity. In this case, unit is eV, or multiple of that.
    :kword xlog:        Plot x-axis with log instead of lin scale
   
    :kword cbulk:       Center plot on position of total bulk velocity (or if not available,
                        bulk velocity for this population)
    :kword center:      Center plot on provided 3-element velocity vector position (in m/s)
    :kword wflux:       Plot flux instead of distribution function
    :kword keepfmin:    Also draw buffer cells with values below fMin
    :kword differential: Calculate differential flux (dE)

    :kword draw:        Draw image on-screen instead of saving to file (requires x-windowing)

    :kword noborder:    Plot figure edge-to-edge without borders (default off)
    :kword noxlabels:   Suppress x-axis labels and title
    :kword noylabels:   Suppress y-axis labels and title
    :kword scale:       Scale text size (default=1.0)
    :kword thick:       line and axis thickness, default=1.0
    
    :kword nBins:       How many bins to use in histogram (default: 50)
    :kword vmin:        min edge for histogram
    :kword vmax:        max edge for histogram

    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:

    '''

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
        vlsvReader=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        vlsvReader=vlsvobj
    elif ((filedir!=None) and (step!=None)):
        filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
        vlsvReader=pt.vlsvfile.VlsvReader(filename)
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=5*scale # Colour bar ticks
    fontsize4=16*scale # Big label

    # Plot title with time
    timeval=vlsvReader.read_parameter("time")
    if timeval==None:
        timeval=vlsvReader.read_parameter("t")

    if title==None:        
        if timeval == None:    
            plot_title = ''
            print "Unknown time format encountered"
        else:
            plot_title = "t="+str(np.int(timeval))+' s'
    else:
        plot_title = title

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
               sys.exit()
    else:
        if not vlsvReader.check_population(pop):
            print("Unable to detect population "+pop+" in .vlsv file!")
            sys.exit()       

    [xsize, ysize, zsize] = vlsvReader.get_spatial_mesh_size()

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
    if energy!=None:
        velUnitStr = '[eV]'
        velUnit = 1

    if unit!=None:
        velUnit = unit
        if energy==None:
            if unit==1:
                velUnitStr = r'[m/s]'
            else:
                velUnitStr = ''
        else:
            if unit==1000:
                velUnitStr = r'[keV]'
            else:
                velUnitStr = ''

    # Select ploitting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw!=None:
        plt.switch_backend('TkAgg')
    else:
        plt.switch_backend('Agg')  


    if (cellids==None and coordinates==None and coordre==None):
        print("Error: must provide either cell id's or coordinates")
        return -1

    if coordre!=None:
        # Transform to metres
        coordinates = Re*np.asarray(coordre)

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
            sys.exit()
        cidsTemp = []
        for ii in range(xReq.shape[0]):
            cidRequest = (np.int64)(vlsvReader.get_cellid(np.array([xReq[ii],yReq[ii],zReq[ii]])))
            cidNearestVspace = -1
            if cidRequest > 0:
                cidNearestVspace = getNearestCellWithVspace(vlsvReader,cidRequest)
            else:
                print('ERROR: cell not found')
                sys.exit()
            if (cidNearestVspace <= 0):
                print('ERROR: cell with vspace not found')
                sys.exit()
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
        savefigname = outputdir+outputprefix+run+"_vdf_"+pop+stepstr+"_cellid_"+str(cellid)+".png"

        # Read velocity data into histogram
        (checkOk,binsX,edgesX) = vSpaceReduceToSpectrum(vlsvReader,cellid,nBins,vmin,vmax,pop=pop,
                                                       wflux=wflux, cbulk=cbulk, energy=energy,
                                                       center=center,keepfmin=keepfmin)

        # Check that data is ok and not empty
        if checkOk == False:
            print('ERROR: error from velocity space reducer')
            continue

        # Create 300 dpi image of suitable size
        figsize = [4,4]
        fig = plt.figure(figsize=figsize,dpi=300)
    
        # Plot the spectrum
        ax1 = plt.gca() # get current axes
        # find centre values
        centresX = 0.5*(edgesX[1:]+edgesX[:-1])
        # X-axis unit
        centresX = centresX / velUnit
        ax1.plot(centresX, binsX)

        ax1.set_xlabel(velUnitStr, labelpad=1)
        ax1.set_yscale('log')
        if xlog!=None:
            ax1.set_xscale('log')

        # for axiss in ['top','bottom','left','right']:
        #     ax1.spines[axiss].set_linewidth(thick)

        # ax1.xaxis.set_tick_params(width=thick,length=4)
        # ax1.yaxis.set_tick_params(width=thick,length=4)
        # ax1.xaxis.set_tick_params(which='minor',width=thick*0.8,length=2)
        # ax1.yaxis.set_tick_params(which='minor',width=thick*0.8,length=2)

        ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')

        # if noxlabels==None:
        #     plt.xlabel(pltxstr,fontsize=fontsize,weight='black')
        #     plt.xticks(fontsize=fontsize,fontweight='black')
        #     ax1.xaxis.offsetText.set_fontsize(fontsize)
        # else:
        #     for label in ax1.xaxis.get_ticklabels():
        #         label.set_visible(False)
            
        # if noylabels==None:
        #     plt.ylabel(pltystr,fontsize=fontsize,weight='black')
        #     plt.yticks(fontsize=fontsize,fontweight='black')
        #     ax1.yaxis.offsetText.set_fontsize(fontsize)
        # else:
        #     for label in ax1.yaxis.get_ticklabels():
        #         label.set_visible(False)       

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
            plt.savefig(savefigname,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)
        else:
            plt.draw()
            plt.show()
        plt.close()
        plt.clf()

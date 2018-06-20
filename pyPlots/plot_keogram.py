import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import math
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


def line_cut(filename=None,
             coord1=None, coord2=None,offset=None,
             var=None, op=None,
             expression=None,
             pass_vars=None, pass_times=None
             ):

    ''' Performs a line cut for a keogram plot.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword coord1:      Coordinates [x,y,z] (in Re) of the first end of the cut (for line cut only).
    :kword coord2:      Coordinates [x,y,z] (in Re) of the second end of the cut (for line cut only).
    :kword offset:      Offset to add to cut abscissa axis.

    :kword var:         variable to plot, e.g. rho, rhoBeam, beta, temperature, MA, Mms, va, vms,
                        E, B, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc.
    :kword op:          Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude.
 
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

    :returns:           data_cut: A 1D array of values of var along the line cut.
                        cut_abscissa: a 1D array of coordinate along the line from the origin.
                        time_cut: the time of the slice.

    .. code-block:: python

    '''

    Re = 6.371e+6 # Earth radius in m

    # Open vlsv file
    f=pt.vlsvfile.VlsvReader(filename)

    # Get time stamp
    time_cut=f.read_parameter("time")

    # Read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize/Re
    cellids = f.read_variable("CellID")

    # Unit vector giving the direction of the cut
    u_cut = np.array([coord2[0]-coord1[0], coord2[1]-coord1[1], coord2[2]-coord1[2]])
    length_cut = np.sqrt(u_cut[0]**2 + u_cut[1]**2 + u_cut[2]**2)
    u_cut = u_cut / length_cut

    # Construct list of cellIDs to be included in the cut
    cellid1 = (np.int64)(f.get_cellid(np.array([coord1[0]*Re,coord1[1]*Re,coord1[2]*Re])))
    cellid2 = (np.int64)(f.get_cellid(np.array([coord2[0]*Re,coord2[1]*Re,coord2[2]*Re])))

    ncell_cut = (int)(math.floor(length_cut / cellsize))

    cellid_cut = []
    (x_cut,y_cut,z_cut) = ([],[],[])
    
    for icell in range(0,ncell_cut):
        coord_cell = coord1 + icell*cellsize*u_cut
        cellid_cell = (np.int64)(f.get_cellid(np.array([coord_cell[0]*Re,coord_cell[1]*Re,coord_cell[2]*Re])))
        cellid_cut.append(cellid_cell)
        x_cut.append(coord_cell[0])
        y_cut.append(coord_cell[1])
        z_cut.append(coord_cell[2])
        
    cellid_cut.append(cellid2)
    x_cut.append(coord2[0])
    y_cut.append(coord2[1])
    z_cut.append(coord2[2])

    # Building the cut abscissa
    (x0,y0,z0) = (coord1[0],coord1[1],coord1[2])
    cut_abscissa = []
    for ind in range(0,len(x_cut)):
        cut_abscissa.append(np.sqrt((x_cut[ind]-x0)**2+(y_cut[ind]-y0)**2+(z_cut[ind]-z0)**2))


    ##########
    # Read data and calculate required variables in cellIDS belonging to the cut
    ##########
    if expression==None:
        if var == 'rho':
            data_cut = f.read_variable("rho",cellid_cut)

        elif var == 'rhoBeam':
            data_cut = f.read_variable("RhoBackstream",cellid_cut)

        elif var == 'beta':
            data_cut = f.read_variable("beta",cellid_cut)

        elif var == 'temperature':
            data_cut = f.read_variable("Temperature",cellid_cut)

        elif var == 'MA':
            Vmag = f.read_variable("v",operator='magnitude',cellids=cellid_cut)
            va = f.read_variable("va")
            data_cut = Vmag/va

        elif var == 'Mms':
            Vmag = f.read_variable("v",operator='magnitude',cellids=cellid_cut)
            vms = f.read_variable("vms")
            data_cut = Vmag/vms

        elif var == 'va':
            data_cut = f.read_variable("va",cellid_cut)

        elif var == 'vms':
            data_cut = f.read_variable("vms",cellid_cut)

        elif var == 'B':
            if op==None:
                data_cut = f.read_variable("B",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("B",cellid_cut,operator=op)
                # data_cut = data_cut*1e+9 # could be used to ouptut nanotesla instead of tesla

        elif var == 'E':
            if op==None:
                data_cut = f.read_variable("E",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("E",cellid_cut,operator=op)

        elif var == 'V':
            if op==None:
                data_cut = f.read_variable("V",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("V",cellid_cut,operator=op)
                # data_cut = data_cut*1e-3 # Plot this as km/s instead of m/s

        else:
            # Pipe all other vars directly to analysator
            if op==None:
                cb_title = var
                data_cut = f.read_variable(var,cellid_cut)
                # If value was vector value, take magnitude
                if np.ndim(data_cut) != 1:
                    data_cut = np.linalg.norm(np.asarray(data_cut),axis=-1)
            else:
                data_cut = f.read_variable(var,cellid_cut,operator=op)            
            
        if np.ndim(data_cut)!=1:
            print("Error reading variable "+var+"! Exiting.")
            return -1
        

    # If expression or external routine need variables, read them from the file.
    pass_maps=[]
    if pass_vars!=None:
        if pass_times==None:
            # Gather the required variable maps for a single time step
            for mapval in pass_vars:
                pass_map = f.read_variable(mapval,cellid_cut)
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
                step_cellids = fstep.read_variable("CellID",cellid_cut)
                pass_maps.append([])
                for mapval in pass_vars:
                    pass_map = fstep.read_variable(mapval,cellid_cut)
                    pass_maps[tavg_step_i].append(np.ma.asarray(pass_map))

    # Optional user-defined expression used for color panel instead of a single pre-existing var
    if expression!=None:
        data_cut = expression(pass_maps)

    return (data_cut,np.asarray(cut_abscissa)+offset,time_cut)



def circle_cut(filename=None,
               radius=None, origin=None,
               angle1=None, angle2=None, rot_dir=None,
               offset=None, inv_anglax=None,
               var=None, op=None,
               expression=None,
               pass_vars=None, pass_times=None
               ):

    ''' Performs a line cut for a keogram plot.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword radius:      Radius (in Re) of the cut.
    :kword origin:      Coordinates [x,y,z] (in Re) of centre of the circle used for the cut.
    :kword angle1:      Angle (in deg) of the first end of the cut, relative to the +x axis.
    :kword angle2:      Angle (in deg) of the second end of the cut, relative to the +x axis.
    :kword rot_dir:     Direction of rotation [+1 for trigonometric, -1 for clockwise].
    :kword inv_anglax:  Boolean to inverse the angle axis.

    :kword var:         variable to plot, e.g. rho, rhoBeam, beta, temperature, MA, Mms, va, vms,
                        E, B, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc.
    :kword op:          Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude.
 
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

    :returns:           data_cut: A 1D array of values of var along the circle cut.
                        cut_abscissa: a 1D array of angles along the arc from the origin (in deg).
                        time_cut: the time of the slice.

    .. code-block:: python

    '''

    Re = 6.371e+6 # Earth radius in m

    # Open vlsv file
    f=pt.vlsvfile.VlsvReader(filename)

    # Get time stamp
    time_cut=f.read_parameter("time")

    # Read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize/Re
    cellids = f.read_variable("CellID")

    # Make sure angle1 and angle2 are within [-180,180] before converting them
    angle1 = (angle1+180)%360-180.
    angle2 = (angle2+180)%360-180.

    # Convert angle1 and angle2 for suitable construction of the cut
    conversion = 0.
    if rot_dir==1:
        if angle2<angle1:
            conversion = 360
    else:
        if angle2>angle1:
            conversion = -360
    angle2 = angle2 + conversion
        
    
    # Angle variation corresponding to one cell along the cut
    dangle = np.arctan(cellsize/radius)*180/np.pi
    ncell_cut = (int)(math.floor(math.fabs(angle2-angle1)/cellsize))


    # Construct list of cellIDs to be included in the cut
    cellid_cut = []
    angle_cut = []
    
    for icell in range(0,ncell_cut):
        if rot_dir==1:
            angle = angle1 + icell*dangle
        else:
            angle = angle1 - icell*dangle

        # Polar run
        if ysize==1:
            x_pt = origin[0] + radius*np.cos(angle*np.pi/180.)
            y_pt = origin[1]
            z_pt = origin[2] + radius*np.sin(angle*np.pi/180.)
        # Ecliptic run
        else:
            x_pt = origin[0] + radius*np.cos(angle*np.pi/180.)
            y_pt = origin[1] + radius*np.sin(angle*np.pi/180.)
            z_pt = origin[2]

        cellid_cell = (np.int64)(f.get_cellid(np.array([x_pt*Re,y_pt*Re,z_pt*Re])))
        cellid_cut.append(cellid_cell)
        angle_cut.append(angle)

    # Adding the last point with exact angle value
    # Polar run
    if ysize==1:
        x_pt = origin[0] + radius*np.cos(angle2*np.pi/180.)
        y_pt = origin[1]
        z_pt = origin[2] + radius*np.sin(angle2*np.pi/180.)
    # Ecliptic run
    else:
        x_pt = origin[0] + radius*np.cos(angle2*np.pi/180.)
        y_pt = origin[1] + radius*np.sin(angle2*np.pi/180.)
        z_pt = origin[2]

    cellid_cell = (np.int64)(f.get_cellid(np.array([x_pt*Re,y_pt*Re,z_pt*Re])))
    cellid_cut.append(cellid_cell)
    angle_cut.append(angle2)

    cut_abscissa = np.asarray(angle_cut)-conversion+offset
    if inv_anglax:
        cut_abscissa = -cut_abscissa


    ##########
    # Read data and calculate required variables in cellIDS belonging to the cut
    ##########
    if expression==None:
        if var == 'rho':
            data_cut = f.read_variable("rho",cellid_cut)

        elif var == 'rhoBeam':
            data_cut = f.read_variable("RhoBackstream",cellid_cut)

        elif var == 'beta':
            data_cut = f.read_variable("beta",cellid_cut)

        elif var == 'temperature':
            data_cut = f.read_variable("Temperature",cellid_cut)

        elif var == 'MA':
            Vmag = f.read_variable("v",operator='magnitude',cellids=cellid_cut)
            va = f.read_variable("va")
            data_cut = Vmag/va

        elif var == 'Mms':
            Vmag = f.read_variable("v",operator='magnitude',cellids=cellid_cut)
            vms = f.read_variable("vms")
            data_cut = Vmag/vms

        elif var == 'va':
            data_cut = f.read_variable("va",cellid_cut)

        elif var == 'vms':
            data_cut = f.read_variable("vms",cellid_cut)

        elif var == 'B':
            if op==None:
                data_cut = f.read_variable("B",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("B",cellid_cut,operator=op)
                # data_cut = data_cut*1e+9 # could be used to ouptut nanotesla instead of tesla

        elif var == 'E':
            if op==None:
                data_cut = f.read_variable("E",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("E",cellid_cut,operator=op)

        elif var == 'V':
            if op==None:
                data_cut = f.read_variable("V",cellid_cut,operator='magnitude')
            else:
                data_cut = f.read_variable("V",cellid_cut,operator=op)
                # data_cut = data_cut*1e-3 # Plot this as km/s instead of m/s

        else:
            # Pipe all other vars directly to analysator
            if op==None:
                cb_title = var
                data_cut = f.read_variable(var,cellid_cut)
                # If value was vector value, take magnitude
                if np.ndim(data_cut) != 1:
                    data_cut = np.linalg.norm(np.asarray(data_cut),axis=-1)
            else:
                data_cut = f.read_variable(var,cellid_cut,operator=op)            
            
        if np.ndim(data_cut)!=1:
            print("Error reading variable "+var+"! Exiting.")
            return -1
        

    # If expression or external routine need variables, read them from the file.
    pass_maps=[]
    if pass_vars!=None:
        if pass_times==None:
            # Gather the required variable maps for a single time step
            for mapval in pass_vars:
                pass_map = f.read_variable(mapval,cellid_cut)
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
                step_cellids = fstep.read_variable("CellID",cellid_cut)
                pass_maps.append([])
                for mapval in pass_vars:
                    pass_map = fstep.read_variable(mapval,cellid_cut)
                    pass_maps[tavg_step_i].append(np.ma.asarray(pass_map))

    # Optional user-defined expression used for color panel instead of a single pre-existing var
    if expression!=None:
        data_cut = expression(pass_maps)

    return (data_cut,cut_abscissa,time_cut)







def plot_keogram(filedir=None,
                  start=None, stop=None,
                  cut_type=None,
                  coord1=None, coord2=None, offset=0.,
                  radius=None, origin=None,
                  angle1=None, angle2=None, rot_dir=None,
                  inv_anglax=False,
                  outputdir=None, nooverwrite=None,
                  var=None, op=None,
                  title=None, cbtitle=None, draw=None, usesci=None,
                  symlog=None,
                  colormap=None,
                  ylabel=None,
                  run=None,wmark=None, nocb=None,
                  unit=None, thick=1.0,scale=1.0,
                  noborder=None, noxlabels=None, noylabels=None,
                  vmin=None, vmax=None, lin=None,
                  external=None, expression=None, 
                  pass_vars=None, pass_times=None
                  ):

    ''' Plots a keogram and its colourbar.

    :kword filedir:     Directory where bulk files are located.
    :kword start:       Keogram starting index.
    :kword stop:        Keogram stopping index.
    :kword cut_type:    "circle" or "line".
    :kword coord1:      Coordinates [x,y,z] (in Re) of the first end of the cut (for line cut only).
    :kword coord2:      Coordinates [x,y,z] (in Re) of the second end of the cut (for line cut only).
    :kword offset:      Offset to add to cut abscissa axis.
    :kword radius:      Radius (in Re) of the cut (for circular cut only).
    :kword origin:      Coordinates [x,y,z] (in Re) of centre of the circle used for the cut (for circular cut only).
    :kword angle1:      Angle (in deg) of the first end of the cut, relative to the +x axis (for circular cut only).
    :kword angle2:      Angle (in deg) of the second end of the cut, relative to the +x axis (for circular cut only).
    :kword rot_dir:     Direction of rotation [+1 for trigonometric, -1 for clockwise] (for circular cut only).
    :kword inv_anglax:  Boolean to inverse the angle axis (for circular cut only).
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
    :kword nooverwrite: Set to only perform actions if the target output file does not yet exist.
     
    :kword var:         variable to plot, e.g. rho, rhoBeam, beta, temperature, MA, Mms, va, vms,
                        E, B, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc.
    :kword op:          Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude. 
           
    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, nipy_spectral, RdBu, bwr
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
                            
    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:
    plot_keogram(filedir=fileLocation, var="rho", run="BCH",
                  colormap='nipy_spectral', start=2000, stop=3000, outputdir=outputLocation,
                  cut_type='line', coord1=[-40,0,-20], coord2=[-40,0,20])

    plot_keogram(filedir=fileLocation, var="E", op="y", run="BCH",
                  colormap='nipy_spectral', start=2000, stop=3000, outputdir=outputLocation,
                  cut_type='circle', radius=12., angle1=-135., angle2=135., rot_dir=-1, inv_anglax=True, offset=-180.)

    '''

    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    # watermarkimage=os.path.expandvars('$HOME/appl_taito/analysator/pyPlot/logo_color.png')

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
        # Default values
        colormap="hot_desaturated"
        if op!=None:
            colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=5*scale # Colour bar ticks

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

    # Plot title with time #TODO
    if title==None:
        plot_title = "t="+str(np.int(tstart))+' s'
    else:
        plot_title = title

    # step, used for file name
    if start!=None:
        stepstr = '_'+str(start).rjust(7,'0')
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
    savefigname = outputdir+run+"_keogram_"+varstr+opstr+stepstr+".png"

    # Check if target file already exists and overwriting is disabled
    if (nooverwrite!=None and os.path.exists(savefigname)):
        return

    Re = 6.371e+6 # Earth radius in m
    #read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")

    # Select window to draw
    if cut_type=='line':
        if ((coord1!=None) and (coord2!=None)):
            length_line = np.sqrt((coord2[0]-coord1[0])**2 + (coord2[1]-coord1[1])**2 + (coord2[2]-coord1[2])**2)
            boxcoords = [tstart,tstop,0.,length_line]
        else:
            print("ERROR: needs to specify coordinates of line cut end points")
            return
    elif cut_type=='circle':
        # If no rotation direction is specified, default to trigonometric
        if rot_dir==None:
            rot_dir = 1
        if ((radius!=None) and (angle1!=None) and (angle2!=None)):
            boxcoords = [tstart,tstop,angle1,angle2]
        else:
            print("ERROR: needs to specify radius of circle cut and angles of end points relative to +x axis")
            return
    else:
        print("ERROR: needs to specify cut type")
        return


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
       


    # TODO Building the datamap corresponding to the keogram    
    datamap = np.array([])
    time_keogram = np.array([])

    for step in range(start,stop):
        filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'

        if cut_type=='line':
            (column,cut_abscissa,time_cut) = line_cut(filename=filename,coord1=coord1,coord2=coord2,offset=offset,var=var,op=op, 
                                                           expression=expression,pass_vars=pass_vars, pass_times=pass_times)
        else:
            (column,cut_abscissa,time_cut) = circle_cut(filename=filename,radius=radius,origin=origin,angle1=angle1,angle2=angle2,
                                                           rot_dir=rot_dir,offset=offset,inv_anglax=inv_anglax,
                                                           var=var,op=op,expression=expression,pass_vars=pass_vars,pass_times=pass_times)

        datamap = np.append(datamap,column)
        time_keogram = np.append(time_keogram,time_cut)


    # Determine sizes of the keogram array
    sizes=[time_keogram.size,cut_abscissa.size]


    # Reshape data to an ordered 2D array that can be plotted
    if np.ndim(datamap) != 2:
        datamap = datamap.reshape([sizes[0],sizes[1]])

    datamap = np.transpose(datamap)
        

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


    # Scale data extent and plot box
    boxcoords[2] = cut_abscissa[0]
    boxcoords[3] = cut_abscissa[-1]


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
    figsize = [4.0,3.15]
    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=figsize,dpi=300)
    
    # Generates the mesh to map the data to.
    # Note, datamap is still of shape [ysize,xsize] (?)
    [XmeshXY,YmeshXY] = scipy.meshgrid(time_keogram,cut_abscissa)

    fig1 = plt.pcolormesh(XmeshXY,YmeshXY,datamap, cmap=colormap,norm=norm)
    ax1 = plt.gca() # get current axes

    # Save the keogram to file in case further processing is needed
    datamap.dump(outputdir+'datamap')
    np.save(outputdir+'time_keogram',time_keogram)
    np.save(outputdir+'cut_abscissa',cut_abscissa)

    # Title and plot limits
    if len(plot_title)!=0:
        ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')

    plt.xlim([boxcoords[0],boxcoords[1]])
    plt.ylim([boxcoords[2],boxcoords[3]])

    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(thick)
    ax1.xaxis.set_tick_params(width=thick,length=3)
    ax1.yaxis.set_tick_params(width=thick,length=3)
    #ax1.xaxis.set_tick_params(which='minor',width=3,length=5)
    #ax1.yaxis.set_tick_params(which='minor',width=3,length=5)

    if noxlabels==None:
        plt.xlabel('Time [s]',fontsize=fontsize,weight='black')
        plt.xticks(fontsize=fontsize,fontweight='black')
        ax1.xaxis.offsetText.set_fontsize(fontsize)# set axis exponent offset font sizes
    if noylabels==None:
        if ylabel!=None:
            plt.ylabel(ylabel)
        elif cut_type=='line': # Line cut
            plt.ylabel('s ['+unitstr+']',fontsize=fontsize,weight='black')
        else: # Circular cut
            plt.ylabel('Angle [deg]',fontsize=fontsize,weight='black')
        plt.yticks(fontsize=fontsize,fontweight='black')
        ax1.yaxis.offsetText.set_fontsize(fontsize)# set axis exponent offset font sizes


    # Optional external additional plotting routine overlayed on color plot
    # Uses the same pass_maps variable as expressions
    if external!=None:
        extresult=external(ax1, XmeshXY,YmeshXY, pass_maps)

    # Colour bar title construction
    if expression==None:
        if var == 'rho':
            cb_title = r"$n_\mathrm{p} [\mathrm{m}^{-3}]$"
        elif var == 'rhoBeam':
            cb_title = r"$\rho_{\mathrm{beam}} [\mathrm{m}^{-3}]$"
        elif var == 'beta':
            cb_title = r"$\beta$"
        elif var == 'temperature':
            cb_title = r"$T$ [K]"
        elif var == 'MA':
            cb_title = r"$\mathrm{M}_\mathrm{A}$"
        elif var == 'Mms':
            cb_title = r"$\mathrm{M}_\mathrm{ms}$"
        elif var == 'va':
            cb_title = r"$v_\mathrm{A}$"
        elif var == 'vms':
            cb_title = r"$v_\mathrm{ms}$"
        elif var == 'B':
            if op==None:
                cb_title = r"$|B|$ [T]"
            else:
                cb_title = r"$B_"+op+"$ [T]"
        elif var == 'E':
            if op==None:
                cb_title = r"$|E|$ [V/m]"
            else:
                cb_title = r"$E_"+op+"$ [V/m]"
        elif var == 'V':
            if op==None:
                cb_title = r"$|V|\,[\mathrm{m}\,\mathrm{s}^{-1}]$"
            else:
                cb_title = r"$V_"+op+"\,[\mathrm{m}\,\mathrm{s}^{-1}]$"
        else:
            # Pipe all other vars directly to analysator
            if op==None:
                cb_title = var
                # If value was vector value, take magnitude
                if np.ndim(datamap) != 1:
                    cb_title = r"$|"+var+"|$"
            else:
                cb_title = r+""+var+"$_"+op+"$"

    if cbtitle==None:
        if expression!=None:
            cb_title_use = expression.__name__.replace("_","\_") # replaces underscores so math mode subscript mode isn't activated
        else:
            cb_title_use = cb_title
    else:
        cb_title_use = cbtitle

    # Colourbar title
    if len(cb_title_use)!=0:
        cb_title_locy = 1.0 + 0.08#/ratio
        plt.text(1.0, 1.03, cb_title_use, fontsize=fontsize,weight='black', transform=ax1.transAxes)

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

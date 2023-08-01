# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from skimage import measure
import scipy
import os, sys
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import LinearLocator
import matplotlib.ticker as mtick
import colormaps as cmaps
from matplotlib.cbook import get_sample_data
from distutils.version import LooseVersion, StrictVersion
import ids3d


def plot_isosurface(filename=None,
                    vlsvobj=None,
                    filedir=None, step=None,
                    outputdir=None, nooverwrite=None,
                    #
                    surf_var=None, surf_op=None, surf_level=None,
                    color_var=None, color_op=None,
                    surf_step=1,
                    #
                    title=None, cbtitle=None,
                    draw=None, usesci=None,
                    symmetric=False,
                    highres=None,
                    boxm=[],boxre=[],
                    colormap=None,
                    run=None,wmark=None, nocb=False,
                    unit=None, thick=1.0,scale=1.0,
                    vscale=1.0,
                    vmin=None, vmax=None, lin=None, symlog=None,
                    angle = [30.,90.]
                    ):

    ''' Plots a coloured plot with axes and a colour bar.

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/ or override with PTOUTPUTDIR)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
    :kword nooverwrite: Set to only perform actions if the target output file does not yet exist                    
     
    :kword surf_var:    Variable to read for defining surface
    :kword surf_op:     Operator to use for variable to read surface
    :kword surf_level:  Level at which to define surface
    :kword surf_step:   Vertex stepping for surface generation: larger value returns coarser surface

    :kword color_var:   Variable to read for coloring surface
    :kword color_op:    Operator to use for variable to color surface ('x', 'y', 'z'; if None and color_var is a vector, the magnitude is taken)
           
    :kword boxm:        zoom box extents [x0,x1,y0,y1] in metres (default and truncate to: whole simulation box)
    :kword boxre:       zoom box extents [x0,x1,y0,y1] in Earth radii (default and truncate to: whole simulation box)
    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, parula, nipy_spectral, RdBu, bwr
    :kword run:         run identifier, used for constructing output filename
    :kword title:       string to use as plot title instead of time
    :kword cbtitle:     string to use as colorbar title instead of map name
    :kword unit:        Plot axes using 10^{unit} m (default: Earth radius R_E)

    :kwird usesci:      Use scientific notation for colorbar ticks? (default: 1)
    :kword vscale:      Scale all values with this before plotting. Useful for going from e.g. m^-3 to cm^-3
                        or from tesla to nanotesla. Guesses correct units for colourbar for some known
                        variables. Set to None to seek for a default scaling.
    :kword vmin,vmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword symmetric:   Set the absolute value of vmin and vmax to the greater of the two
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2.
    :kword wmark:       If set to non-zero, will plot a Vlasiator watermark in the top left corner.
    :kword highres:     Creates the image in high resolution, scaled up by this value (suitable for print). 
    :kword draw:        Set to nonzero in order to draw image on-screen instead of saving to file (requires x-windowing)

    :kword scale:       Scale text size (default=1.0)
    :kword thick:       line and axis thickness, default=1.0
    :kword nocb:        Set to suppress drawing of colourbar
    :kword angle:       Viewing elevation and azimuthal angles in degrees

    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

    # Example usage:

    '''
    # Verify the location of this watermark image
    watermarkimage=os.path.join(os.path.dirname(__file__), 'logo_color.png')
    # watermarkimage=os.path.expandvars('$HOME/appl_taito/analysator/pyPlot/logo_color.png')

    outputprefix = ''
    if outputdir==None:
        outputdir=pt.plot.defaultoutputdir
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
        filename = glob.glob(filedir+'bulk*'+str(step).rjust(7,'0')+'.vlsv')[0]
        #filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
        f=pt.vlsvfile.VlsvReader(filename)
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return
                
    # Scientific notation for colorbar ticks?
    if usesci==None:
        usesci=1
    
    if colormap==None:
        # Default values
        colormap="hot_desaturated"
        if color_op!=None:
            colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=7*scale # Colour bar ticks

    # Plot title with time
    timeval=None
    timeval=f.read_parameter("time")
    if timeval==None:
        timeval=f.read_parameter("t")
    if timeval==None:
        print("Unknown time format encountered")

    # Plot title with time
    if title==None:        
        if timeval == None:    
            print("Unknown time format encountered")
            plot_title = ''
        else:
            #plot_title = "t="+str(np.int(timeval))+' s'
            plot_title = "t="+'{:4.2f}'.format(timeval)+' s'
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
    if run==None:
        run='plot'

    # Verify validity of operator
    surf_opstr=''
    color_opstr=''
    if color_op!=None:
        if color_op!='x' and color_op!='y' and color_op!='z':
            print("Unknown operator "+color_op+", defaulting to None/magnitude for a vector.")
            color_op=None            
        else:
            # For components, always use linear scale, unless symlog is set
            color_opstr='_'+color_op
            if symlog==None:
                lin=9
    # Verify validity of operator
    if surf_op!=None:
        if surf_op!='x' and surf_op!='y' and surf_op!='z':
            print("Unknown operator "+surf_op)
            surf_op=None            
        else:
            surf_opstr='_'+surf_op

    # Output file name
    surf_varstr=surf_var.replace("/","_")
    if color_var!=None:
        color_varstr=color_var.replace("/","_")
    else:
        color_varstr="solid"
    savefigname = outputdir+outputprefix+run+"_isosurface_"+surf_varstr+surf_opstr+"-"+color_varstr+color_opstr+stepstr+".png"

    # Check if target file already exists and overwriting is disabled
    if (nooverwrite!=None and os.path.exists(savefigname)):
        # Also check that file is not empty
        if os.stat(savefigname).st_size > 0:
            return
        else:
            print("Found existing file "+savefigname+" of size zero. Re-rendering.")


    Re = 6.371e+6 # Earth radius in m
    # read in mesh size and cells in ordinary space
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    xsize = int(xsize)
    ysize = int(ysize)
    zsize = int(zsize)    
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    cellsize = (xmax-xmin)/xsize
    cellids = f.read_variable("CellID")

    # Read the FSgrid mesh
    try:
        [xsizefg, ysizefg, zsizefg] = f.get_fsgrid_mesh_size()
        xsizefg = int(xsizefg)
        ysizefg = int(ysizefg)
        zsizefg = int(zsizefg)
        [xminfg, yminfg, zminfg, xmaxfg, ymaxfg, zmaxfg] = f.get_fsgrid_mesh_extent()
        cellsizefg = (xmaxfg-xminfg)/xsizefg
        pt.plot.plot_helpers.CELLSIZE = cellsizefg
    except:
        if xsize!=1 and ysize!=1 and zsize!=1:
            print("Did not find fsgrid data, but found 3D DCCRG mesh. Attempting to adapt.")
            [xsizefg, ysizefg, zsizefg] = [xsize * 2**f.get_max_refinement_level(), ysize * 2**f.get_max_refinement_level(), zsize * 2**f.get_max_refinement_level()]
            [xminfg, yminfg, zminfg, xmaxfg, ymaxfg, zmaxfg] = [xmin, ymin, zmin, xmax, ymax, zmax]
            cellsizefg = cellsize
            pt.plot.plot_helpers.CELLSIZE = cellsize
        else:
            print("Found 2D DCCRG mesh without FSgrid data. Exiting.")
            return -1

    # sort the cellid and the datamap list
    indexids = cellids.argsort()
    cellids = cellids[indexids]

    # find the highest refiment level
    reflevel = ids3d.refinement_level(xsize, ysize, zsize, cellids[-1])
    for i in range(5): # Check if Vlasov grid doesn't reach maximum (fsgrid) refinement
        if xsize*(2**(reflevel + i)) == xsizefg:
            reflevel += i
            break


    if (xsize==1) or (ysize==1) or (zsize==1):
        print("Error: isosurface plotting requires 3D spatial domain!")
        return

    simext=[xmin,xmax,ymin,ymax,zmin,zmax]
    sizes=[xsize,ysize,zsize]

    # Select window to draw
    if len(boxm)==6:
        boxcoords=boxm
    elif len(boxre)==6:
        boxcoords=[i*Re for i in boxre]
    else:
        boxcoords=simext

    # If box extents were provided manually, truncate to simulation extents
    boxcoords[0] = max(boxcoords[0],simext[0])
    boxcoords[1] = min(boxcoords[1],simext[1])
    boxcoords[2] = max(boxcoords[2],simext[2])
    boxcoords[3] = min(boxcoords[3],simext[3])
    boxcoords[4] = max(boxcoords[4],simext[4])
    boxcoords[5] = min(boxcoords[5],simext[5])

    # Axes and units (default R_E)
    if unit!=None: # Use m or km or other
        if unit==0:
            unitstr = pt.plot.rmstring('m')
        if unit==3:
            unitstr = pt.plot.rmstring('km')
        else:
            unitstr = r'10^{'+str(int(unit))+'} '+pt.plot.rmstring('m')
        unit = np.power(10,int(unit))
    else:
        unitstr = pt.plot.rmstring('R')+'_'+pt.plot.rmstring('E')
        unit = Re
    unitstr = pt.plot.mathmode(unitstr)
        
    # Scale data extent and plot box
    simext_org = simext
    simext=[i/unit for i in simext]
    boxcoords=[i/unit for i in boxcoords]

    if color_var != None:
        if color_op==None:
            color_op="pass"
        datamap_info = f.read_variable_info(color_var, operator=color_op)

        cb_title_use = datamap_info.latex
        # Check if vscale results in standard unit
        vscale, _, datamap_unit_latex = datamap_info.get_scaled_units(vscale=vscale)

        # Add unit to colorbar title
        if datamap_unit_latex:
            cb_title_use = cb_title_use + "\,["+datamap_unit_latex+"]"
    else: # color_var==None
        cb_title_use = ""
        nocb=1

    if cbtitle is not None:
        # Here allow underscores for manual math mode
        cb_title_use = cbtitle      


    if surf_op is None:
            surf_op="pass"
    surf_datamap_info = f.read_variable_info(surf_var, operator=surf_op)

    surf_datamap = surf_datamap_info.data

    # Verify data shape
    if np.ndim(surf_datamap)==0:
        print("Error, read only single surface variable value from vlsv file!",surf_datamap.shape)
        return -1
    

    #Define the box limits
    low = np.array([boxcoords[0], boxcoords[2], boxcoords[4]])*unit
    up = np.array([boxcoords[1], boxcoords[3], boxcoords[5]])*unit

    # Choose CellIDs inside the box
    ids, idx = ids3d.ids3d_box(cellids, low, up, reflevel, xsize, ysize, zsize, [xmin, ymin, zmin, xmax, ymax, zmax])

    surf_datamap = surf_datamap[indexids]
    surf_datamap = surf_datamap[idx]

    # Reshape surface data
    if np.ndim(surf_datamap)==1:
        surf_data_in_box = ids3d.idmesh3d2(ids, surf_datamap, reflevel, xsize, ysize, zsize,  None)
    elif np.ndim(surf_datamap)==2:
        surf_data_in_box = ids3d.idmesh3d2(ids, surf_datamap, reflevel, xsize, ysize, zsize, surf_datamap.shape[1])
    elif np.ndim(surf_datamap)==3:
        surf_data_in_box = ids3d.idmesh3d2(ids, surf_datamap, reflevel, xsize, ysize, zsize, (surf_datamap.shape[1],surf_datamap.shape[2]))
    else:
        print("Dimension error in constructing 2D AMR slice!")
        return -1

    # Remove empty rows, columns and tubes
    empty_0 = np.where(~np.all(surf_data_in_box==0, axis=0))
    empty_1 = np.where(~np.all(surf_data_in_box==0, axis=1))

    MaskX = empty_1[0]
    MaskY = empty_0[0]
    MaskZ = empty_1[1]

    cropped_surf_data = surf_data_in_box[min(MaskX):max(MaskX)+1,:,:]
    cropped_surf_data = cropped_surf_data[:,min(MaskY):max(MaskY)+1,:]
    cropped_surf_data = cropped_surf_data[:,:,min(MaskZ):max(MaskZ)+1]


    # Keep removing the face of the rectangle with the largest number of invalid values until the surface has only proper values
    counter = 0
    while True:
        zero_counts = np.zeros(6)
        zero_counts[0] = np.count_nonzero(cropped_surf_data[0, :, :]  == 0) # Face 0 (front)
        zero_counts[1] = np.count_nonzero(cropped_surf_data[-1, :, :] == 0)  # Face 1 (back)
        zero_counts[2] = np.count_nonzero(cropped_surf_data[:, 0, :]  == 0) # Face 2 (left)
        zero_counts[3] = np.count_nonzero(cropped_surf_data[:, -1, :] == 0)  # Face 3 (right)
        zero_counts[4] = np.count_nonzero(cropped_surf_data[:, :, 0]  == 0) # Face 4 (top)
        zero_counts[5] = np.count_nonzero(cropped_surf_data[:, :, -1] == 0)  # Face 5 (bottom)

        face_to_be_removed = np.argmax(zero_counts)

        if face_to_be_removed==0:
            cropped_surf_data = cropped_surf_data[1:,:,:]
        elif face_to_be_removed==1:
            cropped_surf_data = cropped_surf_data[:-1,:,:]
        elif face_to_be_removed==2:
            cropped_surf_data = cropped_surf_data[:,1:,:]
        elif face_to_be_removed==3:
            cropped_surf_data = cropped_surf_data[:,:-1,:]
        elif face_to_be_removed==4:
            cropped_surf_data = cropped_surf_data[:,:,1:]
        elif face_to_be_removed==5:
            cropped_surf_data = cropped_surf_data[:,:,:-1]

        counter+=1
        if counter > 50:    # Failsafe
            print("Error with boundaries! Exiting.")
            return -1
        
        if np.all(zero_counts==0):
            break 
    

    if surf_level==None:
        surf_level = 0.5*(np.amin(cropped_surf_data)+np.amax(cropped_surf_data))
    print("Minimum found surface value "+str(np.amin(cropped_surf_data))+" surface level "+str(surf_level)+" max "+str(np.amax(cropped_surf_data)))

    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw is not None:
        if str(matplotlib.get_backend()) is not pt.backend_interactive: #'TkAgg': 
            plt.switch_backend(pt.backend_interactive)
    else:
        if str(matplotlib.get_backend()) is not pt.backend_noninteractive: #'Agg':
            plt.switch_backend(pt.backend_noninteractive)  


    spacing = ((xmaxfg-xminfg)/(xsizefg*unit), (ymaxfg-yminfg)/(ysizefg*unit), (zmaxfg-zminfg)/(zsizefg*unit))

    verts, faces, normals, arrays = measure.marching_cubes(cropped_surf_data, level=surf_level,
                                                                   spacing=spacing,
                                                                   gradient_direction='descent',
                                                                   step_size=surf_step,
                                                                   allow_degenerate=True, method='lewiner', mask=cropped_surf_data!=0)
    


    #offset with respect to simulation domain corner
    verts[:,0] = verts[:,0] + low[0]/unit
    verts[:,1] = verts[:,1] + low[1]/unit
    verts[:,2] = verts[:,2] + low[2]/unit

    # A box that contains the isosurface
    contour_box_low = np.array([min(verts[:,0])-1,min(verts[:,1])-1,min(verts[:,2])-1])*unit
    contour_box_up = np.array([max(verts[:,0])+1,max(verts[:,1])+1,max(verts[:,2])+1])*unit



    # Next find color variable values at vertices
    if color_var != None:
        nverts = len(verts[:,0])
        print("Extracting color values for "+str(nverts)+" vertices and "+str(len(faces[:,0]))+" faces.")
        all_coords = np.empty((nverts, 3))
        for i in np.arange(nverts):            
            # # due to mesh generation, some coordinates may be outside simulation domain
            # WARNING this means it might be doing wrong things in the periodic dimension of 2.9D runs.
            coords = verts[i,:]*unit 
            coords[0] = max(coords[0],simext_org[0]+0.1*cellsize)
            coords[0] = min(coords[0],simext_org[1]-cellsize)
            coords[1] = max(coords[1],simext_org[2]+0.1*cellsize)
            coords[1] = min(coords[1],simext_org[3]-cellsize)
            coords[2] = max(coords[2],simext_org[4]+0.1*cellsize)
            coords[2] = min(coords[2],simext_org[5]-cellsize)
            all_coords[i] = coords


        # Choose CellIDs inside the isosurface box
        color_ids, color_idx = ids3d.ids3d_box(cellids, contour_box_low, contour_box_up, reflevel, xsize, ysize, zsize, [xmin, ymin, zmin, xmax, ymax, zmax])

        # Read the variables to be plotted inside the box
        if f.check_variable(color_var)!=True:
            print("Error, color variable "+color_var+" not found!")
            return -1
        if color_op=="pass":
            vg_colors = f.read_variable(color_var, color_ids)
            # If value was vector value, take magnitude
            if np.ndim(vg_colors) != 1:
                vg_colors = np.linalg.norm(np.asarray(vg_colors),axis=-1)
        else:
            vg_colors = f.read_variable(color_var, color_ids, operator=color_op)            
        if np.ndim(vg_colors)!=1:
            print("Error reading color variable "+color_var+"! Exiting.")
            return -1


        vg_coords = f.read_variable("vg_coordinates", color_ids)

        # Interpolate the values to the surface
        interpolation = scipy.interpolate.RBFInterpolator(vg_coords, vg_colors, neighbors = 27)

        color_data = interpolation(all_coords)*vscale

        # Make sure color data is 1-dimensional (e.g. magnitude of E instead of 3 components)
        if np.ndim(color_data)!=1:
            color_data=np.linalg.norm(color_data, axis=-1)

    if color_var==None:
        # dummy norm
        print("No surface color given, using dummy setup")
        norm = BoundaryNorm([0,1], ncolors=cmapuse.N, clip=True)
        vminuse=0
        vmaxuse=1
    else:
        # If automatic range finding is required, find min and max of array
        # Performs range-finding on a masked array to work even if array contains invalid values
        color_data = np.ma.masked_invalid(color_data)
        if vmin!=None:
            vminuse=vmin
        else: 
            vminuse=np.ma.amin(color_data)
        if vmax!=None:
            vmaxuse=vmax
        else:
            vmaxuse=np.ma.amax(color_data)       

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

        # Check that lower bound is valid for logarithmic plots
        if (vminuse <= 0) and (lin==None) and (symlog==None):
            # Drop negative and zero values
            vminuse = np.ma.amin(np.ma.masked_less_equal(color_data,0))

        # Special case of very small vminuse values
        if ((vmin is None) or (vmax is None)) and (vminuse > 0) and (vminuse < vmaxuse*1.e-5):
            vminuse = vmaxuse*1e-5
            if lin is not None:
                vminuse = 0

        # If symlog scaling is set:
        if symlog!=None:
            if symlog>0:
                linthresh = symlog 
            else:
                linthresh = max(abs(vminuse),abs(vmaxuse))*1.e-2

        # Lin or log colour scaling, defaults to log
        if lin is None:
            # Special SymLogNorm case
            if symlog is not None:
                if LooseVersion(matplotlib.__version__) < LooseVersion("3.2.0"):
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
            ticks = np.linspace(vminuse,vmaxuse,num=lin)            

        print("Selected color range: "+str(vminuse)+" to "+str(vmaxuse))

    # Create 300 dpi image of suitable size
    fig = plt.figure(figsize=[5,5],dpi=450)
    #ax1 = fig.gca(projection='3d')
    ax1 = fig.add_subplot(111, projection='3d')

    # If requested high res image
    if highres:
        highresscale = 2
        if ((type(highres) is float) or (type(highres) is int)):
            highresscale = float(highres)
            if np.isclose(highresscale, 1.0):
                highresscale = 2
        figsize= [x * highresscale for x in figsize]
        fontsize=fontsize*highresscale
        fontsize2=fontsize2*highresscale
        fontsize3=fontsize3*highresscale
        scale=scale*highresscale
        thick=thick*highresscale
        streamlinethick=streamlinethick*highresscale
        vectorsize=vectorsize*highresscale

    # Generate virtual bounding box to get equal aspect
    maxrange = np.array([boxcoords[1]-boxcoords[0], boxcoords[3]-boxcoords[2], boxcoords[5]-boxcoords[4]]).max() / 2.0
    midvals = np.array([boxcoords[1]+boxcoords[0], boxcoords[3]+boxcoords[2], boxcoords[5]+boxcoords[4]]) / 2.0


    # Three options:
    # 2.9D ecliptic, 2.9D polar, or 3D
    if ysize < 0.2*xsize: # 2.9D polar, perform rotation
        generatedsurface = ax1.plot_trisurf(verts[:,2], verts[:,0], verts[:,1], triangles=faces,
                                            cmap=cmapuse, norm=norm, vmin=vminuse, vmax=vmaxuse, 
                                            lw=0, shade=False, edgecolors=None, antialiased=False)
        ax1.set_xlabel("z ["+unitstr+"]", fontsize=fontsize3)
        ax1.set_ylabel("x ["+unitstr+"]", fontsize=fontsize3)
        ax1.set_zlabel("y ["+unitstr+"]", fontsize=fontsize3)


        # Set camera angle
        ax1.view_init(elev=angle[0], azim=angle[1])
        # Set virtual bounding box
        ax1.set_xlim([midvals[2]-maxrange, midvals[2]+maxrange])
        ax1.set_ylim([midvals[0]-maxrange, midvals[0]+maxrange])
        ax1.set_zlim([midvals[1]-maxrange, midvals[1]+maxrange])
        ax1.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
        
    else: # 3D or 2.9D ecliptic, leave as is
        generatedsurface = ax1.plot_trisurf(verts[:,0], verts[:,1], verts[:,2], triangles=faces,
                                            cmap=cmapuse, norm=norm, vmin=vminuse, vmax=vmaxuse, 
                                            lw=0.2, shade=False, edgecolors=None)
        ax1.set_xlabel("x ["+unitstr+"]", fontsize=fontsize3)
        ax1.set_ylabel("y ["+unitstr+"]", fontsize=fontsize3)
        ax1.set_zlabel("z ["+unitstr+"]", fontsize=fontsize3)
        # Set camera angle
        ax1.view_init(elev=angle[0], azim=angle[1])
        # Set virtual bounding box
        ax1.set_xlim([midvals[0]-maxrange, midvals[0]+maxrange])
        ax1.set_ylim([midvals[1]-maxrange, midvals[1]+maxrange])
        ax1.set_zlim([midvals[2]-maxrange, midvals[2]+maxrange])
        ax1.tick_params(labelsize=fontsize3)#,width=1.5,length=3)


    # Setting per-triangle colours for plot_trisurf needs to be done
    # as a separate set_array call.
    if color_var != None:
        # Find face-averaged colors
        # (simply setting the array to color_data failed for some reason)
        colors = np.mean(color_data[faces], axis=1)
        generatedsurface.set_array(colors)
        
    # Title and plot limits
    if len(plot_title)!=0:
        ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')


    # ax1.set_aspect('equal') #<<<--- this does not work for 3D plots!

    # for axis in ['top','bottom','left','right']:
    #     ax1.spines[axis].set_linewidth(thick)
    # ax1.xaxis.set_tick_params(width=thick,length=3)
    # ax1.yaxis.set_tick_params(width=thick,length=3)
    # #ax1.xaxis.set_tick_params(which='minor',width=3,length=5)
    # #ax1.yaxis.set_tick_params(which='minor',width=3,length=5)    


    if not nocb:
        # First draw colorbar
        if usesci==0:        
            cb = fig.colorbar(generatedsurface,ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmt), drawedges=False, fraction=0.023, pad=0.02)
        else:
            cb = fig.colorbar(generatedsurface,ticks=ticks,format=mtick.FuncFormatter(pt.plot.cbfmtsci),drawedges=False, fraction=0.046, pad=0.04)

        if len(cb_title_use)!=0:
            cb_title_use = pt.plot.mathmode(pt.plot.bfstring(cb_title_use))

        if lin is None:
            pt.plot.cb_linear = False
        else:
            pt.plot.cb_linear = True

        cb.ax.tick_params(labelsize=fontsize3)#,width=1.5,length=3)
        cb.outline.set_linewidth(thick)
        cb.ax.set_title(cb_title_use)
        cb.ax.title.set_horizontalalignment('center')

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

        # if too many subticks in logarithmic colorbar:
        if lin is None and symlog is None:
            nlabels = len(cb.ax.yaxis.get_ticklabels()) #/ ratio
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
    if wmark!=None:        
        wm = plt.imread(get_sample_data(watermarkimage))
        newax = fig.add_axes([0.01, 0.90, 0.3, 0.08], anchor='NW', zorder=-1)
        newax.imshow(wm)
        newax.axis('off')

    plt.tight_layout()
    savefig_pad=0.1 # The default is 0.1
    bbox_inches=None


    # Save output or draw on-screen
    if draw==None:
        # Note: generated title can cause strange PNG header problems
        # in rare cases. This problem is under investigation, but is related to the exact generated
        # title string. This try-catch attempts to simplify the time string until output succedes.
        try:
            plt.savefig(savefigname,dpi=450, bbox_inches=bbox_inches, pad_inches=savefig_pad)
            savechange=0
        except:
            savechange=1
            plot_title = "t="+'{:4.1f}'.format(timeval)+' s '
            ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')                
            try:
                plt.savefig(savefigname,dpi=450, bbox_inches=bbox_inches, pad_inches=savefig_pad)
            except:
                plot_title = "t="+str(np.int(timeval))+' s   '
                ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')                
                try:
                    plt.savefig(savefigname,dpi=450, bbox_inches=bbox_inches, pad_inches=savefig_pad)
                except:
                    plot_title = ""
                    ax1.set_title(plot_title,fontsize=fontsize2,fontweight='bold')                
                    try:
                        plt.savefig(savefigname,dpi=450, bbox_inches=bbox_inches, pad_inches=savefig_pad)
                    except:
                        print("Error:", sys.exc_info())
                        print("Error with attempting to save figure, sometimes due to matplotlib LaTeX integration.")
                        print("Usually removing the title should work, but this time even that failed.")                        
                        savechange = -1
        if savechange>0:
            print("Due to rendering error, replaced image title with "+plot_title)
        if savechange>=0:
            print(savefigname+"\n")
    else:
        plt.draw()
        plt.show()

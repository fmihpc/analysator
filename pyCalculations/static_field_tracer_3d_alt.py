import numpy as np
import warnings
from scipy import interpolate


def static_field_tracer_3d_alt( vlsvReader, coord_list, max_iterations, dx, direction='+', fg_b = None ):
    ''' trace_static_field() integrates along the (static) magnetic field to calculate a final position 
        based on static_field_tracer3d, which in turn is based on Analysator's static_field_tracer

        :param vlsvReader:         An open vlsv file
        :param coord_list:         3-element list of 1D numpy arrays representing starting positions [x, y, z]
        :param max_iterations:     The maximum amount of iterations before the algorithm stops
        :param dx:                 One iteration step length
        :param direction:          '+' or '-' or '+-' Follow field in the plus direction or minus direction
        :param bvar:               String, variable name to trace [default 'B']
        :returns:                  List of coordinates    
    '''
    #    if (bvar is not 'B'):
    #        warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")
    if direction == '+-':
        backward = static_field_tracer_3d_alt(vlsvReader, coord_list, max_iterations, dx, direction='-', fg_b=fg_b)
        backward.reverse()
        forward = static_field_tracer_3d_alt(vlsvReader, coord_list, max_iterations, dx, direction='+', fg_b=fg_b)
        return backward + forward
    f = vlsvReader
    # Read face_B (denoted 'fg_b' in .vlsv files):
    if (fg_b is None):
        fg_b = f.read_variable('fg_b')        # EGL: fg_b.shape = (1024, 736, 736, 3)  
    # Create x, y, and z coordinates:
    # Read cellids in order to sort variables
    #cellids = vlsvReader.read_variable("CellID")
    xsize = fg_b.shape[0]
    ysize = fg_b.shape[1]
    zsize = fg_b.shape[2]
    xmin = f.read_parameter('xmin')
    xmax = f.read_parameter('xmax')
    ymin = f.read_parameter('ymin')
    ymax = f.read_parameter('ymax')
    zmin = f.read_parameter('zmin')
    zmax = f.read_parameter('zmax')
    sizes = np.array([xsize, ysize, zsize])
    maxs = np.array([xmax, ymax, zmax])
    mins = np.array([xmin, ymin, zmin])
    dcell = (maxs - mins)/(sizes.astype('float'))
    x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
    y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
    z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
    coordinates = np.array([x,y,z])
    # Debug:
    #print("SIZE WRONG: " + str(len(x)) + " " + str(sizes[0]))
    # Create grid interpolation
    interpolator_face_B_0 = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), fg_b[:,:,:,0], bounds_error = False, fill_value = np.nan)
    interpolator_face_B_1 = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), fg_b[:,:,:,1], bounds_error = False, fill_value = np.nan)
    interpolator_face_B_2 = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), fg_b[:,:,:,2], bounds_error = False, fill_value = np.nan)
    interpolators = [interpolator_face_B_0, interpolator_face_B_1, interpolator_face_B_2]
   #######################################################
    if direction == '-':
        multiplier = -1
    else:
        multiplier = 1
    points = [np.array(coord_list).T]           # list of 2d arrays each with shape [N, 3]
    point = points[0]
    N = coord_list[0].size
    B_unit = np.zeros([3, N])                     # B_unit has shape [3,N] for speed considerations
    for i in range(max_iterations):
#        point = points[j]
        B_unit[0, :] = interpolators[0](point)
        B_unit[1, :] = interpolators[1](point)
        B_unit[2, :] = interpolators[2](point)
        B_mag = np.linalg.norm(B_unit, axis=(0))
        B_unit[0, :] = B_unit[0, :] / B_mag
        B_unit[1, :] = B_unit[1, :] / B_mag
        B_unit[2, :] = B_unit[2, :] / B_mag
        new_point = point + multiplier*B_unit.T * dx
        point = new_point
        points.append( point )
   #######################################################
    return points





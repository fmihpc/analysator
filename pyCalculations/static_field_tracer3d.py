import numpy as np
import scipy as sp
import pytools as pt
import warnings
from scipy import interpolate


def static_field_tracer3d( vlsvReader, x0, max_iterations, dx, direction='+', bvar='fg_b' ):
   ''' Field tracer in a static frame
       :param vlsvReader:         An open vlsv file
       :param x:                  Starting point for the field trace, 3-element numpy array
       :param max_iterations:     The maximum amount of iteractions before the algorithm stops
       :param dx:                 One iteration step length
       :param direction:          '+' or '-' or '+-' Follow field in the plus direction or minus direction
       :param bvar:               String, variable name to trace [default 'fg_b'].
                                  OR an ndarray, i.e. obtained vlsvReader.reader_parameter() externally and passed into the function
       :returns:                  List of coordinates
   '''

   if(bvar is not 'fg_b'):
     warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")

   if direction == '+-':
     backward = static_field_tracer3d(vlsvReader, x0, max_iterations, dx, direction='-', bvar=bvar)
     backward.reverse()
     forward = static_field_tracer3d(vlsvReader, x0, max_iterations, dx, direction='+', bvar=bvar)
     return backward + forward

   f = vlsvReader
   # Read cellids in order to sort variables
   cellids = vlsvReader.read_variable("CellID")
   xsize = f.read_parameter("xcells_ini")
   ysize = f.read_parameter("ycells_ini")
   zsize = f.read_parameter("zcells_ini")
   xmin = f.read_parameter('xmin')
   xmax = f.read_parameter('xmax')
   ymin = f.read_parameter('ymin')
   ymax = f.read_parameter('ymax')
   zmin = f.read_parameter('zmin')
   zmax = f.read_parameter('zmax')
   
   # Read face_B:
   if type(bvar) == np.ndarray:
      face_B = bvar
   else:
      face_B = f.read_variable(bvar)
   
   xsize = face_B.shape[0]
   ysize = face_B.shape[1]
   zsize = face_B.shape[2]
   
   sizes = np.array([xsize, ysize, zsize])
   maxs = np.array([xmax, ymax, zmax])
   mins = np.array([xmin, ymin, zmin])
   dcell = (maxs - mins)/(sizes.astype('float'))

   # Create x, y, and z coordinates:
   x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
   y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
   z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
   # Debug:
   if( len(x) != sizes[0] ):
      print("SIZE WRONG: " + str(len(x)) + " " + str(sizes[0]))
   
   # Create grid interpolation
   
   interpolator_face_B_0 = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), face_B[:,:,:,0])
   interpolator_face_B_1 = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), face_B[:,:,:,1])
   interpolator_face_B_2 = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), face_B[:,:,:,2])
   interpolators = [interpolator_face_B_0, interpolator_face_B_1, interpolator_face_B_2]
   
   #######################################################
   if direction == '-':
      multiplier = -1
   else:
      multiplier = 1
   
   points = [x0]
   for i in range(max_iterations):
      previous_point = points[-1]
      B_unit = np.zeros(3)
      B_unit[0] = interpolators[0]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit[1] = interpolators[1]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit[2] = interpolators[2]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit = B_unit / float(np.linalg.norm(B_unit))
      points.append( previous_point + multiplier*B_unit * dx )
   #######################################################
   
   return points

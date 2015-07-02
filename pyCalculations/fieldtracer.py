import numpy as np
import scipy as sp
from scipy import interpolate



def field_tracer( vlsvReader, x0, max_iterations, dx, direction='+' ):
   ''' Field tracer in a static frame

       :param vlsvReader:         An open vlsv file
       :param x:                  Starting point for the field trace
       :param max_iterations:     The maximum amount of iteractions before the algorithm stops
       :param dx:                 One iteration step length
       :param direction:          '+' or '-' Follow field in the plus direction or minus direction
       :returns:                  Field points in array format [x0,x1,x2,x3]
   '''
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

   sizes = np.array([xsize, ysize, zsize])
   maxs = np.array([xmax, ymax, zmax])
   mins = np.array([xmin, ymin, zmin])
   dcell = (maxs - mins)/(sizes.astype('float'))

   # Neglect the indices where there is only 1 cell per direction
   if xsize <= 1:
      indices = [2,1]
   if ysize <= 1:
      indices = [2,0]
   if zsize <= 1:
      indices = [1,0]

   # Read face_B:
   face_B = f.read_variable("B")
   face_Bx = face_B[:,0]
   face_By = face_B[:,1]
   face_Bz = face_B[:,2]

   face_Bx = face_Bx[cellids.argsort()].reshape(sizes[indices])
   face_By = face_By[cellids.argsort()].reshape(sizes[indices])
   face_Bz = face_Bz[cellids.argsort()].reshape(sizes[indices])

   face_B = np.array([face_Bx, face_By, face_Bz])

   # Create x, y, and z coordinates:
   x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
   y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
   z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
   coordinates = np.array([x,y,z])
   # Debug:
   if( len(x) != sizes[0] ):
      print "SIZE WRONG: " + str(len(x)) + " " + str(sizes[0])


   # Create grid interpolation
   interpolator_face_B_0 = interpolate.RectBivariateSpline(coordinates[indices[0]] - 0.5*dcell[indices[0]], coordinates[indices[1]], face_B[indices[0]], kx=2, ky=2, s=0)
   interpolator_face_B_1 = interpolate.RectBivariateSpline(coordinates[indices[0]], coordinates[indices[1]] - 0.5*dcell[indices[1]], face_B[indices[1]], kx=2, ky=2, s=0)
   #interpolator_face_B_2= interpolate.RectBivariateSpline(coordinates[indices[0]], coordinates[indices[1]] - 0.5*dcell[indices[1]], face_B[indices[2]], kx=2, ky=2, s=0)
   interpolators = [interpolator_face_B_0, interpolator_face_B_1]#, interpolator_face_B_2]

   # Print the point data for face_B0 and face_B1
   #print interpolators[0]([x0[indices[0]]], [x0[indices[1]]])
   #print interpolators[1]([x0[indices[0]]], [x0[indices[1]]])
   #print interpolators[2]([x[indices[1]]], [x[indices[0]]])

   #TODO: Cython
   #######################################################
   if direction == '-':
      multiplier = -1
   else:
      multiplier = 1

   points = [x0]
   for i in range(max_iterations):
      previous_point = points[len(points)-2]
      B_unit = np.zeros(3)
      B_unit[indices[0]] = interpolators[0](previous_point[indices[0]], previous_point[indices[1]])
      B_unit[indices[1]] = interpolators[1](previous_point[indices[0]], previous_point[indices[1]])
      B_unit = B_unit / float(np.linalg.norm(B_unit))
      points.append( previous_point + multiplier*B_unit * dx )
   #######################################################

   from vtkwriter import write_vtk_file

   write_vtk_file( "test.vtk", points )




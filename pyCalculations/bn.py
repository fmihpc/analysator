# File which includes calculation for B_n for bowshock

import numpy as np

def rotation_matrix_2d( angle ):
   ''' Creates a rotation matrix that can be used to rotate a 2d vector by an angle in the counter-clockwise direction

       :param angle: Rotation angle
       :returns: The rotation matrix
   '''
   return np.array([[np.cos(angle), -1*np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])


def B_n_bowshock( vlsvReader, points ):
   ''' Calculates B_n vector(s) for given points

       :param vlsvReader: Some VlsvFile with a file opened
       :param points: Some set of points that define the bowshock border

       .. code-block::

          # Example:
          f = VlsvFile("testfile.vlsv")
          B_n_vectors = B_n_bowshock( vlsvReader=f, points=np.array([[2,3,4], [2,3,5], [5, 7, 7]])
          print B_n_vectors

       .. note:: B_n vector will be of length points-1

       .. note:: THIS IS ONLY IN THE X-Y PLANE! THIS WILL NOT WORK WITH SIMULATIONS THAT HAVE FOR EX. Y-Z PLANE
   '''
   from cutthrough import cut_through
   if len(points) < 2:
      print "Invalid point length, points array must be more than 1"
      return

   B_n_vectors = []

   # Keep file open for fast reading
   vlsvReader.optimize_open_file()

   for i in xrange(len(points)-1):
      # Get point 1 and point 2
      point1 = np.array(points[i], copy=False)
      point2 = np.array(points[i+1], copy=False)

      # Calculate the normal unit vector
      n_unit = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
      n_unit = n_unit * np.array([1,1,0])

      # Get the points shifted by unit vector times 3e6
      point1_shifted = point1 + n_unit * 3.0e6
      point2_shifted = point2 + n_unit * 3.0e6
      
     # get the cell ids
      cutthrough = cut_through( vlsvReader, point1_shifted, point2_shifted )
      cellids = cutthrough[0].data
      print len(cellids)
      # Read B values for cellids:
      B_values = vlsvReader.read_variable( name="B", operator="pass", cellids=cellids )

      # Go through each B (if more than one value):
      if len(np.shape(B_values)) == 2:
         for B in B_values:
            B = B * np.array([1,1,0])
            B_n_vectors.append( np.arccos( B.dot( n_unit ) / np.linalg.norm(B)) )
      else:
         # Only one B value:
         B_values = B_values * np.array([1,1,0])
         B_n_vectors.append( np.arccos( B_values.dot( n_unit ) / np.linalg.norm(B_values) ) )

   # Close file
   vlsvReader.optimize_close_file()

   B_n_vectors = np.array( B_n_vectors, copy=False )

   return B_n_vectors














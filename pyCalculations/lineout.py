# Cut-throughs from vlsv files

import numpy as np
import sys

def lineout( vlsvReader, point1, point2, variable, operator="pass",interpolation_order=1, points=100 ):
   ''' Returns cell ids and distances from point 1 for every cell in a line between given point1 and point2

       :param vlsvReader:       Some open VlsvReader
       :type vlsvReader:        :class:`vlsvfile.VlsvReader`
       :param point1:           The starting point of a cut-through line
       :param point2:           The ending point of a cut-through line
       :returns: an 2D array containing distances and values in the following format: [ [distances...] [values] ]

       .. code-block:: python

          Example:
          vlsvReader = VlsvReader(\"testfile.vlsv\")

   '''
   # Transform point1 and point2 into numpy array:
   point1 = np.array(point1)
   point2 = np.array(point2)
   # Get parameters from the file to determine a good length between points (step length):

   # Make sure point1 and point2 are inside bounds
   if vlsvReader.get_cellid(point1) == 0:
      print "ERROR, POINT1 IN CUT-THROUGH OUT OF BOUNDS!"
   if vlsvReader.get_cellid(point2) == 0:
      print "ERROR, POINT2 IN CUT-THROUGH OUT OF BOUNDS!"

   value_len=len(np.atleast_1d(vlsvReader.read_interpolated_variable( variable, point1, operator)))
   
   if value_len==1:
      values=np.zeros(points)
   else:
      values=np.zeros((points,value_len))

   distance=np.zeros(points)
   coordinates=np.zeros((points,3))
   
   for i in range(points):
      relative_coordinate=(point2 - point1) * i / (points-1)
      if interpolation_order==1:
         values[i]=vlsvReader.read_interpolated_variable( variable, point1 + relative_coordinate, operator)
      elif interpolation_order==0:
         values[i]=vlsvReader.read_variable(variable, vlsvReader.get_cellid(point1 + relative_coordinate), operator)

      distance[i]=np.sqrt(sum(relative_coordinate**2))
      coordinates[i]=point1 +  relative_coordinate

   return (distance,coordinates,values)






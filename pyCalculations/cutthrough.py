# Cut-throughs from vlsv files

import numpy as np
import sys

def get_cellids_coordinates_distances( vlsvReader, xmax, xmin, xcells, ymax, ymin, ycells, zmax, zmin, zcells, cell_lengths, point1, point2 ):
   ''' Calculates coordinates to be used in the cut_through. The coordinates are calculated so that every cell gets picked in the coordinates.
       :param vlsvReader:       Some open VlsvReader
       :type vlsvReader:        :class:`vlsvfile.VlsvReader`
       :param point1:           The starting point of a cut-through line
       :param point2:           The ending point of a cut-through line
       :returns: Coordinates of for the cell cut-through as well as cell ids and distances
   '''
   point1 = np.array(point1, copy=False)
   point2 = np.array(point2, copy=False)


   distances = []

   cellids = []

   coordinates = []

   epsilon = sys.float_info.epsilon*2

   iterator = point1
   unit_vector = (point2 - point1) / np.linalg.norm(point2 - point1 + epsilon)
   while True:
      # Get the cell id
      cellid = vlsvReader.get_cellid(iterator)
      if cellid == 0:
         print "ERROR, invalid cell id!"
         return
      # Get the max and min boundaries:
      min_bounds = vlsvReader.get_cell_coordinates(cellid) - 0.5 * cell_lengths
      max_bounds = np.array([min_bounds[i] + cell_lengths[i] for i in range(0,3)])
      # Check which boundary we hit first when we move in the unit_vector direction:

      coefficients_min = np.divide((min_bounds - iterator), unit_vector)
      coefficients_max = np.divide((max_bounds - iterator) , unit_vector)
      # Remove negative coefficients:
      for i in xrange(3):
         if coefficients_min[i] <= 0:
            coefficients_min[i] = sys.float_info.max
         if coefficients_max[i] <= 0:
            coefficients_max[i] = sys.float_info.max



      # Find the minimum coefficient for calculating the minimum distance from a boundary
      coefficient = min([min(coefficients_min), min(coefficients_max)]) * 1.01
      # Append the coordinate:
      coordinates.append( iterator + coefficient * unit_vector )

      # Append the distance:
      distances.append( np.linalg.norm( coordinates[len(coordinates)-1] - point1 ) )

      # Append the cell id:
      cellids.append( vlsvReader.get_cellid( coordinates[len(coordinates)-1] ) )

      # Move the iterator to the next cell. Note: Epsilon is here to ensure there are no bugs with float rounding
      iterator = iterator + coefficient * unit_vector
      # Check to make sure the iterator is not moving past point2:
      if min((point2 - iterator)* unit_vector) < 0:
         break
   # Return the coordinates, cellids and distances for processing
   from output import output_1d
   return output_1d( [np.array(cellids, copy=False), np.array(distances, copy=False), np.array(coordinates, copy=False)], ["CellID", "distances", "coordinates"], ["", "m", "m"] )


def cut_through( vlsvReader, point1, point2 ):
   ''' Returns cell ids and distances from point 1 for every cell in a line between given point1 and point2

       :param vlsvReader:       Some open VlsvReader
       :type vlsvReader:        :class:`vlsvfile.VlsvReader`
       :param point1:           The starting point of a cut-through line
       :param point2:           The ending point of a cut-through line
       :returns: an array containing cell ids, coordinates and distances in the following format: [cell ids, distances, coordinates]

       .. code-block:: python

          Example:
          vlsvReader = VlsvReader(\"testfile.vlsv\")
          cut_through = cut_through(vlsvReader, [0,0,0], [2,5e6,0])
          cellids = cut_through[0]
          distances = cut_through[1]
          print \"Cell ids: \" + str(cellids)
          print \"Distance from point 1 for every cell: \" + str(distances)
   '''
   # Transform point1 and point2 into numpy array:
   point1 = np.array(point1)
   point2 = np.array(point2)
   # Get parameters from the file to determine a good length between points (step length):
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = vlsvReader.read_parameter(name="xcells_ini")
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = vlsvReader.read_parameter(name="ycells_ini")
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = vlsvReader.read_parameter(name="zcells_ini")

   # Make sure point1 and point2 are inside bounds
   if vlsvReader.get_cellid(point1) == 0:
      print "ERROR, POINT1 IN CUT-THROUGH OUT OF BOUNDS!"
   if vlsvReader.get_cellid(point2) == 0:
      print "ERROR, POINT2 IN CUT-THROUGH OUT OF BOUNDS!"

   #Calculate cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   # Get list of coordinates for the cells:
   return get_cellids_coordinates_distances( vlsvReader, xmax, xmin, xcells, ymax, ymin, ycells, zmax, zmin, zcells, cell_lengths, point1, point2 )






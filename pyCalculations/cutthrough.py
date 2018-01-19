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


   distances = [0]

   cellids = [vlsvReader.get_cellid(point1)]

   coordinates = [point1]

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

      # Get the cell id in the new coordinates
      newcoordinate = iterator + coefficient * unit_vector
      newcellid = vlsvReader.get_cellid( newcoordinate )
      # Check if valid cell id:
      if newcellid == 0:
         break
      # Append the cell id:
      cellids.append( newcellid )


      # Append the coordinate:
      coordinates.append( newcoordinate )

      # Append the distance:
      distances.append( np.linalg.norm( newcoordinate - point1 ) )


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
   mesh_limits = vlsvReader.get_spatial_mesh_extent()
   mesh_size = vlsvReader.get_spatial_mesh_size()
   xmax = mesh_limits[3]
   xmin = mesh_limits[0]
   xcells = mesh_size[0]
   # Do the same for y
   ymax = mesh_limits[4]
   ymin = mesh_limits[1]
   ycells = mesh_size[1]
   # And for z
   zmax = mesh_limits[5]
   zmin = mesh_limits[2]
   zcells = mesh_size[2]

   # Make sure point1 and point2 are inside bounds
   if vlsvReader.get_cellid(point1) == 0:
      print "ERROR, POINT1 IN CUT-THROUGH OUT OF BOUNDS!"
   if vlsvReader.get_cellid(point2) == 0:
      print "ERROR, POINT2 IN CUT-THROUGH OUT OF BOUNDS!"

   #Calculate cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   # Get list of coordinates for the cells:
   return get_cellids_coordinates_distances( vlsvReader, xmax, xmin, xcells, ymax, ymin, ycells, zmax, zmin, zcells, cell_lengths, point1, point2 )






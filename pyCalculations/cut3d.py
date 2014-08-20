# This file contains a function for retrieving a 2d cut from a 2d plane ( Retrieves the cell ids in that area )

#NOT IMPLEMENTED YET

import numpy as np

def cut3d( vlsvReader, xmin, xmax, ymin, ymax, zmin, zmax, variable, operator="pass", trim_array=False ):
   ''' Retrieves variables for the given 3d cut

       :param vlsvReader:         Some VlsvReader with a file open
       :type vlsvReader:          :class:`vlsvfile.VlsvReader`
       :param xmin:               The minimum x coordinate of the 2d cut
       :param xmax:               The maximum x coordinate of the 2d cut
       :param ymin:               The minimum y coordinate of the 2d cut
       :param ymax:               The maximum y coordinate of the 2d cut
       :param zmin:               The minimum z coordinate of the 2d cut
       :param zmax:               The maximum z coordinate of the 2d cut
       :param variable:           Some variable to read from the vlsv file
       :param operator:           The variable operator
       :param trim_array:         If true, shapes the array into an array with minimum amount of dimensions, e.g. if the cut is zmax-zmin=0 then this will return a 2d array

       .. code-block:: python

          Example:
          import pytools as pt
          f = pt.vlsvfile.VlsvReader('example.vlsv')
          three_cut = pt.calculations.cut3d( vlsvReader=f, xmin=1e6, xmax=4e6, ymin=1e6, xmax=4e6, zmin=0, zmax=0, variable="rho" )
          import numpy as np
          # Now three_cut is a three-dimensional array (x,y,z), but we can transform it into a 2-d array (x,y) with:
          dimensions = np.shape( three_cut )
          two_cut = np.reshape( three_cut, dimensions[0:2] )

   '''
   # Get min and max coordinates
   min_coordinates = [xmin, ymin, zmin]
   max_coordinates = [xmax, ymax, zmax]
   coordinates = np.array(min_coordinates)

   # Read the cell lengths:
   ##################################################
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
   ##################################################

   cell_lengths = np.array([
                           (xmax-xmin)/(float)(xcells),
                           (ymax-ymin)/(float)(ycells),
                           (zmax-zmin)/(float)(zcells)
                           ])

   # Create 3d array (This array will be returned)
   array_dimensions = np.array([
                               (int)((max_coordinates[0] - min_coordinates[0]) / cell_lengths[0] + 1),
                               (int)((max_coordinates[1] - min_coordinates[1]) / cell_lengths[1] + 1),
                               (int)((max_coordinates[2] - min_coordinates[2]) / cell_lengths[2] + 1)
                               ])
   array = [[np.zeros(array_dimensions[0]) for i in xrange(array_dimensions[1])] for j in xrange(array_dimensions[2])]

   # Optimize file read:
   vlsvReader.optimize_open_file()
   # Loop through the array
   for k in xrange(array_dimensions[2]):
      for j in xrange(array_dimensions[1]):
         for i in xrange(array_dimensions[0]):
            coordinates = np.array([
                                   min_coordinates[0] + i*cell_lengths[0],
                                   min_coordinates[1] + j*cell_lengths[1],
                                   min_coordinates[2] + k*cell_lengths[2]
                                   ])
            #print str(k) + " " + str(j) + " " + str(i) + " " + str(np.shape(array))
            array[k][j][i] = vlsvReader.read_variable(variable, cellids=vlsvReader.get_cellid(coordinates), operator=operator)

   # Close optimization
   vlsvReader.optimize_close_file()

   # Check for trimming:
   if trim_array == True:
      shape = np.shape(array)
      new_shape = []
      # Loop through shape:
      for i in xrange(len(shape)):
         # If the array has no dimensions in some directions, dont add it to the new_shape (e.g. if the cut3d has only one layer of cells in z-direction):
         if shape[i] > 1:
            new_shape.append( shape[i] )
      # Reshape:
      array = np.reshape(array, new_shape)
   return array


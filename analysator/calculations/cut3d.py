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

# This file contains a function for retrieving a 2d cut from a 2d plane ( Retrieves the cell ids in that area )

#NOT IMPLEMENTED YET

import numpy as np
import logging
from analysator.calculations.ids3d import ids3d_box

def cut3d( vlsvReader, xmin, xmax, ymin, ymax, zmin, zmax, variable, operator="pass", trim_array=False ):
   ''' Retrieves variables for the given 3d cut

       :param vlsvReader:         Some VlsvReader with a file open
       :type vlsvReader:          :class:`vlsvfile.VlsvReader`
       :param xmin:               The minimum x coordinate of the 3d cut
       :param xmax:               The maximum x coordinate of the 3d cut
       :param ymin:               The minimum y coordinate of the 3d cut
       :param ymax:               The maximum y coordinate of the 3d cut
       :param zmin:               The minimum z coordinate of the 3d cut
       :param zmax:               The maximum z coordinate of the 3d cut
       :param variable:           Some variable to read from the vlsv file
       :param operator:           The variable operator
       :param trim_array:         If true, shapes the array into an array with minimum amount of dimensions, e.g. if the cut is zmax-zmin=0 then this will return a 2d array

       .. code-block:: python

          Example:
          import analysator as pt
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
   ##################################################

# use this instead...
#def ids3d_box(cellids, low, up, reflevel,
#          xsize, ysize, zsize,
#          spatial_mesh_extent):
#    """Returns lists of CellIDs and the corresponding indices inside a rectangular 3D box.
#    :param cellids:             List of cellids in the simulation box
#    :param low:                 An array holding the lower boundaries of the requested box
#    :param up:                  An array holding the upper boundaries of the requested box
#    :param reflevel:            Highest refinement level in the simulation
#    :param xsize:               Number of cells in x direction
#    :param ysize:               Number of cells in y direction
#    :param zsize:               Number of cells in z direction
#    :param spatial_mesh_extent: Smallest and largest values of each spatial component
#    """



   cids = vlsvReader.read_variable("CellID")
   coords = vlsvReader.get_cell_coordinates(cids)
   
   # start by finding the contained cellids
   ids,indices = ids3d_box(cids, min_coordinates, max_coordinates, vlsvReader.get_max_refinement_level(),
                           xcells, ycells, zcells, mesh_limits)

   # adjust cell_lengths to smallest dxs found
   cell_lengths = np.min(vlsvReader.read_variable("vg_dx", cellids=ids),axis=0)

   # Create 3d array (This array will be returned)
   array_dimensions = np.array([
                               (int)((max_coordinates[0] - min_coordinates[0]) / cell_lengths[0] + 1),
                               (int)((max_coordinates[1] - min_coordinates[1]) / cell_lengths[1] + 1),
                               (int)((max_coordinates[2] - min_coordinates[2]) / cell_lengths[2] + 1)
                               ])

   # Create the uniform coordinate array:
   max_coords_adjusted = min_coordinates + array_dimensions*cell_lengths
   Xs,Ys,Zs = np.meshgrid(np.linspace(min_coordinates[0],max_coords_adjusted[0], array_dimensions[0], endpoint=False),
                          np.linspace(min_coordinates[1],max_coords_adjusted[1], array_dimensions[1], endpoint=False),
                          np.linspace(min_coordinates[2],max_coords_adjusted[2], array_dimensions[2], endpoint=False),
                              indexing="ij")

   coords_lin = np.reshape(np.stack((Xs,Ys,Zs),axis=-1), (-1,3))

   cellids_array = vlsvReader.get_cellid(coords_lin)
   data_arr = vlsvReader.read_variable(variable, cellids=cellids_array, operator=operator)
   test_var = vlsvReader.read_variable(variable, cellids=[1], operator=operator)
   var_dims = test_var.shape

   data_arr = np.reshape(data_arr, (*array_dimensions,*var_dims))
   # Transpose to fit the previous implementation
   data_arr = data_arr.transpose((2,1,0, *[i+3 for i in range(len(var_dims))]))

   # Check for trimming:
   if trim_array == True:
      shape = np.shape(data_arr)
      new_shape = []
      # Loop through shape:
      for i in range(len(shape)):
         # If the array has no dimensions in some directions, dont add it to the new_shape (e.g. if the cut3d has only one layer of cells in z-direction):
         if shape[i] > 1:
            new_shape.append( shape[i] )
      # Reshape:
      data_arr = np.reshape(data_arr, new_shape)
   return data_arr


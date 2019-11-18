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

import numpy as np
import pylab as pl
from getcellid import *

def get_cell_indices( bounds, cellid ):
   cellid = (int)(cellid - 1)
   # Get xmin, xmax, etc
   xmin = bounds[0]
   ymin = bounds[0 + 3]
   zmin = bounds[0 + 3 + 3]
   
   xmax = bounds[1]
   ymax = bounds[1 + 3]
   zmax = bounds[1 + 3 + 3]

   xcells = bounds[2]
   ycells = bounds[2 + 3]
   zcells = bounds[2 + 3 + 3]

   # Get cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   cellindices = np.zeros(3)

   cellindices[0] = (int)(cellid)%(int)(xcells)
   cellindices[1] = ((int)(cellid)/(int)(xcells))%(int)(ycells)
   cellindices[2] = (int)(cellid)/(int)(xcells*ycells)

   return cellindices.astype(int)

#def get_cell_id_from_indices( bounds, cellindices ):
#   # Get xmin, xmax, etc
#   xmin = bounds[0]
#   ymin = bounds[0 + 3]
#   zmin = bounds[0 + 3 + 3]
#   
#   xmax = bounds[1]
#   ymax = bounds[1 + 3]
#   zmax = bounds[1 + 3 + 3]
#
#   xcells = bounds[2]
#   ycells = bounds[2 + 3]
#   zcells = bounds[2 + 3 + 3]
#
#   return cellindices[0] + cellindices[1]
   

#def get_3d_cellids( bounds, BBOX ):
#   # Get xmin, xmax, etc
#   xmin = bounds[0]
#   ymin = bounds[0 + 3]
#   zmin = bounds[0 + 3 + 3]
#   
#   xmax = bounds[1]
#   ymax = bounds[1 + 3]
#   zmax = bounds[1 + 3 + 3]
#
#   xcells = bounds[2]
#   ycells = bounds[2 + 3]
#   zcells = bounds[2 + 3 + 3]
#
#   # Get cell lengths:
#   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
#
#   minCellId = get_cell_id(bounds, [BBOX[0], BBOX[2], BBOX[4]])
#   maxCellId = 

def get_2d_array( fileNames, variables, BBOX ):
   '''Fetches the cell's data from filenames
      :param fileNames         List of the files
      :param variables         list of variables
      :param BBOX              boundary box, (=[xmin, xmax, ymin, ymax, zmin, zmax])
      :returns array of variables
   '''
   cellids = []

   fileNames = np.atleast_1d(fileNames)
   variables = np.atleast_1d(variables)

   # Get the cell ids inside the boundary box:
   vlsvReader = VlsvReader(fileNames[0])
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

   # Get bounds for reading cell ids from coordinates
   bounds = [xmin, xmax, xcells, ymin, ymax, ycells, zmin, zmax, zcells]

   # Get cell lengths
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   array = {}
   #Get cell ids:
   #Iterate through the BBOX coordinates:
   minCoordinates = [BBOX[0], BBOX[2], BBOX[4]]
   maxCoordinates = [BBOX[1], BBOX[3], BBOX[5]]
   coordinates = np.array(minCoordinates)

   xindex = 0
   yindex = 0
   zindex = 0

   while coordinates[2] <= maxCoordinates[2]:
      while coordinates[1] <= maxCoordinates[1]:
         while coordinates[0] <= maxCoordinates[0]:
            cellid = get_cell_id( bounds, coordinates )
            cellids.append(cellid)
            # Move to the next cell in x direction
            coordinates[0] = coordinates[0] + cell_lengths[0]
         # Move to the next cell in y direction
         coordinates[1] = coordinates[1] + cell_lengths[1]
         coordinates[0] = minCoordinates[0]
      # Move to the next cell in z direction
      coordinates[2] = coordinates[2] + cell_lengths[2]
      coordinates[1] = minCoordinates[1]
      coordinates[0] = minCoordinates[0]
   print("")
   # Create an array for holding variables:
   maxIndices = get_cell_indices( bounds, max(cellids) )
   minIndices = get_cell_indices( bounds, min(cellids) )

   # Create an array for holding variables
   two_d_variables = []
   for i in range(maxIndices[1] - minIndices[1] + 1):
      two_d_variables.append(np.zeros((maxIndices[0] - minIndices[0] + 1)))
   two_d_variables = np.array(two_d_variables)

   # Input positions of the cell ids
   positions = {}
   for i in cellids:
      positions[i] = get_cell_indices( bounds, i ) - minIndices

   variable_dict = {}
   for i in variables:
      variable_dict[i] = []

   print(cellids)
   # Read variables:
   for f in fileNames:
      # Open file
      vlsvReader = VlsvReader(f)
      # Read variables:
      for i in variables:
         # Read in the variable
         variableArray = vlsvReader.read_variables_for_cellids(i, cellids)
         # Input into correct positions:
         for j in range(len(cellids)):
            position = positions[cellids[j]]
            two_d_variables[position[1]][position[0]] = variableArray[j]
         # Input variables:
         variable_dict[i].append(np.array(two_d_variables))
   # Return variables:
   return variable_dict

def fourier_2d_array( variables, showplots=False, saveplots="none", showraw=False, kaiserwindowparameter=0, t_start=0, dt=1.0 ):
   # Go through all the variables:
   for i in variables.items():
      variable = i[0]
      arrays = i[1]
      index = 0
      for array in arrays:
         pl.figure()
         # Plot the image
         if showraw == True:
            pl.subplot(2,1,1)
         pl.imshow(np.abs(np.fft.fftshift(np.fft.fft2(np.outer(np.kaiser(len(array), kaiserwindowparameter),np.kaiser(len(array[0]), kaiserwindowparameter))*(array-np.mean(array))))), interpolation="Nearest")
         if showraw == True:
            pl.subplot(2,1,2)
            pl.imshow(array, interpolation="Nearest")
         # Save the image:
         pl.title(variable + " " + str(t_start + index*dt))
         if saveplots != "none":
            pl.savefig(saveplots + "_" + variable + "_" + str(index) + ".png")
         if showplots == False:
            pl.close()
         index = index + 1
   if showplots == True:
      pl.show()
   


















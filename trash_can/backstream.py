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


def get_rho_nonbackstream( vlsvReader, cellid, radius ):
   ''' Calculates the non backstream population's (solar wind's) contributions on rho
       :param vlsvReader          Some VlsvReader with a file open
       :param cellid              The cellid whose rho to calculate
       Returns the value of rho
   '''
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   # Get cells:
   vcellids = velocity_cell_data.keys()
   # Get avgs data:
   avgs = velocity_cell_data.values()
   # Get a list of velocity coordinates shifted by the solar wind bulk velocity:
   V_sw = np.array([-500000,0,0])
   v = vlsvReader.get_velocity_cell_coordinates(vcellids) - V_sw
   # Get sum of radiuses (to the power of second):
   radiuses = np.sum(np.abs(v)**2,axis=-1)
   # Go through velocity cells and calculate the data:
   #radius2 = 468621**2
   radius2 = radius**2
   rho_nonbackstream = 0
   for i in xrange(len(radiuses)):
      if radiuses[i] <= radius2:
         rho_nonbackstream = rho_nonbackstream + avgs[i]

   # Get the volume of a velocity cell:
   vxblocks = vlsvReader.read_parameter("vxblocks_ini")
   vxmin = vlsvReader.read_parameter("vxmin")
   vxmax = vlsvReader.read_parameter("vxmax")
   DV3 = ((vxmax-vxmin)/((float)(vxblocks)*4))**3
   # Return the value
   return rho_nonbackstream*DV3

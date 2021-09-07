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
import scipy as sp
import pytools as pt
import warnings
from scipy import interpolate

def fg_PoyntingFlux(bulkReader):
   b = bulkReader.read_fg_variable_as_volumetric('fg_b')
   print('b.shape=',b.shape)
   e = bulkReader.read_fg_variable_as_volumetric('fg_e')
   print('e.shape=',e.shape)

   mu_0 = 1.25663706144e-6
   S = np.cross(e,b)/mu_0
   return S


def fg_divPoynting(bulkReader, dx=1e6):
   S = fg_PoyntingFlux(bulkReader)

   #divS = pt.plot.plot_helpers.numdiv(pt.plot.plot_helpers.TransposeVectorArray(np.squeeze(S))).T
   divS = (np.roll(S[:,:,:,0],-1, 0) - np.roll(S[:,:,:,0], 1, 0) +
         np.roll(S[:,:,:,1],-1, 1) - np.roll(S[:,:,:,1], 1, 1) +
         np.roll(S[:,:,:,2],-1, 2) - np.roll(S[:,:,:,2], 1, 2))/(2*dx)
   print('divS.shape', divS.shape)
   #divS = np.div

#   print(divS)
   return divS

def fg_vol_jacobian(array, dx=1e6):
   dFx_dx, dFx_dy, dFx_dz = np.gradient(array[:,:,:,0], dx)
   dFy_dx, dFy_dy, dFy_dz = np.gradient(array[:,:,:,1], dx)
   dFz_dx, dFz_dy, dFz_dz = np.gradient(array[:,:,:,2], dx)

def fg_vol_curl(array, dx=1e6):
   dummy,  dFx_dy, dFx_dz = np.gradient(array[:,:,:,0], dx)
   dFy_dx, dummy,  dFy_dz = np.gradient(array[:,:,:,1], dx)
   dFz_dx, dFz_dy, dummy  = np.gradient(array[:,:,:,2], dx)

   rotx = dFz_dy - dFy_dz
   roty = dFx_dz - dFz_dx
   rotz = dFy_dx - dFx_dy

   return np.stack([rotx, roty, rotz], axis=-1)
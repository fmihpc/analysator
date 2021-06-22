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

def get_fg_divPoynting:
   b = bulkReader.read_fg_variable_as_volumetric('fg_b')
   print('b.shape=',b.shape)
   e = bulkReader.read_fg_variable_as_volumetric('fg_e')
   print('e.shape=',e.shape)

   mu_0 = 1.25663706144e-6
   S = np.cross(e,b)/mu_0
   print('S.shape', S.shape)

   #divS = pt.plot.plot_helpers.numdiv(pt.plot.plot_helpers.TransposeVectorArray(np.squeeze(S))).T
   divS = (np.roll(S[:,:,:,0],-1, 0) - np.roll(S[:,:,:,0], 1, 0) +
         np.roll(S[:,:,:,1],-1, 1) - np.roll(S[:,:,:,1], 1, 1) +
         np.roll(S[:,:,:,2],-1, 2) - np.roll(S[:,:,:,2], 1, 2))/(2*1e6)
   print('divS.shape', divS.shape)
   #divS = np.div

   fig=plt.pyplot.figure()
   ax=plt.pyplot.axes()
   im=ax.imshow(divS[:,int(divS.shape[1]/2),:], origin='lower', norm=plt.colors.SymLogNorm(base=10, linthresh=1e-12, linscale = 1.0, vmin=-1e-7, vmax=1e-7, clip=True), cmap='vik')
   plt.pyplot.colorbar(im, ax=ax)
   plt.pyplot.savefig('./divS_centered.png')
   print(divS)
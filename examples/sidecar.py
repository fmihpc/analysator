# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2022 University of Helsinki
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


# This file contains an example on how to use the vlsvReader and vlsvWriter
# functionalities to derive quantities (especially differential operations 
# on an AMR mesh via uniform grid resampling) from .vlsv variables and store
# them to a sidecar file AMR grid.
#
# The file first defines a reader for input data and a writer for the output,
# and copies over a list of variables desired for the output (CellID is 
# definitely worth copying over). A definition for curl on fsgrid is given,
# after which fg_b (a face-centered field) is re-centered as volumetric. J is
# then calculated on fsgrid from the volumetric B, which is then downsampled 
# onto SpatialGrid via averaging and written into the sidecar file.
#
# A Script like this can be adapter for personal use, and
# passed to a job submission system for generating a multitude of sidecars
# via e.g. array jobs.

import pytools as pt
import numpy as np

q = 1.60217662e-19
me = 9.10938356e-31
mp = 1.6726219e-27
eps0 = 8.854e-12 

import time

# Change filenames and paths to your liking
fn = '/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/bulk1.0001450.vlsv'
f = pt.vlsvfile.VlsvReader(fn)
writer = pt.vlsvfile.VlsvWriter(f, 'pysidecar_bulk1.0001450.vlsv')

print(f.get_all_variables())
print(f)
cellIds=f.read_variable("CellID")
argsorti=f.read_variable("CellID").argsort()


fgsize=f.get_fsgrid_mesh_size()
fgext=f.get_fsgrid_mesh_extent()
t = time.time()
print("copy_var_list")
writer.copy_variables_list(f, ["CellID"])
print('writer init and var copy, elapsed:', time.time()-t)


def fg_vol_curl(reader, array):
   dx = reader.get_fsgrid_cell_size()
   dummy,  dFx_dy, dFx_dz = np.gradient(array[:,:,:,0], *dx)
   dFy_dx, dummy,  dFy_dz = np.gradient(array[:,:,:,1], *dx)
   dFz_dx, dFz_dy, dummy  = np.gradient(array[:,:,:,2], *dx)

   rotx = dFz_dy - dFy_dz
   roty = dFx_dz - dFz_dx
   rotz = dFy_dx - dFx_dy

   return np.stack([rotx, roty, rotz], axis=-1)

print('Processing for J, elapsed', time.time()-t)
fg_b_vol=f.read_fg_variable_as_volumetric('fg_b')
print('fgbvol read, shape ', fg_b_vol.shape, 'elapsed:', time.time()-t)

fg_J = fg_vol_curl(f, fg_b_vol)/1.25663706144e-6
print('fgbvol curled, elapsed:', time.time()-t)

f.map_vg_onto_fg()
print('vg mapped to fg, elapsed:', time.time()-t)

vg_J = f.fsgrid_array_to_vg(fg_J)
print('fg_J mapped to vg, elapsed:', time.time()-t)

writer.write(vg_J,'vg_J','VARIABLE','SpatialGrid')
print('J written, elapsed:', time.time()-t)
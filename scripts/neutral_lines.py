# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2024 University of Helsinki
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

import pytools as pt
import numpy as np
import sys,time

file_id = int(sys.argv[1])



path = "/wrk-vakka/users/mjalho/xo-paper/KomarRepro/theta150/"
fn = "jacobs.vlsv"
f = pt.vlsvfile.VlsvReader(path+fn)

#FHA
fnout = path+"xo.vlsv"

cids = f.read_variable("CellID")



t = time.time()

fw = pt.vlsvfile.VlsvWriter(f, fnout, copy_meshes="SpatialGrid")

#fw.copy_variables(f,["CellID","vg_b_vol","LMN_magnetopause/vg_L","LMN_magnetopause/vg_N","LMN_magnetopause/vg_jacobian"])
fw.copy_variables(f,["CellID","vg_b_vol"])
fw.write(f.read_variable("vg_j"),"vg_j","VARIABLE", "SpatialGrid")
fw.write(f.read_variable("vg_dxs"), "vg_dxs", "VARIABLE", "SpatialGrid")
LMNs = f.read_variable("vg_lmn",cids)
fw.write(LMNs.reshape((-1,9)), "vg_LMN", "VARIABLE","SpatialGrid")
fw.write(f.read_variable("vg_mdd_dimensionality"), "vg_MDD_dimensionality", "VARIABLE", "SpatialGrid")
fw.write(f.read_variable("vg_lmn_neutral_line_distance"), "vg_LN_null_line_distance", "VARIABLE", "SpatialGrid")
fw.write(f.read_variable("vg_lmn_L_flip_distance"), "vg_L_flip_distance", "VARIABLE", "SpatialGrid")

# print(f.read_variable("vg_lmn_neutral_line_distance",cids))
LMN_jacob = f.read_variable("vg_jacobian_B", cids)
LMN_jacob = np.reshape(LMN_jacob,(LMN_jacob.shape[0],3,3))
LMN_jacob = np.transpose(LMNs,(0, 2, 1)) @ LMN_jacob @ LMNs

dBNdL = LMN_jacob[:,2,0]
fw.write(dBNdL, "vg_dBNdL", "VARIABLE", "SpatialGrid")

dBLdN = LMN_jacob[:,0,2]
fw.write(dBLdN, "vg_dBLdN", "VARIABLE", "SpatialGrid")
print('things written, elapsed:', time.time()-t)


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

import pytools as pt
import numpy as np
import sys,time

file_id = int(sys.argv[1])



#path = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/visualizations/lmn/"
path = "/wrk-vakka/users/mjalho/xo-paper/KomarRepro/theta150/"
#fn = "jlsidecar_mva_bulk1.{:07d}.vlsv".format(file_id)
fn = "jacobs.vlsv"
f = pt.vlsvfile.VlsvReader(path+fn)

#FHA
#fn = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/pysidecar_bulk1.{:07d}.vlsv".format(file_id)
#f = pt.vlsvfile.VlsvReader(fn)

#fnout = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/visualizations/lmn/pyXO4/pyXO_bulk1.{:07d}.vlsv".format(file_id)
fnout = path+"xo.vlsv"
# fnout = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/visualizations/lmn/pyXO5/pyXO_bulk1.{:07d}.vlsv".format(file_id)

#fnout = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/XO/pyXO_bulk1.{:07d}.vlsv".format(file_id)

cids = f.read_variable("CellID")
# ci = np.argwhere(cids==355282689)

#cids = [355282688,355282689,355282690,355282691]
#cids =-1
#cid = 355282689
# print(f.read_variable_info("vg_jacobian_B",cid))
# print(f.read_variable("vg_gtg",cid))
# print(f.read_variable("vg_ggt",cid))
# # print(f.read_variable("vg_b_vol",cid))


# print(f.read_variable("vg_J",cid))


# print(f.get_cell_coordinates(cid)/6371e3)

# print(f.read_variable("LMN_magnetopause/vg_L",cid[0]))
# print(f.read_variable("LMN_magnetopause/vg_L",cid[1]))


# print("lmn output",f.read_variable("vg_dxs",cid))
# # print("vg_L",f.read_variable("LMN_magnetopause/vg_L",cid))
# # print("vg_M",f.read_variable("LMN_magnetopause/vg_M",cid))
# # print("vg_N",f.read_variable("LMN_magnetopause/vg_N",cid))
# # print("Ds",f.read_variable("vg_mdd_dimensionality",cid))
# print("XO_d",f.read_variable("vg_lmn_neutral_line_distance",cid))

#sys.exit()


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

sys.exit()

LMN_jacob = f.read_variable("LMN_magnetopause/vg_jacobian")
L = f.read_variable("LMN_magnetopause/vg_L")
M = f.read_variable("LMN_magnetopause/vg_M")
N = f.read_variable("LMN_magnetopause/vg_N")
B = f.read_variable("vg_b_vol")
exts = f.get_spatial_mesh_extent()
exts = np.array([exts[3]-exts[0],exts[4]-exts[1],exts[5]-exts[2]])

BL = np.sum(L*B, axis=-1)
BM = np.sum(M*B, axis=-1)
BN = np.sum(N*B, axis=-1)
dotp = np.sum(L*N, axis=-1)

gradBL = LMN_jacob[:,0:3]
gradBM = LMN_jacob[:,3:6]
gradBN = LMN_jacob[:,6:9]

gradBLn = np.linalg.norm(gradBL,axis=-1)
gradBMn = np.linalg.norm(gradBM,axis=-1)
gradBNn = np.linalg.norm(gradBN,axis=-1)

dBNdL = LMN_jacob[:,6]
det = np.zeros_like(BL)

# Distance to zero plane for BL an BN
sL = BL/(gradBLn*1e-9)
sM = BM/(gradBMn*1e-9)
sN = BN/(gradBNn*1e-9)

ci = np.argwhere(cids==355282689)

#cids = [355282689]
dx = np.zeros_like(B)
for i, cid in enumerate(cids):
   dx[i]= f.get_cell_dx(cid)


# Vectors to zero intercept plane
L_zero_intercept = gradBL/np.broadcast_to(gradBLn,(3,cids.size)).transpose()
L_zero_intercept = L_zero_intercept*np.broadcast_to(sL,(3,cids.size)).transpose()/dx

M_zero_intercept = gradBM/np.broadcast_to(gradBMn,(3,cids.size)).transpose()
M_zero_intercept = M_zero_intercept*np.broadcast_to(sM,(3,cids.size)).transpose()/dx

N_zero_intercept = gradBN/np.broadcast_to(gradBNn,(3,cids.size)).transpose()
N_zero_intercept = N_zero_intercept*np.broadcast_to(sN,(3,cids.size)).transpose()/dx

# the intercept vectors are also normals to the planes of zero value; 
n_line = np.cross(L_zero_intercept,N_zero_intercept, axis=-1)

#Find a line intercept point and its norm
n_line_intercept = (L_zero_intercept + N_zero_intercept)
s_line = np.linalg.norm(n_line_intercept,axis=-1)


det = (np.all(np.abs(L_zero_intercept) < 0.5,axis=-1) & np.all(np.abs(N_zero_intercept) < 0.5,axis=-1)) # not quite accurate, but let's try

# for i, cid in enumerate(cids):
#    det[i] = dBNdL[i]*float(np.all(np.abs(L_zero_intercept[i]) < 0.5) and np.all(np.abs(N_zero_intercept[i]) < 0.5)) # not quite accurate, but let's try

print(det[ci], s_line[ci])
#fw.write(det, "XO_determinant", "VARIABLE", "SpatialGrid")

fw.write(s_line, "XO_lined", "VARIABLE", "SpatialGrid")
fw.write(L_zero_intercept, "L0_distance", "VARIABLE", "SpatialGrid")
fw.write(M_zero_intercept, "M0_distance", "VARIABLE", "SpatialGrid")
fw.write(N_zero_intercept, "N0_distance", "VARIABLE", "SpatialGrid")
fw.write(dBNdL, "dBNdL", "VARIABLE", "SpatialGrid")


def test():
   cid = 355282689

   LMN_jacob = np.reshape(f.read_variable("LMN_magnetopause/vg_jacobian", cid),(3,3))
   # dBNdl is -8.54339e-07
   dBNdL = LMN_jacob[2,0] # as in VisIt exprs!

   exts = f.get_spatial_mesh_extent()
   exts = np.array([exts[3]-exts[0],exts[4]-exts[1],exts[5]-exts[2]])
   dx = (exts/f.get_spatial_mesh_size())*2**-f.get_amr_level(cid)
   coords = f.get_cell_coordinates(cid)

   L = f.read_variable("LMN_magnetopause/vg_L", cid)
   M = f.read_variable("LMN_magnetopause/vg_M", cid)
   N = f.read_variable("LMN_magnetopause/vg_N", cid)

   B = f.read_variable("vg_b_vol", cid)

   BL = np.dot(L,B)
   BM = np.dot(M,B)
   BN = np.dot(N,B)

   print(BL,BM,BN)

   gradBL = LMN_jacob[0,:]
   print("test gradBL", gradBL)
   sL = np.linalg.norm(gradBL)/BL
   print("test sL", sL)

   gradBN = LMN_jacob[2,:]
   print("test gradBN", gradBN)

   sN = np.linalg.norm(gradBN)/BN
   print("test sN", sN)


   print(dx)
   print(sL)

   L_zero_intercept = sL*gradBL/np.linalg.norm(gradBL)/dx
   print("test L_zero_intercept", L_zero_intercept)
   L_zero_in_cell = np.all(np.abs(L_zero_intercept) < 0.5)

   N_zero_intercept = sN*gradBN/np.linalg.norm(gradBN)/dx
   N_zero_in_cell = np.all(np.abs(N_zero_intercept) < 0.5)

   print(L_zero_in_cell,N_zero_in_cell, dBNdL/np.abs(dBNdL))

#test()
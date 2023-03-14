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

''' A file for finding null lines and miscellaneous data
'''
import numpy as np
import sys

#def fBox
   
def LMN_null_lines_FOTE(LMNs, jacobs, Bs, dxs, coords):
   ''' Function that uses a linear approximation of B from B and its Jacobian
   in the LMN coordinates to get an approximate distance to the neutral line.
   inputs:
   * LMN basis vectors stack
   * Jacobian of B in 9-element vector stack
   * vg_b_vol

   return a measure for closeness to the cell center for a neutral line and other data.
   For the point d closest to the cell center on the line:
      case a: d within the cube, return -1/d
      case b: d outside the cube, return d
      This keeps a monotonous measure with return values < 0 having certainly the line within 
      the cells, while values > 0 showing how far off a linear model places the line.
      Todo:
      Proper Ursian SDF/SDB
   '''
   stack = True
   if(Bs.shape == (3,)):  # I want these to be stacks
      stack = False
      LMNs = np.array([LMNs])
      jacobs = np.array([jacobs])
      dxs = np.array(dxs)
      Bs = np.array([Bs])

   jacobs = np.reshape(jacobs,(jacobs.shape[0],3,3))
   Ls = LMNs[:,:,0]
   Ms = LMNs[:,:,1]
   Ns = LMNs[:,:,2]
   n_cells = Bs.shape[0]
   #rot = [basis_L[1:3,i] basis_M[1:3,i] basis_N[1:3,i]]
   #so basically just LMN
   # rotate the jacobians
   jacobs = np.transpose(LMNs,(0, 2,1)) @ jacobs @ LMNs

   BL = np.sum(Ls*Bs, axis=-1)
   BM = np.sum(Ms*Bs, axis=-1)
   BN = np.sum(Ns*Bs, axis=-1)
   unit = 1 #1e-9
   gradBL = jacobs[:, 0,:]*unit
   #gradBM = jacobs[1,:]
   gradBN = jacobs[:, 2,:]*unit

   gradBLn = np.linalg.norm(gradBL,axis=-1)
   gradBNn = np.linalg.norm(gradBN,axis=-1)

   gradBL = gradBL/np.broadcast_to(gradBLn,(3,n_cells)).transpose()
   gradBN = gradBN/np.broadcast_to(gradBNn,(3,n_cells)).transpose()

   # Distance to zero plane for BL an BN; nb. jacobian in nT/m..
   sL = -BL/(gradBLn)
   #sM = BM/(gradBMn*1e-9)
   sN = -BN/(gradBNn)

   n_line = np.cross(gradBL,gradBN, axis=-1)
   n_line = n_line/np.broadcast_to(np.linalg.norm(n_line,axis=-1),(3,n_cells)).transpose()

   #Find a line intercept point and its norm
   n_line_intercept = np.zeros_like(gradBLn)#(L_zero_intercept + N_zero_intercept) # well this is just wrong when not perp
   n_line_intercept.fill(np.nan)

   dots = np.sum(gradBL*gradBN, axis=-1)
   c1 = sL - sN*(dots)/(1 - dots**2)
   c2 = sN - sL*(dots)/(1 - dots**2)
   n_line_intercept = (gradBL*np.broadcast_to(c1,(3,n_cells)).transpose() + gradBN*np.broadcast_to(c2,(3,n_cells)).transpose())/dxs

   #rotate the n_line_intercept to xyz instead of LMN
   nn = n_line_intercept[:,np.newaxis,:]

   n_line_intercept = np.matmul(nn, LMNs).squeeze() # inverse of rotation for row vectors
   n_line = np.matmul(n_line[:,np.newaxis,:], LMNs).squeeze()

   # Make certain that we have a point along the line that is closest to the cell center
   par_dist = np.sum(n_line_intercept*n_line, axis=-1)[:,np.newaxis]
   n_line_intercept = n_line_intercept - n_line*np.broadcast_to(np.sum(n_line_intercept*n_line, axis=-1),(3,n_cells)).transpose() # is it not nearest to 0 already?
   s_line = np.linalg.norm(n_line_intercept,axis=-1)

   #Let's go back a dx or two so we need only to look at positive intercepts
   n_line_intercept = n_line_intercept - n_line*2*dxs


# This is horrible, I sort of wish we had spherical shells
   # mirror across the origin on all axes if on negative hemisphere -
   # easier to check for intersection only with xyz = 0.5 planes
   nsz = (3, n_line_intercept.shape[0],n_line_intercept.shape[1])
   positive_intercepts_0 = np.ones(nsz)
   positive_intercepts_1 = np.ones(nsz)
   positive_intercepts_0.fill(np.nan)
   positive_intercepts_1.fill(np.nan)
   lineintercept = n_line_intercept.copy()
   linevec = n_line/np.broadcast_to(np.linalg.norm(n_line,axis=-1),(3,n_cells)).transpose()
   
   # find distances to the planes delimiting the cell along the line
   
   d_fw_planes_0 = np.divide(-0.5-lineintercept,linevec)
   d_fw_planes_1 = np.divide(0.5-lineintercept,linevec)

   hits = np.all(np.abs(n_line_intercept) < 0.5, axis=-1) # these are certain, also init hits array
   eps = 1e-12

   # fw intercepts
   for d in [0,1,2]:
      ds = (np.repeat(np.array([d_fw_planes_0[:,d]]).transpose(),3,axis=1))
      
      ds = linevec*ds
      
      new_intercepts = lineintercept + ds
      mask = np.all(np.abs(new_intercepts) <= 0.5+eps, axis=-1)
      positive_intercepts_0[d,mask,:] = new_intercepts[mask,:]
      hits = hits | mask

      ds = (np.repeat(np.array([d_fw_planes_1[:,d]]).transpose(),3,axis=1))
      ds = linevec*ds
      new_intercepts = lineintercept + ds
      mask = np.all(np.abs(new_intercepts) <= 0.5+eps, axis=-1)
      positive_intercepts_1[d,mask,:] = new_intercepts[mask,:]
      hits = hits | mask
   
   np.divide(-1, s_line,out=s_line, where=hits)
   # np.add(s_line, 1,out=s_line, where=hits) # so that barely hitting approaches 0-
   # np.subtract(s_line, 0.5, out=s_line, where=np.logical_not(hits)) # so that barely missing approaches 0+

# these are a bit fiddly yet?

   tmpout="/proj/mjalho/analysator/scripts/tmp.dat"
   # print(n_line_intercept)
   # print(n_line_intercept*dxs)
   a = coords+(n_line_intercept)*dxs
   b = n_line
   # print("a",a)
   # print("coords",coords)
   # print("b",b)
   # print(s_line)
   # print(b)
   # print(n_line_intercept)
   # print(dxs)
   par_dist = np.sum((n_line_intercept*dxs)*n_line, axis=-1)[:,np.newaxis]
   # print(par_dist)
   stck = np.hstack((a,b,s_line[:,np.newaxis],coords,-(n_line_intercept)*dxs,
                  par_dist))
      #print(stck)
      #np.save(tmpout, stck[np.all(np.isfinite(stck),axis=1) & (s_line < 1),:]) # get close hits
   if stack:
      return s_line, stck
   else:
      return s_line[0], stck[0]

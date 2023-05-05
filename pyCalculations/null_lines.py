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

# Box: correct distance to corners hg_sdf
def fBox(pp, bb=np.array([0.5,0.5,0.5])):
   p = np.atleast_2d(pp)
   b = np.atleast_2d(bb)
   b = np.broadcast_to(b, p.shape)
   d = np.abs(p) - b
   vec30 = np.zeros_like(d)
   return np.linalg.norm(np.fmax(d, vec30),axis=-1) + np.amax(np.fmin(d, vec30),axis=-1)


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
   unit = 1#1e-9
   gradBL = jacobs[:, 0,:]*unit
   #gradBM = jacobs[1,:]
   gradBN = jacobs[:, 2,:]*unit

   gradBLn = np.linalg.norm(gradBL,axis=-1)
   gradBNn = np.linalg.norm(gradBN,axis=-1)
   mask = np.isfinite(gradBLn) & np.isfinite(gradBNn) # require sane gradients
   mask = mask & (gradBLn >0) & (gradBNn >0)          # require nonzero gradients
   mask = mask & np.all(np.isfinite(Ls), axis=-1)     # these should be covered by above..
   mask = mask & np.all(np.isfinite(Ms), axis=-1)
   mask = mask & np.all(np.isfinite(Ns), axis=-1)


   gradBL = gradBL/np.broadcast_to(gradBLn,(3,n_cells)).transpose()
   gradBN = gradBN/np.broadcast_to(gradBNn,(3,n_cells)).transpose()
   dots = np.sum(gradBL*gradBN, axis=-1)
   mask = mask & (dots != 1)
   # Distance to zero plane for BL an BN; nb. jacobian in nT/m..
   sL = -BL/(gradBLn)
   #sM = BM/(gradBMn*1e-9)
   sN = -BN/(gradBNn)

   n_line = np.full(gradBL.shape,np.nan)

   n_line[mask,:] = np.cross(gradBL,gradBN, axis=-1)[mask,:]
   n_line = n_line/np.broadcast_to(np.linalg.norm(n_line,axis=-1),(3,n_cells)).transpose()

   #Find a line intercept point and its norm
   n_line_intercept = np.ones_like(gradBL)#(L_zero_intercept + N_zero_intercept) # well this is just wrong when not perp
   print("Shape of n_line_intercept", n_line_intercept.shape)
   n_line_intercept.fill(np.inf)

   c1 = sL - sN*(dots)/(1 - dots**2)
   c2 = sN - sL*(dots)/(1 - dots**2)
   print(len(mask) - np.sum(mask), "bad points")
   n_line_intercept[mask,:] = (((gradBL*np.broadcast_to(c1,(3,n_cells)).transpose() + gradBN*np.broadcast_to(c2,(3,n_cells)).transpose())/dxs))[mask,:]

   #rotate the n_line_intercept to xyz instead of LMN
   nn = n_line_intercept[:,np.newaxis,:]

   n_line_intercept[mask,:] = np.matmul(nn, LMNs).squeeze()[mask,:] # inverse of rotation for row vectors
   n_line[mask,:] = np.matmul(n_line[:,np.newaxis,:], LMNs).squeeze()[mask,:]

   # Make certain that we have a point along the line that is closest to the cell center
   par_dist = np.sum(n_line_intercept*n_line, axis=-1)[:,np.newaxis]
   n_line_intercept = n_line_intercept - n_line*np.broadcast_to(np.sum(n_line_intercept*n_line, axis=-1),(3,n_cells)).transpose() # is it not nearest to 0 already?
   #s_line = np.linalg.norm(n_line_intercept,axis=-1)
   s_line = np.full(gradBLn.shape,np.inf)
   s_line[mask] = fBox(n_line_intercept, np.array([0.5,0.5,0.5]))[mask] # Proper SDF, but this is just the initial point

   #Let's go back a dx or two so we need only to look at forward intercepts
   init_interval = 4
   n_line_intercept = n_line_intercept - n_line*init_interval/2


# This is horrible, I sort of wish we had spherical shells
   # mirror across the origin on all axes if on negative hemisphere -
   # easier to check for intersection only with xyz = 0.5 planes
   # nsz = (3, n_line_intercept.shape[0],n_line_intercept.shape[1])
   # positive_intercepts_0 = np.ones(nsz)
   # positive_intercepts_1 = np.ones(nsz)
   # positive_intercepts_0.fill(np.nan)
   # positive_intercepts_1.fill(np.nan)
   lineintercept = n_line_intercept
   linevec = n_line/np.broadcast_to(np.linalg.norm(n_line,axis=-1),(3,n_cells)).transpose()
   
   # # find distances to the planes delimiting the cell along the line
   
   # d_fw_planes_0 = np.divide(-0.5-lineintercept,linevec)
   # d_fw_planes_1 = np.divide(0.5-lineintercept,linevec)

   #hits = np.all(np.abs(n_line_intercept) < 0.5, axis=-1) # these are certain, also init hits array
   #eps = 1e-12

   # fw intercepts
   if False:
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

   # numsteps = 41 # Sample this many points with SDF - brute force for now; 41 points is 0.1dx for +-2dx
   # for d in np.linspace(0,4,num=numsteps):
   #    sdfs = fBox(lineintercept + linevec*d)
   #    np.fmin(s_line, sdfs, out=s_line, where = mask)

   # Golden section search. The SDF contains absolute values and it's not too bad when vectorized.
   in_stepping = mask.copy()
   bracket_dir = np.ones_like(s_line)

   gratio = 0.38197
   lower = gratio*init_interval
   a = np.zeros_like(s_line)
   b = np.ones_like(s_line)*lower
   c = np.ones_like(s_line)*init_interval
   x = np.ones_like(s_line)*b+(c-b)*gratio
   #print(a,b,c,x)
   in_stepping = (
                  in_stepping & 
                  (fBox(lineintercept + a[:,np.newaxis]*linevec) > fBox(lineintercept + b[:,np.newaxis]*linevec)) & 
                  (fBox(lineintercept + c[:,np.newaxis]*linevec) > fBox(lineintercept + b[:,np.newaxis]*linevec))
                  )

   tolerance = 0.01 # interval to reach in dx
   niter = 0
   maxiters = 20
   db = fBox(lineintercept + b[:,np.newaxis]*linevec)
   dx = fBox(lineintercept + x[:,np.newaxis]*linevec)
   while(np.any(in_stepping[mask])): 

      db[in_stepping] = fBox(lineintercept[in_stepping,:] + b[in_stepping,np.newaxis]*linevec[in_stepping,:])
      dx[in_stepping] = fBox(lineintercept[in_stepping,:] + x[in_stepping,np.newaxis]*linevec[in_stepping,:])
      a[(bracket_dir == 1) & (db < dx) & in_stepping] = a[(bracket_dir == 1) & (db < dx) & in_stepping]
      b[(bracket_dir == 1) & (db < dx) & in_stepping] = b[(bracket_dir == 1) & (db < dx) & in_stepping]
      c[(bracket_dir == 1) & (db < dx) & in_stepping] = x[(bracket_dir == 1) & (db < dx) & in_stepping]

      a[(bracket_dir == 1) & (db >= dx) & in_stepping] = b[(bracket_dir == 1) & (db >= dx) & in_stepping]
      b[(bracket_dir == 1) & (db >= dx) & in_stepping] = x[(bracket_dir == 1) & (db >= dx) & in_stepping]
      c[(bracket_dir == 1) & (db >= dx) & in_stepping] = c[(bracket_dir == 1) & (db >= dx) & in_stepping]

      a[(bracket_dir == -1) & (db < dx) & in_stepping] = x[(bracket_dir == -1) & (db < dx) & in_stepping]
      b[(bracket_dir == -1) & (db < dx) & in_stepping] = b[(bracket_dir == -1) & (db < dx) & in_stepping]
      c[(bracket_dir == -1) & (db < dx) & in_stepping] = c[(bracket_dir == -1) & (db < dx) & in_stepping]

      a[(bracket_dir == -1) & (db >= dx) & in_stepping] = b[(bracket_dir == -1) & (db >= dx) & in_stepping]
      b[(bracket_dir == -1) & (db >= dx) & in_stepping] = x[(bracket_dir == -1) & (db >= dx) & in_stepping]
      c[(bracket_dir == -1) & (db >= dx) & in_stepping] = c[(bracket_dir == -1) & (db >= dx) & in_stepping]

      bracket_dir[(b - a) < (c - b)] = 1
      bracket_dir[(b - a) >= (c - b)] = -1
      x[in_stepping & (bracket_dir == -1)] = b[in_stepping & (bracket_dir == -1)] - ((b-a)*gratio)[in_stepping & (bracket_dir == -1)]
      x[in_stepping & (bracket_dir ==  1)] = b[in_stepping & (bracket_dir ==  1)] + ((c-b)*gratio)[in_stepping & (bracket_dir ==  1)]

      in_stepping = in_stepping & ((c - a) > tolerance)
      niter = niter+1
      print(niter,"iterations, remaining: ", np.sum(in_stepping), 
            "right ", np.sum(bracket_dir[in_stepping] ==1), "left ", np.sum(bracket_dir[in_stepping]==-1))
      if(niter >= maxiters):
         break
      
   s_line[mask] = np.fmin(db,dx)[mask]




   # np.divide(-1, s_line,out=s_line, where=hits) # SDF takes care of this now!
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
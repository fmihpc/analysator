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

   divS = (np.roll(S[:,:,:,0],-1, 0) - np.roll(S[:,:,:,0], 1, 0) +
         np.roll(S[:,:,:,1],-1, 1) - np.roll(S[:,:,:,1], 1, 1) +
         np.roll(S[:,:,:,2],-1, 2) - np.roll(S[:,:,:,2], 1, 2))/(2*dx)
   print('divS.shape', divS.shape)

   return divS

def fg_vol_jacobian(b, dx=1e6):
   dFx_dx, dFx_dy, dFx_dz = np.gradient(b[:,:,:,0], dx)
   dFy_dx, dFy_dy, dFy_dz = np.gradient(b[:,:,:,1], dx)
   dFz_dx, dFz_dy, dFz_dz = np.gradient(b[:,:,:,2], dx)

   return np.stack()

def fg_vol_curl(array, dx=1e6):
   dummy,  dFx_dy, dFx_dz = np.gradient(array[:,:,:,0], dx)
   dFy_dx, dummy,  dFy_dz = np.gradient(array[:,:,:,1], dx)
   dFz_dx, dFz_dy, dummy  = np.gradient(array[:,:,:,2], dx)

   rotx = dFz_dy - dFy_dz
   roty = dFx_dz - dFz_dx
   rotz = dFy_dx - dFx_dy

   return np.stack([rotx, roty, rotz], axis=-1)

def fg_vol_div(array, dx=1e6):
   dFx_dx, dummy, dummy = np.gradient(array[:,:,:,0], dx)
   dummy, dFy_dy, dummy = np.gradient(array[:,:,:,1], dx)
   dummy, dummy, dFz_dz = np.gradient(array[:,:,:,2], dx)

   return dFx_dx+dFy_dy+dFz_dz

def vfield3_dot(a, b):
    """Calculates dot product of vectors a and b in 3D vector field"""

    return np.sum(np.multiply(a,b), axis=-1)
    #return (
    #    a[:, :, :, 0] * b[:, :, :, 0]
    #    + a[:, :, :, 1] * b[:, :, :, 1]
    #    + a[:, :, :, 2] * b[:, :, :, 2]
    #)

def vfield3_matder(a, b, dr):
    """Calculates material derivative of 3D vector fields a and b"""

    bx = b[:, :, :, 0]
    by = b[:, :, :, 1]
    bz = b[:, :, :, 2]

    grad_bx = np.gradient(bx, dr)
    grad_by = np.gradient(by, dr)
    grad_bz = np.gradient(bz, dr)

    resx = vfield3_dot(a, grad_bx)
    resy = vfield3_dot(a, grad_by)
    resz = vfield3_dot(a, grad_bz)

    return np.stack((resx, resy, resz), axis=-1)

def vfield3_normalise(a):

    amag = np.linalg.norm(a, axis=-1)

    res=a/amag
    return res
    #resx = a[:, :, :, 0] / amag
    #resy = a[:, :, :, 1] / amag
    #resz = a[:, :, :, 2] / amag

    #return np.stack((resx, resy, resz), axis=-1)

def vfield3_curvature(a, dr=1000e3):
   a = vfield3_normalise(a)
   return vfield3_matder(a, a, dr)

def ballooning_crit(B, P, beta, dr=1000e3):

    n = vfield3_curvature(B, dr)

    nnorm = vfield3_normalise(n)

    gradP = np.gradient(P,dr)

    kappaP = vfield3_dot(nnorm, gradP) / P
    
    kappaC = vfield3_dot(nnorm, n)

    balloon = (2 + beta) / 4.0 * kappaP / (kappaC + 1e-27)

    return (balloon, kappaC, gradP)
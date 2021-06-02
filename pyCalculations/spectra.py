# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2021 University of Helsinki
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

# Function to reduce the velocity space in a spatial cell to an omnidirectional energy spectrum
def get_spectrum_energy(vlsvReader,
                  cid,
                  population="proton",
                  fMin = 1e-21,
                  EMin=100,
                  EMax=80e3,
                  nBins=66,
                  mass=1.6726219e-27, # default: mp
                  q=1.60217662e-19,   # default: qe
                  frame=None,
                  restart=True):
   import numpy as np
   import pytools as pt
   EkinBinEdges = np.logspace(np.log10(EMin),np.log10(EMax),nBins)
   vlsvReader = pt.vlsvfile.VlsvReader(vlsvReader)
   # check if velocity space exists in this cell
   if not restart and vlsvReader.read_variable('fSaved',cid) != 1.0:
      return (False,np.zeros(nBins), EkinBinEdges)
   if vlsvReader.check_variable('MinValue') == True:
      fMin = vlsvReader.read_variable('MinValue',cid)
   #print('Cell ' + str(cid).zfill(9))
   velcells = vlsvReader.read_velocity_cells(cid, population)
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()))
   V2 = np.sum(np.square(V),1)
   Ekin = 0.5*mass*V2/q
   f = list(zip(*velcells.items()))


   # check that velocity space has cells - still return a zero histogram in the same shape
   if(len(f) > 0):
      f = np.asarray(f[1])
   else:
      return (False,np.zeros(nBins), EkinBinEdges)
   ii_f = np.where(f >= fMin)
   if len(ii_f) < 1:
      return (False,np.zeros(nBins), EkinBinEdges)
   f = f[ii_f]
   Ekin = Ekin[ii_f]
   V2 = V2[ii_f]

   Ekin[Ekin < min(EkinBinEdges)] = min(EkinBinEdges)
   Ekin[Ekin > max(EkinBinEdges)] = max(EkinBinEdges)
   # normalization
   fv = f*np.sqrt(V2) # use particle flux as weighting
   # compute histogram
   (nhist,edges) = np.histogram(Ekin,bins=EkinBinEdges,weights=fv,normed=0)
   # normalization
   dE = EkinBinEdges[1:] - EkinBinEdges[0:-1]
   nhist = np.divide(nhist,(dE*4*np.pi))
   return (True,nhist,edges)


# Function to reduce the velocity space in a spatial cell to an omnidirectional energy spectrum
def get_spectrum_modvelocity(vlsvReader,
                  cid,
                  population="proton",
                  fMin = 1e-21,
                  VMin=100,
                  VMax=2e6,
                  nBins=66,
                  frame=None,
                  restart=True):
   import numpy as np
   import pytools as pt
   VBinEdges = np.logspace(np.log10(VMin),np.log10(VMax),nBins)
   vlsvReader = pt.vlsvfile.VlsvReader(vlsvReader)
   # check if velocity space exists in this cell
   if not restart and vlsvReader.read_variable('fSaved',cid) != 1.0:
      return (False,np.zeros(nBins), VBinEdges)
   if vlsvReader.check_variable('MinValue') == True:
      fMin = vlsvReader.read_variable('MinValue',cid)
   #print('Cell ' + str(cid).zfill(9))
   velcells = vlsvReader.read_velocity_cells(cid, population)
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()))
   V2 = np.sum(np.square(V),1)
   Vmod = np.sqrt(V2)
   f = list(zip(*velcells.items()))


   # check that velocity space has cells - still return a zero histogram in the same shape
   if(len(f) > 0):
      f = np.asarray(f[1])
   else:
      return (False,np.zeros(nBins), VBinEdges)
   ii_f = np.where(f >= fMin)
   if len(ii_f) < 1:
      return (False,np.zeros(nBins), VBinEdges)
   f = f[ii_f]
   Vmod = Vmod[ii_f]
   V2 = V2[ii_f]

   Vmod[Vmod < min(VBinEdges)] = min(VBinEdges)
   Vmod[Vmod > max(VBinEdges)] = max(VBinEdges)
   # normalization
   fv = f*np.sqrt(V2) # use particle flux as weighting
   # compute histogram
   (nhist,edges) = np.histogram(Vmod,bins=VBinEdges,weights=fv,normed=0)
   # normalization
   dV = VBinEdges[1:] - VBinEdges[0:-1]
   nhist = np.divide(nhist,(dV*4*np.pi))
   return (True,nhist,edges)

# Function to reduce the velocity space in a spatial cell to an omnidirectional energy spectrum
def get_spectrum_alongaxis_vel(vlsvReader,
                  cid,
                  population="proton",
                  vector=None,
                  vectorVar="vg_b_vol",
                  fMin=1e-21,
                  VMin=-2e6,
                  VMax=2e6,
                  nBins=200,
                  frame=None,
                  restart=True):
   import numpy as np
   import pytools as pt
   vlsvReader = pt.vlsvfile.VlsvReader(vlsvReader)

   if vectorVar is not None and vector is None:
      vector=vlsvReader.read_variable(vectorVar, cid)
   vector = vector/np.linalg.norm(vector)
   VBinEdges = np.linspace(VMin,VMax,nBins)
   # check if velocity space exists in this cell
   if not restart and vlsvReader.read_variable('fSaved',cid) != 1.0:
      return (False,np.zeros(nBins), VBinEdges)
   if vlsvReader.check_variable('MinValue') == True:
      fMin = vlsvReader.read_variable('MinValue',cid)
   #print('Cell ' + str(cid).zfill(9))
   velcells = vlsvReader.read_velocity_cells(cid, population)
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()))
   V2 = np.sum(np.square(V),1)
   Vproj = np.dot(V,vector)
   f = list(zip(*velcells.items()))


   # check that velocity space has cells - still return a zero histogram in the same shape
   if(len(f) > 0):
      f = np.asarray(f[1])
   else:
      return (False,np.zeros(nBins), VBinEdges)
   ii_f = np.where(f >= fMin)
   if len(ii_f) < 1:
      return (False,np.zeros(nBins), VBinEdges)
   f = f[ii_f]
   Vproj = Vproj[ii_f]
   V2 = V2[ii_f]

   Vproj[Vproj < min(VBinEdges)] = min(VBinEdges)
   Vproj[Vproj > max(VBinEdges)] = max(VBinEdges)
   # normalization
   fv = f*np.sqrt(V2) # use particle flux as weighting
   # compute histogram
   (nhist,edges) = np.histogram(Vproj,bins=VBinEdges,weights=fv,normed=0)
   # normalization
   dV = VBinEdges[1:] - VBinEdges[0:-1]
   nhist = np.divide(nhist,(dV*4*np.pi))
   return (True,nhist,edges)


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
# Weighted by particle flux/none
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
                  weight='flux',
                  restart=True):
   import numpy as np
   import pytools
   EkinBinEdges = np.logspace(np.log10(EMin),np.log10(EMax),nBins)
   vlsvReader = pytools.vlsvfile.VlsvReader(vlsvReader)
   # check if velocity space exists in this cell
   if not restart and vlsvReader.read_variable('fSaved',cid) != 1.0:
      return (False,np.zeros(nBins), EkinBinEdges)
   if vlsvReader.check_variable('MinValue') == True:
      fMin = vlsvReader.read_variable('MinValue',cid)
   #print('Cell ' + str(cid).zfill(9))
   velcells = vlsvReader.read_velocity_cells(cid, population)
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()), pop=population)
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
   if(weight == 'flux'):
      fw = f*np.sqrt(V2) # use particle flux as weighting
      units = "s/(m^2 4pi eV)"
      latexunits = '$\mathrm{s}\,\mathrm{m}^{-2}\,(4\pi\,\mathrm{eV})^{-1}$'
      latex='$f(\vec{r})v\,\DeltaE\,sr^-1$'
   else:
      fw = f
      units = "1/(m^3 4pi eV)"
      latexunits = '$\mathrm{m}^{-3}\,(4\pi\,\mathrm{eV})^{-1}$'
      latex='$f(\vec{r})\,\DeltaE^-1\,sr^-1$'
      weight = 'particles'
   dV = np.prod(vlsvReader.get_velocity_mesh_dv(population))
   # compute histogram
   (nhist,edges) = np.histogram(Ekin,bins=EkinBinEdges,weights=fw*dV,normed=0)
   # normalization
   ftotal = np.sum(f)*dV
   #print('ftotal', ftotal, 'sum hist', np.sum(nhist))
   dE = EkinBinEdges[1:] - EkinBinEdges[0:-1]
   nhist = np.divide(nhist,(dE*4*np.pi))
   vari = pytools.calculations.VariableInfo(nhist,
                                       name="Omnidirectional energy spectrum "+population+' ('+weight+')',
                                       units=units,
                                       latex=latex,
                                       latexunits=latexunits)
   return (True,vari,edges)


# Function to reduce the velocity space in a spatial cell to an omnidirectional speed spectrum
def get_spectrum_modvelocity(vlsvReader,
                  cid,
                  population="proton",
                  fMin = 1e-21,
                  VMin=100,
                  VMax=2e6,
                  nBins=66,
                  frame=None,
                  weight='flux',
                  restart=True):
   import numpy as np
   import pytools as pt
   VBinEdges = np.logspace(np.log10(VMin),np.log10(VMax),nBins)
   vlsvReader = pytools.vlsvfile.VlsvReader(vlsvReader)
   # check if velocity space exists in this cell
   if not restart and vlsvReader.read_variable('fSaved',cid) != 1.0:
      return (False,np.zeros(nBins), VBinEdges)
   if vlsvReader.check_variable('MinValue') == True:
      fMin = vlsvReader.read_variable('MinValue',cid)
   #print('Cell ' + str(cid).zfill(9))
   velcells = vlsvReader.read_velocity_cells(cid, population)
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()), pop=population)
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
   if(weight == 'flux'):
      fw = f*np.sqrt(V2) # use particle flux as weighting
   else:
      fw = f
   # compute histogram
   (nhist,edges) = np.histogram(Vmod,bins=VBinEdges,weights=fw,normed=0)
   # normalization
   dv = VBinEdges[1:] - VBinEdges[0:-1]
   nhist = np.divide(nhist,(dv*4*np.pi))
   return (True,nhist,edges)

# Function to reduce the velocity space in a spatial cell to an spectrum in velocity along a vector
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
                  weight='flux',
                  restart=True):
   import numpy as np
   import pytools as pt
   vlsvReader = pytools.vlsvfile.VlsvReader(vlsvReader)

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
   V = vlsvReader.get_velocity_cell_coordinates(list(velcells.keys()), pop=population)
   V2 = np.sum(np.square(V),1)
   Vproj = np.dot(V,vector)
   f = list(velcells.values())
   

   # check that velocity space has cells - still return a zero histogram in the same shape
   if(len(f) > 0):
      f = np.asarray(f)
   else:
      return (False,np.zeros(nBins), VBinEdges)
   ii_f = np.where(f >= fMin)
   #ii_f = np.where(f >= fMin or True)
   if len(ii_f) < 1:
      return (False,np.zeros(nBins), VBinEdges)
   #f = f[ii_f]
   #Vproj = Vproj[ii_f]
   #V2 = V2[ii_f]

   Vproj[Vproj < min(VBinEdges)] = min(VBinEdges)
   Vproj[Vproj > max(VBinEdges)] = max(VBinEdges)
   # normalization
   # normalization
   if(weight == 'flux'):
      fw = f*np.sqrt(V2) # use particle flux as weighting
      units = "s^2/(m^3)"
      latexunits = '$\mathrm{s}\,\mathrm{m}^{-4}$'
      latex='$f(\vec{r})\,\Deltav^-1$'
   else:
      units = "s/(m^4 4pi)"
      latexunits = '$\mathrm{s}\mathrm{m}^{-4}$'
      latex='$f(\vec{r})\,\Deltav^-1\$'
      weight = 'particles'
      fw = f
   
   # compute histogram
   #print(fw)
   #print(len(Vproj), len(fw))
   dV = np.prod(vlsvReader.get_velocity_mesh_dv(population))
   (nhist,edges) = np.histogram(Vproj,bins=VBinEdges,weights=fw*dV,normed=0)
   # normalization
   dv = VBinEdges[1:] - VBinEdges[0:-1]
   nhist = np.divide(nhist,dv)
   
   vari = pytools.calculations.VariableInfo(nhist,
                                    name="Parallel spectrum of "+population+' ('+weight+')',
                                    units=units,
                                    latex=latex,
                                    latexunits=latexunits)

   return (True,vari,edges)


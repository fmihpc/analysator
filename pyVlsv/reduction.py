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

''' A file for doing data reduction on variables
'''
import numpy as np
import pylab as pl
from reducer import DataReducerVariable
from rotation import rotateTensorToVector, rotateArrayTensorToVector
from gyrophaseangle import gyrophase_angles
import vlsvvariables
import sys

mp = 1.672622e-27
elementalcharge = 1.6021773e-19

def pass_op( variable ):
   # do nothing
   return variable

def x_component( variable ):
   if np.ndim(variable) == 1:
      return variable[0]
   return variable[:,0]

def y_component( variable ):
   if np.ndim(variable) == 1:
      return variable[1]
   return variable[:,1]

def z_component( variable ):
   if np.ndim(variable) == 1:
      return variable[2]
   return variable[:,2]

def magnitude( variable ):
   return np.linalg.norm(np.asarray(variable),axis=-1)

def absolute( variable ):
   return abs(variable)

def sumv( variable ):
   # Note: this is used to sum over multipops, thus the summing axis is zero
   if np.ndim(variable) > 3:
      print('Error: Number of dimensions is too large')
      return
   else:
      # First dimension: populations
      # Second dimension: cells
      # Third dimension: components
      return np.sum(np.array(variable),axis=0)

def condition_matrix_array( condition, matrices ):
   # This routine is still very slow due to for-loops
   ''' Repeats condition n times and forms an array of it
       :param: condition    Some matrix of conditions
       :param: matrices     An array of matrices on which to apply the condition
       :returns: An array of arrays which only have the conditioned elements from each matrix

       .. code-block:: python
          Example: Pick the diagonal components
          diagonal_condition = np.array([
                      [True,  False, False],
                      [False, True,  False],
                      [False, False, True]
                              ])
          matrices = np.array([
                      [[3,  1, 2],
                       [0,  3, 0],
                       [5,  1, 3]],
                      [[5,  1, 2],
                       [0,  5, 2],
                       [1,  0, 5]]
                     ])
          print condition_matrix_array( diagonal_condition, matrices )
          # Output:
          # array([[3,3,3], [5,5,5]])

   '''
   # Make sure the matrices is in the correct shape
   if np.ndim(matrices) == np.ndim(condition):
      matrices = np.reshape(matrices, tuple(np.concatenate(([1],[x for x in matrices.shape]))) )
   # Repeat the condition array n times
   condition_array = condition*np.reshape(np.ones(np.product(matrices.shape)), matrices.shape)
   # Get the number of true hits in the condition:
   num_of_true = 0
   for i in np.ravel(condition):
      if i == True:
         num_of_true = num_of_true + 1
   # Extract the elements that fill the condition:
   extracted = np.extract( condition_array, matrices )
   # Reshape the matrices:
   extracted = np.reshape(extracted, (len(matrices), num_of_true))
   if len(extracted) == 1:
      extracted = extracted[0]
   # Return the extracted elements
   return extracted



def restart_B( variables ):
   ''' Data reducer function for summing perturbed_B and background_B
   '''
   return variables[0]+variables[1]

def restart_V( variables ):
   ''' Data reducer function for finding bulk V for a restart fike
   '''
   moments = variables[0]

   if np.ndim(moments)==1: # single cell
      if len(moments)==4:  # pre-multipop restart
         V = moments[1:4]/moments[0]
      elif len(moments)==5: # multipop restart
         V = moments[1:4]
      else:
         print("Unrecognized length for moments!")
         return None
   else: # array of cells
      if len(moments[0,:])==4:  # pre-multipop restart
         rho = np.ma.masked_less_equal(moments[:,0],0)
         V = np.ma.divide(moments[:,1:4],rho[:,np.newaxis])
      elif len(moments[0,:])==5: # multipop restart
         V = moments[:,1:4]
      else:
         print("Unrecognized length for moments!")
         return None
   return V

def restart_rho( variables ):
   ''' Data reducer function for calculating proton number density rho from restart file
   '''
   moments = variables[0]

   if np.ndim(moments)==1: # single cell
      if len(moments)==4:  # pre-multipop restart
         rho = moments[0]
      else:
         print("Unable to determine rho from moments!")
         return None
   else: # array of cells
      if len(moments[0,:])==4:  # pre-multipop restart
         rho = moments[:,0]
      else:
         print("Unable to determine rho from moments!")
         return None
   return rho

def restart_rhom( variables ):
   ''' Data reducer function for calculating rhom from restart file
   '''
   moments = variables[0]
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp

   if np.ndim(moments)==1: # single cell
      if len(moments)==4:  # pre-multipop restart
         rhom = moments[0]*mass
      elif len(moments)==5: # multipop restart
         rhom = moments[0]
      else:
         print("Unrecognized length for moments!")
         return None
   else: # array of cells
      if len(moments[0,:])==4:  # pre-multipop restart
         rhom = moments[:,0]*mass
      elif len(moments[0,:])==5: # multipop restart
         rhom = moments[:,0]
      else:
         print("Unrecognized length for moments!")
         return None
   return rhom

def restart_rhoq( variables ):
   ''' Data reducer function for calculating rhoq from restart file
   '''
   moments = variables[0]
   charge = vlsvvariables.speciescharge[vlsvvariables.activepopulation]*elementalcharge

   if np.ndim(moments)==1: # single cell
      if len(moments)==4:  # pre-multipop restart
         rhoq = moments[0]*charge
      elif len(moments)==5: # multipop restart
         rhoq = moments[4]
      else:
         print("Unrecognized length for moments!")
         return None
   else: # array of cells
      if len(moments[0,:])==4:  # pre-multipop restart
         rhoq = moments[:,0]*charge
      elif len(moments[0,:])==5: # multipop restart
         rhoq = moments[:,4]
      else:
         print("Unrecognized length for moments!")
         return None
   return rhoq


def rhom( variables ):
   ''' Data reducer function for calculating rhom from pre-multipop file
   '''
   rho = variables[0]
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   return rho*mass

def rhoq( variables ):
   ''' Data reducer function for calculating rhoq from pre-multipop file
   '''
   rho = variables[0]
   charge = vlsvvariables.speciescharge[vlsvvariables.activepopulation]*elementalcharge
   return rho*charge

def v( variables ):
   ''' Data reducer function for getting velocity from rho and rho_v
       variables[0] = rho_v and variables[1] = rho
   '''
   epsilon = sys.float_info.epsilon
   rho_v = np.array(variables[0])
   rho = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[1])),0)
   if np.ndim(rho) == 0:
      return np.ma.divide(rho_v,rho)
   else:
      return np.ma.divide(rho_v,rho[:,np.newaxis])

def vms( variables ):
   ''' Data reducer function for getting magnetosonic velocity 
       input: P, rhom, B
   '''
   epsilon = sys.float_info.epsilon
   mu_0 = 1.25663706144e-6
   P = variables[0]
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[1])),0)
   B = variables[2]
   if np.ndim(B) == 1:
      Btot = np.sqrt( np.sum( np.square(B) ) )
   else:
      Btot = np.sqrt( np.sum(np.square(B),axis=1) )
   vs = np.sqrt( np.ma.divide( P*5.0/3.0, rho_m ) )
   vA = np.ma.divide( Btot,np.sqrt( mu_0*rho_m ) )
   vms = np.sqrt( np.square(vs) + np.square(vA) )
   return vms

def vs( variables ):
   ''' Data reducer function for getting the sound speed
       input: P, rhom
   '''
   epsilon = sys.float_info.epsilon
   P = variables[0]
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[1])),0)
   vs = np.sqrt( np.ma.divide( P*5.0/3.0, rho_m ) )
   return vs

def va( variables ):
   ''' Data reducer function for getting the Alfven velocity
       input: rhom, B
   '''
   epsilon = sys.float_info.epsilon
   mu_0 = 1.25663706144e-6
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[0])),0)
   B = variables[1]
   if np.ndim(B) == 1:
      Btot = np.sqrt( np.sum( np.square(B) ) )
   else:
      Btot = np.sqrt( np.sum(np.square(B),axis=1) )
   vA = np.ma.divide( Btot,np.sqrt( rho_m*mu_0 ) )
   return vA

def MA( variables ):
   ''' Data reducer function for getting the Alfvenic Mach number
   '''
   bulkv = np.linalg.norm(variables[0],axis=-1)
   Alfvenspeed = np.ma.masked_less_equal(variables[1],0)
   MA = np.ma.divide(bulkv, Alfvenspeed)
   return MA

def Mms( variables ):
   ''' Data reducer function for getting the magnetosonic Mach number
   '''
   bulkv = np.linalg.norm(variables[0],axis=-1)
   magnetosonicspeed = np.ma.masked_less_equal(variables[1],0)
   Mms = np.ma.divide(bulkv, magnetosonicspeed)
   return Mms

def ParallelVectorComponent( variables ):
   ''' Data reducer function for vector component parallel to the magnetic field (or another vector)
   '''
   inputvector = variables[0]
   bgvector = variables[1]
   if( np.ndim(inputvector)==1 ):
      bgnorm = bgvector/np.linalg.norm(bgvector)
      return (inputvector*bgnorm).sum()
   else:
      bgnorm = np.ma.divide(bgvector, np.linalg.norm(bgvector, axis=-1)[:,np.newaxis])
      return (inputvector*bgnorm).sum(-1)

def PerpendicularVectorComponent( variables ):
   ''' Data reducer function for vector component perpendicular to the magnetic field (or another vector)
   '''
   inputvector = variables[0]
   bgvector = variables[1]
   if( np.ndim(inputvector)==1 ):
      bgnorm = bgvector / np.linalg.norm(bgvector)
      vpara = (inputvector*bgnorm).sum()
      vmag = np.linalg.norm(inputvector)
      return np.sqrt(vmag*vmag - vpara*vpara)
   else:
      bgnorm = np.ma.divide(bgvector, np.linalg.norm(bgvector, axis=-1)[:,np.newaxis])
      vpara = (inputvector*bgnorm).sum(-1)
      vmag = np.linalg.norm(inputvector, axis=-1)
      return np.sqrt(vmag*vmag - vpara*vpara)
   
def FullTensor( variables ):
   ''' Data reducer function to reconstruct a full tensor from diagonal and off-diagonal
       components (e.g. the pressure tensor)
   '''
   TensorDiagonal = np.array(variables[0])
   TensorOffDiagonal = np.array(variables[1])
   if(np.ndim(TensorDiagonal)==1 ):      
      TensorDiagonal = TensorDiagonal.reshape(1,3)[0]
      TensorOffDiagonal = TensorOffDiagonal.reshape(1,3)[0]
      return np.array([[TensorDiagonal[0], TensorOffDiagonal[2], TensorOffDiagonal[1]],
                       [TensorOffDiagonal[2], TensorDiagonal[1], TensorOffDiagonal[0]],
                       [TensorOffDiagonal[1], TensorOffDiagonal[0], TensorDiagonal[2]]])
   else:
      result = np.empty([len(TensorDiagonal[:,0]),3,3]) # Warning, unitialized!
      result[:,0,0] = TensorDiagonal[:,0]
      result[:,0,1] = TensorOffDiagonal[:,2]
      result[:,0,2] = TensorOffDiagonal[:,1]
      result[:,1,0] = TensorOffDiagonal[:,2]
      result[:,1,1] = TensorDiagonal[:,1]
      result[:,1,2] = TensorOffDiagonal[:,0]
      result[:,2,0] = TensorOffDiagonal[:,1]
      result[:,2,1] = TensorOffDiagonal[:,0]
      result[:,2,2] = TensorDiagonal[:,2]
      return result

def RotatedTensor( variables ):
   ''' Data reducer for rotating e.g. the pressure tensor to align the z-component 
       with a vector, e.g. the magnetic field
   '''
   Tensor = variables[0]
   B = variables[1]
   if( np.ndim(B)==1 ):
      B = B.reshape(1,3)[0]
      Tensor = Tensor.reshape(3,3)
      return rotateTensorToVector(Tensor, B)
   else:
      return rotateArrayTensorToVector(Tensor, B)

def ParallelTensorComponent( variables ):
   ''' Data reducer for finding the parallel component from a rotated field-aligned tensor
   '''
   RotatedTensor = variables[0]
   if( np.ndim(RotatedTensor)==2 ):
      return RotatedTensor[2,2]
   else:
      return RotatedTensor[:,2,2]

def PerpendicularTensorComponent( variables ):
   ''' Data reducer for finding the perpendicular component from a rotated field-aligned tensor
       e.g. perpendicular pressure
   '''
   RotatedTensor = variables[0]
   if( np.ndim(RotatedTensor)==2 ):
      return 0.5*(RotatedTensor[0,0] + RotatedTensor[1,1])
   else:
      return 0.5*(RotatedTensor[:,0,0] + RotatedTensor[:,1,1])

def Anisotropy( variables ):
   ''' Data reducer for finding the ratio of perpendicular to parallel components of a tensor
   '''
   RotatedTensor = variables[0]
   if( np.ndim(RotatedTensor)==2 ):
      divisor = np.ma.masked_equal(np.ma.masked_invalid(RotatedTensor[2,2]),0)
      return 0.5*np.ma.divide(RotatedTensor[0,0] + RotatedTensor[1,1], divisor)
   else:
      divisor = np.ma.masked_equal(np.ma.masked_invalid(RotatedTensor[:,2,2]),0)
      return 0.5*np.ma.divide(RotatedTensor[:,0,0] + RotatedTensor[:,1,1], divisor)

def Pressure( variables ):
   ''' Data reducer for finding the scalar pressure
   '''
   PTensorDiagonal = variables[0]
   return 1.0/3.0 * np.ma.sum(np.ma.masked_invalid(PTensorDiagonal),axis=-1)

def Pdyn( variables ):
   ''' Data reducer function for dynamic pressure
       input: V, rhom
   '''
   Vmag = np.linalg.norm(np.array(variables[0]), axis=-1)
   rhom = np.array(variables[1])
   return Vmag*Vmag*rhom

def Pdynx( variables ):
   ''' Data reducer function for dynamic pressure with just V_x
       input: V, rhom
   '''
   rhom = np.array(variables[1])
   V = np.array(variables[0])
   if( np.ndim(V)==2 ):
      Vx = V[:,0]
   else:
      Vx = V[0]
   return Vx*Vx*rhom

def Poynting( variables ):
   ''' Data reducer for the Poynting vector
   '''
   E=np.array(variables[0])
   B=np.array(variables[1])
   return np.cross(E, B) / 1.25663706144e-6

def Temperature( variables ):
   ''' Data reducer for converting pressure to temperature
   '''
   kb = 1.38065e-23
   Pressure = variables[0] # either a tensor, vector, array, or value
   rho = variables[1] # eithern array or a value
   divisor = np.ma.masked_less_equal( np.ma.masked_invalid(np.array(rho)),0) * kb
   # assumes first dimension is either single cell or a single-dimension array
   if np.ndim(divisor)==0: # single cell
      if np.ndim(Pressure)==0:
         return np.ma.divide(Pressure, divisor)
      elif np.ndim(Pressure)==1:
         return np.ma.divide(Pressure, divisor[np.newaxis])
      elif np.ndim(Pressure)==2:
         return np.ma.divide(Pressure, divisor[np.newaxis,np.newaxis])
   else: # array of cellids
      if np.ndim(Pressure)==1:
         return np.ma.divide(Pressure, divisor)
      elif np.ndim(Pressure)==2:
         return np.ma.divide(Pressure, divisor[:,np.newaxis])
      elif np.ndim(Pressure)==3:
         return np.ma.divide(Pressure, divisor[:,np.newaxis,np.newaxis])
   # Should not reach here...
   print("Error finding dimensions in calculating temperature!")
   return -1

def aGyrotropy( variables ):
   ''' Data reducer function to evaluate agyrotropy
   from equation (6) of M. Swisdak, Quantifying gyrotropy in magnetic reconnection
   https://doi.org/10.1002/2015GL066980
   (read also the appendix)
   '''
   PTensorRot = np.array(variables[0])
   # np.array([[PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]],
   #           [PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]],
   #           [PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]]])
   if(np.ndim(PTensorRot)==2):
      PPerp = 0.5*(PTensorRot[0,0]+PTensorRot[1,1])
      numerator = PTensorRot[1,0]*PTensorRot[1,0] + PTensorRot[2,0]*PTensorRot[2,0] + PTensorRot[2,1]*PTensorRot[2,1]
      denominator = PPerp*PPerp + 2.0*PPerp*PTensorRot[2,2]
      if denominator > 1.e-20:
         aGyro = numerator/denominator
      else:
         aGyro = 0.0
      return np.array(aGyro)
   else:
      PPerp = 0.5*(PTensorRot[:,0,0]+PTensorRot[:,1,1])
      numerator = np.ma.masked_less_equal( PTensorRot[:,1,0]*PTensorRot[:,1,0] + PTensorRot[:,2,0]*PTensorRot[:,2,0] + PTensorRot[:,2,1]*PTensorRot[:,2,1], 0)
      denominator = np.ma.masked_less_equal( PPerp*PPerp + 2.0*PPerp*PTensorRot[:,2,2], 0)
      aGyro = np.ma.divide(numerator, denominator)
      return np.array(aGyro)

def beta( variables ):
   ''' Data reducer for finding the plasma beta
   '''
   Pressure = variables[0]
   Magneticfield = variables[1]   
   return 2.0 * 1.25663706144e-6 * np.ma.divide(Pressure, np.sum(np.asarray(Magneticfield)**2,axis=-1))

def rMirror( variables ):
   # More efficient file access, now just takes PTensor and B
   PT = variables[0]
   B = variables[1]
   PTRotated =  RotatedTensor([PT,B])
   TAniso = Anisotropy([PTRotated]) # PAniso == TAniso
   PPerp = PerpendicularTensorComponent([PTRotated])
   betaPerp = beta([PPerp,B])
   return betaPerp * (TAniso - 1)   

def thermalvelocity( variables ):
   Temperature = variables[0]
   kb = 1.38065e-23
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   # Corrected to calculate the mean speed sqrt(8kT/pi m)
   thermalvelocity = np.sqrt(Temperature*(kb*8./(mass*3.14159)))
   return thermalvelocity

def Vstream( variables ):
   rhoVstream = variables[0]
   rhostream = variables[1]
   rhoVNonBackstream = variables[2]
   rhoNonBackstream = variables[3]
   # get velocity of both populations:
   vBackstream = v( [rhoVBackstream, rhoBackstream] )
   vNonBackstream = v( [rhoVNonBackstream, rhoNonBackstream] )
   vBeam = vBackstream - vNonBackstream
   return vBeam # <- is a vector quantity

def v_beam( variables ):
   vBackstream = variables[0]
   vNonBackstream = variables[1]
   vBeam = vBackstream - vNonBackstream
   return vBeam # <- is a vector quantity

def v_beam_ratio( variables ):
   vBackstream = magnitude(variables[0])
   vNonBackstream = magnitude(variables[1])
   divisor = np.ma.masked_less_equal( np.ma.masked_invalid(vNonBackstream),0)
   return np.ma.divide(vBackstream,divisor)

def Bz_linedipole_avg( variables ):
   # This reducer needs to be verified
   x = variables[0]
   y = variables[1]
   z = variables[2]
   dx = variables[3]
   dy = variables[4]
   dz = variables[5]
   return -126.2e6*((dx+x)/(dx*(z**2+(dx+x)**2)) - x/(dx*(z**2+x**2)))

def Bz_linedipole_diff( variables ):
   # This reducer needs to be verified
   Bb = variables[0]
   Bzldp = variables[1]
   print(Bzldp.shape)
   divisor = np.ma.masked_less_equal( np.ma.masked_invalid(magnitude(Bb)),0)
   return np.ma.divide(np.abs(Bb[:,2] - Bzldp), divisor)

def gyrophase_relstddev( variables, velocity_cell_data, velocity_coordinates ):
   # This reducer needs to be verified
   bulk_velocity = variables[0]
   B = variables[1]
   B_unit = B / np.linalg.norm(B)
   
   gyrophase_data = gyrophase_angles(bulk_velocity, B_unit, velocity_cell_data, velocity_coordinates)
   histo = pl.hist(gyrophase_data[0].data, weights=gyrophase_data[1].data, bins=36, range=[-180.0,180.0], log=False, normed=1)
   return np.std(histo[0])/np.mean(histo[0])

def Dng( variables ):
   # This reducer needs to be verified
   # This routine is still very slow due to for-loops
   PTensor = variables[0]
   PParallel = variables[1]
   PPerpendicular = variables[2]
   B = variables[3]
   # following Aunai et al 2103 syntax
   if( np.ndim(PTensor)==2 ):
      B_normalized = np.divide(B, np.sqrt(np.sum(B[:]**2)))
      idMatrix = np.diag(np.ones(3))
      G = PPerpendicular*idMatrix + (PParallel - PPerpendicular) * np.outer(B_normalized, B_normalized)
      N = PTensor - G
      return 2.0 * np.linalg.norm(N, 'fro') / PTensor.trace()
   else:
      B_normalized = np.divide(B, np.sqrt(np.sum(B[:]**2, axis=1))[:,None])
      idMatrix = np.diag(np.ones(3))
      G = [PPerpendicular[i]*idMatrix + (PParallel[i] - PPerpendicular[i]) * np.outer(B_normalized[i], B_normalized[i]) for i in np.arange(len(PParallel))]
      N = PTensor - G
      return [np.divide(2.0*np.linalg.norm(N[i], 'fro'), PTensor[i].trace()) for i in np.arange(len(PParallel))]

def ion_inertial( variables ):
   rho = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[0])),0)
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   charge = vlsvvariables.speciescharge[vlsvvariables.activepopulation]*elementalcharge
   c = 2.9979e8
   epsilonnought = 8.8542e-12
   omegapi = np.sqrt(rho/(mass*epsilonnought))*charge
   di = np.ma.divide(c,omegapi)
   return di

def firstadiabatic( variables ):
   Tperp = variables[0]
   bvector = variables[1]
   B = np.linalg.norm(bvector, axis=-1)
   B = np.ma.masked_less_equal(np.ma.masked_invalid(B),0)
   return np.ma.divide(Tperp,B)

#list of operators. The user can apply these to any variable,
#including more general datareducers. Can only be used to reduce one
#variable at a time
data_operators = {}
data_operators["pass"] = pass_op
data_operators["magnitude"] = magnitude
data_operators["absolute"] = absolute
data_operators["sum"] = sumv
data_operators["x"] = x_component
data_operators["y"] = y_component
data_operators["z"] = z_component

# Hacky way to add vector access
def makelambda(index):
   return lambda data: data[index] if np.ndim(data)==1 else data[:,index]
for i in range(50):
   data_operators[i] = makelambda(i)
   data_operators[str(i)] = data_operators[i]

   
# When vlsvreader tries to read data, it will check using the following order:
# 1) Is the variable directly in the file?
# 2) Is the name something that exists in the file, but only per-population? (answer: sum over populations)
# 3) Is the name a regular datareducer?
# 4) Is there a multipop datareducer for this variable
# The same logic is used for the variables required by datareducers as well.

# if the requested variable starts with "vg_", use the v5 datareducers, otherwise the originals.

#datareducers with more complex, case dependent structure.
datareducers = {}
datareducers["v"] =                      DataReducerVariable(["rho_v", "rho"], v, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vms"] =                    DataReducerVariable(["pressure", "rhom", "b"], vms, "m/s", 1, latex=r"$v_\mathrm{ms}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vs"] =                     DataReducerVariable(["pressure", "rhom"], vs, "m/s", 1, latex=r"$v_\mathrm{s}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["va"] =                     DataReducerVariable(["rhom", "b"], va, "m/s", 1, latex=r"$v_\mathrm{A}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["ma"] =                     DataReducerVariable(["v", "va"], MA, "", 1, latex=r"$M_\mathrm{A}$",latexunits=r"")
datareducers["mms"] =                    DataReducerVariable(["v", "vms"], Mms, "", 1, latex=r"$M_\mathrm{ms}$",latexunits=r"")

datareducers["vparallel"] =              DataReducerVariable(["v", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_\parallel$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vperpendicular"] =         DataReducerVariable(["v", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_\perp$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vparallelbackstream"] =    DataReducerVariable(["vbackstream", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{st}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vperpendicularbackstream"]=DataReducerVariable(["vbackstream", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{st}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vparallelnonbackstream"] =    DataReducerVariable(["vnonbackstream", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{th}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vperpendicularnonbackstream"]=DataReducerVariable(["vnonbackstream", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{th}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

datareducers["eparallel"] =              DataReducerVariable(["e", "b"], ParallelVectorComponent, "V/m", 1, latex=r"$E_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["eperpendicular"] =         DataReducerVariable(["e", "b"], PerpendicularVectorComponent, "V/m", 1, latex=r"$E_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["ejeparallel"] =              DataReducerVariable(["eje", "b"], ParallelVectorComponent, "V/m", 1, latex=r"$EJE_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["ejeperpendicular"] =         DataReducerVariable(["eje", "b"], PerpendicularVectorComponent, "V/m", 1, latex=r"$EJE_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")

datareducers["pdyn"] =            DataReducerVariable(["v", "rhom"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn}$",latexunits=r"Pa")
datareducers["pdynx"] =            DataReducerVariable(["v", "rhom"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,x}$",latexunits=r"Pa")

datareducers["poynting"] = DataReducerVariable(["e", "b"], Poynting, "W/m2", 3, latex=r"$S$", latexunits=r"\mathrm{W}\,\mathrm{m}^{-2}$")
datareducers["di"] =              DataReducerVariable(["proton/rho"], ion_inertial, "m", 1, latex=r"$d_\mathrm{i}$",latexunits=r"$\mathrm{m}$")
datareducers["firstadiabatic"] =    DataReducerVariable(["tperpendicular","b"], firstadiabatic, "K/T", 1, latex=r"$T_\perp B^{-1}$",latexunits=r"$\mathrm{K}\,\mathrm{T}^{-1}$")

# Reducers for simplifying access calls for old and/or new output data versions
datareducers["vbackstream"] =            DataReducerVariable(["rhovbackstream", "rhobackstream"], v, "m/s", 3, latex=r"$V_\mathrm{st}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vnonbackstream"] =         DataReducerVariable(["rhovnonbackstream", "rhononbackstream"], v, "m/s", 3, latex=r"$V_\mathrm{th}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["rhom"] =                   DataReducerVariable(["rho"], rhom, "kg/m3", 1, latex=r"$\rho_m$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
datareducers["rhoq"] =                   DataReducerVariable(["rho"], rhoq, "C/m3", 1, latex=r"$\rho_q$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")
# Reducers for restart files
datareducers["b"] =                      DataReducerVariable(["background_b", "perturbed_b"], restart_B, "T", 3, latex=r"$B$",latexunits=r"T")
datareducers["restart_v"] =              DataReducerVariable(["moments"], restart_V, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["restart_rho"] =            DataReducerVariable(["moments"], restart_rho, "1/m3", 1, latex=r"$n_\mathrm{p}$",latexunits=r"$\mathrm{m}^{-3}$")
datareducers["restart_rhom"] =           DataReducerVariable(["moments"], restart_rhom, "kg/m3", 1, latex=r"$\rho_m$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
datareducers["restart_rhoq"] =           DataReducerVariable(["moments"], restart_rhoq, "C/m3", 1, latex=r"$\rho_q$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")

datareducers["pressure"] =               DataReducerVariable(["ptensordiagonal"], Pressure, "Pa", 1, latex=r"$P$", latexunits=r"Pa")
datareducers["ptensor"] =                DataReducerVariable(["ptensordiagonal", "ptensoroffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}$", latexunits=r"Pa")
datareducers["ptensorrotated"] =         DataReducerVariable(["ptensor", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}$", latexunits=r"Pa")
datareducers["pparallel"] =              DataReducerVariable(["ptensorrotated"], ParallelTensorComponent, "Pa", 1, latex=r"$P_\parallel$", latexunits=r"Pa")
datareducers["pperpendicular"] =         DataReducerVariable(["ptensorrotated"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_\perp$", latexunits=r"Pa")
datareducers["pperpoverpar"] =           DataReducerVariable(["ptensorrotated"], Anisotropy, "", 1, latex=r"$P_\perp P_\parallel^{-1}$", latexunits=r"")
datareducers["panisotropy"] =           DataReducerVariable(["ptensorrotated"], Anisotropy, "", 1, latex=r"$P_\perp P_\parallel^{-1}$", latexunits=r"")
datareducers["agyrotropy"] =             DataReducerVariable(["ptensorrotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag}$", latexunits=r"")

datareducers["pbackstream"] =                 DataReducerVariable(["ptensorbackstreamdiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{st}$", latexunits=r"Pa")
datareducers["ptensorbackstream"] =           DataReducerVariable(["ptensorbackstreamdiagonal", "ptensorbackstreamoffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}$", latexunits=r"Pa")
datareducers["ptensorrotatedbackstream"] =    DataReducerVariable(["ptensorbackstream", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}^\mathrm{R}$", latexunits=r"Pa")
datareducers["pparallelbackstream"] =         DataReducerVariable(["ptensorrotatedbackstream"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{st}}$", latexunits=r"Pa")
datareducers["pperpendicularbackstream"] =    DataReducerVariable(["ptensorrotatedbackstream"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{st}}$", latexunits=r"Pa")
datareducers["pperpoverparbackstream"] =      DataReducerVariable(["ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{st}} P_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["agyrotropybackstream"] =        DataReducerVariable(["ptensorrotatedbackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,st}$", latexunits=r"")

datareducers["pnonbackstream"] =              DataReducerVariable(["ptensornonbackstreamdiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{th}$", latexunits=r"Pa")
datareducers["ptensornonbackstream"] =        DataReducerVariable(["ptensornonbackstreamdiagonal", "ptensornonbackstreamoffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}$", latexunits=r"Pa")
datareducers["ptensorrotatednonbackstream"] = DataReducerVariable(["ptensornonbackstream", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}^\mathrm{R}$", latexunits=r"Pa")
datareducers["pparallelnonbackstream"] =      DataReducerVariable(["ptensorrotatednonbackstream"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{th}}$", latexunits=r"Pa")
datareducers["pperpendicularnonbackstream"] = DataReducerVariable(["ptensorrotatednonbackstream"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{th}}$", latexunits=r"Pa")
datareducers["pperpoverparnonbackstream"] =   DataReducerVariable(["ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{th}} P_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")
datareducers["agyrotropynonbackstream"] =     DataReducerVariable(["ptensorrotatednonbackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,th}$", latexunits=r"")

# Note: Temperature summing over multipop works only if only one population.
datareducers["temperature"] =            DataReducerVariable(["pressure", "rho"], Temperature, "K", 1, latex=r"$T$", latexunits=r"K")
datareducers["ttensor"] =                DataReducerVariable(["ptensor", "rho"], Temperature, "K", 9, latex=r"$\mathcal{T}$", latexunits=r"K")
datareducers["ttensorrotated"] =         DataReducerVariable(["ptensorrotated", "rho"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}$", latexunits=r"K")
datareducers["tparallel"] =              DataReducerVariable(["pparallel", "rho"], Temperature, "K", 1, latex=r"$T_\parallel$", latexunits=r"K")
datareducers["tperpendicular"] =         DataReducerVariable(["pperpendicular", "rho"], Temperature, "K", 1, latex=r"$T_\perp$", latexunits=r"K")

datareducers["tbackstream"] =            DataReducerVariable(["pbackstream", "rhobackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{st}$", latexunits=r"K")
datareducers["ttensorbackstream"] =      DataReducerVariable(["ptensorbackstream", "rhobackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{st}$", latexunits=r"K")
datareducers["ttensorrotatedbackstream"]=DataReducerVariable(["ptensorrotatedbackstream", "rhobackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{st}^\mathrm{R}$", latexunits=r"K")
datareducers["tparallelbackstream"] =    DataReducerVariable(["pparallelbackstream", "rhobackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{st}}$", latexunits=r"K")
datareducers["tperpendicularbackstream"]=DataReducerVariable(["pperpendicularbackstream", "rhobackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{st}}$", latexunits=r"K")

datareducers["tnonbackstream"] =              DataReducerVariable(["pnonbackstream", "rhononbackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{th}$", latexunits=r"K")
datareducers["ttensornonbackstream"] =        DataReducerVariable(["ptensornonbackstream", "rhononbackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{th}$", latexunits=r"K")
datareducers["ttensorrotatednonbackstream"] = DataReducerVariable(["ptensorrotatednonbackstream", "rhononbackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{th}^\mathrm{R}$", latexunits=r"K")
datareducers["tparallelnonbackstream"] =      DataReducerVariable(["pparallelnonbackstream", "rhononbackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{th}}$", latexunits=r"K")
datareducers["tperpendicularnonbackstream"] = DataReducerVariable(["pperpendicularnonbackstream", "rhononbackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{th}}$", latexunits=r"K")

# These ratios are identical to the pressure ratios
datareducers["tperpoverpar"] =                DataReducerVariable(["ptensorrotated"], Anisotropy, "", 1, latex=r"$T_\perp T_\parallel^{-1}$", latexunits=r"")
datareducers["tperpoverparbackstream"] =      DataReducerVariable(["ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{st}} T_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["tperpoverparnonbackstream"] =   DataReducerVariable(["ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{th}} T_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")

# These ratios are identical to the pressure ratiosd
datareducers["betaperpoverpar"] =             DataReducerVariable(["ptensorrotated"], Anisotropy, "", 1, latex=r"$\beta_\perp \beta_\parallel^{-1}$", latexunits=r"")
datareducers["betaperpoverparbackstream"] =   DataReducerVariable(["ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{st}} \beta_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["betaperpoverparnonbackstream"] =DataReducerVariable(["ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{th}} \beta_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")

datareducers["beta"] =                   DataReducerVariable(["pressure", "b"], beta ,"", 1, latex=r"$\beta$", latexunits=r"")
datareducers["betaparallel"] =           DataReducerVariable(["pparallel", "b"], beta ,"", 1, latex=r"$\beta_\parallel$", latexunits=r"")
datareducers["betaperpendicular"] =      DataReducerVariable(["pperpendicular", "b"], beta ,"", 1, latex=r"$\beta_\perp$", latexunits=r"")

datareducers["rmirror"] =                DataReducerVariable(["ptensor", "b"], rMirror, "", 1, latex=r"$R_\mathrm{m}$")
datareducers["dng"] =                    DataReducerVariable(["ptensor", "pparallel", "pperpendicular", "b"], Dng, "", 1, latex=r"$\mathrm{Dng}$") # I think this has vector length 1?
datareducers["vbeam"] =                  DataReducerVariable(["vbackstream", "vnonbackstream"], v_beam, "m/s", 3, latex=r"$V_\mathrm{st}-V$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vbeamratio"] =             DataReducerVariable(["vbackstream", "vnonbackstream"], v_beam_ratio, "", 1, latex=r"$V_\mathrm{st} V^{-1}$", latexunits=r"")
datareducers["thermalvelocity"] =               DataReducerVariable(["temperature"], thermalvelocity, "m/s", 1, latex=r"$v_\mathrm{th}$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["bz_linedipole_avg"] =      DataReducerVariable(["x", "y", "z", "dx", "dy", "dz"], Bz_linedipole_avg, "T", 1, latex=r"$\langle B_{z,\mathrm{ld}}\rangle$")
datareducers["bz_linedipole_diff"] =     DataReducerVariable(["b", "bz_linedipole_avg"], Bz_linedipole_diff, "", 1, latex=r"$\Delta B_{z,\mathrm{ld}}$")

#reducers with useVspace
datareducers["gyrophase_relstddev"] =    DataReducerVariable(["v", "b"], gyrophase_relstddev, "", 1, useVspace=True) # I think this has vector length 1?





#multipopdatareducers
multipopdatareducers = {}
multipopdatareducers["pop/rhom"] =                   DataReducerVariable(["pop/rho"], rhom, "kg/m3", 1, latex=r"$\rho_{m,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
multipopdatareducers["pop/rhoq"] =                   DataReducerVariable(["pop/rho"], rhoq, "C/m3", 1, latex=r"$\rho_{q,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")

multipopdatareducers["pop/pdyn"] =            DataReducerVariable(["pop/v", "pop/rhom"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP}}$",latexunits=r"Pa")
multipopdatareducers["pop/pdynx"] =            DataReducerVariable(["pop/v", "pop/rhom"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP},x}$",latexunits=r"Pa")

multipopdatareducers["pop/vparallel"] =              DataReducerVariable(["pop/v", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopdatareducers["pop/vperpendicular"] =         DataReducerVariable(["pop/v", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopdatareducers["pop/vparallelbackstream"] =         DataReducerVariable(["pop/vbackstream", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{st},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopdatareducers["pop/vparallelnonbackstream"] =         DataReducerVariable(["pop/vnonbackstream", "b"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{th},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopdatareducers["pop/vperpendicularbackstream"] =         DataReducerVariable(["pop/vbackstream", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{st},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopdatareducers["pop/vperpendicularnonbackstream"] =         DataReducerVariable(["pop/vnonbackstream", "b"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{th},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")


multipopdatareducers["pop/pressure"] =               DataReducerVariable(["pop/ptensordiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensor"] =                DataReducerVariable(["pop/ptensordiagonal", "pop/ptensoroffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensorrotated"] =         DataReducerVariable(["pop/ptensor", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/pparallel"] =              DataReducerVariable(["pop/ptensorrotated"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpendicular"] =         DataReducerVariable(["pop/ptensorrotated"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpoverpar"] =           DataReducerVariable(["pop/ptensorrotated"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}} P_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/agyrotropy"] =             DataReducerVariable(["pop/ptensorrotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP}$", latexunits=r"")

multipopdatareducers["pop/pbackstream"] =            DataReducerVariable(["pop/ptensorbackstreamdiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensorbackstream"] =      DataReducerVariable(["pop/ptensorbackstreamdiagonal", "pop/ptensorbackstreamoffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensorrotatedbackstream"]=DataReducerVariable(["pop/ptensorbackstream", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/pparallelbackstream"] =    DataReducerVariable(["pop/ptensorrotatedbackstream"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpendicularbackstream"]=DataReducerVariable(["pop/ptensorrotatedbackstream"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpoverparbackstream"] = DataReducerVariable(["pop/ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}} P_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/agyrotropybackstream"] =   DataReducerVariable(["pop/ptensorrotatedbackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,st}$", latexunits=r"")

multipopdatareducers["pop/pnonbackstream"] =              DataReducerVariable(["pop/ptensornonbackstreamdiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensornonbackstream"] =        DataReducerVariable(["pop/ptensornonbackstreamdiagonal", "pop/ptensornonbackstreamoffdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/ptensorrotatednonbackstream"] = DataReducerVariable(["pop/ptensornonbackstream", "b"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/pparallelnonbackstream"] =      DataReducerVariable(["pop/ptensorrotatednonbackstream"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpendicularnonbackstream"] = DataReducerVariable(["pop/ptensorrotatednonbackstream"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopdatareducers["pop/pperpoverparnonbackstream"] =   DataReducerVariable(["pop/ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}} P_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/agyrotropynonbackstream"] =     DataReducerVariable(["pop/ptensorrotatednonbackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,th}$", latexunits=r"")

multipopdatareducers["pop/temperature"] =            DataReducerVariable(["pop/pressure", "pop/rho"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/ttensor"] =                DataReducerVariable(["pop/ptensor", "pop/rho"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/ttensorrotated"] =         DataReducerVariable(["pop/ptensorrotated", "pop/rho"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/tparallel"] =              DataReducerVariable(["pop/pparallel", "pop/rho"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"K")
multipopdatareducers["pop/tperpendicular"] =         DataReducerVariable(["pop/pperpendicular", "pop/rho"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"K")

multipopdatareducers["pop/tbackstream"] =            DataReducerVariable(["pop/pbackstream", "pop/rhobackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/ttensorbackstream"] =      DataReducerVariable(["pop/ptensorbackstream", "pop/rhobackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/ttensorrotatedbackstream"]=DataReducerVariable(["pop/ptensorrotatedbackstream", "pop/rhobackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/tparallelbackstream"] =    DataReducerVariable(["pop/pparallelbackstream", "pop/rhobackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")
multipopdatareducers["pop/tperpendicularbackstream"]=DataReducerVariable(["pop/pperpendicularbackstream", "pop/rhobackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")

multipopdatareducers["pop/tnonbackstream"] =              DataReducerVariable(["pop/pnonbackstream", "pop/rhononbackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/ttensornonbackstream"] =        DataReducerVariable(["pop/ptensornonbackstream", "pop/rhononbackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/ttensorrotatednonbackstream"] = DataReducerVariable(["pop/ptensorrotatednonbackstream", "pop/rhononbackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/tparallelnonbackstream"] =      DataReducerVariable(["pop/pparallelnonbackstream", "pop/rhononbackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")
multipopdatareducers["pop/tperpendicularnonbackstream"] = DataReducerVariable(["pop/pperpendicularnonbackstream", "pop/rhononbackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")

# These ratios are identical to the pressure ratios
multipopdatareducers["pop/tperpoverpar"] =                 DataReducerVariable(["pop/ptensorrotated"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} T_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/tperpoverparbackstream"] =       DataReducerVariable(["pop/ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}} T_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/tperpoverparnonbackstream"] =    DataReducerVariable(["pop/ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}} T_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")

# These ratios are identical to the pressure ratios
multipopdatareducers["pop/betaperpoverpar"] =              DataReducerVariable(["pop/ptensorrotated"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}} \beta_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/betaperpoverparbackstream"] =    DataReducerVariable(["pop/ptensorrotatedbackstream"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,st}} \beta_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/betaperpoverparnonbackstream"] = DataReducerVariable(["pop/ptensorrotatednonbackstream"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,th}} \beta_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")

multipopdatareducers["pop/thermalvelocity"] =               DataReducerVariable(["temperature"], thermalvelocity, "m/s", 1, latex=r"$v_\mathrm{th,REPLACEPOP}$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopdatareducers["pop/firstadiabatic"] =    DataReducerVariable(["pop/tperpendicular","b"], firstadiabatic, "K/T", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} B^{-1}$",latexunits=r"$\mathrm{K}\,\mathrm{T}^{-1}$")

# Do these betas make sense per-population?
multipopdatareducers["pop/beta"] =                   DataReducerVariable(["pop/pressure", "b"], beta ,"", 1, latex=r"$\beta_\mathrm{REPLACEPOP}$", latexunits=r"")
multipopdatareducers["pop/betaparallel"] =           DataReducerVariable(["pop/pparallel", "b"], beta ,"", 1, latex=r"$\beta_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"")
multipopdatareducers["pop/betaperpendicular"] =      DataReducerVariable(["pop/pperpendicular", "b"], beta ,"", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"")

multipopdatareducers["pop/rmirror"] =                DataReducerVariable(["pop/ptensor", "b"], rMirror, "", 1, latex=r"$R_\mathrm{m,REPLACEPOP}$")
multipopdatareducers["pop/dng"] =                    DataReducerVariable(["pop/ptensor", "pop/pparallel", "pop/pperpendicular", "b"], Dng, "", 1, latex=r"$\mathrm{Dng}_\mathrm{REPLACEPOP}$")


##########################################
#  v5reducers for use with Vlasiator5 data
##########################################

v5reducers = {}
v5reducers["vg_vms"] =                    DataReducerVariable(["vg_pressure", "vg_rhom", "vg_b_vol"], vms, "m/s", 1, latex=r"$v_\mathrm{ms}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
v5reducers["vg_vs"] =                     DataReducerVariable(["vg_pressure", "vg_rhom"], vs, "m/s", 1, latex=r"$v_\mathrm{s}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
v5reducers["vg_va"] =                     DataReducerVariable(["vg_rhom", "vg_b_vol"], va, "m/s", 1, latex=r"$v_\mathrm{A}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
v5reducers["vg_ma"] =                     DataReducerVariable(["vg_v", "vg_va"], MA, "", 1, latex=r"$M_\mathrm{A}$",latexunits=r"")
v5reducers["vg_mms"] =                    DataReducerVariable(["vg_v", "vg_vms"], Mms, "", 1, latex=r"$M_\mathrm{ms}$",latexunits=r"")

v5reducers["vg_v_parallel"] =              DataReducerVariable(["vg_v", "vg_b_vol"], ParallelVectorComponent, "m/s", 1, latex=r"$V_\parallel$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
v5reducers["vg_v_perpendicular"] =         DataReducerVariable(["vg_v", "vg_b_vol"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_\perp$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

#v5reducers["vg_e_parallel"] =              DataReducerVariable(["vg_e_vol", "vg_b_vol"], ParallelVectorComponent, "V/m", 1, latex=r"$E_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
#v5reducers["vg_e_perpendicular"] =         DataReducerVariable(["vg_e_vol", "vg_b_vol"], PerpendicularVectorComponent, "V/m", 1, latex=r"$E_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
#v5reducers["vg_poynting"] = DataReducerVariable(["vg_e_vol", "vg_b_vol"], Poynting, "W/m2", 3, latex=r"$S$", latexunits=r"\mathrm{W}\,\mathrm{m}^{-2}$")

v5reducers["vg_eje_parallel"] =              DataReducerVariable(["vg_eje", "vg_b_vol"], ParallelVectorComponent, "V/m", 1, latex=r"$EJE_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
v5reducers["vg_eje_perpendicular"] =         DataReducerVariable(["vg_eje", "vg_b_vol"], PerpendicularVectorComponent, "V/m", 1, latex=r"$EJE_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")

v5reducers["vg_pdyn"] =            DataReducerVariable(["vg_v", "vg_rhom"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn}$",latexunits=r"Pa")
v5reducers["vg_pdynx"] =            DataReducerVariable(["vg_v", "vg_rhom"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,x}$",latexunits=r"Pa")

v5reducers["vg_di"] =              DataReducerVariable(["proton/vg_rho"], ion_inertial, "m", 1, latex=r"$d_\mathrm{i}$",latexunits=r"$\mathrm{m}$")

v5reducers["vg_pressure"] =               DataReducerVariable(["vg_ptensor_diagonal"], Pressure, "Pa", 1, latex=r"$P$", latexunits=r"Pa")
v5reducers["vg_ptensor"] =                DataReducerVariable(["vg_ptensor_diagonal", "vg_ptensor_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}$", latexunits=r"Pa")
v5reducers["vg_ptensor_rotated"] =         DataReducerVariable(["vg_ptensor", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}$", latexunits=r"Pa")
v5reducers["vg_p_parallel"] =              DataReducerVariable(["vg_ptensor_rotated"], ParallelTensorComponent, "Pa", 1, latex=r"$P_\parallel$", latexunits=r"Pa")
v5reducers["vg_p_perpendicular"] =         DataReducerVariable(["vg_ptensor_rotated"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_\perp$", latexunits=r"Pa")
v5reducers["vg_p_anisotropy"] =           DataReducerVariable(["vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$P_\perp P_\parallel^{-1}$", latexunits=r"")
v5reducers["vg_agyrotropy"] =             DataReducerVariable(["vg_ptensor_rotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag}$", latexunits=r"")

v5reducers["vg_p_nonthermal"] =                 DataReducerVariable(["vg_ptensor_nonthermal_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{st}$", latexunits=r"Pa")
v5reducers["vg_ptensor_nonthermal"] =           DataReducerVariable(["vg_ptensor_nonthermal_diagonal", "vg_ptensor_nonthermal_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}$", latexunits=r"Pa")
v5reducers["vg_ptensor_rotated_nonthermal"] =    DataReducerVariable(["vg_ptensor_nonthermal", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}^\mathrm{R}$", latexunits=r"Pa")
v5reducers["vg_p_parallel_nonthermal"] =         DataReducerVariable(["vg_ptensor_rotated_nonthermal"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{st}}$", latexunits=r"Pa")
v5reducers["vg_p_perpendicular_nonthermal"] =    DataReducerVariable(["vg_ptensor_rotated_nonthermal"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{st}}$", latexunits=r"Pa")
v5reducers["vg_p_anisotropy_nonthermal"] =      DataReducerVariable(["vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{st}} P_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
v5reducers["vg_agyrotropy_nonthermal"] =        DataReducerVariable(["vg_ptensor_rotated_nonthermal"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,st}$", latexunits=r"")

v5reducers["vg_p_thermal"] =              DataReducerVariable(["vg_ptensor_thermal_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{th}$", latexunits=r"Pa")
v5reducers["vg_ptensor_thermal"] =        DataReducerVariable(["vg_ptensor_thermal_diagonal", "vg_ptensor_thermal_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}$", latexunits=r"Pa")
v5reducers["vg_ptensor_rotated_thermal"] = DataReducerVariable(["vg_ptensor_thermal", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}^\mathrm{R}$", latexunits=r"Pa")
v5reducers["vg_p_parallel_thermal"] =      DataReducerVariable(["vg_ptensor_rotated_thermal"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{th}}$", latexunits=r"Pa")
v5reducers["vg_p_perpendicular_thermal"] = DataReducerVariable(["vg_ptensor_rotated_thermal"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{th}}$", latexunits=r"Pa")
v5reducers["vg_p_anisotropy_thermal"] =   DataReducerVariable(["vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{th}} P_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")
v5reducers["vg_agyrotropy_thermal"] =     DataReducerVariable(["vg_ptensor_rotated_thermal"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,th}$", latexunits=r"")

# Note: Temperature summing over multipop works only if only one population.
v5reducers["vg_temperature"] =            DataReducerVariable(["vg_pressure", "vg_rhom"], Temperature, "K", 1, latex=r"$T$", latexunits=r"K")
v5reducers["vg_ttensor"] =                DataReducerVariable(["vg_ptensor", "vg_rhom"], Temperature, "K", 9, latex=r"$\mathcal{T}$", latexunits=r"K")
v5reducers["vg_ttensor_rotated"] =         DataReducerVariable(["vg_ptensor_rotated", "vg_rhom"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}$", latexunits=r"K")
v5reducers["vg_t_parallel"] =              DataReducerVariable(["vg_p_parallel", "vg_rhom"], Temperature, "K", 1, latex=r"$T_\parallel$", latexunits=r"K")
v5reducers["vg_t_perpendicular"] =         DataReducerVariable(["vg_p_perpendicular", "vg_rhom"], Temperature, "K", 1, latex=r"$T_\perp$", latexunits=r"K")

 # These ratios are identical to the pressure ratios
v5reducers["vg_t_anisotropy"] =                DataReducerVariable(["vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$T_\perp T_\parallel^{-1}$", latexunits=r"")
v5reducers["vg_t_anisotropy_nonthermal"] =      DataReducerVariable(["vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{st}} T_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
v5reducers["vg_t_anisotropy_thermal"] =   DataReducerVariable(["vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{th}} T_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")
 # These ratios are identical to the pressure ratios
v5reducers["vg_beta_anisotropy"] =             DataReducerVariable(["vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$\beta_\perp \beta_\parallel^{-1}$", latexunits=r"")
v5reducers["vg_beta_anisotropy_nonthermal"] =   DataReducerVariable(["vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{st}} \beta_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
v5reducers["vg_beta_anisotropy_thermal"] =DataReducerVariable(["vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{th}} \beta_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")

v5reducers["vg_beta"] =                   DataReducerVariable(["vg_pressure", "vg_b_vol"], beta ,"", 1, latex=r"$\beta$", latexunits=r"")
v5reducers["vg_beta_parallel"] =           DataReducerVariable(["vg_p_parallel", "vg_b_vol"], beta ,"", 1, latex=r"$\beta_\parallel$", latexunits=r"")
v5reducers["vg_beta_perpendicular"] =      DataReducerVariable(["vg_p_perpendicular", "vg_b_vol"], beta ,"", 1, latex=r"$\beta_\perp$", latexunits=r"")

v5reducers["vg_rmirror"] =                DataReducerVariable(["vg_ptensor", "vg_b_vol"], rMirror, "", 1, latex=r"$R_\mathrm{m}$")

v5reducers["vg_restart_v"] =              DataReducerVariable(["moments"], restart_V, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
v5reducers["vg_restart_rho"] =            DataReducerVariable(["moments"], restart_rho, "1/m3", 1, latex=r"$n_\mathrm{p}$",latexunits=r"$\mathrm{m}^{-3}$")
v5reducers["vg_restart_rhom"] =           DataReducerVariable(["moments"], restart_rhom, "kg/m3", 1, latex=r"$\rho_m$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
v5reducers["vg_restart_rhoq"] =           DataReducerVariable(["moments"], restart_rhoq, "C/m3", 1, latex=r"$\rho_q$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")


#multipopv5reducers
multipopv5reducers = {}
multipopv5reducers["pop/vg_rhom"] =                   DataReducerVariable(["pop/vg_rho"], rhom, "kg/m3", 1, latex=r"$\rho_{m,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
multipopv5reducers["pop/vg_rhoq"] =                   DataReducerVariable(["pop/vg_rho"], rhoq, "C/m3", 1, latex=r"$\rho_{q,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")
multipopv5reducers["pop/vg_pdyn"] =            DataReducerVariable(["pop/vg_v", "pop/vg_rhom"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP}}$",latexunits=r"Pa")
multipopv5reducers["pop/vg_pdynx"] =            DataReducerVariable(["pop/vg_v", "pop/vg_rhom"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP},x}$",latexunits=r"Pa")

multipopv5reducers["pop/vg_v_parallel"] =              DataReducerVariable(["pop/vg_v", "vg_b_vol"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopv5reducers["pop/vg_v_perpendicular"] =         DataReducerVariable(["pop/vg_v", "vg_b_vol"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopv5reducers["pop/vg_v_parallel_thermal"] =              DataReducerVariable(["pop/vg_v_thermal", "vg_b_vol"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{th},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopv5reducers["pop/vg_v_perpendicular_thermal"] =         DataReducerVariable(["pop/vg_v_thermal", "vg_b_vol"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{th},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopv5reducers["pop/vg_v_parallel_nonthermal"] =              DataReducerVariable(["pop/vg_v_nonthermal", "vg_b_vol"], ParallelVectorComponent, "m/s", 1, latex=r"$V_{\parallel,\mathrm{st},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopv5reducers["pop/vg_v_perpendicular_nonthermal"] =         DataReducerVariable(["pop/vg_v_nonthermal", "vg_b_vol"], PerpendicularVectorComponent, "m/s", 1, latex=r"$V_{\perp,\mathrm{st},\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopv5reducers["pop/vg_pressure"] =               DataReducerVariable(["pop/vg_ptensor_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor"] =                DataReducerVariable(["pop/vg_ptensor_diagonal", "pop/vg_ptensor_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor_rotated"] =         DataReducerVariable(["pop/vg_ptensor", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_parallel"] =              DataReducerVariable(["pop/vg_ptensor_rotated"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_perpendicular"] =         DataReducerVariable(["pop/vg_ptensor_rotated"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_anisotropy"] =           DataReducerVariable(["pop/vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}} P_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_agyrotropy"] =             DataReducerVariable(["pop/vg_ptensor_rotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP}$", latexunits=r"")

multipopv5reducers["pop/vg_p_nonthermal"] =            DataReducerVariable(["pop/vg_ptensor_nonthermal_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor_nonthermal"] =      DataReducerVariable(["pop/vg_ptensor_nonthermal_diagonal", "pop/vg_ptensor_nonthermal_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor_rotated_nonthermal"]=DataReducerVariable(["pop/vg_ptensor_nonthermal", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_parallel_nonthermal"] =    DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_perpendicular_nonthermal"]=DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_anisotropy_nonthermal"] = DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}} P_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_agyrotropy_nonthermal"] =   DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,st}$", latexunits=r"")

multipopv5reducers["pop/vg_p_thermal"] =              DataReducerVariable(["pop/vg_ptensor_thermal_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor_thermal"] =        DataReducerVariable(["pop/vg_ptensor_thermal_diagonal", "pop/vg_ptensor_thermal_offdiagonal"], FullTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_ptensor_rotated_thermal"] = DataReducerVariable(["pop/vg_ptensor_thermal", "vg_b_vol"], RotatedTensor, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_parallel_thermal"] =      DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], ParallelTensorComponent, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_perpendicular_thermal"] = DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], PerpendicularTensorComponent, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopv5reducers["pop/vg_p_anisotropy_thermal"] =   DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}} P_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_agyrotropy_thermal"] =     DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,th}$", latexunits=r"")

multipopv5reducers["pop/vg_temperature"] =            DataReducerVariable(["pop/vg_pressure", "pop/vg_rho"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor"] =                DataReducerVariable(["pop/vg_ptensor", "pop/vg_rho"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor_rotated"] =         DataReducerVariable(["pop/vg_ptensor_rotated", "pop/vg_rho"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_parallel"] =              DataReducerVariable(["pop/vg_p_parallel", "pop/vg_rho"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_perpendicular"] =         DataReducerVariable(["pop/vg_p_perpendicular", "pop/vg_rho"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"K")

multipopv5reducers["pop/vg_t_nonthermal"] =            DataReducerVariable(["pop/vg_p_nonthermal", "pop/vg_rho_nonthermal"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor_nonthermal"] =      DataReducerVariable(["pop/vg_ptensor_nonthermal", "pop/vg_rho_nonthermal"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor_rotated_nonthermal"]=DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal", "pop/vg_rho_nonthermal"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_parallel_nonthermal"] =    DataReducerVariable(["pop/vg_p_parallel_nonthermal", "pop/vg_rho_nonthermal"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_perpendicular_nonthermal"]=DataReducerVariable(["pop/vg_p_perpendicular_nonthermal", "pop/vg_rho_nonthermal"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")

multipopv5reducers["pop/vg_t_thermal"] =              DataReducerVariable(["pop/vg_p_thermal", "pop/vg_rho_thermal"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor_thermal"] =        DataReducerVariable(["pop/vg_ptensor_thermal", "pop/vg_rho_thermal"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopv5reducers["pop/vg_ttensor_rotated_thermal"] = DataReducerVariable(["pop/vg_ptensor_rotated_thermal", "pop/vg_rho_thermal"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_parallel_thermal"] =      DataReducerVariable(["pop/vg_p_parallel_thermal", "pop/vg_rho_thermal"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")
multipopv5reducers["pop/vg_t_perpendicular_thermal"] = DataReducerVariable(["pop/vg_p_perpendicular_thermal", "pop/vg_rho_thermal"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")

 # These ratios are identical to the pressure ratios
multipopv5reducers["pop/vg_t_anisotropy"] =                 DataReducerVariable(["pop/vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} T_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_t_anisotropy_nonthermal"] =       DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}} T_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_t_anisotropy_thermal"] =    DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}} T_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_beta_anisotropy"] =              DataReducerVariable(["pop/vg_ptensor_rotated"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}} \beta_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_beta_anisotropy_nonthermal"] =    DataReducerVariable(["pop/vg_ptensor_rotated_nonthermal"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,st}} \beta_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopv5reducers["pop/vg_beta_anisotropy_thermal"] = DataReducerVariable(["pop/vg_ptensor_rotated_thermal"], Anisotropy, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,th}} \beta_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")

multipopv5reducers["pop/vg_thermalvelocity"] =               DataReducerVariable(["vg_temperature"], thermalvelocity, "m/s", 1, latex=r"$v_\mathrm{th,REPLACEPOP}$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopv5reducers["pop/vg_firstadiabatic"] =    DataReducerVariable(["pop/vg_t_perpendicular","vg_b_vol"], firstadiabatic, "K/T", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} B^{-1}$",latexunits=r"$\mathrm{K}\,\mathrm{T}^{-1}$")

# Do these betas make sense per-population?
multipopv5reducers["pop/vg_beta"] =                   DataReducerVariable(["pop/vg_pressure", "vg_b_vol"], beta ,"", 1, latex=r"$\beta_\mathrm{REPLACEPOP}$", latexunits=r"")
multipopv5reducers["pop/vg_beta_parallel"] =           DataReducerVariable(["pop/vg_p_parallel", "vg_b_vol"], beta ,"", 1, latex=r"$\beta_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"")
multipopv5reducers["pop/vg_beta_perpendicular"] =      DataReducerVariable(["pop/vg_p_perpendicular", "vg_b_vol"], beta ,"", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"")

multipopv5reducers["pop/vg_rmirror"] =                DataReducerVariable(["pop/vg_ptensor", "vg_b_vol"], rMirror, "", 1, latex=r"$R_\mathrm{m,REPLACEPOP}$")
multipopv5reducers["pop/vg_dng"] =                    DataReducerVariable(["pop/vg_ptensor", "pop/vg_p_parallel", "pop/vg_p_perpendicular", "vg_b_vol"], Dng, "", 1, latex=r"$\mathrm{Dng}_\mathrm{REPLACEPOP}$")


multipopv5reducers["pop/vg_temperature"] =            DataReducerVariable(["pop/vg_pressure", "pop/vg_rho"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopv5reducers["pop/vg_pressure"] =               DataReducerVariable(["pop/vg_ptensor_diagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP}$", latexunits=r"Pa")

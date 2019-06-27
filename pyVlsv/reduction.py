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
   mp = 1.672622e-27

   if np.ndim(moments)==1: # single cell
      if len(moments)==4:  # pre-multipop restart
         rhom = moments[0]*mp
      elif len(moments)==5: # multipop restart
         rhom = moments[0]
      else:
         print("Unrecognized length for moments!")
         return None
   else: # array of cells
      if len(moments[0,:])==4:  # pre-multipop restart
         rhom = moments[:,0]*mp
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
   charge = 1.6021773e-19

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
   mp = 1.672622e-27
   return rho*mp

def rhoq( variables ):
   ''' Data reducer function for calculating rhoq from pre-multipop file
   '''
   rho = variables[0]
   charge = 1.6021773e-19
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
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   mu_0 = 1.25663706144e-6
   P = variables[0]
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[1])),0) *mp
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
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   P = variables[0]
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[1])),0) *mp
   vs = np.sqrt( np.ma.divide( P*5.0/3.0, rho_m ) )
   return vs

def va( variables ):
   ''' Data reducer function for getting the Alfven velocity
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   mu_0 = 1.25663706144e-6
   rho_m = np.ma.masked_less_equal(np.ma.masked_invalid(np.array(variables[0])),0)*mp
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

def VParallel( variables ):
   ''' Data reducer function for bulk flow parallel to the magnetic field
   '''
   bulkv = variables[0]
   bvector = variables[1]
   if( np.ndim(bulkv)==1 ):
      bnorm = bvector / np.linalg.norm(bvector)
      return (bulkv*bnorm).sum()
   else:
      bnorm = bvector / np.linalg.norm(bvector, axis=-1)[:,np.newaxis]
      return (bulkv*bnorm).sum(-1)

def VPerpendicular( variables ):
   ''' Data reducer function for bulk flow perpendicular to the magnetic field
   '''
   bulkv = variables[0]
   bvector = variables[1]
   if( np.ndim(bulkv)==1 ):
      bnorm = bvector / np.linalg.norm(bvector)
      vpara = (bulkv*bnorm).sum()
      vmag = np.linalg.norm(bulkv)
      return np.sqrt(vmag*vmag - vpara*vpara)
   else:
      bnorm = bvector / np.linalg.norm(bvector, axis=-1)[:,np.newaxis]
      vpara = (bulkv*bnorm).sum(-1)
      vmag = np.linalg.norm(bulkv, axis=-1)
      return np.sqrt(vmag*vmag - vpara*vpara)
   
def PTensor( variables ):
   ''' Data reducer function to reconstruct the pressure tensor from
       the vlsv diagonal and off-diagonal components.
   '''
   PTensorDiagonal = np.array(variables[0])
   PTensorOffDiagonal = np.array(variables[1])
   if(np.ndim(PTensorDiagonal)==1 ):      
      PTensorDiagonal = PTensorDiagonal.reshape(1,3)[0]
      PTensorOffDiagonal = PTensorOffDiagonal.reshape(1,3)[0]
      return np.array([[PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]],
                       [PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]],
                       [PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]]])
   else:
      result = np.empty([len(PTensorDiagonal[:,0]),3,3]) # Warning, unitialized!
      result[:,0,0] = PTensorDiagonal[:,0]
      result[:,0,1] = PTensorOffDiagonal[:,2]
      result[:,0,2] = PTensorOffDiagonal[:,1]
      result[:,1,0] = PTensorOffDiagonal[:,2]
      result[:,1,1] = PTensorDiagonal[:,1]
      result[:,1,2] = PTensorOffDiagonal[:,0]
      result[:,2,0] = PTensorOffDiagonal[:,1]
      result[:,2,1] = PTensorOffDiagonal[:,0]
      result[:,2,2] = PTensorDiagonal[:,2]
      return result

def PTensorRotated( variables ):
   ''' Data reducer for rotating the pressure tensor to align the z-component with the magnetic field
   '''
   PTensor = variables[0]
   B = variables[1]
   if( np.ndim(B)==1 ):
      B = B.reshape(1,3)[0]
      PTensor = PTensor.reshape(3,3)
      return rotateTensorToVector(PTensor, B)
   else:
      return rotateArrayTensorToVector(PTensor, B)

def PParallel( variables ):
   ''' Data reducer for finding the parallel pressure
   '''
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return PTensorRotated[2,2]
   else:
      return PTensorRotated[:,2,2]

def PPerpendicular( variables ):
   ''' Data reducer for finding the perpendicular pressure
   '''
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return 0.5*(PTensorRotated[0,0] + PTensorRotated[1,1])
   else:
      return 0.5*(PTensorRotated[:,0,0] + PTensorRotated[:,1,1])

def PPerpOverPar( variables ):
   ''' Data reducer for finding the ratio of perpendicular to parallel pressure
   '''
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      divisor = np.ma.masked_equal(np.ma.masked_invalid(PTensorRotated[2,2]),0)
      return 0.5*np.ma.divide(PTensorRotated[0,0] + PTensorRotated[1,1], divisor)
   else:
      divisor = np.ma.masked_equal(np.ma.masked_invalid(PTensorRotated[:,2,2]),0)
      return 0.5*np.ma.divide(PTensorRotated[:,0,0] + PTensorRotated[:,1,1], divisor)

def Pressure( variables ):
   ''' Data reducer for finding the scalar pressure
   '''
   PTensorDiagonal = variables[0]
   return 1.0/3.0 * np.ma.sum(np.ma.masked_invalid(PTensorDiagonal),axis=-1)

def Pdyn( variables ):
   ''' Data reducer function for dynamic pressure
   '''
   mp = 1.672622e-27
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   Vmag = np.linalg.norm(np.array(variables[0]), axis=-1)
   rhom = np.array(variables[1])*mass
   return Vmag*Vmag*rhom

def Pdynx( variables ):
   ''' Data reducer function for dynamic pressure with just V_x
   '''
   mp = 1.672622e-27
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   rhom = np.array(variables[1])*mass
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
   Pressure = variables[0] # either a tensor, vector, array, or value
   rho = variables[1] # eithern array or a value
   divisor = np.ma.masked_less_equal( np.ma.masked_invalid(np.array(rho)),0) * 1.38065e-23
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
   return 2.0 * 1.25663706144e-6 * Pressure / np.sum(np.asarray(Magneticfield)**2,axis=-1)

def rMirror( variables ):
   # More efficient file access, now just takes PTensor and B
   PT = variables[0]
   B = variables[1]
   PTRotated =  PTensorRotated([PT,B])
   TAniso = PPerpOverPar([PTRotated]) # PAniso == TAniso
   PPerp = PPerpendicular([PTRotated])
   betaPerp = beta([PPerp,B])
   return betaPerp * (TAniso - 1)   

def v_thermal( variables ):
   Temperature = variables[0]
   k = 1.38065e-23
   mp = 1.672622e-27
   mass = vlsvvariables.speciesamu[vlsvvariables.activepopulation]*mp
   # Corrected to calculate the mean speed sqrt(8kT/pi m)
   vThermal = np.sqrt(Temperature*(k*8./(mass*3.14159)))
   return vThermal

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
   print Bzldp.shape
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
   charge = 1.6021773e-19
   mp = 1.672622e-27
   c = 2.9979e8
   epsilonnaught = 8.8542e-12
   omegapi = np.sqrt(rho/(mp*epsilonnaught))*charge
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

#datareducers with more complex, case dependent structure.
datareducers = {}
datareducers["V"] =                      DataReducerVariable(["rho_v", "rho"], v, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$") # Scripts should transition to using capital V
datareducers["v"] =                      DataReducerVariable(["rho_v", "rho"], v, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vms"] =                    DataReducerVariable(["Pressure", "rho", "B"], vms, "m/s", 1, latex=r"$v_\mathrm{ms}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vs"] =                     DataReducerVariable(["Pressure", "rho"], vs, "m/s", 1, latex=r"$v_\mathrm{s}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["va"] =                     DataReducerVariable(["rho", "B"], va, "m/s", 1, latex=r"$v_\mathrm{A}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["MA"] =                     DataReducerVariable(["V", "va"], MA, "", 1, latex=r"$M_\mathrm{A}$",latexunits=r"")
datareducers["Mms"] =                    DataReducerVariable(["V", "vms"], Mms, "", 1, latex=r"$M_\mathrm{ms}$",latexunits=r"")

datareducers["VParallel"] =              DataReducerVariable(["V", "B"], VParallel, "m/s", 1, latex=r"$V_\parallel$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VPerpendicular"] =         DataReducerVariable(["V", "B"], VPerpendicular, "m/s", 1, latex=r"$V_\perp$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VParallelBackstream"] =    DataReducerVariable(["VBackstream", "B"], VParallel, "m/s", 1, latex=r"$V_{\parallel,\mathrm{st}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VPerpendicularBackstream"]=DataReducerVariable(["VBackstream", "B"], VPerpendicular, "m/s", 1, latex=r"$V_{\perp,\mathrm{st}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VParallelNonBackstream"] =    DataReducerVariable(["VNonBackstream", "B"], VParallel, "m/s", 1, latex=r"$V_{\parallel,\mathrm{th}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VPerpendicularNonBackstream"]=DataReducerVariable(["VNonBackstream", "B"], VPerpendicular, "m/s", 1, latex=r"$V_{\perp,\mathrm{th}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

datareducers["EParallel"] =              DataReducerVariable(["E", "B"], VParallel, "V/m", 1, latex=r"$E_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["EPerpendicular"] =         DataReducerVariable(["E", "B"], VPerpendicular, "V/m", 1, latex=r"$E_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["EJEParallel"] =              DataReducerVariable(["EJE", "B"], VParallel, "V/m", 1, latex=r"$EJE_\parallel$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")
datareducers["EJEPerpendicular"] =         DataReducerVariable(["EJE", "B"], VPerpendicular, "V/m", 1, latex=r"$EJE_\perp$",latexunits=r"$\mathrm{V}\,\mathrm{m}^{-1}$")

datareducers["Pdyn"] =            DataReducerVariable(["V", "rho"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn}$",latexunits=r"Pa")
datareducers["Pdynx"] =            DataReducerVariable(["V", "rho"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,x}$",latexunits=r"Pa")

datareducers["Poynting"] = DataReducerVariable(["E", "B"], Poynting, "W/m2", 3, latex=r"$S$", latexunits=r"\mathrm{W}\,\mathrm{m}^{-2}$")
datareducers["di"] =              DataReducerVariable(["rho"], ion_inertial, "m", 1, latex=r"$d_\mathrm{i}$",latexunits=r"$\mathrm{m}$")
datareducers["dimp"] =              DataReducerVariable(["proton_rho"], ion_inertial, "m", 1, latex=r"$d_\mathrm{i}$",latexunits=r"$\mathrm{m}$")
datareducers["firstadiabatic"] =    DataReducerVariable(["TPerpendicular","B"], firstadiabatic, "K/T", 1, latex=r"$T_\perp B^{-1}$",latexunits=r"$\mathrm{K}\,\mathrm{T}^{-1}$")

# Reducers for simplifying access calls for old and/or new output data versions
datareducers["VBackstream"] =            DataReducerVariable(["RhoVBackstream", "RhoBackstream"], v, "m/s", 3, latex=r"$V_\mathrm{st}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["VNonBackstream"] =         DataReducerVariable(["RhoVNonBackstream", "RhoNonBackstream"], v, "m/s", 3, latex=r"$V_\mathrm{th}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["rhom"] =                   DataReducerVariable(["rho"], rhom, "kg/m3", 1, latex=r"$\rho_m$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
datareducers["rhoq"] =                   DataReducerVariable(["rho"], rhoq, "C/m3", 1, latex=r"$\rho_q$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")
# Reducers for restart files
datareducers["B"] =                      DataReducerVariable(["background_B", "perturbed_B"], restart_B, "T", 3, latex=r"$B$",latexunits=r"T")
datareducers["restart_V"] =              DataReducerVariable(["moments"], restart_V, "m/s", 3, latex=r"$V$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["restart_rho"] =            DataReducerVariable(["moments"], restart_rho, "1/m3", 1, latex=r"$n_\mathrm{p}$",latexunits=r"$\mathrm{m}^{-3}$")
datareducers["restart_rhom"] =           DataReducerVariable(["moments"], restart_rhom, "kg/m3", 1, latex=r"$\rho_m$",latexunits=r"$\mathrm{kg}\,\mathrm{m}^{-3}$")
datareducers["restart_rhoq"] =           DataReducerVariable(["moments"], restart_rhoq, "C/m3", 1, latex=r"$\rho_q$",latexunits=r"$\mathrm{C}\,\mathrm{m}^{-3}$")

datareducers["Pressure"] =               DataReducerVariable(["PTensorDiagonal"], Pressure, "Pa", 1, latex=r"$P$", latexunits=r"Pa")
datareducers["PTensor"] =                DataReducerVariable(["PTensorDiagonal", "PTensorOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}$", latexunits=r"Pa")
datareducers["PTensorRotated"] =         DataReducerVariable(["PTensor", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}$", latexunits=r"Pa")
datareducers["PParallel"] =              DataReducerVariable(["PTensorRotated"], PParallel, "Pa", 1, latex=r"$P_\parallel$", latexunits=r"Pa")
datareducers["PPerpendicular"] =         DataReducerVariable(["PTensorRotated"], PPerpendicular, "Pa", 1, latex=r"$P_\perp$", latexunits=r"Pa")
datareducers["PPerpOverPar"] =           DataReducerVariable(["PTensorRotated"], PPerpOverPar, "", 1, latex=r"$P_\perp P_\parallel^{-1}$", latexunits=r"")
datareducers["aGyrotropy"] =             DataReducerVariable(["PTensorRotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag}$", latexunits=r"")

datareducers["PBackstream"] =                 DataReducerVariable(["PTensorBackstreamDiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{st}$", latexunits=r"Pa")
datareducers["PTensorBackstream"] =           DataReducerVariable(["PTensorBackstreamDiagonal", "PTensorBackstreamOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}$", latexunits=r"Pa")
datareducers["PTensorRotatedBackstream"] =    DataReducerVariable(["PTensorBackstream", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{st}^\mathrm{R}$", latexunits=r"Pa")
datareducers["PParallelBackstream"] =         DataReducerVariable(["PTensorRotatedBackstream"], PParallel, "Pa", 1, latex=r"$P_{\parallel,\mathrm{st}}$", latexunits=r"Pa")
datareducers["PPerpendicularBackstream"] =    DataReducerVariable(["PTensorRotatedBackstream"], PPerpendicular, "Pa", 1, latex=r"$P_{\perp,\mathrm{st}}$", latexunits=r"Pa")
datareducers["PPerpOverParBackstream"] =      DataReducerVariable(["PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$P_{\perp,\mathrm{st}} P_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["aGyrotropyBackstream"] =        DataReducerVariable(["PTensorRotatedBackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,st}$", latexunits=r"")

datareducers["PNonBackstream"] =              DataReducerVariable(["PTensorNonBackstreamDiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{th}$", latexunits=r"Pa")
datareducers["PTensorNonBackstream"] =        DataReducerVariable(["PTensorNonBackstreamDiagonal", "PTensorNonBackstreamOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}$", latexunits=r"Pa")
datareducers["PTensorRotatedNonBackstream"] = DataReducerVariable(["PTensorNonBackstream", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{th}^\mathrm{R}$", latexunits=r"Pa")
datareducers["PParallelNonBackstream"] =      DataReducerVariable(["PTensorRotatedNonBackstream"], PParallel, "Pa", 1, latex=r"$P_{\parallel,\mathrm{th}}$", latexunits=r"Pa")
datareducers["PPerpendicularNonBackstream"] = DataReducerVariable(["PTensorRotatedNonBackstream"], PPerpendicular, "Pa", 1, latex=r"$P_{\perp,\mathrm{th}}$", latexunits=r"Pa")
datareducers["PPerpOverParNonBackstream"] =   DataReducerVariable(["PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$P_{\perp,\mathrm{th}} P_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")
datareducers["aGyrotropyNonBackstream"] =     DataReducerVariable(["PTensorRotatedNonBackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,th}$", latexunits=r"")

datareducers["Temperature"] =            DataReducerVariable(["Pressure", "rho"], Temperature, "K", 1, latex=r"$T$", latexunits=r"K")
datareducers["TTensor"] =                DataReducerVariable(["PTensor", "rho"], Temperature, "K", 9, latex=r"$\mathcal{T}$", latexunits=r"K")
datareducers["TTensorRotated"] =         DataReducerVariable(["PTensorRotated", "rho"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}$", latexunits=r"K")
datareducers["TParallel"] =              DataReducerVariable(["PParallel", "rho"], Temperature, "K", 1, latex=r"$T_\parallel$", latexunits=r"K")
datareducers["TPerpendicular"] =         DataReducerVariable(["PPerpendicular", "rho"], Temperature, "K", 1, latex=r"$T_\perp$", latexunits=r"K")

datareducers["TBackstream"] =            DataReducerVariable(["PBackstream", "RhoBackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{st}$", latexunits=r"K")
datareducers["TTensorBackstream"] =      DataReducerVariable(["PTensorBackstream", "RhoBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{st}$", latexunits=r"K")
datareducers["TTensorRotatedBackstream"]=DataReducerVariable(["PTensorRotatedBackstream", "RhoBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{st}^\mathrm{R}$", latexunits=r"K")
datareducers["TParallelBackstream"] =    DataReducerVariable(["PParallelBackstream", "RhoBackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{st}}$", latexunits=r"K")
datareducers["TPerpendicularBackstream"]=DataReducerVariable(["PPerpendicularBackstream", "RhoBackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{st}}$", latexunits=r"K")

datareducers["TNonBackstream"] =              DataReducerVariable(["PNonBackstream", "RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{th}$", latexunits=r"K")
datareducers["TTensorNonBackstream"] =        DataReducerVariable(["PTensorNonBackstream", "RhoNonBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{th}$", latexunits=r"K")
datareducers["TTensorRotatedNonBackstream"] = DataReducerVariable(["PTensorRotatedNonBackstream", "RhoNonBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{th}^\mathrm{R}$", latexunits=r"K")
datareducers["TParallelNonBackstream"] =      DataReducerVariable(["PParallelNonBackstream", "RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{th}}$", latexunits=r"K")
datareducers["TPerpendicularNonBackstream"] = DataReducerVariable(["PPerpendicularNonBackstream", "RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{th}}$", latexunits=r"K")
 # These ratios are identical to the pressure ratios
datareducers["TPerpOverPar"] =                DataReducerVariable(["PTensorRotated"], PPerpOverPar, "", 1, latex=r"$T_\perp T_\parallel^{-1}$", latexunits=r"")
datareducers["TPerpOverParBackstream"] =      DataReducerVariable(["PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$T_{\perp,\mathrm{st}} T_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["TPerpOverParNonBackstream"] =   DataReducerVariable(["PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$T_{\perp,\mathrm{th}} T_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")
datareducers["betaPerpOverPar"] =             DataReducerVariable(["PTensorRotated"], PPerpOverPar, "", 1, latex=r"$\beta_\perp \beta_\parallel^{-1}$", latexunits=r"")
datareducers["betaPerpOverParBackstream"] =   DataReducerVariable(["PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$\beta_{\perp,\mathrm{st}} \beta_{\parallel,\mathrm{st}}^{-1}$", latexunits=r"")
datareducers["betaPerpOverParNonBackstream"] =DataReducerVariable(["PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$\beta_{\perp,\mathrm{th}} \beta_{\parallel,\mathrm{th}}^{-1}$", latexunits=r"")

datareducers["beta"] =                   DataReducerVariable(["Pressure", "B"], beta ,"", 1, latex=r"$\beta$", latexunits=r"")
datareducers["betaParallel"] =           DataReducerVariable(["PParallel", "B"], beta ,"", 1, latex=r"$\beta_\parallel$", latexunits=r"")
datareducers["betaPerpendicular"] =      DataReducerVariable(["PPerpendicular", "B"], beta ,"", 1, latex=r"$\beta_\perp$", latexunits=r"")

datareducers["Rmirror"] =                DataReducerVariable(["PTensor", "B"], rMirror, "", 1, latex=r"$R_\mathrm{m}$")
datareducers["Dng"] =                    DataReducerVariable(["PTensor", "PParallel", "PPerpendicular", "B"], Dng, "", 1, latex=r"$\mathrm{Dng}$") # I think this has vector length 1?
datareducers["vBeam"] =                  DataReducerVariable(["VBackstream", "VNonBackstream"], v_beam, "m/s", 3, latex=r"$V_\mathrm{st}-V$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["vBeamRatio"] =             DataReducerVariable(["VBackstream", "VNonBackstream"], v_beam_ratio, "", 1, latex=r"$V_\mathrm{st} V^{-1}$", latexunits=r"")
datareducers["vThermal"] =               DataReducerVariable(["Temperature"], v_thermal, "m/s", 1, latex=r"$v_\mathrm{th}$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
datareducers["Bz_linedipole_avg"] =      DataReducerVariable(["X", "Y", "Z", "DX", "DY", "DZ"], Bz_linedipole_avg, "T", 1, latex=r"$\langle B_{z,\mathrm{ld}}\rangle$")
datareducers["Bz_linedipole_diff"] =     DataReducerVariable(["B", "Bz_linedipole_avg"], Bz_linedipole_diff, "", 1, latex=r"$\Delta B_{z,\mathrm{ld}}$")

#reducers with useVspace
datareducers["gyrophase_relstddev"] =    DataReducerVariable(["v", "B"], gyrophase_relstddev, "", 1, useVspace=True) # I think this has vector length 1?





#multipopdatareducers
multipopdatareducers = {}
multipopdatareducers["pop/Pdyn"] =            DataReducerVariable(["pop/V", "pop/rho"], Pdyn, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP}}$",latexunits=r"Pa")
multipopdatareducers["pop/Pdynx"] =            DataReducerVariable(["pop/V", "pop/rho"], Pdynx, "Pa", 1, latex=r"$P_\mathrm{dyn,\mathrm{REPLACEPOP},x}$",latexunits=r"Pa")

multipopdatareducers["pop/VParallel"] =              DataReducerVariable(["pop/V", "B"], VParallel, "m/s", 1, latex=r"$V_{\parallel,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")
multipopdatareducers["pop/VPerpendicular"] =         DataReducerVariable(["pop/V", "B"], VPerpendicular, "m/s", 1, latex=r"$V_{\perp,\mathrm{REPLACEPOP}}$",latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

multipopdatareducers["pop/Pressure"] =               DataReducerVariable(["pop/PTensorDiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensor"] =                DataReducerVariable(["pop/PTensorDiagonal", "pop/PTensorOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensorRotated"] =         DataReducerVariable(["pop/PTensor", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"Pa")
multipopdatareducers["pop/PParallel"] =              DataReducerVariable(["pop/PTensorRotated"], PParallel, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpendicular"] =         DataReducerVariable(["pop/PTensorRotated"], PPerpendicular, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpOverPar"] =           DataReducerVariable(["pop/PTensorRotated"], PPerpOverPar, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP}} P_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/aGyrotropy"] =             DataReducerVariable(["pop/PTensorRotated"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP}$", latexunits=r"")

multipopdatareducers["pop/PBackstream"] =            DataReducerVariable(["pop/PTensorBackstreamDiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensorBackstream"] =      DataReducerVariable(["pop/PTensorBackstreamDiagonal", "pop/PTensorBackstreamOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensorRotatedBackstream"]=DataReducerVariable(["pop/PTensorBackstream", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"Pa")
multipopdatareducers["pop/PParallelBackstream"] =    DataReducerVariable(["pop/PTensorRotatedBackstream"], PParallel, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpendicularBackstream"]=DataReducerVariable(["pop/PTensorRotatedBackstream"], PPerpendicular, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpOverParBackstream"] = DataReducerVariable(["pop/PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,st}} P_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/aGyrotropyBackstream"] =   DataReducerVariable(["pop/PTensorRotatedBackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,st}$", latexunits=r"")

multipopdatareducers["pop/PNonBackstream"] =              DataReducerVariable(["pop/PTensorNonBackstreamDiagonal"], Pressure, "Pa", 1, latex=r"$P_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensorNonBackstream"] =        DataReducerVariable(["pop/PTensorNonBackstreamDiagonal", "pop/PTensorNonBackstreamOffDiagonal"], PTensor, "Pa", 9, latex=r"$\mathcal{P}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/PTensorRotatedNonBackstream"] = DataReducerVariable(["pop/PTensorNonBackstream", "B"], PTensorRotated, "Pa", 9, latex=r"$\mathcal{P}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"Pa")
multipopdatareducers["pop/PParallelNonBackstream"] =      DataReducerVariable(["pop/PTensorRotatedNonBackstream"], PParallel, "Pa", 1, latex=r"$P_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpendicularNonBackstream"] = DataReducerVariable(["pop/PTensorRotatedNonBackstream"], PPerpendicular, "Pa", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"Pa")
multipopdatareducers["pop/PPerpOverParNonBackstream"] =   DataReducerVariable(["pop/PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$P_{\perp,\mathrm{REPLACEPOP,th}} P_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/aGyrotropyNonBackstream"] =     DataReducerVariable(["pop/PTensorRotatedNonBackstream"], aGyrotropy, "", 1, latex=r"$Q_\mathrm{ag,REPLACEPOP,th}$", latexunits=r"")

multipopdatareducers["pop/Temperature"] =            DataReducerVariable(["pop/Pressure", "pop/rho"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/TTensor"] =                DataReducerVariable(["pop/PTensor", "pop/rho"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/TTensorRotated"] =         DataReducerVariable(["pop/PTensorRotated", "pop/rho"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP}$", latexunits=r"K")
multipopdatareducers["pop/TParallel"] =              DataReducerVariable(["pop/PParallel", "pop/rho"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"K")
multipopdatareducers["pop/TPerpendicular"] =         DataReducerVariable(["pop/PPerpendicular", "pop/rho"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"K")

multipopdatareducers["pop/TBackstream"] =            DataReducerVariable(["pop/PBackstream", "pop/RhoBackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/TTensorBackstream"] =      DataReducerVariable(["pop/PTensorBackstream", "pop/RhoBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/TTensorRotatedBackstream"]=DataReducerVariable(["pop/PTensorRotatedBackstream", "pop/RhoBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,st}$", latexunits=r"K")
multipopdatareducers["pop/TParallelBackstream"] =    DataReducerVariable(["pop/PParallelBackstream", "pop/RhoBackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")
multipopdatareducers["pop/TPerpendicularBackstream"]=DataReducerVariable(["pop/PPerpendicularBackstream", "pop/RhoBackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}}$", latexunits=r"K")

multipopdatareducers["pop/TNonBackstream"] =              DataReducerVariable(["pop/PNonBackstream", "pop/RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/TTensorNonBackstream"] =        DataReducerVariable(["pop/PTensorNonBackstream", "pop/RhoNonBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/TTensorRotatedNonBackstream"] = DataReducerVariable(["pop/PTensorRotatedNonBackstream", "pop/RhoNonBackstream"], Temperature, "K", 9, latex=r"$\mathcal{T}^\mathrm{R}_\mathrm{REPLACEPOP,th}$", latexunits=r"K")
multipopdatareducers["pop/TParallelNonBackstream"] =      DataReducerVariable(["pop/PParallelNonBackstream", "pop/RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_{\parallel,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")
multipopdatareducers["pop/TPerpendicularNonBackstream"] = DataReducerVariable(["pop/PPerpendicularNonBackstream", "pop/RhoNonBackstream"], Temperature, "K", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}}$", latexunits=r"K")

 # These ratios are identical to the pressure ratios
multipopdatareducers["pop/TPerpOverPar"] =                 DataReducerVariable(["pop/PTensorRotated"], PPerpOverPar, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} T_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/TPerpOverParBackstream"] =       DataReducerVariable(["pop/PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,st}} T_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/TPerpOverParNonBackstream"] =    DataReducerVariable(["pop/PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP,th}} T_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/betaPerpOverPar"] =              DataReducerVariable(["pop/PTensorRotated"], PPerpOverPar, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}} \beta_{\parallel,\mathrm{REPLACEPOP}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/betaPerpOverParBackstream"] =    DataReducerVariable(["pop/PTensorRotatedBackstream"], PPerpOverPar, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,st}} \beta_{\parallel,\mathrm{REPLACEPOP,st}}^{-1}$", latexunits=r"")
multipopdatareducers["pop/betaPerpOverParNonBackstream"] = DataReducerVariable(["pop/PTensorRotatedNonBackstream"], PPerpOverPar, "", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP,th}} \beta_{\parallel,\mathrm{REPLACEPOP,th}}^{-1}$", latexunits=r"")

multipopdatareducers["pop/vThermal"] =               DataReducerVariable(["Temperature"], v_thermal, "m/s", 1, latex=r"$v_\mathrm{th,REPLACEPOP}$", latexunits=r"$\mathrm{m}\,\mathrm{s}^{-1}$")

datareducers["pop/firstadiabatic"] =    DataReducerVariable(["pop/TPerpendicular","B"], firstadiabatic, "K/T", 1, latex=r"$T_{\perp,\mathrm{REPLACEPOP}} B^{-1}$",latexunits=r"$\mathrm{K}\,\mathrm{T}^{-1}$")

# multipopdatareducers["pop/TxRotated"] =              DataReducerVariable(["pop/TTensorRotated"], TxRotated, "K")
# multipopdatareducers["pop/TyRotated"] =              DataReducerVariable(["pop/TTensorRotated"], TyRotated, "K")

# Do these betas make sense per-population?
multipopdatareducers["pop/beta"] =                   DataReducerVariable(["pop/Pressure", "B"], beta ,"", 1, latex=r"$\beta_\mathrm{REPLACEPOP}$", latexunits=r"")
multipopdatareducers["pop/betaParallel"] =           DataReducerVariable(["pop/PParallel", "B"], beta ,"", 1, latex=r"$\beta_{\parallel,\mathrm{REPLACEPOP}}$", latexunits=r"")
multipopdatareducers["pop/betaPerpendicular"] =      DataReducerVariable(["pop/PPerpendicular", "B"], beta ,"", 1, latex=r"$\beta_{\perp,\mathrm{REPLACEPOP}}$", latexunits=r"")

multipopdatareducers["pop/Rmirror"] =                DataReducerVariable(["pop/PTensor", "B"], rMirror, "", 1, latex=r"$R_\mathrm{m,REPLACEPOP}$")
multipopdatareducers["pop/Dng"] =                    DataReducerVariable(["pop/PTensor", "pop/PParallel", "pop/PPerpendicular", "B"], Dng, "", 1, latex=r"$\mathrm{Dng}_\mathrm{REPLACEPOP}$")


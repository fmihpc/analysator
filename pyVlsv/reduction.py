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
import filemanagement
# Input paths:
fullPath = filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__))
# Input current folder's path
filemanagement.sys.path.insert(0, fullPath)
# Input folder paths
filemanagement.sys.path.insert(0, fullPath + "/" + "pyCalculations")
from reducer import DataReducerVariable
from rotation import rotateTensorToVector
from gyrophaseangle import gyrophase_angles
import sys



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
   return np.sum(np.asarray(variable)**2,axis=-1)**(0.5)

def sumv( variable ):
   if np.ndim(variable) > 3:
      print('Error: Number of dimensions is too large')
      return
   else:
      res = []
      if np.size(variable[0][0]) == 1:
         for i in range(len(variable[0])):
            summation = 0
            for j in range(len(variable)):
               summation += variable[j][i]
            res.append(summation)
      else: 
         for i in range(len(variable[0])):
            summation = np.zeros_like(variable[0][0])
            for j in range(len(variable)):
               summation += np.array(variable[j][i])
            if i==0:
               res = summation
            else:
               res = np.vstack((res, summation)) 
      return np.array(res)

def condition_matrix_array( condition, matrices ):
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

#do nothing
def pass_op( variable ):
   return variable


def v( variables ):
   ''' Data reducer function for getting velocity from rho and rho_v
       variables[0] = rho_v and variables[1] = rho
   '''
   epsilon = sys.float_info.epsilon
   rho_v = variables[0]
   rho = variables[1] + epsilon
   if np.ndim(rho) == 0:
      return np.divide(rho_v,rho)
   else:
      rho = np.reshape(rho,(len(rho),1))
      return np.divide(rho_v,rho)

def vms( variables ):
   ''' Data reducer function for getting magnetosonic velocity
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   mu_0 = 1.25663706144e-6
   P = variables[0]
   rho_m = (variables[1] + epsilon)*mp
   B = variables[2]
   if np.ndim(B) == 1:
      Btot = np.sqrt( np.sum( np.square(B) ) )
   else:
      Btot = np.sqrt( np.sum(np.square(B),axis=1) )
   vs = np.sqrt( np.divide( P*5.0/3.0, rho_m ) )
   vA = np.divide( Btot,np.sqrt( mu_0*rho_m ) )
   vms = np.sqrt( np.square(vs) + np.square(vA) )
   return vms

def vs( variables ):
   ''' Data reducer function for getting the sound speed
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   P = variables[0]
   rho_m = (variables[1] + epsilon)*mp
   vs = np.sqrt( np.divide( P*5.0/3.0, rho_m ) )
   return vs

def va( variables ):
   ''' Data reducer function for getting the Alfven velocity
   '''
   epsilon = sys.float_info.epsilon
   mp = 1.672622e-27
   mu_0 = 1.25663706144e-6
   rho_m = (variables[0] + epsilon)*mp
   B = variables[1]
   if np.ndim(B) == 1:
      Btot = np.sqrt( np.sum( np.square(B) ) )
   else:
      Btot = np.sqrt( np.sum(np.square(B),axis=1) )
   vA = np.divide( Btot,np.sqrt( mu_0*rho_m ) )
   return vA

def VParallel( variables ):
   bulkv = variables[0]
   bvector = variables[1]
   if( np.ndim(bulkv)==1 ):
      bnorm = bvector / np.linalg.norm(bvector)
      return (bulkv*bnorm).sum()
   else:
      bnorm = bvector / np.linalg.norm(bvector, axis=-1)[:,np.newaxis]
      return (bulkv*bnorm).sum(-1)

def VPerpendicular( variables ):
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
   PTensorDiagonal = variables[0]
   PTensorOffDiagonal = variables[1]
   if(np.ndim(PTensorDiagonal)==1 ):
      PTensorDiagonal = PTensorDiagonal.reshape(1,3)[0]
      PTensorOffDiagonal = PTensorOffDiagonal.reshape(1,3)[0]
      return np.array([[PTensorDiagonal[0], PTensorOffDiagonal[2], PTensorOffDiagonal[1]],[PTensorOffDiagonal[2], PTensorDiagonal[1], PTensorOffDiagonal[0]],[PTensorOffDiagonal[1], PTensorOffDiagonal[0], PTensorDiagonal[2]]])
   else:
      result = []
      for index in np.arange(len(PTensorDiagonal)):
         result.append(np.array([[PTensorDiagonal[index,0], PTensorOffDiagonal[index,2], PTensorOffDiagonal[index,1]],[PTensorOffDiagonal[index,2], PTensorDiagonal[index,1], PTensorOffDiagonal[index,0]],[PTensorOffDiagonal[index,1], PTensorOffDiagonal[index,0], PTensorDiagonal[index,2]]]))
      return np.asarray(result)

def PTensorRotated( variables ):
   PTensor = variables[0]
   B = variables[1]
   if( np.ndim(B)==1 ):
      B = B.reshape(1,3)[0]
      PTensor = PTensor.reshape(3,3)
      return rotateTensorToVector(PTensor, B)
   else:
      result = []
      for index in np.arange(len(PTensor)):
         result.append(rotateTensorToVector(PTensor[index], B[index]))
      return np.asarray(result)

def Pressure( variables ):
   PTensorDiagonal = variables[0]
   if(np.ndim(PTensorDiagonal)==1 ):
      return 1.0/3.0 * np.sum(PTensorDiagonal)
   return 1.0/3.0 * (PTensorDiagonal[:,0]+PTensorDiagonal[:,1]+PTensorDiagonal[:,2])

def TTensor( variables ):
   PTensor = variables[0]
   rho = variables[1]
   if(rho.size == 1):
      return np.divide(PTensor, rho * 1.38065e-23)
   else:
      result = []
      for index in np.arange(len(rho)):
         result.append(np.divide(PTensor[index], (rho[index] + 1.0) * 1.38065e-23))
      return np.asarray(result)

def TTensorRotated( variables):
   TTensor = variables[0]
   B = variables[1]
   if( np.ndim(B)==1 ):
      B = B.reshape(1,3)[0]
      TTensor = TTensor.reshape(3,3)
      return rotateTensorToVector(TTensor, B)
   else:
      result = []
      for index in np.arange(len(TTensor)):
         result.append(rotateTensorToVector(TTensor[index], B[index]))
      return np.asarray(result)

def Temperature( variables ):
   Pressure = variables[0]
   rho = variables[1]
   return np.divide(Pressure, (rho + 1.0) * 1.38065e-23)

def TPerpendicular( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return 0.5*(TTensorRotated[0,0] + TTensorRotated[1,1])
   else:
      return 0.5*(TTensorRotated[:,0,0] + TTensorRotated[:,1,1])

def TParallel( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return TTensorRotated[2,2]
   else:
      return TTensorRotated[:,2,2]

def TPerpendicularBackstream( variables ):
   TTensorRotatedBackstream = variables[0]
   if( np.ndim(TTensorRotatedBackstream)==2 ):
      return 0.5*(TTensorRotatedBackstream[0,0] + TTensorRotatedBackstream[1,1])
   else:
      return 0.5*(TTensorRotatedBackstream[:,0,0] + TTensorRotatedBackstream[:,1,1])

def TParallelBackstream( variables ):
   TTensorRotatedBackstream = variables[0]
   if( np.ndim(TTensorRotatedBackstream)==2 ):
      return TTensorRotatedBackstream[2,2]
   else:
      return TTensorRotatedBackstream[:,2,2]


def TxRotated( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return TTensorRotated[0,0]
   else:
      return TTensorRotated[:,0,0]

def TyRotated( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return TTensorRotated[1,1]
   else:
      return TTensorRotated[:,1,1]

def PPerpendicular( variables ):
   PTensor = variables[0]
   PParallel = variables[1]
   if( np.ndim(PTensorRotated)==2 ):
      return 0.5*(PTensor.trace() - PParallel)
   else:
      return 0.5*(PTensor.transpose().trace() - PParallel)

def PParallel( variables ):
   PTensor = variables[0]
   B = variables[1]
   if( np.ndim(PTensor)==2 ):
      B_normalized = np.divide(B, np.sqrt(np.sum(B[:]**2)))
      return B_normalized.dot(PTensor.dot(B_normalized))
   else:
      B_normalized = np.divide(B, np.sqrt(np.sum(B[:]**2, axis=1))[:,None])
      return np.asarray([B_normalized[i].dot(PTensor[i].dot(B_normalized[i])) for i in np.arange(len(PTensor))])

def TPerpOverPar( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return 0.5*np.divide(TTensorRotated[0,0] + TTensorRotated[1,1], TTensorRotated[2,2])
   else:
      return 0.5*np.divide(TTensorRotated[:,0,0] + TTensorRotated[:,1,1], TTensorRotated[:,2,2])

def PPerpOverPar( variables ):
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return 0.5*np.divide(PTensorRotated[0,0] + PTensorRotated[1,1], PTensorRotated[2,2])
   else:
      return 0.5*np.divide(PTensorRotated[:,0,0] + PTensorRotated[:,1,1], PTensorRotated[:,2,2])

def beta( variables ):
   Pressure = variables[0]
   B_magnitude = magnitude(variables[1])
   return 2.0 * 1.25663706144e-6 * Pressure / (B_magnitude*B_magnitude)

def betaParallel( variables ):
   Pressure = variables[0]
   B_magnitude = magnitude(variables[1])
   return 2.0 * 1.25663706144e-6 * Pressure / (B_magnitude*B_magnitude)

def betaPerpendicular( variables ):
   Pressure = variables[0]
   B_magnitude = magnitude(variables[1])
   return 2.0 * 1.25663706144e-6 * Pressure / (B_magnitude*B_magnitude)

def betaPerpOverPar( variables ):
   betaPerp = variables[0]
   betaPar = variables[1]
   return np.divide(betaPerp, betaPar)

def rMirror( variables ):
   TPerpOverPar = variables[0]
   betaPerpendicular = variables[1]
   return betaPerpendicular * (TPerpOverPar - 1)

def v_beam( variables ):
   rhoVBackstream = variables[0]
   rhoBackstream = variables[1]
   rhoVNonBackstream = variables[2]
   rhoNonBackstream = variables[3]
   # get velocity of both populations:
   vBackstream = v( [rhoVBackstream, rhoBackstream] )
   vNonBackstream = v( [rhoVNonBackstream, rhoNonBackstream] )
   vBeam = vBackstream - vNonBackstream
   return vBeam

def v_beam_ratio( variables ):
   rhoVBackstream = variables[0]
   rhoBackstream = variables[1]
   rhoVNonBackstream = variables[2]
   rhoNonBackstream = variables[3]
   # get velocity of both populations:
   vBackstream = magnitude(v( [rhoVBackstream, rhoBackstream] ))
   vNonBackstream = magnitude(v( [rhoVNonBackstream, rhoNonBackstream] ))
   vBeamRatio = vBackstream / vNonBackstream
   return vBeamRatio

def v_thermal( variables ):
   temperatureBackstream = variables[0]
   k = 1.38065e-23
   ion_mass = 1.672622e-27
   vThermal = np.sqrt(k*temperatureBackstream/ion_mass)
   return vThermal

def v_thermal_vector( variables ):
   TTensorRotatedBackstream = variables[0]
   # Take diagonal components of each:
   condition = [
               [True,  False, False],
               [False, True,  False],
               [False, False, True]
               ]
   diagonal_elements = condition_matrix_array( condition, TTensorRotatedBackstream )
   # We have the diagonal elements of the Temperature tensor now
   # Do the math..
   k = 1.38065e-23
   ion_mass = 1.672622e-27
   v_thermal_vector = np.sqrt( diagonal_elements*k/ion_mass )
   return v_thermal_vector

def Bz_linedipole_avg( variables ):
   x = variables[0]
   y = variables[1]
   z = variables[2]
   dx = variables[3]
   dy = variables[4]
   dz = variables[5]
   return -126.2e6*((dx+x)/(dx*(z**2+(dx+x)**2)) - x/(dx*(z**2+x**2)))

def Bz_linedipole_diff( variables ):
   Bb = variables[0]
   Bzldp = variables[1]
   print Bzldp.shape
   return np.divide(np.abs(Bb[:,2] - Bzldp), magnitude(Bb))

def gyrophase_relstddev( variables, velocity_cell_data, velocity_coordinates ):
   bulk_velocity = variables[0]
   B = variables[1]
   B_unit = B / np.linalg.norm(B)
   
   gyrophase_data = gyrophase_angles(bulk_velocity, B_unit, velocity_cell_data, velocity_coordinates)
   histo = pl.hist(gyrophase_data[0].data, weights=gyrophase_data[1].data, bins=36, range=[-180.0,180.0], log=False, normed=1)
   return np.std(histo[0])/np.mean(histo[0])

def Dng( variables ):
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

#datareducers with more complex, case dependent structure.
datareducers = {}
datareducers["v"] =                      DataReducerVariable(["rho_v", "rho"], v, "m/s")
datareducers["vms"] =                    DataReducerVariable(["Pressure", "rho", "B"], vms, "m/s")
datareducers["vs"] =                     DataReducerVariable(["Pressure", "rho"], vs, "m/s")
datareducers["va"] =                     DataReducerVariable(["rho", "B"], va, "m/s")
datareducers["VParallel"] =              DataReducerVariable(["v", "B"], VParallel, "m/s")
datareducers["VPerpendicular"] =         DataReducerVariable(["v", "B"], VPerpendicular, "m/s")
datareducers["PTensor"] =                DataReducerVariable(["PTensorDiagonal", "PTensorOffDiagonal"], PTensor, "Pa")
datareducers["PTensorBackstream"] =      DataReducerVariable(["PTensorBackstreamDiagonal", "PTensorBackstreamOffDiagonal"], PTensor, "Pa")
datareducers["PTensorRotated"] =         DataReducerVariable(["PTensor", "B"], PTensorRotated, "Pa")
datareducers["Pressure"] =               DataReducerVariable(["PTensorDiagonal"], Pressure, "Pa")
datareducers["PBackstream"] =            DataReducerVariable(["PTensorBackstreamDiagonal"], Pressure, "Pa")
datareducers["TTensor"] =                DataReducerVariable(["PTensor", "rho"], TTensor, "K")
datareducers["TTensorRotated"] =         DataReducerVariable(["TTensor", "B"], TTensorRotated, "K")
datareducers["TTensorBackstream"] =      DataReducerVariable(["PTensorBackstream", "RhoBackstream"], TTensor, "K")
datareducers["TTensorRotatedBackstream"]=DataReducerVariable(["TTensorBackstream", "B"], TTensorRotated, "K")
datareducers["Temperature"] =            DataReducerVariable(["Pressure", "rho"], Temperature, "K")
datareducers["TBackstream"] =            DataReducerVariable(["PBackstream", "RhoBackstream"], Temperature, "K")
datareducers["TParallel"] =              DataReducerVariable(["TTensorRotated"], TParallel, "K")
datareducers["TParallelBackstream"] =    DataReducerVariable(["TTensorRotatedBackstream"], TParallelBackstream, "K")
datareducers["TPerpendicularBackstream"]=DataReducerVariable(["TTensorRotatedBackstream"], TPerpendicularBackstream, "K")
datareducers["TPerpendicular"] =         DataReducerVariable(["TTensorRotated"], TPerpendicular, "K")
datareducers["TxRotated"] =              DataReducerVariable(["TTensorRotated"], TxRotated, "K")
datareducers["TyRotated"] =              DataReducerVariable(["TTensorRotated"], TyRotated, "K")
datareducers["TPerpOverPar"] =           DataReducerVariable(["TTensorRotated"], TPerpOverPar, "K")
datareducers["PParallel"] =              DataReducerVariable(["PTensor", "B"], PParallel, "Pa")
datareducers["PPerpendicular"] =         DataReducerVariable(["PTensor", "PParallel"], PPerpendicular, "Pa")
datareducers["PPerpOverPar"] =           DataReducerVariable(["PTensorRotated"], PPerpOverPar, "K")
datareducers["beta"] =                   DataReducerVariable(["Pressure", "B"], beta ,"")
datareducers["betaParallel"] =           DataReducerVariable(["PParallel", "B"], beta ,"")
datareducers["betaPerpendicular"] =      DataReducerVariable(["PPerpendicular", "B"], beta ,"")
datareducers["betaPerpOverPar"] =        DataReducerVariable(["betaPerpendicular", "betaParallel"], betaPerpOverPar, "")
datareducers["Rmirror"] =                DataReducerVariable(["TPerpOverPar", "betaPerpendicular"], rMirror, "")
datareducers["Dng"] =                    DataReducerVariable(["PTensor", "PParallel", "PPerpendicular", "B"], Dng, "")

datareducers["vBeam"] =                  DataReducerVariable(["RhoVBackstream", "RhoBackstream", "RhoVNonBackstream", "RhoNonBackstream"], v_beam, "m/s")
datareducers["vBeamRatio"] =             DataReducerVariable(["RhoVBackstream", "RhoBackstream", "RhoVNonBackstream", "RhoNonBackstream"], v_beam_ratio, "")
datareducers["vThermal"] =               DataReducerVariable(["TBackstream"], v_thermal, "m/s")
datareducers["vThermalVector"] =         DataReducerVariable(["TTensorRotatedBackstream"], v_thermal_vector, "m/s")
datareducers["Bz_linedipole_avg"] =      DataReducerVariable(["X", "Y", "Z", "DX", "DY", "DZ"], Bz_linedipole_avg, "T")
datareducers["Bz_linedipole_diff"] =     DataReducerVariable(["B", "Bz_linedipole_avg"], Bz_linedipole_diff, "")

#reducers with useVspace
datareducers["gyrophase_relstddev"] =    DataReducerVariable(["v", "B"], gyrophase_relstddev, "", useVspace=True)


#list of operators. The user can apply these to any variable,
#including more general datareducers. Can only be used to reduce one
#variable at a time
data_operators = {}
data_operators["pass"] = pass_op
data_operators["magnitude"] = magnitude
data_operators["sum"] = sumv
data_operators["x"] = x_component
data_operators["y"] = y_component
data_operators["z"] = z_component



#multipopdatareducers
multipopdatareducers = {}
#multipopdatareducers["pop/vms"] =                    DataReducerVariable(["pop/Pressure", "pop/rho", "B"], vms, "m/s")
#multipopdatareducers["pop/vs"] =                     DataReducerVariable(["pop/Pressure", "pop/rho"], vs, "m/s")
#multipopdatareducers["pop/va"] =                     DataReducerVariable(["pop/rho", "B"], va, "m/s")
multipopdatareducers["pop/VParallel"] =              DataReducerVariable(["pop/V", "B"], VParallel, "m/s")
multipopdatareducers["pop/VPerpendicular"] =         DataReducerVariable(["pop/V", "B"], VPerpendicular, "m/s")
multipopdatareducers["pop/PTensor"] =                DataReducerVariable(["pop/PTensorDiagonal", "pop/PTensorOffDiagonal"], PTensor, "Pa")
multipopdatareducers["pop/PTensorBackstream"] =      DataReducerVariable(["pop/PTensorBackstreamDiagonal", "pop/PTensorBackstreamOffDiagonal"], PTensor, "Pa")
multipopdatareducers["pop/PTensorRotated"] =         DataReducerVariable(["pop/PTensor", "B"], PTensorRotated, "Pa")
multipopdatareducers["pop/Pressure"] =               DataReducerVariable(["pop/PTensorDiagonal"], Pressure, "Pa")
multipopdatareducers["pop/PBackstream"] =            DataReducerVariable(["pop/PTensorBackstreamDiagonal"], Pressure, "Pa")
multipopdatareducers["pop/TTensor"] =                DataReducerVariable(["pop/PTensor", "pop/rho"], TTensor, "K")
multipopdatareducers["pop/TTensorRotated"] =         DataReducerVariable(["pop/TTensor", "B"], TTensorRotated, "K")
multipopdatareducers["pop/TTensorBackstream"] =      DataReducerVariable(["pop/PTensorBackstream", "pop/RhoBackstream"], TTensor, "K")
multipopdatareducers["pop/TTensorRotatedBackstream"]=DataReducerVariable(["pop/TTensorBackstream", "B"], TTensorRotated, "K")
multipopdatareducers["pop/Temperature"] =            DataReducerVariable(["pop/Pressure", "pop/rho"], Temperature, "K")
multipopdatareducers["pop/TBackstream"] =            DataReducerVariable(["pop/PBackstream", "pop/RhoBackstream"], Temperature, "K")
multipopdatareducers["pop/TParallel"] =              DataReducerVariable(["pop/TTensorRotated"], TParallel, "K")
multipopdatareducers["pop/TPerpendicular"] =         DataReducerVariable(["pop/TTensorRotated"], TPerpendicular, "K")
multipopdatareducers["pop/TxRotated"] =              DataReducerVariable(["pop/TTensorRotated"], TxRotated, "K")
multipopdatareducers["pop/TyRotated"] =              DataReducerVariable(["pop/TTensorRotated"], TyRotated, "K")
multipopdatareducers["pop/TPerpOverPar"] =           DataReducerVariable(["pop/TTensorRotated"], TPerpOverPar, "K")
multipopdatareducers["pop/PParallel"] =              DataReducerVariable(["pop/PTensor", "B"], PParallel, "Pa")      
multipopdatareducers["pop/PPerpendicular"] =         DataReducerVariable(["pop/PTensor", "pop/PParallel"], PPerpendicular, "Pa")
multipopdatareducers["pop/PPerpOverPar"] =           DataReducerVariable(["pop/PTensorRotated"], PPerpOverPar, "K")
multipopdatareducers["pop/beta"] =                   DataReducerVariable(["pop/Pressure", "B"], beta ,"")
multipopdatareducers["pop/betaParallel"] =           DataReducerVariable(["pop/PParallel", "B"], beta ,"")
multipopdatareducers["pop/betaPerpendicular"] =      DataReducerVariable(["pop/PPerpendicular", "B"], beta ,"")
multipopdatareducers["pop/betaPerpOverPar"] =        DataReducerVariable(["pop/betaPerpendicular", "pop/betaParallel"], betaPerpOverPar, "")
multipopdatareducers["pop/Rmirror"] =                DataReducerVariable(["pop/TPerpOverPar", "pop/betaPerpendicular"], rMirror, "")
multipopdatareducers["pop/Dng"] =                    DataReducerVariable(["pop/PTensor", "pop/PParallel", "pop/PPerpendicular", "B"], Dng, "")
#multipopdatareducers["pop/vBeam"] =                  DataReducerVariable(["pop/RhoVBackstream", "pop/RhoBackstream", "pop/RhoVNonBackstream", "pop/RhoNonBackstream"], v_beam, "m/s")
#multipopdatareducers["pop/vBeamRatio"] =             DataReducerVariable(["pop/RhoVBackstream", "pop/RhoBackstream", "pop/RhoVNonBackstream", "pop/RhoNonBackstream"], v_beam_ratio, "")
multipopdatareducers["pop/vThermal"] =               DataReducerVariable(["pop/TBackstream"], v_thermal, "m/s")
multipopdatareducers["pop/vThermalVector"] =         DataReducerVariable(["pop/TTensorRotatedBackstream"], v_thermal_vector, "m/s")


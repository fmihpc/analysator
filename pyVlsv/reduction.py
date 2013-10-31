''' A file for doing data reduction on variables
'''
import numpy as np
import filemanagement
# Input paths:
fullPath = filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__))
# Input current folder's path
filemanagement.sys.path.insert(0, fullPath)
# Input folder paths
filemanagement.sys.path.insert(0, fullPath + "/" + "pyCalculations")
from reducer import DataReducerVariable
from rotation import rotateTensorToVector

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

#do nothing
def pass_op( variable ):
   return variable


def v( variables ):
   ''' Data reducer function for getting velocity from rho and rho_v
       variables[0] = rho_v and variables[1] = rho
   '''
   rho_v = variables[0]
   rho = np.reshape(variables[1],(len(variables[1]),1))
   return np.divide(rho_v,rho)

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
      PTensorDiagonal = PTensorDiagonal.reshape(1,3)
   return PTensorDiagonal[:,0]+PTensorDiagonal[:,1]+PTensorDiagonal[:,2]

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

def PPerpendicular( variables ):
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return 0.5*(PTensorRotated[0,0] + PTensorRotated[1,1])
   else:
      return 0.5*(PTensorRotated[:,0,0] + PTensorRotated[:,1,1])

def PParallel( variables ):
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return PTensorRotated[2,2]
   else:
      return PTensorRotated[:,2,2]

def TPerpOverPar( variables ):
   TTensorRotated = variables[0]
   if( np.ndim(TTensorRotated)==2 ):
      return 0.5*(TTensorRotated[0,0] + TTensorRotated[1,1]) / TTensorRotated[2,2]
   else:
      return 0.5*(TTensorRotated[:,0,0] + TTensorRotated[:,1,1]) / TTensorRotated[:,2,2]

def PPerpOverPar( variables ):
   PTensorRotated = variables[0]
   if( np.ndim(PTensorRotated)==2 ):
      return 0.5*(PTensorRotated[0,0] + PTensorRotated[1,1]) / PTensorRotated[2,2]
   else:
      return 0.5*(PTensorRotated[:,0,0] + PTensorRotated[:,1,1]) / PTensorRotated[:,2,2]

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

def rMirror( variables ):
   TPerpOverPar = variables[0]
   betaPerpendicular = variables[1]
   return betaPerpendicular * (TPerpOverPar - 1)

#datareducers with more complex, case dependent structure. 
datareducers = {}
datareducers["v"] =                 DataReducerVariable(["rho_v", "rho"], v, "m/s")
datareducers["PTensor"] =           DataReducerVariable(["PTensorDiagonal", "PTensorOffDiagonal"], PTensor, "Pa")
datareducers["PTensorRotated"] =    DataReducerVariable(["PTensor", "B"], PTensorRotated, "Pa")
datareducers["Pressure"] =          DataReducerVariable(["PTensorDiagonal"], Pressure, "Pa")
datareducers["TTensor"] =           DataReducerVariable(["PTensor", "rho"], TTensor, "K")
datareducers["TTensorRotated"] =    DataReducerVariable(["TTensor", "B"], TTensorRotated, "K")
datareducers["Temperature"] =       DataReducerVariable(["Pressure", "rho"], Temperature, "K")
datareducers["TParallel"] =         DataReducerVariable(["TTensorRotated"], TParallel, "K")
datareducers["TPerpendicular"] =    DataReducerVariable(["TTensorRotated"], TPerpendicular, "K")
datareducers["TPerpOverPar"] =      DataReducerVariable(["TTensorRotated"], TPerpOverPar, "K")
datareducers["PParallel"] =         DataReducerVariable(["PTensorRotated"], PParallel, "Pa")
datareducers["PPerpendicular"] =    DataReducerVariable(["PTensorRotated"], PPerpendicular, "Pa")
datareducers["PPerpOverPar"] =      DataReducerVariable(["PTensorRotated"], PPerpOverPar, "K")
datareducers["beta"] =              DataReducerVariable(["Pressure", "B"], beta ,"ND")
datareducers["betaParallel"] =      DataReducerVariable(["PParallel", "B"], beta ,"ND")
datareducers["betaPerpendicular"] = DataReducerVariable(["PPerpendicular", "B"], beta ,"ND")
datareducers["Rmirror"]           = DataReducerVariable(["TPerpOverPar", "betaPerpendicular"], rMirror, "ND")


#list of operators. The user can apply these to any variable,
#including more general datareducers. Can only be used to reduce one
#variable at a time
datareduction_operators = {}
datareduction_operators["pass"] = pass_op
datareduction_operators["magnitude"] = magnitude
datareduction_operators["x"] = x_component
datareduction_operators["y"] = y_component
datareduction_operators["z"] = z_component










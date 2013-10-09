''' A file for doing data reduction on variables
'''
import numpy as np
from reducer import DataReducerVariable

def B_x( B ):
   if len(B.shape) == 1:
      return B[0]
   return B[:,0]

def B_y( B ):
   if len(B.shape) == 1:
      return B[1]
   return B[:,1]

def B_z( B ):
   if len(B.shape) == 1:
      return B[2]
   return B[:,2]

def B_mag( B ):
   return np.sum(np.abs(B)**2,axis=-1)**(1./2)

def E_x( E ):
   if len(B.shape) == 1:
      return E[0]
   return E[:,0]

def E_y( E ):
   if len(E.shape) == 1:
      return E[1]
   return E[:,1]

def E_z( E ):
   if len(E.shape) == 1:
      return E[2]
   return E[:,2]

def E_mag( E ):
   return np.sum(np.abs(E)**2,axis=-1)**(1./2)

def v( variables ):
   ''' Data reducer function for getting velocity from rho and rho_v
       variables[0] = rho_v and variables[1] = rho
   '''
   rho_v = variables[0]
   rho = variables[1]
   return rho_v / rho

def v_x( variables ):
   result = v(variables)
   if len(result.shape) == 1:
      return result[0]
   return result[:,0]

def v_y( variables ):
   result = v(variables)
   if len(result.shape) == 1:
      return result[1]
   return result[:,1]

def v_z( variables ):
   result = v(variables)
   if len(result.shape) == 1:
      return result[2]
   return result[:,2]

def v_mag( variables ):
   return np.sum(np.abs(v(variables))**2,axis=-1)**(1./2)


datareducers = {}
datareducers["B_x"] = DataReducerVariable("B", B_x, "B_x", "T")
datareducers["B_y"] = DataReducerVariable("B", B_y, "B_y", "T")
datareducers["B_z"] = DataReducerVariable("B", B_z, "B_z", "T")
datareducers["B_mag"] = DataReducerVariable("B", B_mag, "B_mag", "T")

datareducers["E_x"] = DataReducerVariable("E", E_x, "E_x", "V/m")
datareducers["E_y"] = DataReducerVariable("E", E_y, "E_y", "V/m")
datareducers["E_z"] = DataReducerVariable("E", E_z, "E_z", "V/m")
datareducers["E_mag"] = DataReducerVariable("E", E_mag, "E_mag", "V/m")

datareducers["v"] = DataReducerVariable(["rho_v", "rho"], v, "v", "m/s")
datareducers["v_x"] = DataReducerVariable(["rho_v", "rho"], v_x, "v_x", "m/s")
datareducers["v_y"] = DataReducerVariable(["rho_v", "rho"], v_y, "v_y", "m/s")
datareducers["v_z"] = DataReducerVariable(["rho_v", "rho"], v_z, "v_z", "m/s")
datareducers["v_mag"] = DataReducerVariable(["rho_v", "rho"], v_mag, "v_mag", "m/s")





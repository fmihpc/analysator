''' A file for doing data reduction on variables
'''
import numpy as np
from reducer import DataReducerVariable

def x_component( variable ):
   if len(variable.shape) == 1:
      return variable[0]
   return variable[:,0]

def y_component( variable ):
   if len(variable.shape) == 1:
      return variable[1]
   return variable[:,1]

def z_component( variable ):
   if len(variable.shape) == 1:
      return variable[2]
   return variable[:,2]

def magnitude( variable ):
   return np.sum(variable**2,axis=-1)**(1./2)

#do nothing
def pass_op( variable ):
   return variable


def v( variables ):
   ''' Data reducer function for getting velocity from rho and rho_v
       variables[0] = rho_v and variables[1] = rho
   '''
   rho_v = variables[0]
   rho = variables[1]
   return rho_v / rho


#datareducers with more complex, case dependent structure. 
datareducers = {}
datareducers["v"] = DataReducerVariable(["rho_v", "rho"], v, "m/s")


#list of operators. The user can apply these to any variable,
#including more general datareducers. Can only be used to reduce one
#variable at a time
datareduction_operators = {}
datareduction_operators["pass"] = pass_op
datareduction_operators["magnitude"] = magnitude
datareduction_operators["x"] = x_component
datareduction_operators["y"] = y_component
datareduction_operators["z"] = z_component










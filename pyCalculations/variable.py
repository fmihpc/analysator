# This file has a class "Variable" that holds all the important data for variables e.g. variable's name, the units and the data on the variable
import numpy as np

class VariableInfo:
   ''' A class/struct for holding variable info. This includes the variable data in array form, the name and the units.
       Variable info is in: data, name, units.
       NOTE:
       LIST VARIABLES IN VARIBLE INFO:
       data              Data of the variable (array list)
       name              Name of the variable
       units             Units of the variable
   '''
   def __init__(self, data_array, name="", units=""):
      self.data = np.asarray(data_array)
      self.name = name
      self.units = units

   def __repr__(self):
      return "VariableInfo(Name: \'" + str(self.name) + "\' Units: \'" + str(self.units) + "\')"

   def __str__(self):
      return self.__repr__()

   def get_variable(self, index ):
      ''' Creates a new variable with identical data but only certain index is included
          E.g. if the data is an array of 3d vectors, get_variable(0) would return the variable with data[:,0] as the data
          :param index         Vector index
          :returns an edited version of the variable
      '''
      import numpy as np
      if len(self.data) <= 0:
         print "BAD DATA LENGTH"
         return []
      if len(np.atleast_1d(self.data[0])) <= index:
         print "BAD INDEX, THE INDEX IS LARGER THAN VECTOR SIZE!"
         return []
      return VariableInfo( self.data[:,index], self.name, self.units )



def get_data( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance
       :param variable           The variable as a VariableInfo instance or a list
       returns data of the variable
   '''
   if isinstance(variable, VariableInfo):
      return variable.data
   else:
      return variable

def get_name( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance
       :param variable           The variable as a VariableInfo instance or a list
       returns the name of the variable or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.name
   else:
      return ""

def get_units( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance
       :param variable           The variable as a VariableInfo instance or a list
       returns the units of the variable or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.units
   else:
      return ""


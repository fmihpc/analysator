# This file has a class "Variable" that holds all the important data for variables e.g. variable's name, the units and the data on the variable

class VariableInfo:
   ''' A class/struct for holding variable info. This includes the variable data in array form, the name and the units
   '''
   def __init__(self, data_array, name="", units=""):
      self.data = data_array
      self.name = name
      self.units = units



def get_data( variable ):
   if isinstance(variable, VariableInfo):
      return variable.data
   else:
      return variable

def get_name( variable ):
   if isinstance(variable, VariableInfo):
      return variable.name
   else:
      return ""

def get_units( variable ):
   if isinstance(variable, VariableInfo):
      return variable.units
   else:
      return ""


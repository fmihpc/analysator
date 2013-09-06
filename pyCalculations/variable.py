# This file has a class "Variable" that holds all the important data for variables e.g. variable's name, the units and the data on the variable

class VariableInfo:
   def __init__(self, data_array):
      self.data = data_array
      self.name = ""
      self.units = ""



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


# This file has a class "Variable" that holds all the important data for variables e.g. variable's name, the units and the data on the variable

class VariableInfo:
   self.name
   self.data
   self.units
   def __init__(self, data_array):
      self.data = data_array
      self.name = ""
      self.units = ""

# This file has a function for reading multiple cellids in and storing the values into an array

def cellids( vlsvReader, cellids, variable, units="" ):
   ''' Reads the given variable's values for given cell ids
       :param vlsvReader           Some VlsvFile with a file open
       :param cellids              Some cell ids
       :param variables            Some list of variable names E.g. ["rho", "B"]
       :returns the variable values for given cell ids
   '''
   import numpy as np
   from variable import get_data
   # Read variables in:
   variable_values = vlsvReader.read_variables_for_cellids( variable, get_data(cellids) )
   from output import output_1d
   return output_1d( [variable_values], [variable], [units] )

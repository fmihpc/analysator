# This file has a function for reading multiple cellids in and storing the values into an array

def cellids( vlsvReader, cellids, variable, units="" ):
   ''' Reads the given variable's values for given cell ids

       :param vlsvReader:           Some VlsvFile with a file open
       :type vlsvReader:            :class:`vlsvreader.VlsvFile`
       :param cellids:              Some cell ids
       :param variable:             Some variable name
       :returns: the variable values for given cell ids

       .. seealso:: modules :py:mod::'vlsvreader'
   '''
   import numpy as np
   from variable import get_data
   # Read variables in:
   variable_values = vlsvReader.read_variables_for_cellids( variable, get_data(cellids) )
   from output import output_1d
   return output_1d( [variable_values], [variable], [units] )

# This file contains static spacecraft function for fetching variable data from multiple vlsv files and looking at the data in time

import numpy as np

def cell_time_evolution( vlsvReader_list, variables, cellids, units="" ):
   ''' Returns variable data from a time evolution of some certain cell ids

       :param vlsvReader_list:         List containing VlsvReaders with a file open
       :type vlsvReader_list:          :class:`vlsvfile.VlsvReader`
       :param variables:               Name of the variables
       :param cellids:                 List of cell ids
       :param units:                   List of units for the variables (OPTIONAL)
       :returns: an array containing the data for the time evolution for every cell id

       .. code-block:: python

          # Example of the return list with 3 variable names and 2 cell ids:
          [cellid1variable1, cellid1variable2, cellid1variable3, cellid2variable1, cellid2variable2, cellid3variable3]
   
          # Example of usage:
          time_data = cell_time_evolution( vlsvReader_list=[VlsvReader("bulk.000.vlsv"), VlsvReader("bulk.001.vlsv"), VlsvReader("bulk.002.vlsv")], variables=["rho", "Pressure", "B"], cellids=[2,4], units=["N", "Pascal", "T"] )
   '''
   vlsvReader_list = np.atleast_1d(vlsvReader_list)
   variables = np.atleast_1d(variables)
   cellids = np.atleast_1d(cellids)
   parameters = ["t","tstep","fileIndex"]
   parameter_units=["s","",""]
   #construct empty units, if none are given
   if (units == "") or (len(units) != len(variables)):
      units=[ "" for i in xrange(len(variables))]
   #construct data 
   data = [[] for i in xrange(len(parameters)+len(cellids)*len(variables))]
   for t in xrange(len(vlsvReader_list)):
      # Get the vlsv reader
      vlsvReader = vlsvReader_list[t]
      # Open the vlsv reader's file:
      vlsvReader.optimize_open_file()
      #go through parameters
      for j in xrange(len(parameters)):
         # Read the parameter
         # Save the data into the right slot in the data array:
         data[j].append(vlsvReader.read_parameter( parameters[j]))

      # Go through variables:
      for j in xrange(len(variables)):
         variable = variables[j]
         # Read the variable for all cell ids
         #variables_for_cellids = vlsvReader.read_variables_for_cellids( variable, cellids )
         # Save the data into the right slot in the data array:
         for i in xrange(len(cellids)):
            data[len(parameters)+i*len(variables)+j].append(vlsvReader.read_variable( variable, cellids[i] ))
      # For optimization purposes we are now freeing vlsvReader's memory
      # Note: Upon reading data vlsvReader created an internal hash map that takes a lot of memory
      vlsvReader.optimize_clear_fileindex_for_cellid()
      # Close the vlsv reader's file:
      vlsvReader.optimize_close_file()
   from output import output_1d
   return output_1d( data, 
                     parameters +  [variables[(int)(i)%(int)(len(variables))] for i in xrange(len(data)-len(parameters))], 
                     parameter_units + [units[(int)(i)%(int)(len(units))] for i in xrange(len(data)-len(parameters))] )



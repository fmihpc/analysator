# This file contains static spacecraft function for fetching variable data from multiple vlsv files and looking at the data in time

def virtual_spacecraft( vlsvReader_list, variable_names, cellids )
   ''' Returns variable data from a time evolution of some certain cell ids
       :param vlsvReader_list         List containing VlsvFiles with a file open
       :param variable_names          Name of the variables
       :param cellids                 List of cell ids
       :returns an array containing the data for the time evolution for every cell id
       Example of the return list with 3 variable names and 2 cell ids:
       [cellid1variable1, cellid1variable2, cellid1variable3, cellid2variable1, cellid2variable2, cellid3variable3]
   '''
   vlsvReader_list = np.atleast_1d(vlsvReader_list)
   variable_names = np.atleast_1d(variable_names)
   cellids = np.atleast_1d(cellids)
   data = [[[] for j in xrange(len(variable_names))] for i in xrange(len(cellids))]
   for t in xrange(len(vlsvReader_list)):
      # Get the vlsv reader
      vlsvReader = vlsvReader_list[t]
      # Go through variables:
      for j in xrange(len(variable_names)):
         variable = variable_names[i]
         # Read the variable for all cell ids
         variables_for_cellids = vlsvReader.read_variables_for_cellids( variable, cellids )
         # Save the data into the right slot in the data array:
         for i in xrange(len(cellids)):
            data[i][j].append(variables_for_cellids[i])
      # For optimization purposes we are now freeing vlsvReader's memory
      # Note: Upon reading data vlsvReader created an internal hash map that takes a lot of memory
      vlsvReader.clear_fileindex_for_cellid()
   return np.ravel(data)



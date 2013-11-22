# File that creates the output format for arrays

def output_1d( arrays, names, units="" ):
   ''' Creates an output out of 1d arrays

       :param arrays:    Some arrays with data in them
       :param names:     Names for the array elements
       :param units:     Units for the arrays (optional)
       :returns: the arrays in a new format (dictionary currently)

       .. note::

          arrays and names must be of same length

       .. code-block:: python

          #Example usage:
          output_1d( [[2,1,23], [5,78,4], [2,3,2]], ["rho", "B", "Pressure"] )
          This would interpret [2,1,23] as an array called \"rho\"
   '''
   if units == "":
      units = ["" for i in xrange(len(arrays))]
   if( (len(arrays) != len(names)) or (len(arrays) != len(units)) ):
      print "BAD ARRAY AND NAME LENGTH IN OUTPUT_1D (pyCalculations/output.py)"
      return []
   new_format = []
   from variable import VariableInfo
   for i in xrange(len(arrays)):
      variable = VariableInfo(arrays[i])
      variable.name = names[i]
      variable.units = units[i]
      new_format.append(variable)
   if len(new_format) == 1:
      return new_format[0]
   else:
      return new_format

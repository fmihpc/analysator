# A file for testing python tools
import pytools as pt
import numpy as np

def vlsv_test( filename, datafilename ):
   ''' A function for testing that the current trunk revision is working. The way this functions is that we have a test vlsv file and a test data file name. There should be a data file made by a working trunk revision and then this function should be run to create a new data file. They two data files should look indentical if nothing is wrong.
       :param filename:           Name of the test vlsv file
       :param datafilename:       Name of the data file to write
   '''
   # Open a file for reading:
   vlsvReader = pt.vlsvfile.VlsvReader(filename)

   vlsvReader.list()

   # Read some data.
   data = []
   data.append( str( vlsvReader.read_variable("rho") ) )
   data.append( str( vlsvReader.read_variable_info("rho") ) )

   data.append( str( vlsvReader.check_variable("rho") ) )
   data.append( str( vlsvReader.check_variable("B_vol") ) )
   data.append( str( vlsvReader.check_variable("false_test") ) )

   data.append( str( vlsvReader.get_cellid_locations() ) )

   cellid = 135291
   coordinates = vlsvReader.get_cell_coordinates( cellid )

   data.append( str(vlsvReader.get_cellid( coordinates )) )

   if int(data[len(data)-1]) != cellid:
      print "ERROR, BAD CELLID/COORDINATES IN VLSV TEST " + str(data[len(data)-1])
      data.append( "ERROR, BAD CELLID/COORDINATES IN VLSV TEST" )

   data.append( str(coordinates) )

   data.append( str( vlsvReader.read_parameter("xcells_ini") ) )

   data.append( str( vlsvReader.read_velocity_cells( cellid ) ) )

   data.append( str( vlsvReader.read_blocks( cellid ) ) )

   # Get a cut through
   point1 = np.array( [0, 0, 0] )
   point2 = np.array( [50e6, 60e6, 0] )
   cutthrough = pt.calculations.cut_through( vlsvReader, point1, point2 )

   data.append( str( cutthrough[0].data ) )
   data.append( str( cutthrough[1].data ) )

   # Get pitch angles:
   pitch_angles = pt.calculations.pitch_angles( vlsvReader, cellid, cosine=True, plasmaframe=False )
   data.append( str( pitch_angles[0].data ) )
   data.append( str( pitch_angles[1].data ) )
   

   # Open a file for writing:
   out = open(datafilename, 'w')
   # Write data into the testfile
   for i in data:
      out.write(i)
      out.write("\n")

#
#def run_tests():
#   vlsv_test()
#   
#
############## MAIN #########################
#run_tests()

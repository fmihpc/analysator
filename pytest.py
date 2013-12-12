# A file for testing python tools
import pytools as pt
import numpy as np
from hashlib import md5

def vlsv_test( filename, datafilename ):
   ''' A function for testing that the current trunk revision is working. The way this functions is that we have a test vlsv file and a test data file name. There should be a data file made by a working trunk revision and then this function should be run to create a new data file. They two data files should look indentical if nothing is wrong.
       :param filename:           Name of the test vlsv file
       :param datafilename:       Name of the data file to write
   '''
   # Open a file for reading:
   vlsvReader = pt.vlsvreader.VlsvFile(filename)

   # Read some data.
   data = []
   data.append( ["rho_val", str( vlsvReader.read_variable("rho") )] )
   data.append( ["rho_variable_info", str( vlsvReader.read_variable_info("rho") )] )

   data.append( ["rho_check", str( vlsvReader.check_variable("rho") )] )
   data.append( ["B_vol_check", str( vlsvReader.check_variable("B_vol") )] )
   data.append( ["False_test_check", str( vlsvReader.check_variable("false_test") )] )

   data.append( ["cellid_locations", str( vlsvReader.get_cellid_locations() )] )

   cellid = 135291
   coordinates = vlsvReader.get_cell_coordinates( cellid )

   data.append( ["cellid", str(vlsvReader.get_cellid( coordinates ))] )

   if int(data[len(data)-1][1]) != cellid:
      print "ERROR, BAD CELLID/COORDINATES IN VLSV TEST " + str(data[len(data)-1])
      data.append( "ERROR, BAD CELLID/COORDINATES IN VLSV TEST", "" )

   data.append( ["coordinates", str(coordinates)] )

   data.append( ["parameter_xcells", str( vlsvReader.read_parameter("xcells_ini") )] )

   data.append( ["velocity_cells", str( vlsvReader.read_velocity_cells( cellid ) )] )

   data.append( ["blocks_check", str( vlsvReader.read_blocks( cellid ) )] )

   # Get a cut through
   point1 = np.array( [0, 0, 0] )
   point2 = np.array( [50e6, 60e6, 0] )
   cutthrough = pt.calculations.cut_through( vlsvReader, point1, point2 )

   data.append( ["cutthrough_cellids", str( cutthrough[0].data )] )
   data.append( ["cutthrough_distances", str( cutthrough[1].data )] )

   # Get pitch angles:
   pitch_angles = pt.calculations.pitch_angles( vlsvReader, cellid, cosine=True, plasmaframe=False )
   data.append( ["pitch_anlges_avgs", str( pitch_angles[0].data )] )
   data.append( ["pitch_angles_data", str( pitch_angles[1].data )] )
   
   # Open a file for writing:
   out = open(datafilename, 'w')
   # Write data into the testfile
   for i in data:
      # Convert string to md5
      outwrite = md5()
      outwrite.update(i[1])
      out.write(i[0] + ": " + str(outwrite.hexdigest()))
      out.write("\n")

#
#def run_tests():
#   vlsv_test()
#   
#
############## MAIN #########################
#run_tests()

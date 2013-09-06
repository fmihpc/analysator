# Cut-throughs from vlsv files

import numpy as np

def cut_through( vlsvReader, point1, point2 ):
   ''' Returns cell ids and distances from point 1 for every cell in a line between given point1 and point2
       :param vlsvReader       Some open VlsvFile
       :param point1           The starting point of a cut-through line
       :param point2           The ending point of a cut-through line
       :returns an array containing cell ids, coordinates and distances in the following format: [cell ids, distances]
   '''
   # Transform point1 and point2 into numpy array:
   point1 = np.array(point1)
   point2 = np.array(point2)
   # Get parameters from the file to determine a good length between points (step length):
   # Get xmax, xmin and xcells_ini
   xmax = vlsvReader.read_parameter(name="xmax")
   xmin = vlsvReader.read_parameter(name="xmin")
   xcells = vlsvReader.read_parameter(name="xcells_ini")
   # Do the same for y
   ymax = vlsvReader.read_parameter(name="ymax")
   ymin = vlsvReader.read_parameter(name="ymin")
   ycells = vlsvReader.read_parameter(name="ycells_ini")
   # And for z
   zmax = vlsvReader.read_parameter(name="zmax")
   zmin = vlsvReader.read_parameter(name="zmin")
   zcells = vlsvReader.read_parameter(name="zcells_ini")
   #Calculate cell lengths:
   cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])

   # Get coordinates between the points:
   coordinateList = []
   unitVec = (point2 - point1)/np.linalg.norm(point2 - point1)
   # Get a vector of cell length's length in the direction of (point2-point1)
   oneVec = unitVec * np.min(cell_lengths)
   # Get the number of points:
   numberOfPoints = (int)(np.linalg.norm(point2 - point1) / np.linalg.norm(oneVec)) + 1
   # Some routine checks:
   if numberOfPoints <= 2:
      print "ERROR, TOO SMALL A DISTANCE BETWEEN POINT1 AND POINT2"
      return
   # Input coordinates:
   for i in xrange((int)(np.linalg.norm(point2 - point1) / np.linalg.norm(oneVec)) + 1):
      coordinateList.append( point1 + oneVec*i )
   cellids = []
   coordinates = []
   distances = []
   curDists = []
   appDist = False
   for coordinate in coordinateList:
      # Get the cell indices
      cellindices = np.array([(int)((coordinate[0] - xmin)/(float)(cell_lengths[0])), (int)((coordinate[1] - ymin)/(float)(cell_lengths[1])), (int)((coordinate[2] - zmin)/(float)(cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * xcells + cellindices[2] * xcells * ycells + 1
      curDists.append(np.linalg.norm(coordinate - point1))
      # If the cell id is already in the list, don't append:
      if len(cellids) > 0:
         if cellid == cellids[len(cellids)-1]:
            continue
      # Append everything:
      cellids.append(cellid)
      coordinates.append(coordinate)
      if appDist == True:
         distances.append(np.median(curDists))
         curDists = []
      else:
         appDist = True
   if len(curDists) != 0:
      distances.append(np.median(curDists))
   else:
      distances.append(distances[len(distances)-1])
   # Return the cut through:
   # Note: Make the output into dictionary
   from output import output_1d
   #return [np.array(cellids[0:len(cellids)-2]), np.array(distances[0:len(distances)-2])]
   return output_1d( [np.array(cellids[0:len(cellids)-2]), np.array(distances[0:len(distances)-2])], ["CellID", "distances"], ["", "m"])













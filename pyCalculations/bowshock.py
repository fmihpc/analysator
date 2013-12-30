# This file has functions concerning bow shock and its detection

import numpy as np

def rankine_hugoniot( V_2, V_1, rho_2, rho_1, B_2, B_1, P_2, P_1, shock_normal ):
   ''' Calculates the rankine hugoniot jump conditions for oblique shock with the given variables

       :param V_2: the speed on the second side of the shock border
       :param V_1: The speed on the first side of the shock border
       :param rho_2: The charge density on the second side of the shock border
       :param rho_1: The charge density on the first side of the shock border
       :param P_2: The pressure on the second side of the shock border
       :param P_1: The pressure on the first side of the shock border
       :param shock_normal: The shock normal vector
       :returns: Vector containing the results of rankine-hugoniot conditions (The closer values are to zero, the closer the rankine-hugoniot condiions are to their theoretical value)
   '''
   # Make sure the shock normal vector is a unit one
   shock_normal = shock_normal / np.linalg.norm(shock_normal)
   # Calculate the V_x and V_y where x is in the direction of the shock normal
   print V_1
   V_1x = np.dot( V_1, shock_normal )
   V_1y = V_1 - V_1x
   V_2x = np.dot( V_2, shock_normal )
   V_2y = V_2 - V_2x
   # Calculate the Xi value
   Xi = rho_2 / rho_1
   # Calculate the B_x and B_y values
   B_1x = np.dot( B_1, shock_normal )
   B_1y = B_1 - B_1x
   B_2x = np.dot( B_2, shock_normal )
   B_2y = B_2 - B_2x
   # Calculate the alfven velocity norm
   vacuum_permeability = 1.257e-6
   V_a1 = np.linalg.norm(B_1) / np.sqrt(vacuum_permeability * rho_1)
   # Calculate the rankine-hugoniot jump conditions for oblique shocks
   gamma = 5.0/3.0
   rankine_hugoniot = [
             V_2x / V_1x - 1/Xi,
             V_2y / V_1y - (np.linalg.norm(V_1)**2 - np.linalg.norm(V_2)**2)/(np.linalg.norm(V_1)**2 - Xi * V_a1**2),
             B_2x / B_1x - 1,
             B_2y / B_1y - ((np.linalg.norm(V_1)**2 - np.linalg.norm(V_2)**2)*Xi / (np.linalg.norm(V_1)**2 - Xi * V_a1**2)),
             P_2 / P_1 - (Xi + (gamma - 1)*Xi*np.linalg.norm(V_1)**2 / (2*V_s1**2)*(1 - np.linalg.norm(V_2)**2 / np.linalg.norm(V_1)**2))
                      ]
   return rankine_hugoniot

def rankine_hugoniot_two_points( vlsvReader, point1, point2 ):
   ''' Calculates the rankine hugoniot jump conditions for oblique shock for two points

       :param vlsvReader: Some VlsvFile with a file open
       :param point1:     The first point of interest
       :param point2:     The second point of interest
       :returns: Vector containing the results of rankine-hugoniot conditions (The closer values are to zero, the closer the rankine-hugoniot condiions are to their theoretical value)

       .. note:: If point1 is on one side of the shock border and the point2 is on the other side of the shock border, the result should be close to [0, 0, 0, 0, 0]

       .. seealso:: :func:`rankine_hugoniot`

       .. code-block: python
          print rankine_hugoniot_two_points( vlsvReader, [0,0,0], [10e4, 10e4, 0] )
   '''
   # Get the cell ids:
   cellid_1 = vlsvReader.get_cellid( point1 )
   cellid_2 = vlsvReader.get_cellid( point2 )
   if cellid_1 == 0 or cellid_2 == 0:
      print "Error! Point1 or point2 cellid not found!"
      return
   # Read the necessary variables:
   V_1 = vlsvReader.read_variable( "v", cellids=cellid_1 )
   rho_1 = vlsvReader.read_variable( "rho", cellids=cellid_1 )
   B_1 = vlsvReader.read_variable( "B", cellids=cellid_1 )
   P_1 = vlsvReader.read_variable( "Pressure", cellids=cellid_1 )

   V_2 = vlsvReader.read_variable( "v", cellids=cellid_2 )
   rho_2 = vlsvReader.read_variable( "rho", cellids=cellid_2 )
   B_2 = vlsvReader.read_variable( "B", cellids=cellid_2 )
   P_2 = vlsvReader.read_variable( "Pressure", cellids=cellid_2 )

   # Define the shock normal as the vector pointing from point1 to point2
   shock_normal = np.array(point2) - np.array(point1)

   return rankine_hugoniot( V_2, V_1, rho_2, rho_1, B_2, B_1, P_2, P_1, shock_normal )











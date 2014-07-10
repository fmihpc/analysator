# Rankine-Hugoniot conditions to determine the right state from the left state of an MHD oblique shock
import numpy as np
import pylab as pl
import sys
from scipy import optimize
from cutthrough import cut_through
from variable import VariableInfo
from plot import plot_variables
from rankine_solver import solve_rankine_SOE


def rotation_matrix_2d( angle ):
   ''' Creates a rotation matrix that can be used to rotate a 2d vector by an angle in the counter-clockwise direction

       :param angle: Rotation angle
       :returns: The rotation matrix
   '''
   return np.array([[np.cos(angle), -1*np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])


def oblique_solver( Vx1, Vy1, Bx1, By1, T1, P1, rho1 ):
   ''' Calculates the rankine hugoniot jump conditions on the other side of the shock for given parameters

       :param Vx1: Velocity component parallel to the shock normal vector on one side of the shock
       :param Vy1: Velocity component perpendicular to the shock normal vector on one side of the shock
       :param Bx1: Magnetic field component parallel to the shock normal vector on one side of the shock
       :param By1: Magnetic field component perpendicular to the shock normal vector on one side of the shock
       :param T1: Temperature on one side of the shock
       :param rho1: Density on one side of the shock

       :returns: Components on the other side plus pressure P2 and compression ratio X in format [Vx2, Vy2, Bx2, By2, T2, rho2, P2, X]
   '''
   import sympy as sp

   # Constants
   y = 5./3.
   mp = 1.67e-27;
   mu0 = 4.0 * np.pi * 1e-7;
   kb = 1.3806505e-23;
   # Calculate other variables
   theta = np.arccos(Bx1/np.sqrt(Bx1**2 + By1**2));
   vA1 = np.cos(theta)*np.sqrt( Bx1**2 / (mp * rho1 * mu0) )
   V1 = np.sqrt( Vx1**2 + Vy1**2 )
   vs1 = np.sqrt( y * kb * T1 / mp )

   # Declare the rest of the variables
   Vx2, Vy2, Bx2, By2, T2, P2, rho2, X = sp.symbols('Vx2, Vy2, Bx2, By2, T2, P2, rho2, X')
   V2 = sp.sqrt(Vx2*Vx2+Vy2*Vy2)

   # Solve X:
   x = sp.solve( (V1**2-X*vA1**2)**2*(X*vs1**2+1/2.*V1**2*np.cos(theta)**2*(X*(y-1)-(y+1))) + 1/2.*vA1**2*V1**2*np.sin(theta)**2*X*((y+X*(2-y))*V1**2-X*vA1**2*((y-1)-X*(y-1))) , X)
   # Pick X:
   for i in x:
      if i.as_real_imag()[1] == 0 and i.as_real_imag()[0] > 0:
         X = i.as_real_imag()[0]

   # Write down the equations
   equation = []
   equation.append( Vx2/Vx1 - 1/X )
   equation.append( rho1/rho2 - 1/X )
   equation.append( Vy2/Vy1 - (V1*V1-vA1*vA1) / (V1*V1-X*vA1*vA1) )
   equation.append( Bx2/Bx1 - 1 )
   equation.append( By2/By1 - (V1*V1-vA1*vA1)*X/(V1*V1-X*vA1*vA1) )
   equation.append( P2/P1 - (X+(y-1)*X*V1*V1/(2*vs1*vs1))*(1-V2*V2/(V1*V1)) )

   return sp.solve(equation, Vx2, Vy2, Bx2, By2, T2, P2, rho2 )


def oblique_shock( Vx1, Vy1, Bx1, By1, T1, rho1 ):
   ''' Calculates the rankine hugoniot jump conditions on the other side of the shock for given parameters

       :param Vx1: Velocity component parallel to the shock normal vector on one side of the shock
       :param Vy1: Velocity component perpendicular to the shock normal vector on one side of the shock
       :param Bx1: Magnetic field component parallel to the shock normal vector on one side of the shock
       :param By1: Magnetic field component perpendicular to the shock normal vector on one side of the shock
       :param T1: Temperature on one side of the shock
       :param rho1: Density on one side of the shock

       :returns: Components on the other side plus pressure P2 and compression ratio X in format [Vx2, Vy2, Bx2, By2, T2, rho2, P2, X]
   '''

   # Constants
   mp = 1.67e-27;
   mu0 = 4.0 * np.pi * 1e-7;
   kb = 1.3806505e-23;
   
   # Input variables
   Gamma = 5./3.;
   
   # Calculate other variables
   theta = np.arccos(Bx1/np.sqrt(Bx1**2 + By1**2));
   vA1 = np.sqrt( Bx1**2 / (mp * rho1 * mu0) )
   V1 = np.sqrt( Vx1**2 + Vy1**2 )
   P1 = rho1 * kb * T1
   vs1 = np.sqrt( Gamma * kb * T1 / mp )
   
   # Function for solving X, the compression ratio
   def solver(x):
      return (V1**2 - x * vA1**2)**2 * (x * vs1**2 + 0.5 * V1**2 * np.cos(theta)**2 * (x * (Gamma - 1) -(Gamma + 1))) + 0.5 * vA1**2 * V1**2 * np.sin(theta)**2 * x * ((Gamma + x * (2 - Gamma)) * V1**2 - x * vA1**2 * ((Gamma + 1) - x * (Gamma -1)))

   solutions = []

   X = optimize.fsolve(solver, 0.1)[0]

   solutions.append(X)

   epsilon = 0.1

   test_values = np.arange(0.1,60,0.1)
   for value in test_values:
      result = optimize.fsolve(solver, value)[0]
      if ((result <= X + epsilon and result >= X - epsilon) == False) and (solver(result) >= -1*epsilon and solver(result) <= epsilon):
         X = result
         solutions.append(X)
   
   if len(solutions) != 1:
      print "MORE THAN ONE SOLUTION: "
      print solutions
      return

   X = solutions[0]
   print X
   rho2 = rho1 * X
   Vx2 = Vx1 / X
   Vy2 = Vy1 * (V1**2 - vA1**2) / (V1**2 - X * vA1**2)
   
   V2 = np.sqrt( Vx2**2 + Vy2**2 );
   
   Bx2 = Bx1
   By2 = By1 * (V1**2 - vA1**2) * X / (V1**2 - X * vA1**2)
   P2 = P1 * (X + (Gamma - 1) * X * V1**2 * (1 - V2**2 / V1**2) / (2.0 * vs1**2));
   T2 = P2 / (rho2 * kb)
   return [Vx2, Vy2, Bx2, By2, T2, rho2, P2, X ]


def evaluate_rankine( vlsvReader, point1, point2 ):
   ''' A function for comparing the values predicted by Rankine-Hugoniot and the values from the data

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest along the bow-shock
   :param point2: The second point of interest along the bow-shock

   :returns: the evaluation values as a list: [Rhorankine/Rhovlasiator, Bxrankine/Bxvlasiator, Byrankine/Byvlasiator Vxrankine/Vxvlasiator, Vyrankine/Vyvlasiator]
   '''
   # Get spatial grid sizes:
   xcells = (int)(vlsvReader.read_parameter("xcells_ini"))
   ycells = (int)(vlsvReader.read_parameter("ycells_ini"))
   zcells = (int)(vlsvReader.read_parameter("zcells_ini"))

   xmin = vlsvReader.read_parameter("xmin")
   ymin = vlsvReader.read_parameter("ymin")
   zmin = vlsvReader.read_parameter("zmin")
   xmax = vlsvReader.read_parameter("xmax")
   ymax = vlsvReader.read_parameter("ymax")
   zmax = vlsvReader.read_parameter("zmax")

   dx = (xmax - xmin) / (float)(xcells)
   dy = (ymax - ymin) / (float)(ycells)
   dz = (zmax - zmin) / (float)(zcells)

   # Get normal vector from point2 and point1
   point1 = np.array(point1)
   point2 = np.array(point2)
   normal_vector = (point2-point1) / np.linalg.norm(point2 - point1)
   normal_vector = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
   normal_vector = normal_vector * np.array([1,1,0])
   point1_shifted = point1 + 0.5*(point2-point1) + normal_vector * (4*dx)
   point2_shifted = point1 + 0.5*(point2-point1) - normal_vector * (4*dx)
   point1 = np.array(point1_shifted)
   point2 = np.array(point2_shifted)

   # Read cut-through
   cutthrough = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
   # Get cell ids and distances separately
   cellids = cutthrough[0].data
   distances = cutthrough[1]
   # Read data from the file:
   V_data = vlsvReader.read_variable( "v", cellids=cellids )
   B_data = vlsvReader.read_variable( "B", cellids=cellids )
   T_data = vlsvReader.read_variable( "Temperature", cellids=cellids )
   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )
   P_data = vlsvReader.read_variable( "Pressure", cellids=cellids )

   # Transform to de hoffmann teller frame:
   ##############################################
   # V_1 x B_1 = 0:
   V_1 = V_data[0]
   B_1 = B_data[0]
   n = normal_vector
   V_HT = np.cross(n, np.cross(V_1, B_1) ) / np.dot(n, B_1)

   # Check whether V_HT is relativistic or not:
   if np.linalg.norm( V_HT ) > 0.05*3.0e8:
      print "WARNING: When transforming to de Hoffmann teller frame encountered speeds over 0.05 speed of light: " + str(np.linalg.norm( V_HT ))

   # Transform to another frame:
   V_data = V_data - V_HT
   # Check the de-hoffmann teller frame:
   print "HOFFMAN_CHECK: " + str(np.cross(B_1, V_1))
   ##############################################

   # Get parallel and perpendicular components from the first vector (on one side of the shock):
   Vx_data = np.dot(V_data, normal_vector)
   Vy_data = np.sqrt(np.sum(np.abs(V_data - np.outer(Vx_data, normal_vector))**2,axis=-1))
   Bx_data = np.dot(B_data, normal_vector) 
   By_data = np.sqrt(np.sum(np.abs(B_data - np.outer(Bx_data, normal_vector) )**2,axis=-1))

   # Get the upstream and downstream area:
   upstream = 0
   downstream = 0
   for i in xrange(len(rho_data)):
      if rho_data[i] > np.average(rho_data):
         upstream = i-3
         downstream = i+2
         break;

   # Input data: (u = upstream)
   V_u = np.mean(V_data[0:upstream], axis=0)
   B_u = np.mean(B_data[0:upstream], axis=0)
   T_u = np.mean(T_data[0:upstream])
   rho_u = np.mean(rho_data[0:upstream])
   P_u = np.mean(P_data[0:upstream])
   # Get parallel and perpendicular components:
   Vx_u = np.dot(V_u, normal_vector)
   Vy_u = np.linalg.norm(V_u - Vx_u * normal_vector)
   Bx_u = np.dot(B_u, normal_vector)
   By_u = np.linalg.norm(B_u - Bx_u * normal_vector)

   # Input data: (d = downstream)
   V_d = np.mean(V_data[downstream:len(V_data)-1], axis=0)
   B_d = np.mean(B_data[downstream:len(V_data)-1], axis=0)
   T_d = np.mean(T_data[downstream:len(V_data)-1])
   rho_d = np.mean(rho_data[downstream:len(V_data)-1])
   P_d = np.mean(P_data[downstream:len(V_data)-1])
   # Get parallel and perpendicular components:
   Vx_d = np.dot(V_d, normal_vector)
   Vy_d = np.linalg.norm(V_d - Vx_d * normal_vector)
   Bx_d = np.dot(B_d, normal_vector)
   By_d = np.linalg.norm(B_d - Bx_d * normal_vector)


   # Calculate rankine hugoniot jump conditions:
   rankine_conditions = oblique_shock( Vx_u, Vy_u, Bx_u, By_u, T_u, rho_u )

   # Input Rankine (d = downstream):
   Vx_rankine_d = rankine_conditions[0]
   Vy_rankine_d = rankine_conditions[1]
   Bx_rankine_d = rankine_conditions[2]
   By_rankine_d = rankine_conditions[3]
   T_rankine_d = rankine_conditions[4]
   rho_rankine_d = rankine_conditions[5]

   # Get the differences:
   # Data
   Vx = Vx_d - Vx_u
   Vy = Vy_d - Vy_u
   Bx = Bx_d - Bx_u
   By = By_d - By_u
   T = T_d - T_u
   rho = rho_d - rho_u

   # Rankine
   Vx_rankine = Vx_rankine_d - Vx_u
   Vy_rankine = Vy_rankine_d - Vy_u
   Bx_rankine = Bx_rankine_d - Bx_u
   By_rankine = By_rankine_d - By_u
   T_rankine = T_rankine_d - T_u
   rho_rankine = rho_rankine_d - rho_u



   # Return the comparisons:
   return [Vx_rankine/Vx, Vy_rankine/Vy, Bx_rankine/Bx, By_rankine/By, T_rankine/T, rho_rankine/rho]

def plot_rankine( vlsvReader, point1, point2 ):
   ''' A function that plots the theoretical rankine-hugoniot jump condition variables for two given points along with the actual values from a given vlsv file

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest along the bow-shock
   :param point2: The second point of interest along the bow-shock

   :returns: pylab figure
   '''
   # Get spatial grid sizes:
   xcells = (int)(vlsvReader.read_parameter("xcells_ini"))
   ycells = (int)(vlsvReader.read_parameter("ycells_ini"))
   zcells = (int)(vlsvReader.read_parameter("zcells_ini"))

   xmin = vlsvReader.read_parameter("xmin")
   ymin = vlsvReader.read_parameter("ymin")
   zmin = vlsvReader.read_parameter("zmin")
   xmax = vlsvReader.read_parameter("xmax")
   ymax = vlsvReader.read_parameter("ymax")
   zmax = vlsvReader.read_parameter("zmax")

   dx = (xmax - xmin) / (float)(xcells)
   dy = (ymax - ymin) / (float)(ycells)
   dz = (zmax - zmin) / (float)(zcells)

   # Get normal vector from point2 and point1
   point1 = np.array(point1)
   point2 = np.array(point2)
   normal_vector = (point2-point1) / np.linalg.norm(point2 - point1)
   normal_vector = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
   normal_vector = normal_vector * np.array([1,1,0])
   point1_shifted = point1 + 0.5*(point2-point1) + normal_vector * (4*dx)
   point2_shifted = point1 + 0.5*(point2-point1) - normal_vector * (4*dx)
   point1 = np.array(point1_shifted)
   point2 = np.array(point2_shifted)

   # Read cut-through
   cutthrough = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
   # Get cell ids and distances separately
   cellids = cutthrough[0].data
   distances = cutthrough[1]
   # Read data from the file:
   V_data = vlsvReader.read_variable( "v", cellids=cellids )
   B_data = vlsvReader.read_variable( "B", cellids=cellids )
   T_data = vlsvReader.read_variable( "Temperature", cellids=cellids )
   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )
   P_data = vlsvReader.read_variable( "Pressure", cellids=cellids )

   # Transform to de hoffmann teller frame:
   ##############################################
   # V_1 x B_1 = 0:
   V_1 = V_data[0]
   B_1 = B_data[0]
   n = normal_vector
   V_HT = np.cross(n, np.cross(V_1, B_1) ) / np.dot(n, B_1)

   # Check whether V_HT is relativistic or not:
   if np.linalg.norm( V_HT ) > 0.05*3.0e8:
      print "WARNING: When transforming to de Hoffmann teller frame encountered speeds over 0.05 speed of light: " + str(np.linalg.norm( V_HT ))

   # Transform to another frame:
   V_data = V_data - V_HT
   # Check the de-hoffmann teller frame:
   print "HOFFMAN_CHECK: " + str(np.cross(B_1, V_1))
   ##############################################

   # Get parallel and perpendicular components from the first vector (on one side of the shock):
   Vx_data = np.dot(V_data, normal_vector)
   Vy_data = np.sqrt(np.sum(np.abs(V_data - np.outer(Vx_data, normal_vector))**2,axis=-1))
   Bx_data = np.dot(B_data, normal_vector) 
   By_data = np.sqrt(np.sum(np.abs(B_data - np.outer(Bx_data, normal_vector) )**2,axis=-1))


   # Read V, B, T and rho for point1
#   V = vlsvReader.read_variable( "v", cellids=cellids[0] )
#   B = vlsvReader.read_variable( "B", cellids=cellids[0] )
#   T = vlsvReader.read_variable( "Temperature", cellids=cellids[0] )
#   rho = vlsvReader.read_variable( "rho", cellids=cellids[0] )
   V = V_data[0]
   B = B_data[0]
   T = T_data[0]
   rho = rho_data[0]
   P = P_data[0]
   # Get parallel and perpendicular components:
   Vx = np.dot(V, normal_vector)
   Vy = np.linalg.norm(V - Vx * normal_vector)
   Bx = np.dot(B, normal_vector)
   By = np.linalg.norm(B - Bx * normal_vector)

#   # Print the R-H
#   ##########
#   mu0 = 4.0 * np.pi * 1e-7
#   y = 5./3.
#   print solve_rankine_SOE( rho, V, P, B, normal_vector, mu0, y )
#   ##########

   # Calculate rankine hugoniot jump conditions:
   rankine_conditions = oblique_shock( Vx, Vy, Bx, By, T, rho )
   # Print:
   print oblique_solver( Vx, Vy, Bx, By, T, P, rho )

   # Get the upstream and downstream area:
   upstream = 0
   downstream = 0
   for i in xrange(len(rho_data)):
      if rho_data[i] > np.average(rho_data):
         upstream = i-3
         downstream = i+2
         break;

   # Input variables
   Vx_rankine = []
   Vy_rankine = []
   Bx_rankine = []
   By_rankine = []
   T_rankine = []
   rho_rankine = []
   for i in xrange( len(cellids) ):
      if i < upstream+2:
         Vx_rankine.append(Vx)
         Vy_rankine.append(Vy)
         Bx_rankine.append(Bx)
         By_rankine.append(By)
         T_rankine.append(T)
         rho_rankine.append(rho)
      else:
         Vx_rankine.append(rankine_conditions[0])
         Vy_rankine.append(rankine_conditions[1])
         Bx_rankine.append(rankine_conditions[2])
         By_rankine.append(rankine_conditions[3])
         T_rankine.append(rankine_conditions[4])
         rho_rankine.append(rankine_conditions[5])

   # Plot the variables:
   variables = []
   #VariableInfo(self, data_array, name="", units="")
   variables.append( VariableInfo(rho_data, "rho", "m^-3" ) )
   variables.append(VariableInfo(rho_rankine, "rho", "m^-3") )
   variables.append(VariableInfo(Vx_data, "Vx", "m/s"))
   variables.append(VariableInfo(Vx_rankine, "Vx", "m/s"))
   variables.append(VariableInfo(Vy_data, "Vy", "m/s"))
   variables.append(VariableInfo(Vy_rankine, "Vy", "m/s"))
#   variables.append(VariableInfo(Bx_data, "Bx", "T"))
#   variables.append(VariableInfo(Bx_rankine, "Bx", "T"))
   variables.append(VariableInfo(By_data,"By", "T"))
   variables.append(VariableInfo(By_rankine, "By", "T"))
   variables.append(VariableInfo(T_data, "T", "K"))
   variables.append(VariableInfo(T_rankine, "T", "K"))
   fig = pl.figure()
   for i in xrange(len(variables)/2):
      fig.add_subplot(len(variables)/2,1,i+1)

   axes = fig.get_axes()
   # Plot
   from variable import get_data, get_name, get_units
   for i in xrange(len(variables)/2):
      # Get the subplot axis
      ax = axes[i]
      x = get_data(distances)
      y = get_data(variables[2*i])
      ax.plot(get_data(x), get_data(y), lw=2, color='black')
      y = get_data(variables[2*i+1])
      ax.plot(get_data(x), get_data(y), '--', lw=2, color='blue')
      ax.set_ylabel(get_name(variables[2*i]) + " (" + get_units(variables[2*i]) + ")")
      ax.get_xaxis().set_visible(False)
      # Highlight upstream and downstream area
      ax.axvline(distances.data[0], color='black', ls='dotted')
      ax.axvline(distances.data[upstream], color='black', ls='dotted')
      ax.axvline(distances.data[len(distances.data)-1], color='black', ls='dotted')
      ax.axvline(distances.data[downstream], color='black', ls='dotted')
      ax.set_xlim([distances.data[0]-2e5,distances.data[len(distances.data)-1]+2e5])
   axes[len(axes)-1].set_xlabel(get_name(distances) + " (" + get_units(distances) + ")")
   axes[len(axes)-1].set_visible(True)

   

#   # Plot the variables (and move back from de hoffmann teller frame)
#   from plot import plot_multiple_variables
#   variables = []
#   #VariableInfo(self, data_array, name="", units="")
#   variables.append( VariableInfo(rho_data, "rho", "m^-3" ) )
#   variables.append(VariableInfo(rho_rankine, "rho", "m^-3") )
#   variables.append(VariableInfo(Vx_data, "Vx", "m/s"))
#   variables.append(VariableInfo(Vx_rankine, "Vx", "m/s"))
#   variables.append(VariableInfo(Vy_data, "Vy", "m/s"))
#   variables.append(VariableInfo(Vy_rankine, "Vy", "m/s"))
##   variables.append(VariableInfo(Bx_data, "Bx", "T"))
##   variables.append(VariableInfo(Bx_rankine, "Bx", "T"))
#   variables.append(VariableInfo(By_data,"By", "T"))
#   variables.append(VariableInfo(By_rankine, "By", "T"))
#   variables.append(VariableInfo(T_data, "T", "K"))
#   variables.append(VariableInfo(T_rankine, "T", "K"))
#
#   numberOfVariables = len(variables)
#
#   fig = plot_variables( distances, variables, figure=[] )

#   variables.append( VariableInfo(rho_data, "rho", "m^-3" ) )
#   variables.append(VariableInfo(Vx_data, "Vx", "m/s"))
#   variables.append(VariableInfo(Vy_data, "Vy", "m/s"))
##   variables.append(VariableInfo(Bx_data, "Bx", "T"))
##   variables.append(VariableInfo(Bx_rankine, "Bx", "T"))
#   variables.append(VariableInfo(By_data,"By", "T"))
#   variables.append(VariableInfo(T_data, "T", "K"))
#
#   numberOfVariables = len(variables)
#
#   fig = plot_variables( distances, variables, figure=[] )
#   variables.append(VariableInfo(rho_rankine, "rho", "m^-3") )
#   variables.append(VariableInfo(Vx_rankine, "Vx", "m/s"))
#   variables.append(VariableInfo(Vy_rankine, "Vy", "m/s"))
#   variables.append(VariableInfo(By_rankine, "By", "T"))
#   variables.append(VariableInfo(T_rankine, "T", "K"))
#   fig = plot_variables( distances, variables, figure=fig )
   # Print rankine conservation laws (should be 1,1,1,1.. if everything ok)
   print compare_rankine(Vx_data[0], Vy_data[0], Bx_data[0], By_data[0], T_data[0], rho_data[0], P_data[0], Vx_data[len(Vx_data)-1], Vy_data[len(Vx_data)-1], Bx_data[len(Vx_data)-1], By_data[len(Vx_data)-1], T_data[len(Vx_data)-1], rho_data[len(Vx_data)-1], P_data[len(Vx_data)-1] )


   pl.show()

   return 

def plot_rankine_nonshifted( vlsvReader, point1, point2 ):
   ''' A function that plots the theoretical rankine-hugoniot jump condition variables for two given points along with the actual values from a given vlsv file

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest along the bow-shock
   :param point2: The second point of interest along the bow-shock

   :returns: pylab figure
   '''
   # Get spatial grid sizes:
   xcells = (int)(vlsvReader.read_parameter("xcells_ini"))
   ycells = (int)(vlsvReader.read_parameter("ycells_ini"))
   zcells = (int)(vlsvReader.read_parameter("zcells_ini"))

   xmin = vlsvReader.read_parameter("xmin")
   ymin = vlsvReader.read_parameter("ymin")
   zmin = vlsvReader.read_parameter("zmin")
   xmax = vlsvReader.read_parameter("xmax")
   ymax = vlsvReader.read_parameter("ymax")
   zmax = vlsvReader.read_parameter("zmax")

   dx = (xmax - xmin) / (float)(xcells)
   dy = (ymax - ymin) / (float)(ycells)
   dz = (zmax - zmin) / (float)(zcells)

   # Get normal vector from point2 and point1
   normal_vector = (point1-point2) / np.linalg.norm(point2 - point1)
   normal_vector = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
   normal_vector = normal_vector * np.array([1,1,0])

   # Read cut-through
   cutthrough = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
   # Get cell ids and distances separately
   cellids = cutthrough[0].data
   distances = cutthrough[1]
   # Read data from the file:
   V_data = vlsvReader.read_variable( "v", cellids=cellids )
   B_data = vlsvReader.read_variable( "B", cellids=cellids )
   T_data = vlsvReader.read_variable( "Temperature", cellids=cellids )
   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )

   # Transform to de hoffmann teller frame:
   ##############################################
   # V_1 x B_1 = 0:
   V_1 = V_data[0]
   B_1 = B_data[0]
   n = normal_vector
   V_HT = np.cross(n, np.cross(V_1, B_1) ) / np.dot(n, B_1)

   # Check whether V_HT is relativistic or not:
   if np.linalg.norm( V_HT ) > 0.05*3.0e8:
      print "WARNING: When transforming to de Hoffmann teller frame encountered speeds over 0.05 speed of light: " + str(np.linalg.norm( V_HT ))

   # Transform to another frame:
   V_data = V_data - V_HT
   ##############################################

   # Get parallel and perpendicular components from the first vector (on one side of the shock):
   Vx_data = np.dot(V_data, normal_vector)
   Vy_data = np.sqrt(np.sum(np.abs(V_data - np.outer(Vx_data, normal_vector))**2,axis=-1))
   Bx_data = np.dot(B_data, normal_vector) 
   By_data = np.sqrt(np.sum(np.abs(B_data - np.outer(Bx_data, normal_vector) )**2,axis=-1))


   # Read V, B, T and rho for point1
#   V = vlsvReader.read_variable( "v", cellids=cellids[0] )
#   B = vlsvReader.read_variable( "B", cellids=cellids[0] )
#   T = vlsvReader.read_variable( "Temperature", cellids=cellids[0] )
#   rho = vlsvReader.read_variable( "rho", cellids=cellids[0] )
   V = V_data[0]
   B = B_data[0]
   T = T_data[0]
   rho = rho_data[0]
   # Get parallel and perpendicular components:
   Vx = np.dot(V, normal_vector)
   Vy = np.linalg.norm(V - Vx * normal_vector)
   Bx = np.dot(B, normal_vector)
   By = np.linalg.norm(B - Bx * normal_vector)

   # Calculate rankine hugoniot jump conditions:
   rankine_conditions = oblique_shock( Vx, Vy, Bx, By, T, rho )

   # Input variables
   Vx_rankine = []
   Vy_rankine = []
   Bx_rankine = []
   By_rankine = []
   T_rankine = []
   rho_rankine = []
   for i in xrange( len(cellids) ):
      if i < 0.5*len(cellids):
         Vx_rankine.append(Vx)
         Vy_rankine.append(Vy)
         Bx_rankine.append(Bx)
         By_rankine.append(By)
         T_rankine.append(T)
         rho_rankine.append(rho)
      else:
         Vx_rankine.append(rankine_conditions[0])
         Vy_rankine.append(rankine_conditions[1])
         Bx_rankine.append(rankine_conditions[2])
         By_rankine.append(rankine_conditions[3])
         T_rankine.append(rankine_conditions[4])
         rho_rankine.append(rankine_conditions[5])

   # Plot the variables:
   from plot import plot_multiple_variables
   variables = []
   #VariableInfo(self, data_array, name="", units="")
   variables.append( VariableInfo(rho_data, "rho", "m^-3" ) )
   variables.append(VariableInfo(rho_rankine, "rho", "m^-3") )
   variables.append(VariableInfo(Vx_data, "Vx", "m/s"))
   variables.append(VariableInfo(Vx_rankine, "Vx", "m/s"))
   variables.append(VariableInfo(Vy_data, "Vy", "m/s"))
   variables.append(VariableInfo(Vy_rankine, "Vy", "m/s"))
#   variables.append(VariableInfo(Bx_data, "Bx", "T"))
#   variables.append(VariableInfo(Bx_rankine, "Bx", "T"))
   variables.append(VariableInfo(By_data,"By", "T"))
   variables.append(VariableInfo(By_rankine, "By", "T"))
   variables.append(VariableInfo(T_data, "T", "K"))
   variables.append(VariableInfo(T_rankine, "T", "K"))

   numberOfVariables = len(variables)

   fig = plot_variables( distances, variables, figure=[] )

   pl.show()

   return fig


def compare_rankine( Vx1, Vy1, Bx1, By1, T1, rho1, P1, Vx2, Vy2, Bx2, By2, T2, rho2, P2 ):
   ''' Function for calculating the Rankine-Hugoniot jump conditions on both sides of the shock. This function is here to verify whether or not given parameters fill the Rankine-Hugoniot jump conditions.

       :param Vx1: Velocity component parallel to the shock normal vector on upstream side of the shock
       :param Vy1: Velocity component perpendicular to the shock normal vector on upstream side of the shock
       :param Bx1: Magnetic field component parallel to the shock normal vector on upstream side of the shock
       :param By1: Magnetic field component perpendicular to the shock normal vector on upstream side of the shock
       :param T1: Temperature on upstream side of the shock
       :param rho1: Density on upstream side of the shock

       :param Vx2: Velocity component parallel to the shock normal vector on downstream side of the shock
       :param Vy2: Velocity component perpendicular to the shock normal vector on downstream side of the shock
       :param Bx2: Magnetic field component parallel to the shock normal vector on downstream side of the shock
       :param By2: Magnetic field component perpendicular to the shock normal vector on downstream side of the shock
       :param T2: Temperature on downstream side of the shock
       :param rho2: Density on downstream side of the shock

       :returns: The Rankine-Hugoniot jump conditions on both sides of the shock in the form of an array
   '''
   # Constants
   mp = 1.67e-27;
   mu0 = 4.0 * np.pi * 1e-7;
   kb = 1.3806505e-23;
   
   # Input variables
   Gamma = 5./3.;
   
   # Calculate other variables
   theta = np.arccos(Bx1/np.sqrt(Bx1**2 + By1**2));
   vA1 = np.sqrt( Bx1**2 / (mp * rho1 * mu0) )
   V1 = np.sqrt( Vx1**2 + Vy1**2 )
   V2 = np.sqrt( Vx2**2 + Vy2**2 )
   B1 = np.sqrt(Bx1**2+By1**2)
   B2 = np.sqrt(Bx2**2+By2**2)
   vs1 = np.sqrt( Gamma * kb * T1 / mp )
   _V1 = np.array([Vx1, Vy1, 0])
   _V2 = np.array([Vx2, Vy2, 0])
   _B1 = np.array([Bx1, By1, 0])
   _B2 = np.array([Bx2, By2, 0])

  

   # Calculate the R-H jump conditions:
   RH = []
   RH.append( rho1*Vx1 / (rho2*Vx2) )
   RH.append( (rho1*Vx1**2+P1+B1**2/(2*mu0)) / (rho2*Vx2**2+P2+B2**2/(2*mu0)) )
   RH.append( ( rho1*Vx1*Vy1-Bx1*By1/mu0 ) / ( rho2*Vx2*Vy2-Bx2*By2/mu0 ) )
   RH.append( ( rho1*Vx1*(0.5*V1**2+Gamma*P1/(rho1*(Gamma-1))) + Vx1*B1**2/mu0 - np.dot(_V1,_B1)*Bx1/mu0 ) / ( rho2*Vx2*(0.5*V2**2+Gamma*P2/(rho2*(Gamma-2))) + Vx2*B2**2/mu0 - np.dot(_V2,_B2)*Bx2/mu0 ) )
   RH.append( Bx1/Bx2 )
   RH.append( (Vx1*By1-Bx1*Vy1) / (Vx2*By2-Bx2*Vy2) )
   return RH

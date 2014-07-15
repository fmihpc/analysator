# Rankine-Hugoniot conditions to determine the right state from the left state of an MHD oblique shock
import numpy as np
import pylab as pl
import sys
from scipy import optimize
from cutthrough import cut_through
from variable import VariableInfo
from plot import plot_variables
from rankine_solver import solve_rankine_SOE

def rankine( vlsvReader, point1, point2 ):
   ''' A function for comparing the values predicted by Rankine-Hugoniot and the values from the data

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest along the bow-shock
   :param point2: The second point of interest along the bow-shock

   :returns: rankine values as a dictionary
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
   n = (point2-point1) / np.linalg.norm(point2 - point1)
   n = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
   n = n * np.array([1,1,0])
   point1_shifted = point1 + 0.5*(point2-point1) + n * (4*dx)
   point2_shifted = point1 + 0.5*(point2-point1) - n * (4*dx)
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
   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )
   P_data = vlsvReader.read_variable( "Pressure", cellids=cellids )

   # Transform to de hoffmann teller frame:
   ##############################################
   # V_1 x B_1 = 0:
   V_1 = V_data[0]
   B_1 = B_data[0]
   V_HT = np.cross(n, np.cross(V_1, B_1) ) / np.dot(n, B_1)

   # Check whether V_HT is relativistic or not:
   if np.linalg.norm( V_HT ) > 0.05*3.0e8:
      print "WARNING: When transforming to de Hoffmann teller frame encountered speeds over 0.05 speed of light: " + str(np.linalg.norm( V_HT ))

   # Transform to another frame:
   V_data = V_data - V_HT
   V_1 = V_data[0]
   B_1 = B_data[0]
   # Check the de-hoffmann teller frame:
   print "HOFFMAN_CHECK: " + str(np.cross(B_1, V_1))
   ##############################################

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

   # Input data: (d = downstream)
   V_d = np.mean(V_data[downstream:len(V_data)-1], axis=0)
   B_d = np.mean(B_data[downstream:len(V_data)-1], axis=0)

   # Get the tangent:
   t = B_u - np.dot(B_u,n) * n
   t = t / np.linalg.norm(t)

   # Get parallel and perpendicular components from the first vector (on one side of the shock):
#   Vx_data = np.dot(V_data, n)
#   Vy_data = np.sqrt(np.sum(np.abs(V_data - np.outer(Vx_data, n))**2,axis=-1))
#   Bx_data = np.dot(B_data, n)
#   By_data = np.sqrt(np.sum(np.abs(B_data - np.outer(Bx_data, n) )**2,axis=-1))
   Vx_data = np.dot(V_data, n)
   Vy_data = np.dot(V_data, t)
   Bx_data = np.dot(B_data, n) 
   By_data = np.dot(B_data, t)


   # Get parallel and perpendicular components:
   Vx_u = np.dot(V_u, n)
   Vy_u = np.dot(V_u, t)
   Bx_u = np.dot(B_u, n)
   By_u = np.dot(B_u, t)

   # Get parallel and perpendicular components:
   Vx_d = np.dot(V_d, n)
   Vy_d = np.dot(V_d, t)
   Bx_d = np.dot(B_d, n)
   By_d = np.dot(B_d, t)


   # Input the rest of the data:
   # Upstream
   rho_u = np.mean(rho_data[0:upstream])
   P_u = np.mean(P_data[0:upstream])
   # Downstream
   rho_d = np.mean(rho_data[downstream:len(V_data)-1])
   P_d = np.mean(P_data[downstream:len(V_data)-1])

   # Calculate rankine hugoniot jump conditions:
   rankine_conditions_dict = oblique_solver( Vx_u, Vy_u, Bx_u, By_u, P_u, rho_u )

   rankine_conditions = []

   for i in rankine_conditions_dict.iteritems():
      rankine_conditions.append(i[1])

   Vx_rankine_d = rankine_conditions[5]
   Vy_rankine_d = rankine_conditions[4]
   Bx_rankine_d = rankine_conditions[1]
   By_rankine_d = rankine_conditions[0]
   rho_rankine_d = rankine_conditions[3]
   P_rankine_d = rankine_conditions[2]


#   # Input Rankine (d = downstream):
#   Vx_rankine_d = rankine_conditions[0]
#   Vy_rankine_d = rankine_conditions[1]
#   Bx_rankine_d = rankine_conditions[2]
#   By_rankine_d = rankine_conditions[3]
#   rho_rankine_d = rankine_conditions[5]
#   P_rankine_d = rankine_conditions[6]

   # Get the differences:
   # Data
   Vx = Vx_d - Vx_u
   Vy = Vy_d - Vy_u
   Bx = Bx_d - Bx_u
   By = By_d - By_u
   rho = rho_d - rho_u
   P = P_d - P_u

   # Rankine
   Vx_rankine = Vx_rankine_d - Vx_u
   Vy_rankine = Vy_rankine_d - Vy_u
   Bx_rankine = Bx_rankine_d - Bx_u
   By_rankine = By_rankine_d - By_u
   rho_rankine = rho_rankine_d - rho_u
   P_rankine = P_rankine_d - P_u


   # Return all the data:
   return {"Vx_rankine_u": Vx_u,"Vy_rankine_u": Vy_u,"Bx_rankine_u": Bx_u,"By_rankine_u": By_u, "rho_rankine_u": rho_u, "P_rankine_u": P_u, "Vx_rankine_d": Vx_rankine_d, "Vy_rankine_d": Vy_rankine_d, "Bx_rankine_d": Bx_rankine_d, "By_rankine_d": By_rankine_d, "rho_rankine_d": rho_rankine_d, "P_rankine_d": P_rankine_d, "Vx_u": Vx_u, "Vy_u": Vy_u, "Bx_u": Bx_u, "By_u": By_u, "rho_u": rho_u, "P_u": P_u,"Vx_d": Vx_d, "Vy_d": Vy_d, "Bx_d": Bx_d, "By_d": By_d, "rho_d": rho_d, "P_d": P_d,"Vx_data": Vx_data, "Vy_data": Vy_data, "Bx_data": Bx_data, "By_data": By_data, "rho_data": rho_data, "P_data": P_data, "cutthrough": cutthrough}


def rotation_matrix_2d( angle ):
   ''' Creates a rotation matrix that can be used to rotate a 2d vector by an angle in the counter-clockwise direction

       :param angle: Rotation angle
       :returns: The rotation matrix
   '''
   return np.array([[np.cos(angle), -1*np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])


def oblique_solver( Vx1, Vy1, Bx1, By1, P1, rho1 ):
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
   vA1 = np.sqrt( (Bx1**2+By1**2) / (mp * rho1 * mu0) )
   V1 = np.sqrt( Vx1**2 + Vy1**2 )
   vs1 = np.sqrt( y * P1 / (mp*rho1) )


   # Print out temperature: T = P/nk
   print "Temperature:" + str(P1/(rho1*kb))

   # Declare the rest of the variables
   Vx2, Vy2, Bx2, By2, T2, P2, rho2, X = sp.symbols('Vx2, Vy2, Bx2, By2, T2, P2, rho2, X')
   V2 = sp.sqrt(Vx2*Vx2+Vy2*Vy2)

   # Solve X:
   x = sp.solve( (V1**2-X*vA1**2)**2*(X*vs1**2+1/2.*V1**2*np.cos(theta)**2*(X*(y-1)-(y+1))) + 1/2.*vA1**2*V1**2*np.sin(theta)**2*X*((y+X*(2-y))*V1**2-X*vA1**2*((y-1)-X*(y-1))) , X)
   print "x: " + str(x)
   # Pick X:
   for i in x:
      if i.as_real_imag()[0] > 0:
         X = i.as_real_imag()[0]

   # Write down the equations
   equation = []
   equation.append( Vx2/Vx1 - 1/X )
   equation.append( rho1/rho2 - 1/X )
   equation.append( Vy2/Vy1 - (V1**2-vA1**2) / (V1**2-X*vA1**2) )
   equation.append( Bx2/Bx1 - 1 )
   equation.append( By2/By1 - (V1**2-vA1**2)*X/(V1**2-X*vA1**2) )
   equation.append( P2/P1 - (X+((y-1)*X*V1**2/(2*vs1**2))*(1-V2**2/(V1**2))) )

   print "alfven:"
   print vA1
   print vs1
   print (X+((y-1)*X*V1**2/(2*vs1**2))*(1-V2**2/(V1**2)))

   return sp.solve(equation, Vx2, Vy2, Bx2, By2, P2, rho2 , dict=False, set=False)

def oblique_shock( Vx1, Vy1, Bx1, By1, T1, P1, rho1 ):
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
   vA1 = np.sqrt( (Bx1**2+By1**2) / (mp * rho1 * mu0) )
   V1 = np.sqrt( Vx1**2 + Vy1**2 )
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

   :returns: the evaluation values as a list: [Vx_rankine/Vx, Vy_rankine/Vy, By_rankine/By, rho_rankine/rho, P_rankine/P]
   '''
   # Get the variables
   rankine_dict = rankine( vlsvReader, point1, point2 )

   Vx_d = rankine_dict["Vx_d"]
   Vx_u = rankine_dict["Vx_u"]
   Vy_d = rankine_dict["Vy_d"]
   Vy_u = rankine_dict["Vy_u"]
   Bx_d = rankine_dict["Bx_d"]
   Bx_u = rankine_dict["Bx_u"]
   By_d = rankine_dict["By_d"]
   By_u = rankine_dict["By_u"]
   rho_d = rankine_dict["rho_d"]
   rho_u = rankine_dict["rho_u"]
   P_d = rankine_dict["P_d"]
   P_u = rankine_dict["P_u"]

   Vx_rankine_d = rankine_dict["Vx_rankine_d"]
   Vy_rankine_d = rankine_dict["Vy_rankine_d"]
   Bx_rankine_d = rankine_dict["Bx_rankine_d"]
   By_rankine_d = rankine_dict["By_rankine_d"]
   rho_rankine_d = rankine_dict["rho_rankine_d"]
   P_rankine_d = rankine_dict["P_rankine_d"]

   # Get the differences:
   # Data
   Vx = Vx_d - Vx_u
   Vy = Vy_d - Vy_u
   Bx = Bx_d - Bx_u
   By = By_d - By_u
   rho = rho_d - rho_u
   P = P_d - P_u

   # Rankine
   Vx_rankine = Vx_rankine_d - Vx_u
   Vy_rankine = Vy_rankine_d - Vy_u
   Bx_rankine = Bx_rankine_d - Bx_u
   By_rankine = By_rankine_d - By_u
   rho_rankine = rho_rankine_d - rho_u
   P_rankine = P_rankine_d - P_u

   # Return the comparisons:
   return [Vx_rankine/Vx, Vy_rankine/Vy, By_rankine/By, rho_rankine/rho, P_rankine/P]

def plot_rankine( vlsvReader, point1, point2, savename="" ):
   ''' A function that plots the theoretical rankine-hugoniot jump condition variables for two given points along with the actual values from a given vlsv file

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest along the bow-shock
   :param point2: The second point of interest along the bow-shock

   :returns: pylab figure
   '''
   # Get the variables
   rankine_dict = rankine( vlsvReader, point1, point2 )

   distances = rankine_dict["cutthrough"][1]

   Vx_d = rankine_dict["Vx_d"]
   Vx_u = rankine_dict["Vx_u"]
   Vy_d = rankine_dict["Vy_d"]
   Vy_u = rankine_dict["Vy_u"]
   Bx_d = rankine_dict["Bx_d"]
   Bx_u = rankine_dict["Bx_u"]
   By_d = rankine_dict["By_d"]
   By_u = rankine_dict["By_u"]
   rho_d = rankine_dict["rho_d"]
   rho_u = rankine_dict["rho_u"]
   P_d = rankine_dict["P_d"]
   P_u = rankine_dict["P_u"]

   Vx_rankine_d = rankine_dict["Vx_rankine_d"]
   Vx_rankine_u = rankine_dict["Vx_rankine_u"]
   Vy_rankine_d = rankine_dict["Vy_rankine_d"]
   Vy_rankine_u = rankine_dict["Vy_rankine_u"]
   Bx_rankine_d = rankine_dict["Bx_rankine_d"]
   Bx_rankine_u = rankine_dict["Bx_rankine_u"]
   By_rankine_d = rankine_dict["By_rankine_d"]
   By_rankine_u = rankine_dict["By_rankine_u"]
   rho_rankine_d = rankine_dict["rho_rankine_d"]
   rho_rankine_u = rankine_dict["rho_rankine_u"]
   P_rankine_d = rankine_dict["P_rankine_d"]
   P_rankine_u = rankine_dict["P_rankine_u"]

   Vx_data = rankine_dict["Vx_data"]
   Vy_data = rankine_dict["Vy_data"]
   Bx_data = rankine_dict["Bx_data"]
   By_data = rankine_dict["By_data"]
   rho_data = rankine_dict["rho_data"]
   P_data = rankine_dict["P_data"]

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
   rho_rankine = []
   P_rankine = []
   for i in xrange( len(rho_data) ):
      if i < upstream+2:
         Vx_rankine.append(Vx_rankine_u)
         Vy_rankine.append(Vy_rankine_u)
         Bx_rankine.append(Bx_rankine_u)
         By_rankine.append(By_rankine_u)
         rho_rankine.append(rho_rankine_u)
         P_rankine.append(P_rankine_u)
      else:
         Vx_rankine.append(Vx_rankine_d)
         Vy_rankine.append(Vy_rankine_d)
         Bx_rankine.append(Bx_rankine_d)
         By_rankine.append(By_rankine_d)
         rho_rankine.append(rho_rankine_d)
         P_rankine.append(P_rankine_d)

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
   variables.append(VariableInfo(P_data, "P", "Pascals"))
   variables.append(VariableInfo(P_rankine, "P", "Pascals"))

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
#
#   numberOfVariables = len(variables)
#
#   fig = plot_variables( distances, variables, figure=[] )
#   variables.append(VariableInfo(rho_rankine, "rho", "m^-3") )
#   variables.append(VariableInfo(Vx_rankine, "Vx", "m/s"))
#   variables.append(VariableInfo(Vy_rankine, "Vy", "m/s"))
#   variables.append(VariableInfo(By_rankine, "By", "T"))
#   fig = plot_variables( distances, variables, figure=fig )

   if savename == "":
      pl.show()
   else:
      pl.savefig(savename)
      pl.close()

   return 

#def plot_rankine_nonshifted( vlsvReader, point1, point2 ):
#   ''' A function that plots the theoretical rankine-hugoniot jump condition variables for two given points along with the actual values from a given vlsv file
#
#   :param vlsvReader: Some open vlsv reader file
#   :param point1: The first point of interest along the bow-shock
#   :param point2: The second point of interest along the bow-shock
#
#   :returns: pylab figure
#   '''
#   # Get spatial grid sizes:
#   xcells = (int)(vlsvReader.read_parameter("xcells_ini"))
#   ycells = (int)(vlsvReader.read_parameter("ycells_ini"))
#   zcells = (int)(vlsvReader.read_parameter("zcells_ini"))
#
#   xmin = vlsvReader.read_parameter("xmin")
#   ymin = vlsvReader.read_parameter("ymin")
#   zmin = vlsvReader.read_parameter("zmin")
#   xmax = vlsvReader.read_parameter("xmax")
#   ymax = vlsvReader.read_parameter("ymax")
#   zmax = vlsvReader.read_parameter("zmax")
#
#   dx = (xmax - xmin) / (float)(xcells)
#   dy = (ymax - ymin) / (float)(ycells)
#   dz = (zmax - zmin) / (float)(zcells)
#
#   # Get normal vector from point2 and point1
#   normal_vector = (point1-point2) / np.linalg.norm(point2 - point1)
#   normal_vector = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
#   normal_vector = normal_vector * np.array([1,1,0])
#
#   # Read cut-through
#   cutthrough = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
#   # Get cell ids and distances separately
#   cellids = cutthrough[0].data
#   distances = cutthrough[1]
#   # Read data from the file:
#   V_data = vlsvReader.read_variable( "v", cellids=cellids )
#   B_data = vlsvReader.read_variable( "B", cellids=cellids )
#   T_data = vlsvReader.read_variable( "Temperature", cellids=cellids )
#   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )
#
#   # Transform to de hoffmann teller frame:
#   ##############################################
#   # V_1 x B_1 = 0:
#   V_1 = V_data[0]
#   B_1 = B_data[0]
#   n = normal_vector
#   V_HT = np.cross(n, np.cross(V_1, B_1) ) / np.dot(n, B_1)
#
#   # Check whether V_HT is relativistic or not:
#   if np.linalg.norm( V_HT ) > 0.05*3.0e8:
#      print "WARNING: When transforming to de Hoffmann teller frame encountered speeds over 0.05 speed of light: " + str(np.linalg.norm( V_HT ))
#
#   # Transform to another frame:
#   V_data = V_data - V_HT
#   ##############################################
#
#   # Get parallel and perpendicular components from the first vector (on one side of the shock):
#   Vx_data = np.dot(V_data, normal_vector)
#   Vy_data = np.sqrt(np.sum(np.abs(V_data - np.outer(Vx_data, normal_vector))**2,axis=-1))
#   Bx_data = np.dot(B_data, normal_vector) 
#   By_data = np.sqrt(np.sum(np.abs(B_data - np.outer(Bx_data, normal_vector) )**2,axis=-1))
#
#
#   # Read V, B, T and rho for point1
##   V = vlsvReader.read_variable( "v", cellids=cellids[0] )
##   B = vlsvReader.read_variable( "B", cellids=cellids[0] )
##   T = vlsvReader.read_variable( "Temperature", cellids=cellids[0] )
##   rho = vlsvReader.read_variable( "rho", cellids=cellids[0] )
#   V = V_data[0]
#   B = B_data[0]
#   T = T_data[0]
#   rho = rho_data[0]
#   # Get parallel and perpendicular components:
#   Vx = np.dot(V, normal_vector)
#   Vy = np.linalg.norm(V - Vx * normal_vector)
#   Bx = np.dot(B, normal_vector)
#   By = np.linalg.norm(B - Bx * normal_vector)
#
#   # Calculate rankine hugoniot jump conditions:
#   rankine_conditions = oblique_shock( Vx, Vy, Bx, By, T, P, rho )
#
#   # Input variables
#   Vx_rankine = []
#   Vy_rankine = []
#   Bx_rankine = []
#   By_rankine = []
#   T_rankine = []
#   rho_rankine = []
#   for i in xrange( len(cellids) ):
#      if i < 0.5*len(cellids):
#         Vx_rankine.append(Vx)
#         Vy_rankine.append(Vy)
#         Bx_rankine.append(Bx)
#         By_rankine.append(By)
#         T_rankine.append(T)
#         rho_rankine.append(rho)
#      else:
#         Vx_rankine.append(rankine_conditions[0])
#         Vy_rankine.append(rankine_conditions[1])
#         Bx_rankine.append(rankine_conditions[2])
#         By_rankine.append(rankine_conditions[3])
#         T_rankine.append(rankine_conditions[4])
#         rho_rankine.append(rankine_conditions[5])
#
#   # Plot the variables:
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
#
#   pl.show()
#
#   return fig



#   # Calculate other variables
#   theta = np.arccos(Bx1/np.sqrt(Bx1**2 + By1**2));
#   vA1 = np.sqrt( Bx1**2 / (mp * rho1 * mu0) )
#   V1 = np.sqrt( Vx1**2 + Vy1**2 )
#   V2 = np.sqrt( Vx2**2 + Vy2**2 )
#   B1 = np.sqrt(Bx1**2+By1**2)
#   B2 = np.sqrt(Bx2**2+By2**2)
#   vs1 = np.sqrt( Gamma * kb * T1 / mp )
#   _V1 = np.array([Vx1, Vy1, 0])
#   _V2 = np.array([Vx2, Vy2, 0])
#   _B1 = np.array([Bx1, By1, 0])
#   _B2 = np.array([Bx2, By2, 0])
#

def compare_rankine( vlsvReader, point1, point2 ):
   ''' Function for calculating the Rankine-Hugoniot jump conditions on both sides of the shock. This function is here to verify whether or not given parameters fill the Rankine-Hugoniot jump conditions.

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
   N = (point2-point1) / np.linalg.norm(point2 - point1)
   N = np.dot(rotation_matrix_2d( 0.5*np.pi ), (point2 - point1)) / np.linalg.norm(point2 - point1)
   N = N * np.array([1,1,0])
   point1_shifted = point1 + 0.5*(point2-point1) + N * (4*dx)
   point2_shifted = point1 + 0.5*(point2-point1) - N * (4*dx)
   point1 = np.array(point1_shifted)
   point2 = np.array(point2_shifted)

   # Read cut-through
   cutthrough = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
   # Get cell ids and distances separately
   cellids = cutthrough[0].data
   distances = cutthrough[1]
   # Read data from the file:
   V = vlsvReader.read_variable( "v", cellids=cellids )
   B = vlsvReader.read_variable( "B", cellids=cellids )
   rho = vlsvReader.read_variable( "rho", cellids=cellids )
   P = vlsvReader.read_variable( "Pressure", cellids=cellids )

   # Get the upstream and downstream area:
   upstream = 0
   downstream = 0
   for i in xrange(len(rho)):
      if rho[i] > np.average(rho):
         upstream = i-3
         downstream = i+2
         break;

   # Input data: (u = upstream)
   V1 = np.mean(V[0:upstream], axis=0)
   B1 = np.mean(B[0:upstream], axis=0)
   rho1 = np.mean(rho[0:upstream])
   P1 = np.mean(P[0:upstream])

   # Input data: (d = downstream)
   V2 = np.mean(V[downstream:len(V)-1], axis=0)
   B2 = np.mean(B[downstream:len(V)-1], axis=0)
   rho2 = np.mean(rho[downstream:len(V)-1])
   P2 = np.mean(P[downstream:len(V)-1])


   # Constants
   mp = 1.67e-27;
   mu0 = 4.0 * np.pi * 1e-7;
   kb = 1.3806505e-23;
   
   # Input variables
   Gamma = 5./3.;


   I1 = P1/((Gamma-1)*rho1)
   I2 = P2/((Gamma-2)*rho2)

   # Calculate the R-H jump conditions:
   RH = []
   RH.append((rho1*np.dot(V1, N)) / (rho2*np.dot(V2, N))) 
   RH.append((rho1*V1*(np.dot(V1, N)) + (P1+np.dot(B1,B1)/(8.*np.pi))*N - (np.dot(B1, N))*B1/(4.*np.pi)) / (rho2*V2*(np.dot(V2, N)) + (P2+np.dot(B2,B2)/(8.*np.pi))*N - (np.dot(B2, N))*B2/(4.*np.pi))) 
   RH.append((np.dot(V1, N)*( (rho1*I1 + rho1*np.dot(V1, V1)/2. + np.dot(B1,B1)/(8.*np.pi)) + (P1+np.dot(B1, B1)/(8.*np.pi)) ) - np.dot(B1,N)*np.dot(B1,V1)/(4.*np.pi)) / (np.dot(V2, N)*( (rho2*I2 + rho2*np.dot(V2, V2)/2. + np.dot(B2,B2)/(8.*np.pi)) + (P2+np.dot(B2, B2)/(8.*np.pi)) ) - np.dot(B2,N)*np.dot(B2,V2)/(4.*np.pi))) 
   RH.append((np.dot(B1, N)) / (np.dot(B2, N))) 
   RH.append((np.cross(N, np.cross(V1, B1) )) / (np.cross(N, np.cross(V2, B2) ))) 

   return RH






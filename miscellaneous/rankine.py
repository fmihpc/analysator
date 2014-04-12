# Rankine-Hugoniot conditions to determine the right state from the left state of an MHD oblique shock



import numpy as np
import pylab as pl
import sys
from scipy import optimize
from cutthrough import cut_through

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
   print "X " + str(X)
   rho2 = rho1 * X
   Vx2 = Vx1 / X
   Vy2 = Vy1 * (V1**2 - vA1**2) / (V1**2 - X * vA1**2)
   
   V2 = Vx2**2 + Vy2**2;
   
   Bx2 = Bx1
   By2 = By1 * (V1**2 - vA1**2) * X / (V1**2 - X * vA1**2)
   P2 = P1 * (X + (Gamma - 1) * X * V1**2 * (1 - V2**2 / V1**2) / (2.0 * vs1**2));
   T2 = P2 / (rho2 * kb)
   return [Vx2, Vy2, Bx2, By2, T2, rho2, P2, X ]


def plot_rankine( vlsvReader, point1, point2 ):
   ''' A function that plots the theoretical rankine-hugoniot jump condition variables for two given points along with the actual values from a given vlsv file

   :param vlsvReader: Some open vlsv reader file
   :param point1: The first point of interest (on one side of the shock)
   :param point2: The second point of interest (on the other side of the shock)

   :returns: pylab figure
   '''
   # Read cut-through
   cut_through = cut_through( vlsvReader=vlsvReader, point1=point1, point2=point2 )
   # Get cell ids and distances separately
   cellids = self.cut_through[0].data
   distances = self.cut_through[1]
   # Read data from the file:
   V_data = vlsvReader.read_variable( "v", cellids=cellids )
   B_data = vlsvReader.read_variable( "B", cellids=cellids )
   T_data = vlsvReader.read_variable( "Temperature", cellids=cellids )
   rho_data = vlsvReader.read_variable( "rho", cellids=cellids )
   # Get normal vector from point2 and point1
   point1 = np.array(point1)
   point2 = np.array(point2)
   normal_vector = (point2-point1) / np.linalg.norm(point2 - point1)

   # Get parallel and perpendicular components:
   Vx_data = np.dot(V_data, normal_vector) 
   Vy_data = np.linalg.norm(V_data - Vx_data * normal_vector) 
   Bx_data = np.dot(B_data, normal_vector) 
   By_data = np.linalg.norm(B_data - Bx_data * normal_vector) 

   # Read V, B, T and rho for point1
   V = self.__vlsvReader.read_variable( "v", cellids=cellids[0] )
   B = self.__vlsvReader.read_variable( "B", cellids=cellids[0] )
   T = self.__vlsvReader.read_variable( "Temperature", cellids=cellids[0] )
   rho = self.__vlsvReader.read_variable( "rho", cellids=cellids[0] )
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
   variables.append(rho_data)
   variables.append(np.array(rho_rankine))
   variables.append(Vx_data)
   variables.append(np.array(Vx_rankine))
   variables.append(Vy_data)
   variables.append(np.array(Vy_rankine))
   variables.append(Bx_data)
   variables.append(np.array(Bx_rankine))
   variables.append(By_data)
   variables.append(np.array(By_rankine))
   variables.append(T_data)
   variables.append(np.array(T_rankine))

   numberOfVariables = len(variables)

   fig = plot_multiple_variables( [distances for i in xrange(numberOfVariables)], variables, figure=[] )
   return fig

















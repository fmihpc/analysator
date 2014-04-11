# Rankine-Hugoniot conditions to determine the right state from the left state of an MHD oblique shock



import numpy as np
import pylab as pl
import sys
from scipy import optimize

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

   X = optimize.fsolve(solver, 0.1)

   solutions.append(X)

   epsilon = 0.1

   test_values = np.arange(0.1,60,0.1)
   for value in test_values:
      result = optimize.fsolve(solver, value)
      if ((result <= X + epsilon and result >= X - epsilon) == False) and (solver(result) >= -1*epsilon and solver(result) <= epsilon):
         X = result
         solutions.append(X)
   
   if len(solutions) != 1:
      print "MORE THAN ONE SOLUTION: "
      print solutions
      return
   
   rho2 = rho1 * X
   Vx2 = Vx1 / X
   Vy2 = Vy1 * (V1**2 - vA1**2) / (V1**2 - X * vA1**2)
   
   V2 = Vx2**2 + Vy2**2;
   
   Bx2 = Bx1
   By2 = By1 * (V1**2 - vA1**2) * X / (V1**2 - X * vA1**2)
   P2 = P1 * (X + (Gamma - 1) * X * V1**2 * (1 - V2**2 / V1**2) / (2.0 * vs1**2));
   T2 = P2 / (rho2 * kb)
   print rho2
   print Vx2
   print Vy2
   print Bx2
   print By2
   print P2
   print T2
   return [Vx2, Vy2, Bx2, By2, T2, rho2, P2, X ]


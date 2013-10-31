'''In this file there are functions that have something to do with fitting data
'''

def subtract_1d_polynomial_fit( x, y ):
   ''' Fits 1d polynomial into the data, subtracts it from the given data and returns the subtracted data.
       :param x           The x-axis
       :param y           Data to manipulate in either variable or raw form
       :returns data from which a 1d polynomial fit has been subtracted
       NOTE: This may be useful for fourier transforms
   '''
   # Fit a polynomial into the data
   from variable import get_data, get_name, get_units
   import numpy as np
   parameters = 2
   def function(args, x, y):
      x = np.array(x)
      y = np.array(y)
      value = args[0] + args[1]*x
      return y - value
   from scipy import optimize
   fit = optimize.leastsq(function, np.ones(parameters), args=(get_data(x), get_data(y)))
   y_fitted = (-1)*function(fit[0], x, 0)
   # Create a new array y2 which has a forced constant amplitude for the (possible) waves:
   y2 = y - y_fitted
   # Return the data
   from output import output_1d
   return output_1d( [y2], [get_name(y)], [get_units(y)] )


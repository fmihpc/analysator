# This class has a function for plotting multiple variables as subplots in one plot figure

import pylab as pl


def plot_multiple_variables( variables_x_list, variables_y_list, figure=[] ):
   ''' Plots multiple variables from the input with pylab
       :param variables_x_list        Some list of variables to be plotted in the x-axis
       :param variables_y_list        Some list of variables to be plotted in the y-axis
       :param figure                  If one wants to plot into an existing figure then the matplotlib figure should be passed as an argument (OPTIONAL)
       :returns a pylab figure with the plot
       Example usage:
       multiplot( [distances, xcoordinates], [rho_values, B_x_values] )
       This would plot rho_values as a function of distance and B_x_values as a function of xcoordinates
       Note:
       multiplot expects variables to be saved in the VariableInfo class
       Notee:
       If for some reason some variable list (x or y) is empty, e.g. variables_x_list = [B_x, [], B_z, rho], then the variable will not be plotted. This can be used if one wants to plot only into certain subplots.
   '''
   import numpy as np
   variables_x_list = np.atleast_1d(variables_x_list)
   variables_y_list = np.atleast_1d(variables_y_list)
   if len(variables_x_list) != len(variables_y_list):
      print "BAD VARIABLE LENGTH"
      return []
   length_of_list = len(variables_x_list)

   if figure != []:
      fig = pl.figure
      if len(fig.get_axes()) < length_of_list:
         for i in (np.arange(length_of_list-len(fig.get_axes())) + len(fig.get_axes())):
            fig.add_subplot(length_of_list,1,i)
   else:
      fig = pl.figure()
      for i in xrange(length_of_list):
         fig.add_subplot(length_of_list,1,i+1)

   axes = fig.get_axes()
   from variable import get_data, get_name, get_units
   for i in xrange(length_of_list):
      
      x = variables_x_list[i]
      y = variables_y_list[i]
      # Check the length of the list
      if (len(np.atleast_1d(x)) == 0) or (len(np.atleast_1d(y)) == 0):
         continue
      ax = axes[i]
      ax.plot(get_data(x), get_data(y))

      if get_units(x) != "":
         ax.set_xlabel(get_name(x) + " [" + get_units(x) + "]")
      else:
         ax.set_xlabel(get_name(x))

      if get_units(y) != "":
         ax.set_ylabel(get_name(y) + " [" + get_units(y) + "]")
      else:
         ax.set_ylabel(get_name(y))

      # Set limits
      xlength = max(get_data(x)) - min(get_data(x))
      ylength = max(get_data(y)) - min(get_data(y))
      ax.set_xlim([min(get_data(x)) - 0.01*xlength, max(get_data(x)) + 0.01*xlength])
      ax.set_ylim([min(get_data(y)) - 0.05*ylength, max(get_data(y)) + 0.05*ylength])

      # Set format
      ax.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))

   return fig




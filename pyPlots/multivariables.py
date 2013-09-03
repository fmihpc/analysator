# This class has a function for plotting multiple variables as subplots in one plot figure

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
   if len(variables_x_list) != len(variables_y_list):
      print "BAD VARIABLE LENGTH"
      return []
   length_of_list = len(variables_x_list)

   if figure != []:
      fig = figure
      if len(fig.get_axes()) < length_of_list:
         for i in (np.arange(length_of_list-len(fig.get_axes())) + len(fig.get_axes())):
            fig.add_subplot(length_of_list,1,i)
   else:
      fig = pl.figure()
      for i in xrange(length_of_list):
         fig.add_subplot(length_of_list,1,i+1)

   axes = fig.get_axes()

   for i in xrange(length_of_list):
      x = variables_x_list[i]
      y = variables_y_list[i]
      # Check the length of the list
      if (len(np.atleast_1d(x)) == 0) or (len(np.atleast_1d(y)) == 0):
         continue
      ax = axes[i]
      ax.plot(x, y)

      if x.units != "":
         ax.set_xlabel(x.name + " [" + x.units + "]")
      else:
         ax.set_xlabel(x.name)

      if y.units != "":
         ax.set_ylabel(y.name + " [" + y.units + "]")
      else:
         ax.set_ylabel(y.name)

      # Set limits
      xlength = max(x) - min(x)
      ylength = max(y) - min(y)
      ax.set_xlim([min(x) - 0.01*xlength, max(x) + 0.01*xlength])
      ax.set_ylim([min(y) - 0.05*ylength, max(y) + 0.05*ylength])

      # Set format
      ax.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))

   return fig




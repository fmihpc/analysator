# This class has a function for plotting multiple variables as subplots in one plot figure


def plot_multiple_variables( variables_x, variables_y ):
   ''' Plots multiple variables from the input with pylab
       :param variables_x        Some list of variables to be plotted in the x-axis
       :param variables_y        Some list of variables to be plotted in the y-axis
       :returns a pylab figure with the plot
       Example usage:
       multiplot( [distances, xcoordinates], [rho_values, B_x_values] )
       This would plot rho_values as a function of distance and B_x_values as a function of xcoordinates
       Note:
       multiplot expects variables to be saved in the VariableInfo class
   '''
   

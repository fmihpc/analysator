import numpy as np
import sys

def vlsv_intpol_points(vlsvReader,points,varlist,operator="pass",interpolation_order=1):
   '''Returns interpolated values of variables at given points
       :param vlsvReader:            Some open VlsvReader
       :type vlsvReader:             :class:`vlsvfile.VlsvReader`
       :param points:                Coordinate points
       :param variable:              Variables
       :param operator:              The operator for the variable, for example "x" for x-component or "magnitude" for magnitude
       :param interpolation_order:   Order of interpolation (0 or 1), defaults to 1
       :returns: A tuple with output: (coordinates,variable_values,header_string)
       .. code-block:: python
          # Example:
          import pytools as pt
          import numpy as np
          f=pt.vlsvfile.VlsvReader(file_name="state00040000.vlsv")
          xmin = f._VlsvReader__xmin
          xmax = f._VlsvReader__xmax
          x = np.linspace(xmin,xmax,100) # points along x-axis
          y = np.zeros(len(x))
          z = np.zeros(len(x))
          points=(np.array([x,y,z])).transpose()
          varlist = ["cellB","n_H+sw_ave"]
          [crd,params,hstr]=pt.calculations.vlsv_intpol_points(f,points,varlist)
          x=crd[:,0]
          y=crd[:,1]
          z=crd[:,2]
          Bx=params[:,0]
          By=params[:,1]
          Bz=params[:,2]
          n=params[:,3]
   '''
   N_points = len(points)
   N_vars = len(varlist)
   if N_vars <= 0:
      print "ERROR: len(varlist) = 0"
      return
   if N_points < 0:
      print "ERROR: len(points) = 0"
      return
   header = "x y z " # header string
   for i in range(N_vars): # loop variable list
      var = varlist[i]
      if vlsvReader.check_variable(var) == False:
         print "ERROR: variable " + var + " does not exist in file " + vlsvReader.file_name
         return
      dim=len(np.atleast_1d(vlsvReader.read_interpolated_variable(var,points[0],operator))) # variable dimensions
      if dim <= 0:
         print "ERROR: bad variable dimension (dim=" + str(dim) + ")"
         return
      values=np.zeros((N_points,dim))
      crds=np.zeros((N_points,3)) # coordinates
      for j in range(N_points):
         if vlsvReader.get_cellid(points[j]) == 0: # coordinates of a point out of domain
            values[j]=np.ones(dim)*np.nan
         elif interpolation_order==1:
            values[j]=vlsvReader.read_interpolated_variable(var,points[j],operator)
         elif interpolation_order==0:
            values[j]=vlsvReader.read_variable(var,vlsvReader.get_cellid(points[j]),operator)
         crds[j]=points[j]
      if i==0:
         res=values
      else:
         res=np.concatenate((res,values),axis=1)
      if dim == 1:
         header = header + var + " "
      else:
         for d in range(dim):
            header = header + var + "(" + str(d) + ") "
   return (crds,res,header)






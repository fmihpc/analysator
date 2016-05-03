import pytools as pt
import numpy as np

def vlsv_intpol_file(file_vlsv,file_orbit,varlist,file_output):
   '''Writes interpolated values of variables in ascii file
       :param file_vlsv:             VLSV file
       :param file_orbit:            Orbit file (columns: x,y,z or t,x,y,z)
       :param varlist:               Variable list
       :param file_output:           Output ascii file (columns: x,y,z,var1,var2,var3,...)
       :returns: none
       .. code-block:: python
          # Example:
          import pytools as pt
          pt.calculations.vlsv_intpol_file("state00040000.vlsv","orbit.dat",["cellB","n_H+sw_ave"],"output.dat")
   '''
   f=pt.vlsvfile.VlsvReader(file_name=file_vlsv)
   points = np.loadtxt(file_orbit,unpack=True).transpose()
   if points.shape[1] == 4:
      points = np.delete(points,0,1) # remove time column
   if points.shape[1] != 3:
      print "ERROR: orbit file must have 3 (x,y,z) or 4 (t,x,y,z) columns"
      return
   [crd,params,hstr]=pt.calculations.vlsv_intpol_points(f,points,varlist)
   d=np.concatenate((crd,params),axis=1)
   np.savetxt(file_output,d,"% 05e",header=hstr,comments="% ")




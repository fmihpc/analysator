from scipy.spatial import Delaunay
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator


class AMRInterpolator(object):
   ''' Class for holding and managing a Delaunay tetrahedralization for cells at refinement interface

   '''

   def __init__(self, reader, method = "RBF", cellids=np.array([1,2,3,4,5],dtype=np.int64)):
      self.__reader = reader
      self.__cellids = np.array(list(set(cellids)),dtype=np.int64)
      # Cannot initialize an empty Delaunay
      #self.__Delaunay = Delaunay(reader.get_cell_coordinates(self.__cellids), incremental = True, qhull_options="QJ Qc Q12")

   def add_cells(self, cells):
      new_cells = [c for c in cells if c not in self.__cellids]
      self.__cellids = np.append(self.__cellids, new_cells)
      
   def get_points(self):
      return self.__reader.get_cell_coordinates(self.__cellids)

   def get_interpolator(self, name, operator, coords, 
                        method="RBF", 
                        methodargs={
                           "RBF":{"neighbors":64},
                           "Delaunay":{"qhull_options":"QJ"}
                           }):
      
      # simplices = self.__Delaunay.find_simplex(coords)
      # pts = self.__Delaunay.points[self.__Delaunay.simplices[simplices].flatten()]

      pts = self.__reader.get_cell_coordinates(self.__cellids)
      vals = self.__reader.read_variable(name, self.__cellids, operator=operator)
      if method == "Delaunay":
         self.__Delaunay = Delaunay(self.reader.get_cell_coordinates(self.__cellids),**methodargs[method])
         return LinearNDInterpolator(self.__Delaunay, vals)
      elif method == "RBF":
         return RBFInterpolator(pts, vals, **methodargs[method]) # Harrison-Stetson number of neighbors


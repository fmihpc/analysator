from scipy.spatial import Delaunay
import numpy as np
from scipy.interpolate import LinearNDInterpolator


class DelaunayWrapper(object):
   ''' Class for holding and managing a Delaunay tetrahedralization for cells at refinement interface

   '''

   def __init__(self, reader):
      self.reader = reader
      self.__cellids = np.array([1,2,3,4,5],dtype=np.int64)
      # Cannot initialize an empty Delaunay
      self.__Delaunay = Delaunay(reader.get_cell_coordinates(self.__cellids), incremental = True, qhull_options="QJ Qc Q12")

   def add_cells(self, cells):
      new_cells = [c for c in cells if c not in self.__cellids]
      self.__Delaunay.add_points(self.reader.get_cell_coordinates(np.array(new_cells)))
      self.__cellids = np.append(self.__cellids, new_cells)

   def get_interpolator(self, name, operator, coords):
      
      simplices = self.__Delaunay.find_simplex(coords)
      pts = self.__Delaunay.points[self.__Delaunay.simplices[simplices].flatten()]

      return LinearNDInterpolator(pts, self.reader.read_variable(name, self.reader.get_cellid(pts), operator=operator))


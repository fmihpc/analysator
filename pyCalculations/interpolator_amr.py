from scipy.spatial import Delaunay
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator
import warnings
from variable import get_data


def f(ksi, fi):
   return (1-ksi[0]) * (1-ksi[1]) * (1-ksi[2]) * fi[0] + \
               ksi[0]  * (1-ksi[1]) * (1-ksi[2]) * fi[1] + \
            (1-ksi[0]) *    ksi[1]  * (1-ksi[2]) * fi[2] + \
               ksi[0]  *    ksi[1]  * (1-ksi[2]) * fi[3] + \
            (1-ksi[0]) * (1-ksi[1]) *    ksi[2]  * fi[4] + \
               ksi[0]  * (1-ksi[1]) *    ksi[2]  * fi[5] + \
            (1-ksi[0]) *    ksi[1]  *    ksi[2]  * fi[6] + \
               ksi[0]  *    ksi[1]  *    ksi[2]  * fi[7]

def df(ksi, fi):
   d0 =  -1 * (1-ksi[1]) * (1-ksi[2]) * fi[0] + \
            1 * (1-ksi[1]) * (1-ksi[2]) * fi[1] + \
         -1 *    ksi[1]  * (1-ksi[2]) * fi[2] + \
            1 *    ksi[1]  * (1-ksi[2]) * fi[3] + \
         -1 * (1-ksi[1]) *    ksi[2]  * fi[4] + \
            1 * (1-ksi[1]) *    ksi[2]  * fi[5] + \
         -1 *    ksi[1]  *    ksi[2]  * fi[6] + \
            1 *    ksi[1]  *    ksi[2]  * fi[7]
   d1 =  (1-ksi[0]) * -1  * (1-ksi[2]) * fi[0] + \
            ksi[0]  * -1  * (1-ksi[2]) * fi[1] + \
         (1-ksi[0]) *  1  * (1-ksi[2]) * fi[2] + \
            ksi[0]  *  1  * (1-ksi[2]) * fi[3] + \
         (1-ksi[0]) * -1  *    ksi[2]  * fi[4] + \
            ksi[0]  * -1  *    ksi[2]  * fi[5] + \
         (1-ksi[0]) *  1  *    ksi[2]  * fi[6] + \
            ksi[0]  *  1  *    ksi[2]  * fi[7]
   d2 =  (1-ksi[0]) * (1-ksi[1]) * -1 * fi[0] + \
            ksi[0]  * (1-ksi[1]) * -1 * fi[1] + \
         (1-ksi[0]) *    ksi[1]  * -1 * fi[2] + \
            ksi[0]  *    ksi[1]  * -1 * fi[3] + \
         (1-ksi[0]) * (1-ksi[1]) *  1 * fi[4] + \
            ksi[0]  * (1-ksi[1]) *  1 * fi[5] + \
         (1-ksi[0]) *    ksi[1]  *  1 * fi[6] + \
            ksi[0]  *    ksi[1]  *  1 * fi[7]
   return np.stack((d0,d1,d2),axis = -1)

def find_ksi(p, verts, tol= 1e-6, maxiters = 200):
   ksi0 = [0.5,0.5,0.5]
   J = df(ksi0, verts)
   print("J", J)
   ksi_n = ksi0
   f_n =  f(ksi_n,verts)-p
   convergence = False
   for i in range(maxiters):
      
      J = df(ksi_n, verts)
      f_n = f(ksi_n,verts)-p
      # print("f(r) = ", f_n)
      # step = np.matmul(np.linalg.inv(J),f_n)
      step = np.linalg.solve(J, -f_n)
      # print("J^-1 f0 = ",step)
      ksi_n1 = step + ksi_n # r_(n+1) 
      ksi_n = ksi_n1
      
      # print(ksi_n1, f(ksi_n1,verts), np.linalg.norm(f(ksi_n1,verts) - p))
      # print("--------------")
      if(np.linalg.norm(f(ksi_n1,verts) - p) < tol):
         convergence = True
         break

   if convergence:
      # print("Converged after ", i, "iterations")
      return ksi_n1
   else:
      warnings.warn("Generalized trilinear interpolation did not converge. Nans inbound.")
      return np.array([np.nan, np.nan, np.nan])
      
class HexahedralTrilinearInterpolator(object):
   ''' Class for doing general hexahedral interpolation, including degenerate hexahedra (...eventually).
   '''

   def __init__(self, pts, vals, **kwargs):
      self.pts = pts
      self.vals = vals
      # Call dual construction from here?
      self.reader = kwargs['reader']
      self.var = kwargs['var']
      self.operator = kwargs['op']

   def __call__(self, pt):
      if(len(pt.shape) == 2):
         vals = []
         for p in pt:
            dual, ksi = self.reader.get_dual(p)
            # print(dual, ksi)
            dual = self.reader._VlsvReader__dual_cells[dual]
            fp = f(ksi, self.reader.read_variable(self.var, np.array(dual), operator=self.operator))
            vals.append(fp)
         return np.array(vals)
      else:
         dual, ksi = self.reader.get_dual(pt)
         dual = self.__dual_cells[dual]
         fp = f(ksi, self.reader.read_variable(self.var, np.array(dual), operator=self.operator))
         return fp


class AMRInterpolator(object):
   ''' Class for holding and managing a Delaunay tetrahedralization for cells at refinement interface

   '''

   def __init__(self, reader, method = "RBF", cellids=np.array([1,2,3,4,5],dtype=np.int64)):
      self.__reader = reader
      self.__cellids = np.array(list(set(cellids)),dtype=np.int64)
      self.duals = {}
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
      methodargs["Trilinear"] = {"reader":self.__reader, "var" : name, "op":operator}

      pts = self.__reader.get_cell_coordinates(self.__cellids)
      vals = self.__reader.read_variable(name, self.__cellids, operator=operator)
      if method == "Delaunay":
         self.__Delaunay = Delaunay(self.reader.get_cell_coordinates(self.__cellids),**methodargs[method])
         return LinearNDInterpolator(self.__Delaunay, vals)
      elif method == "RBF":
         return RBFInterpolator(pts, vals, **methodargs[method]) # Harrison-Stetson number of neighbors
      elif method == "Trilinear":
         return HexahedralTrilinearInterpolator(pts, vals, **methodargs[method])


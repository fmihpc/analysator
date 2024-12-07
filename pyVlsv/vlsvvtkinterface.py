# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2024 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import logging
import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
import os
import sys
import re
import numbers

import vtk.numpy_interface
import vtk.numpy_interface.dataset_adapter
from vlsvreader import VlsvReader

try:
   import vtk
except Exception as e:
   logging.error("VTK import did not succeed")
   raise(e)



'''Extend vtkHyperTreeGrid with a vlsvReader wrapper.

The extension consists of automatic mapping of a vlsv
SpatialGrid to the HTG structure. Functions are provided
for mapping more SpatialGrid data to the HTG structure as
needed. How well does this function as a VTK source already?
'''
class vtkVlsvHyperTreeGrid(vtk.vtkHyperTreeGrid):

   ''' Set the vlsvReader object for the backend,
   otherwise use parent class init, and populate the
   hypertreegrid object with CellID data and the file-
   index information for easier mapping of vlsv data to
   HTG.
   '''
   def __init__(self, vlsvFile, *args, **kwargs):
      self.vlsvreader = vlsvFile #VlsvReader(vlsvFile)
      f = self.vlsvreader
      vtk.vtkHyperTreeGrid.__init__(self, args, kwargs)
      self.Initialize()
      basegridsize = f.get_spatial_mesh_size()
      ext = f.get_fsgrid_mesh_extent()
      nodeCoordinatesX = f.read(tag="MESH_NODE_CRDS_X", mesh='SpatialGrid')
      nodeCoordinatesY = f.read(tag="MESH_NODE_CRDS_Y", mesh='SpatialGrid')
      nodeCoordinatesZ = f.read(tag="MESH_NODE_CRDS_Z", mesh='SpatialGrid')

      self.SetDimensions([len(nodeCoordinatesX),len(nodeCoordinatesY),len(nodeCoordinatesZ)])
      self.SetBranchFactor(2)

      xValues = vtk.vtkDoubleArray()
      xValues.SetNumberOfValues(int(basegridsize[0]+1))
      # print (basegridsizeNodes, nodeCoordinatesX)
      for i,x in enumerate(nodeCoordinatesX):
         xValues.SetValue(i, x)
      self.SetXCoordinates(xValues)

      yValues = vtk.vtkDoubleArray()
      yValues.SetNumberOfValues(int(basegridsize[1]+1))
      for i,x in enumerate(nodeCoordinatesY):
         yValues.SetValue(i, x)
      self.SetYCoordinates(yValues)

      zValues = vtk.vtkDoubleArray()
      zValues.SetNumberOfValues(int(basegridsize[2]+1))
      for i,x in enumerate(nodeCoordinatesZ):
         zValues.SetValue(i, x)
      self.SetZCoordinates(zValues)


      # Here we need to actually init the hypertreegrid,
      # and at least construct the mapping to leaf nodes


      cid_offsets = {0: 0}
      for i in range(1,f.get_max_refinement_level()+1):
         cid_offsets[i] = int(cid_offsets[i-1] + np.prod(np.array([f._VlsvReader__xcells, f._VlsvReader__ycells, f._VlsvReader__zcells]) * 2**(i-1))) # Check next AMR level

      xc = f._VlsvReader__xcells
      yc = f._VlsvReader__ycells
      zc = f._VlsvReader__zcells
      cid_offsets = np.zeros(f.get_max_refinement_level()+1, dtype=np.int64)
      isum = 0
      for i in range(0,f.get_max_refinement_level()):
         isum = isum + 2**(3*i) * xc * yc * zc
         cid_offsets[i+1] = isum

      print(cid_offsets)

      xcells = np.zeros((f.get_max_refinement_level()+1), dtype=np.int64)
      ycells = np.zeros((f.get_max_refinement_level()+1), dtype=np.int64)
      zcells = np.zeros((f.get_max_refinement_level()+1), dtype=np.int64)
      for r in range(f.get_max_refinement_level()+1):
         xcells[r] = f._VlsvReader__xcells*2**(r)
         ycells[r] = f._VlsvReader__ycells*2**(r)
         zcells[r] = f._VlsvReader__zcells*2**(r)
      mins = np.array([f._VlsvReader__xmin,f._VlsvReader__ymin,f._VlsvReader__zmin])

      f._VlsvReader__read_fileindex_for_cellid()

      cidArray = vtk.vtkDoubleArray()
      cidArray.SetName('CellID')
      cidArray.SetNumberOfValues(0)
      self.GetCellData().AddArray(cidArray)

      self.fileIndexArray = vtk.vtkDoubleArray()
      self.fileIndexArray.SetName('fileIndex')
      self.fileIndexArray.SetNumberOfValues(0)
      self.GetCellData().AddArray(self.fileIndexArray)

      f.read_variable_to_cache("CellID")

      cursor = vtk.vtkHyperTreeGridNonOrientedCursor()
      offsetIndex = 0
      numvalues = 0
      self.idxToFileIndex = {}

      def refine_cell(cid, level):
         # ave = np.array([0.0,0.0,0.0,0.0])
         # print(cid, f.get_cell_indices(cid, reflevels = 0), cid in f._VlsvReader__fileindex_for_cellid.keys())
         cell_lengths = np.array([(f._VlsvReader__xmax - f._VlsvReader__xmin)/(xcells[level]),
                                 (f._VlsvReader__ymax - f._VlsvReader__ymin)/(ycells[level]),
                                 (f._VlsvReader__zmax - f._VlsvReader__zmin)/(zcells[level])])
         # print(f.get_cell_coordinates(cid-1))
         cellind = f.get_cell_indices(cid, reflevels = level)

         # print(mins, cellind, cell_lengths)
         cellcoordinates = mins + (cellind + 0.25)*cell_lengths
         # print(cellind,cellcoordinates)

         cell_lengths = np.array([(f._VlsvReader__xmax - f._VlsvReader__xmin)/(xcells[level+1]),
                                 (f._VlsvReader__ymax - f._VlsvReader__ymin)/(ycells[level+1]),
                                 (f._VlsvReader__zmax - f._VlsvReader__zmin)/(zcells[level+1])])
         # print(cell_lengths)
         cellind = np.array((cellcoordinates - mins)/cell_lengths, dtype = np.int32)
         # print(cellind)
         # print(xcells[level],xcells[level+1])
         # print(xcells)
         # sys.exit()
         # delta = [[0,0,0],[0,0,1],[0,1,0],[0,1,1],
         #          [1,0,0],[1,0,1],[1,1,0],[1,1,1]]
         delta = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],
                  [0,0,1],[1,0,1],[0,1,1],[1,1,1]]

         for c in range(8):
               cursor.ToChild(c)  # cell 0/0
               # print(cellind, delta[c])
               ci = cellind + delta[c]
               # print(ci)
               cellid = int(cid_offsets[level+1] + ci[0] + 2**(level+1)*ci[1]*xc + 4**(level+1)*ci[2]*xc*yc + 1)
               # print(f.get_cell_indices(cid, reflevels = level))
               if (cellid in f._VlsvReader__fileindex_for_cellid.keys()):
                  # Assigns an index for the value of the cell pointed to by the current cursor state 
                  idx = cursor.GetGlobalNodeIndex()
                  self.GetCellData().SetActiveScalars('CellID')
                  cidArray.InsertTuple1(idx, cellid)
                  self.GetCellData().SetActiveScalars('fileIndex')
                  self.fileIndexArray.InsertTuple1(idx, f._VlsvReader__fileindex_for_cellid[cellid])
                  # numvalues += 1
                  self.idxToFileIndex[idx] = f._VlsvReader__fileindex_for_cellid[cellid]
                  # ave += np.array([var_x(cellid, 'proton/vg_v'),*var(cellid, 'vg_b_vol')])
               else:
                  # print(cellid)
                  idx = cursor.GetGlobalNodeIndex()
                  cursor.SubdivideLeaf()  # cell 0
                  # print(cellid, "not found, refining again")
                  refine_cell(cellid, level+1)
                  # ave += refval
                  # htg.GetCellData().SetActiveScalars('vg_v')
                  #scalarArray.InsertTuple1(idx, refval[0])
                  # htg.GetCellData().SetActiveVectors('vg_b_vol')
                  #vectorArray.InsertTuple3(idx, *refval[1:])

               cursor.ToParent()  # cell 0
               idx = cursor.GetGlobalNodeIndex()


         return #ave/8


      def handle_cell(cid):
         # Leaf cell, just add
         if (cid in f._VlsvReader__fileindex_for_cellid.keys()):
               # Assigns an index for the value of the cell pointed to by the current cursor state 
               idx = cursor.GetGlobalNodeIndex()
               self.GetCellData().SetActiveScalars('CellID')
               cidArray.InsertTuple1(idx, cid)
               self.GetCellData().SetActiveScalars('fileIndex')
               self.fileIndexArray.InsertTuple1(idx, f._VlsvReader__fileindex_for_cellid[cid])
               # numvalues += 1
               self.idxToFileIndex[idx] = f._VlsvReader__fileindex_for_cellid[cid]
         # Cell not in list -> must be refined (assume cid is in basegrid)
         else:
               idx = cursor.GetGlobalNodeIndex()

               # print(cid , "not found, refining")
               cursor.SubdivideLeaf()  # cell 0
               refine_cell(cid, 0)
               # self.GetCellData().SetActiveScalars('vg_v')
               #scalarArray.InsertTuple1(idx, ave[0])
               # self.GetCellData().SetActiveVectors('vg_b_vol')
               #vectorArray.InsertTuple3(idx, *ave[1:])

               # pass
               # idx = cursor.GetGlobalNodeIndex()
               # scalarArray.InsertTuple1(idx, vardummy(cid))

      nc = 0
      for basegridIdx in range(int(np.prod(basegridsize))):
         cid = basegridIdx +1
         # Set cursor on cell #0 and set the start index
         self.InitializeNonOrientedCursor(cursor, basegridIdx, True)
         cursor.SetGlobalIndexStart(offsetIndex)  # cell 0
         # print("handling ",cid)
         handle_cell(cid)

         offsetIndex += cursor.GetTree().GetNumberOfVertices()
         nc += 1
         if cid > np.prod(basegridsize):
            print("tried to go over")
            break

   '''This function adds one SpatialGrid variable from the reader object and maps
   that to the hypertreegrid object. Variable vector sizes of 1,2,3,4,9 supported.
   '''
   def addArrayFromVlsv(self, varname):
      array = vtk.vtkDoubleArray()
      # cidArray2.DeepCopy(fileIndexArray)
      data = self.vlsvreader.read_variable(varname)
      
      test_var = data[0]
      if np.isscalar(test_var):
         varlen = 1
      else:
         varlen = len(test_var)
      # varlen = data[:,np.newaxis].shape[1]
      print(varlen)
      array.SetNumberOfComponents(varlen)
      array.SetNumberOfTuples(self.fileIndexArray.GetNumberOfTuples())
      array.SetName(varname)
      # print(cidArray2)
      # newdata = np.array([ciddata[self.idxToFileIndex[idx]] for idx in self.idxToFileIndex.keys()])
      if varlen == 1:
         for idx,fileIndex in self.idxToFileIndex.items():
            array.SetTuple1(idx, data[fileIndex])
      elif varlen == 2:
         for idx,fileIndex in self.idxToFileIndex.items():
            array.SetTuple2(idx, *data[fileIndex])
      elif varlen == 3:
         for idx,fileIndex in self.idxToFileIndex.items():
            array.SetTuple3(idx, *data[fileIndex])
      elif varlen == 4:
         for idx,fileIndex in self.idxToFileIndex.items():
            array.SetTuple4(idx, *data[fileIndex])
      elif varlen == 9:
         for idx,fileIndex in self.idxToFileIndex.items():
            array.SetTuple9(idx, *data[fileIndex])
      else:
         raise RuntimeError("No vtk SetTuple wrapper function for varlen = " + str(varlen))
      
      
      self.GetCellData().AddArray(array)
      
      # vtk.numpy_interface.numpy_to_vtk(newdata)
      # cidArray2.SetData(newdata, len(self.idxToFileIndex), 1)

      
      # ... we do not have intermediate node data stored.
      # it could be calculated while constructing, but that
      # defeats the purpose (limit file read/memory use to
      # a threshold). Also, if we can just map directly to
      # the file/a buffer, we do not need to do fancy recursive
      # traversal downsampling.
   

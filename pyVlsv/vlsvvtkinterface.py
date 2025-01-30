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
import numpy as np
import cProfile
import os.path
import pickle

import vtk.numpy_interface
import vtk.numpy_interface.dataset_adapter
import vtk.vtkCommonColor
import vtkmodules.vtkCommonColor
import vtkmodules.vtkCommonColor
import vtkmodules.vtkCommonColor
import pytools as pt
try:
   import vtk
except Exception as e:
   logging.error("VTK import did not succeed")
   raise(e)
import time
# from numba import jit
'''Extend vtkHyperTreeGrid.
This class is not used right now: makes probably more sense to have
most of the features in the reader class below, esp. since the
vtkHyperTreeGridSource is the fast way of constructing the HTG object.

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
      vtk.vtkHyperTreeGrid.__init__(self)#, args, kwargs)
      self.Initialize()
      self.metafile = self.vlsvreader.file_name[:-5]+"_meta"

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

      self.__descriptor = ""

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

      # print(cid_offsets)
      max_ref_level = f.get_max_refinement_level()
      xcells = np.zeros((max_ref_level+1), dtype=np.int64)
      ycells = np.zeros((max_ref_level+1), dtype=np.int64)
      zcells = np.zeros((max_ref_level+1), dtype=np.int64)
      for r in range(max_ref_level+1):
         xcells[r] = xc*2**(r)
         ycells[r] = yc*2**(r)
         zcells[r] = zc*2**(r)
      mins = np.array([f._VlsvReader__xmin,f._VlsvReader__ymin,f._VlsvReader__zmin])

      f._VlsvReader__read_fileindex_for_cellid()

      # cidArray = vtk.vtkDoubleArray()
      # cidArray.SetName('CellID')
      # cidArray.SetNumberOfValues(0)
      # self.GetCellData().AddArray(cidArray)


      f.read_variable_to_cache("CellID")

      cursor = vtk.vtkHyperTreeGridNonOrientedCursor()
      self.idxToFileIndex = {}

      xmin,xmax,ymin,ymax,zmin,zmax = f._VlsvReader__xmin,f._VlsvReader__xmax,f._VlsvReader__ymin,f._VlsvReader__ymax,f._VlsvReader__zmin,f._VlsvReader__zmax,
      
      fileindex_for_cellid = f._VlsvReader__fileindex_for_cellid
      delta = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],
               [0,0,1],[1,0,1],[0,1,1],[1,1,1]],dtype=np.int32)
      
      # @jit(nopython=True)
      def get_cell_indices(cellids, reflevels):
         ''' Single-cellid specialization with external constants.
         Returns a given cell's indices as a numpy array

         :param cellid:            The cell's ID, numpy array
         :param reflevel:          The cell's refinement level in the AMR
         :returns: a numpy array with the coordinates

         .. seealso:: :func:`get_cellid`

         .. note:: The cell ids go from 1 .. max not from 0
         '''

         # # Calculating the index of the first cell at this reflevel
         # index_at_reflevel = np.zeros(max_refinement_level+1, dtype=np.int64)
         # isum = 0
         # for i in range(0,max_refinement_level):
         #    isum = isum + 2**(3*i) * xcells * ycells * zcells
         #    index_at_reflevel[i+1] = isum

         # Get cell indices
         # print(cellids, reflevels, cid_offsets)
         
         cellids = cellids - 1 - cid_offsets[reflevels]
         cellindices = np.full(3, -1,dtype=np.int32)
         cellindices[0] = (cellids)%(xcells[reflevels])
         cellindices[1] = ((cellids)//(xcells[reflevels]))%(ycells[reflevels])
         cellindices[2] = (cellids)//(xcells[reflevels]*ycells[reflevels])

         # print(cellindices, f.get_cell_indices(cellids,reflevels))

         return cellindices

      cell_lengths_levels = []
      for level in range(max_ref_level+1):
          cell_lengths_levels.append(np.array([(xmax - xmin)/(xcells[level]),
                                               (ymax - ymin)/(ycells[level]),
                                               (zmax - zmin)/(zcells[level])],np.int32))
          
         
      # @jit(nopython=True)
      def children(cid, level):
         # cellind = get_cell_indices(cid, level) # get_cellind here for compilation
         cellids = cid - 1 - cid_offsets[level]
         cellind = np.full(3, -1,dtype=np.int32)
         cellind[0] = (cellids)%(xcells[level])
         cellind[1] = ((cellids)//(xcells[level]))%(ycells[level])
         cellind[2] = (cellids)//(xcells[level]*ycells[level])

         # cellcoordinates = mins + (cellind + 0.25)*cell_lengths_levels[level]
         # cellind = np.int64(np.divide((cellcoordinates - mins),cell_lengths_levels[level+1]))
         cellind = cellind*2
         out = np.zeros((8,),dtype=np.int64)
         # c = np.array((0,1,2,3,4,5,6,7))
         out[:] = cid_offsets[level+1] + (cellind[0] + delta[:,0]) + xcells[level+1]*(cellind[1] + delta[:,1]) + (cellind[2] + delta[:,2])*xcells[level+1]*ycells[level+1] + 1
         # out2 = [int(cid_offsets[level+1] + (cellind + delta[c])[0] + xcells[level+1]*(cellind + delta[c])[1] + (cellind + delta[c])[2]*xcells[level+1]*ycells[level+1] + 1) for c in range(8)]
         # print(out, out2)
         return out           


      def refine_cell(cid, level):

         # cellind = get_cell_indices(cid, level)
         # cellcoordinates = mins + (cellind + 0.25)*cell_lengths_levels[level]

         # cellind = np.array(np.divide((cellcoordinates - mins),cell_lengths_levels[level+1]), dtype = np.int32)

         childs = children(cid,level)
         for c,cellid in enumerate(childs):
               # print(c)
               cursor.ToChild(c)  # cell 0/0
               # print(cellind, delta[c])
               # ci = cellind + delta[c]
               # print(ci)
               # cellid = int(cid_offsets[level+1] + ci[0] + xcells[level+1]*ci[1] + ci[2]*xcells[level+1]*ycells[level+1] + 1)
               # print(f.get_cell_indices(cid, reflevels = level))
               if (cellid in fileindex_for_cellid.keys()):
                  # Assigns an index for the value of the cell pointed to by the current cursor state 
                  idx = cursor.GetGlobalNodeIndex()
                  # print(idx, cellid)
                  # self.GetCellData().SetActiveScalars('CellID')
                  # cidArray.InsertTuple1(idx, cellid)
                  # self.GetCellData().SetActiveScalars('fileIndex')
                  # self.fileIndexArray.InsertTuple1(idx, fileindex_for_cellid[cellid])
                  # numvalues += 1
                  self.idxToFileIndex[idx] = fileindex_for_cellid[cellid]
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
               # idx = cursor.GetGlobalNodeIndex()


         return #ave/8

      def handle_cell(cellid):
         # Leaf cell, just add
         if (cellid in fileindex_for_cellid.keys()):
               # Assigns an index for the value of the cell pointed to by the current cursor state 
               idx = cursor.GetGlobalNodeIndex()
               # self.GetCellData().SetActiveScalars('CellID')
               # print(idx, cid)
               # cidArray.InsertTuple1(idx, cellid)
               # self.GetCellData().SetActiveScalars('fileIndex')
               # self.fileIndexArray.InsertTuple1(idx, fileindex_for_cellid[cellid])
               # numvalues += 1
               self.idxToFileIndex[idx] = fileindex_for_cellid[cellid]
         # Cell not in list -> must be refined (assume cid is in basegrid)
         else:
               idx = cursor.GetGlobalNodeIndex()

               # print(cid , "not found, refining")
               cursor.SubdivideLeaf()  # cell 0
               refine_cell(cellid, 0)
               # self.GetCellData().SetActiveScalars('vg_v')
               #scalarArray.InsertTuple1(idx, ave[0])
               # self.GetCellData().SetActiveVectors('vg_b_vol')
               #vectorArray.InsertTuple3(idx, *ave[1:])

               # pass
               # idx = cursor.GetGlobalNodeIndex()
               # scalarArray.InsertTuple1(idx, vardummy(cid))


      def loop():
         offsetIndex = 0

         for basegridIdx in range(int(np.prod(basegridsize))):
            cid = basegridIdx +1
            # Set cursor on cell #0 and set the start index
            self.InitializeNonOrientedCursor(cursor, basegridIdx, True)
            cursor.SetGlobalIndexStart(offsetIndex)  # cell 0
            # print("handling ",cid)
            handle_cell(cid)

            offsetIndex += cursor.GetTree().GetNumberOfVertices()
            if cid > np.prod(basegridsize):
               print("tried to go over")
               break

      import pstats
      from pstats import SortKey
      
      # with cProfile.Profile() as pr:
      reinit = False

      try:
         if reinit:
            raise Exception("reinit = True")
         with open(self.metafile,'rb') as mfile:
            data = pickle.load(mfile)
            # print(data)
            self.idxToFileIndex = data['idxToFileIndexMap']
            self.__descriptor = data['descr']
      except Exception as e:
         print("Re-initializing HTG, no metadata accessible because of ", e)
         reinit = True


      if False:
         t0 = time.time()
         print("Walking HTG")
         loop()
         print("loop done")
         writer = vtk.vtkXMLHyperTreeGridWriter()
         writer.SetFileName("output_FID_1.htg")
         self.fileIndexArray.SetNumberOfValues(1)
         self.fileIndexArray.SetNumberOfTuples(len(self.idxToFileIndex))
         for idx,fileIndex in self.idxToFileIndex.items():
            self.fileIndexArray.InsertTuple1(idx, fileIndex)
         writer.SetInputData(self)
         writer.Write()
         

      # if not reinit:
      #    for idx,fileIndex in self.idxToFileIndex.items():
      #       self.fileIndexArray.InsertTuple1(idx, fileIndex)
      # else:

            
         # ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE).reverse_order()
         # ps.print_stats()
         # self >> writer
         # writer.Write()
         
      

         t1 = time.time()
         
         # ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE).reverse_order()
         # ps.print_stats()
         # print("Actual time elapsed:", time.time()-t0, "of which", time.time()-t1, "spent adding an array")
      
      print("pre-calling the jit function to pre-compile")
      children(1,0)
      
      if reinit:
         self.idxToCellID = {}
         self.idxToFileIndex = {}
         t0 = time.time()
         from io import StringIO
         descr = StringIO()
         print("Building descriptor")
         subdivided = [[] for l in range(max_ref_level+1)]
         idx = 0
         # with cProfile.Profile() as pr:
         for c in range(1,int(np.prod(basegridsize))+1):
            if c in fileindex_for_cellid.keys():
               # self.__descriptor += "."
               descr.write(".")
               self.idxToFileIndex[idx] = fileindex_for_cellid[c]
               self.idxToCellID[idx] = c
            else:
               # self.__descriptor += "R"
               descr.write("R")
               subdivided[0].append(c)
            idx += 1

         # self.__descriptor += "|"
         descr.write("|")
         for l in range(1,max_ref_level+1):
            for c in subdivided[l-1]:
               for child in children(c, l-1):
                  if child in fileindex_for_cellid.keys():
                     # self.__descriptor += "."
                     descr.write(".")
                     self.idxToFileIndex[idx] = fileindex_for_cellid[child]
                     self.idxToCellID[idx] = c
                  else:
                     # self.__descriptor += "R"
                     descr.write("R")
                     subdivided[l].append(child)
                  idx += 1
            if l < max_ref_level:
               # self.__descriptor += "|"
               descr.write("|")
         self.__descriptor = descr.getvalue()
            # ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
            # ps.print_stats()
         print("Actual time spent:",time.time()-t0)
         # self.SetDimensions([len(nodeCoordinatesX),len(nodeCoordinatesY),len(nodeCoordinatesZ)])
         # self.SetBranchFactor(2)
         print("Building with vtkHyperTreeGridSource and descriptor, len ", len(self.__descriptor))

      if reinit:
         with open(self.metafile,'wb') as mfile:
            pickle.dump({"descr":self.__descriptor, "idxToFileIndexMap":self.idxToFileIndex}, mfile)

      if True:
         with cProfile.Profile() as pr:
            # src = vtk.vtkHyperTreeGridSource(max_depth = max_ref_level+1,
            #                                  dimensions=(int(basegridsize[0]+1),int(basegridsize[1]+1),int(basegridsize[2]+1)),
            #                                  grid_scale = (f._VlsvReader__dx,f._VlsvReader__dy,f._VlsvReader__dz),
            #                                  branch_factor = 2,)
            src = vtk.vtkHyperTreeGridSource()
            src.SetMaxDepth(max_ref_level+1)
            src.SetDimensions(int(basegridsize[0]+1),int(basegridsize[1]+1),int(basegridsize[2]+1))
            src.SetGridScale(f._VlsvReader__dx,f._VlsvReader__dy,f._VlsvReader__dz)
            src.SetBranchFactor(2)
            self.__descriptor = src.ConvertDescriptorStringToBitArray(self.__descriptor)
            src.SetDescriptorBits(self.__descriptor)
            # src.SetDescriptor(self.__descriptor)

            src.Update()

            htg = src.GetHyperTreeGridOutput()
            htg.SetXCoordinates(xValues)
            htg.SetYCoordinates(yValues)
            htg.SetZCoordinates(zValues)
   

            # ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE).reverse_order()
            # ps.print_stats()

            # ahtg = src.GetOutput()
            print("foo")
            writer = vtk.vtkXMLHyperTreeGridWriter()
            writer.SetFileName("output_FID_2.htg")
            writer.SetInputData(htg)
            
            print("foois")
            # src >> writer
            writer.Write()
            print("Stats after write")
            ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE).reverse_order()
            ps.print_stats()
            print("foo2")
            self.CopyStructure(htg)
            
            
            
            # #    # Hyper tree grid to unstructured grid filter.
            # htg2ug = vtk.vtkHyperTreeGridToUnstructuredGrid()
            # # # # data = src.RequestData()
            

            # # # shrink = vtk.vtkShrinkFilter(shrink_factor=0.5)

            # mapper = vtk.vtkDataSetMapper()

            # # # # src >> htg2ug >> shrink >> mapper
            # src >> htg2ug >> mapper
            # # htg2ug.SetInputData(src.GetHyperTreeGridOutput())
            # # htg2ug.SetInputData(self)
            # # mapper.SetInputData(htg2ug.GetOutput())

            # actor = vtk.vtkActor()#(mapper=mapper)
            # actor.SetMapper(mapper)
            
            # #actor.property.diffuse_color = vtk.colors.GetColor3d('Burlywood')
            # import vtkmodules.vtkCommonColor as vtkCommonColor
            # colors = vtkCommonColor.vtkNamedColors()
            
            # actor.property.diffuse_color = colors.GetColor3d('Burlywood')

            # renderer = vtk.vtkRenderer()#background=vtk.colors.GetColor3d('SlateGray'))
            # renderer.background = colors.GetColor3d('SlateGray')
            # render_window = vtk.vtkRenderWindow()#(size=(640, 480), window_name='HyperTreeGridSource')
            # render_window.SetSize(640,480)
            # render_window.AddRenderer(renderer)
            # interactor = vtk.vtkRenderWindowInteractor()
            # interactor.render_window = render_window

            # renderer.AddActor(actor)
            # renderer.ResetCamera()
            # renderer.GetActiveCamera().Azimuth(150)
            # renderer.GetActiveCamera().Elevation(30)
            # renderer.ResetCameraClippingRange()

            # render_window.Render()
            # interactor.Start()
            # time.sleep(20)
            # render_window.Close()


      # self.
      # if not reinit:
      self.fileIndexArray = vtk.vtkDoubleArray()
      self.fileIndexArray.SetName('fileIndex')
      
      self.GetCellData().AddArray(self.fileIndexArray)

      # self.fileIndexArray.SetNumberOfValues(1)
      # self.fileIndexArray.SetNumberOfTuples(len(self.idxToFileIndex))

      for idx,fileIndex in self.idxToFileIndex.items():
         self.fileIndexArray.InsertTuple1(idx, fileIndex)
      # self.addArrayFromVlsv("CellID")

      # print(ahtg)


   def findVariablesFromVlsv(self, getReducers = False):

      vars = self.vlsvreader.get_variables()
      if getReducers:
         reducers = self.vlsvreader.get_reducers()
         vars_set = set(vars)
         vars_set.update(reducers)
         vars = list(vars_set)

      # for var in vars:
      #    array = vtk.vtkDoubleArray()
      #    array.SetName(var)
      #    self.GetCellData().AddArray(array)

      return vars

   

   '''This function adds one SpatialGrid variable from the reader object and maps
   that to the hypertreegrid object. Variable vector sizes of 1,2,3,4,9 supported.
   '''
   def addArrayFromVlsv(self, varname):

      # Do not re-add an already existing array
      if self.GetCellData().HasArray(varname):
         print("skipped existing array")
         return True

      array = vtk.vtkDoubleArray()
      # cidArray2.DeepCopy(fileIndexArray)
      data = self.vlsvreader.read_variable(varname)
      
      test_var = data[0]
      if test_var is None:
         logging.warning("varname " + varname + "returned None; skipping")
         return False
      if np.isscalar(test_var):
         varlen = 1
      elif test_var.ndim == 1:
         varlen = len(test_var)         
      elif test_var.ndim > 1:
         varlen = np.prod(test_var.shape)
      else:
         logging.warning("Weird output from " + varname)
         return False
         # varlen = 0
         # print("test var is ", test_var, 'for', varname)

      
      # varlen = data[:,np.newaxis].shape[1]
      # print(varlen)
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
            array.SetTuple9(idx, *np.reshape(data[fileIndex],(9)))
      else:
         raise RuntimeError("No vtk SetTuple wrapper function for varlen = " + str(varlen))
      
      
      self.GetCellData().AddArray(array)
      return True
      
      # vtk.numpy_interface.numpy_to_vtk(newdata)
      # cidArray2.SetData(newdata, len(self.idxToFileIndex), 1)

      
      # ... we do not have intermediate node data stored.
      # it could be calculated while constructing, but that
      # defeats the purpose (limit file read/memory use to
      # a threshold). Also, if we can just map directly to
      # the file/a buffer, we do not need to do fancy recursive
      # traversal downsampling.
   
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase


def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

# Should implement these, more or less:
# https://vtk.org/doc/nightly/html/classvtkXMLHyperTreeGridReader.html

'''Experimental Python wrapper to expose VLSV files as VTK. SpatialGrid supported for now.

This will create an internal VlsvReader object and cache the HTG descriptor and
a corresponding file layout mapping to disk for repeated accesses:
building the descriptor takes ~30s (with Python; maybe the whole loop could be @jit-ted)

This class functions like a VTK source. Calling `addArrayFromVlsv` is necessary
to load and map data from the VLSV to HTG.
'''
class VlsvVtkReader(VTKPythonAlgorithmBase):
   def __init__(self):
      VTKPythonAlgorithmBase.__init__(self, nInputPorts = 0, nOutputPorts=1, outputType='vtkHyperTreeGrid')
      self.__FileName = None
      self.__htg = None
      self.__reader = None
      self.__cellarrays = []
      self._arrayselection = vtk.vtkDataArraySelection()
      self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))

   def SetFileName(self, filename):
      if filename != self.__FileName:
         self.Modified()
         self.__FileName = filename
         if self.__FileName is not None:
            self.__reader = pt.vlsvfile.VlsvReader(self.__FileName)
            fn = os.path.basename(self.__FileName)
            self.__metafile = os.path.join(os.path.dirname(self.__FileName),"vlsvmeta",fn[:-5]+"_meta")

   def GetFileName(self):
      return self.__FileName

   def buildDescriptor(self):
      f = self.__reader
      f._VlsvReader__read_fileindex_for_cellid()
      fileindex_for_cellid = f._VlsvReader__fileindex_for_cellid
      xc = f._VlsvReader__xcells
      yc = f._VlsvReader__ycells
      zc = f._VlsvReader__zcells
      max_ref_level = f.get_max_refinement_level()

      cid_offsets = np.zeros(max_ref_level+1, dtype=np.int64)
      isum = 0
      for i in range(0,max_ref_level):
         isum = isum + 2**(3*i) * xc * yc * zc
         cid_offsets[i+1] = isum


      xcells = np.zeros((max_ref_level+1), dtype=np.int64)
      ycells = np.zeros((max_ref_level+1), dtype=np.int64)
      zcells = np.zeros((max_ref_level+1), dtype=np.int64)
      for r in range(max_ref_level+1):
         xcells[r] = xc*2**(r)
         ycells[r] = yc*2**(r)
         zcells[r] = zc*2**(r)
      delta = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],
               [0,0,1],[1,0,1],[0,1,1],[1,1,1]],dtype=np.int32)


      idxToFileIndex = {}
      from io import StringIO
      descr = StringIO()
      print("Building descriptor")
      subdivided = [[] for l in range(max_ref_level+1)]
      idx = 0
      for c in range(1,int(np.prod(f.get_spatial_mesh_size()))+1):
         if c in fileindex_for_cellid.keys():
            descr.write(".")
            idxToFileIndex[idx] = fileindex_for_cellid[c]
         else:
            descr.write("R")
            subdivided[0].append(c)
         idx += 1

      # @jit(nopython=True)
      def children(cid, level):
         # cellind = get_cell_indices(cid, level) # get_cellind here for compilation
         cellids = cid - 1 - cid_offsets[level]
         cellind = np.full(3, -1,dtype=np.int32)
         cellind[0] = (cellids)%(xcells[level])
         cellind[1] = ((cellids)//(xcells[level]))%(ycells[level])
         cellind[2] = (cellids)//(xcells[level]*ycells[level])

         cellind = cellind*2
         out = np.zeros((8,),dtype=np.int64)
         out[:] = cid_offsets[level+1] + (cellind[0] + delta[:,0]) + xcells[level+1]*(cellind[1] + delta[:,1]) + (cellind[2] + delta[:,2])*xcells[level+1]*ycells[level+1] + 1
         return out           

      
      descr.write("|")
      for l in range(1,max_ref_level+1):
         for c in subdivided[l-1]:
            for child in children(c, l-1):
               if child in fileindex_for_cellid.keys():
                  descr.write(".")
                  idxToFileIndex[idx] = fileindex_for_cellid[child]
               else:
                  descr.write("R")
                  subdivided[l].append(child)
               idx += 1
         if l < max_ref_level:
            descr.write("|")

      return descr.getvalue(), idxToFileIndex


   def getDescriptor(self, reinit=False):

      try:
         if reinit:
            raise Exception("reinit = True")
         with open(self.__metafile,'rb') as mfile:
            data = pickle.load(mfile)
            # print(data)
            self.__idxToFileIndex = data['idxToFileIndexMap']
            self.__descriptor = data['descr']
      except Exception as e:
         print("Re-initializing HTG, no metadata accessible because of ", e)
         self.__descriptor, self.__idxToFileIndex = self.buildDescriptor()
         if not os.path.isdir(os.path.dirname(self.__metafile)):
            os.mkdir(os.path.dirname(self.__metafile))
         with open(self.__metafile,'wb') as mfile:
            pickle.dump({"descr":self.__descriptor, "idxToFileIndexMap":self.__idxToFileIndex}, mfile)

      return self.__descriptor, self.__idxToFileIndex

   def getHTG(self):
      f = self.__reader

      descr, idxToFileIndex = self.getDescriptor()
      basegridsize = f.get_spatial_mesh_size()

      src = vtk.vtkHyperTreeGridSource()
      src.SetMaxDepth(f.get_max_refinement_level()+1)
      src.SetDimensions(int(basegridsize[0]+1),int(basegridsize[1]+1),int(basegridsize[2]+1))
      src.SetGridScale(f._VlsvReader__dx,f._VlsvReader__dy,f._VlsvReader__dz)
      src.SetBranchFactor(2)
      
      #descr = src.ConvertDescriptorStringToBitArray(descr)
      #src.SetDescriptorBits(descr)
      src.SetDescriptor(descr)

      src.Update()
      htg = src.GetHyperTreeGridOutput()

      xValues = vtk.vtkDoubleArray()
      xValues.SetNumberOfValues(int(basegridsize[0]+1))
      nodeCoordinatesX = f.read(tag="MESH_NODE_CRDS_X", mesh='SpatialGrid')
      nodeCoordinatesY = f.read(tag="MESH_NODE_CRDS_Y", mesh='SpatialGrid')
      nodeCoordinatesZ = f.read(tag="MESH_NODE_CRDS_Z", mesh='SpatialGrid')


      for i,x in enumerate(nodeCoordinatesX):
         xValues.SetValue(i, x)
      htg.SetXCoordinates(xValues)

      yValues = vtk.vtkDoubleArray()
      yValues.SetNumberOfValues(int(basegridsize[1]+1))
      for i,x in enumerate(nodeCoordinatesY):
         yValues.SetValue(i, x)
      htg.SetYCoordinates(yValues)

      zValues = vtk.vtkDoubleArray()
      zValues.SetNumberOfValues(int(basegridsize[2]+1))
      for i,x in enumerate(nodeCoordinatesZ):
         zValues.SetValue(i, x)
      htg.SetZCoordinates(zValues)

      htg.fileIndexArray = vtk.vtkDoubleArray()
      htg.fileIndexArray.SetName('fileIndex')

      for idx,fileIndex in self.__idxToFileIndex.items():
         htg.fileIndexArray.InsertTuple1(idx, fileIndex)
      htg.GetCellData().AddArray(htg.fileIndexArray)


      return htg

   '''
   def FillOutputPortInformation(self, port, info):
      pass

   def RequestDataObject(self, request, inInfo, outInfo):
      #This is where you can create output data objects if the output DATA_TYPE_NAME() is not a concrete type.
      pass
   '''
   def findVariablesFromVlsv(self, getReducers = False):

      vars = self.__reader.get_variables()
      if getReducers:
         reducers = self.__reader.get_reducers()
         vars_set = set(vars)
         vars_set.update(reducers)
         vars = list(vars_set)

      return vars

   def RequestInformation(self, request, inInfo, outInfo):

      info = outInfo.GetInformationObject(0)
      # print("VlsvVtkReader RequestInformation:")
      # print(info)
      if self.__htg is None:
         print("Creating htg via RequestInformation")
         self.__htg = self.getHTG()
         vars = self.findVariablesFromVlsv(getReducers=True)
         for name in vars:
            if "vg" in name:
               # print(name)
               self._arrayselection.AddArray(name)
               self._arrayselection.DisableArray(name)
               self.__cellarrays.append(name)
         # self.__htg.addArrayFromVlsv("vg_beta_star")
         # self.__htg.GetCellData().SetActiveScalars("vg_beta_star")

      dims = self.__htg.GetExtent()
      info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(), *dims)

      return 1
   
   '''
   def RequestUpdateExtent(self, request, inInfo, outInfo):
      pass
   '''

   '''This function adds one SpatialGrid variable from the reader object and maps
   that to the hypertreegrid object. Variable vector sizes of 1,2,3,4,9 supported.
   '''
   def addArrayFromVlsv(self, varname):

      htg = self.__htg
      # Do not re-add an already existing array
      if htg.GetCellData().HasArray(varname):
         print("skipped existing array")
         return True

      array = vtk.vtkDoubleArray()
      # cidArray2.DeepCopy(fileIndexArray)
      data = self.__reader.read_variable(varname)
      
      test_var = data[0]
      if test_var is None:
         logging.warning("varname " + varname + "returned None; skipping")
         return False
      if np.isscalar(test_var):
         varlen = 1
      elif test_var.ndim == 1:
         varlen = len(test_var)         
      elif test_var.ndim > 1:
         varlen = np.prod(test_var.shape)
      else:
         logging.warning("Weird output from " + varname)
         return False

      array.SetNumberOfComponents(varlen)
      array.SetNumberOfTuples(htg.fileIndexArray.GetNumberOfTuples())
      array.SetName(varname)
      if varlen == 1:
         for idx,fileIndex in self.__idxToFileIndex.items():
            array.SetTuple1(idx, data[fileIndex])
      elif varlen == 2:
         for idx,fileIndex in self.__idxToFileIndex.items():
            array.SetTuple2(idx, *data[fileIndex])
      elif varlen == 3:
         for idx,fileIndex in self.__idxToFileIndex.items():
            array.SetTuple3(idx, *data[fileIndex])
      elif varlen == 4:
         for idx,fileIndex in self.__idxToFileIndex.items():
            array.SetTuple4(idx, *data[fileIndex])
      elif varlen == 9:
         for idx,fileIndex in self.__idxToFileIndex.items():
            array.SetTuple9(idx, *np.reshape(data[fileIndex],(9)))
      else:
         raise RuntimeError("No vtk SetTuple wrapper function for varlen = " + str(varlen))
      
      
      htg.GetCellData().AddArray(array)
      return True

   def RequestData(self, request, inInfo, outInfo):

      # print("VlsvVtkReader RequestData:")
      # print(outInfo)

      if self.__htg is None:
         print("Creating htg via RequestData")
         self.__htg = self.getHTG()
      
      for name in self.__cellarrays:
         if self._arrayselection.ArrayIsEnabled(name):
            success = self.addArrayFromVlsv(name)
            if not success:
               self._arrayselection.RemoveArrayByName(name)


      output = vtk.vtkHyperTreeGrid.GetData(outInfo)



      output.ShallowCopy(self.__htg)

      return 1
   
   def ReadAllScalarsOn(self):
      self._arrayselection.EnableAllArrays()
      # for name in self.__cellarrays:
         
      #    self.__htg.addArrayFromVlsv(name)


   def GetDataArraySelection(self):
      return self._arrayselection


def __main__():
   import pytools as pt
   # This initializes a hypertreegrid from the given reader.
   reader = VlsvVtkReader()
   reader.SetFileName("/home/mjalho/Downloads/bulk.0002189.vlsv")
   # reader.SetFileName("/home/mjalho/bulk1.0001418.vlsv")
   # reader.ReadAllScalarsOn()
   reader.Update()
   cf = vtk.vtkHyperTreeGridContour()
   cf.SetInputConnection(reader.GetOutputPort())
   cf.SetValue(0, 1)

   m = vtk.vtkPolyDataMapper()
   
   m.SetInputConnection(cf.GetOutputPort())

   a = vtk.vtkActor()
   a.SetMapper(m)

   ren = vtk.vtkRenderer()
   ren.AddActor(a)

   renWin = vtk.vtkRenderWindow()
   renWin.AddRenderer(ren)
   renWin.SetSize(600, 600)

   # renWin.Render()
   
   # time.sleep(10)
   
   htg = reader.GetOutputDataObject(0)
   # print(htg)
   # htg.__vtkVlsvHyperTreeGrid_VlsvAddArray('vg_b_vol')
   # htg = reader.GetOutputData()
   # print(htg)
   # These functions grab one SpatialGrid variable and map that to 
   # the hypertreegrid. Variable vector sizes of 1,2,3,4,9 supported.
   # htg.addArrayFromVlsv("vg_b_vol")
   # htg.addArrayFromVlsv("vg_beta_star")

   writer = vtk.vtkXMLHyperTreeGridWriter()
   writer.SetFileName("output_FID.htg")
   
   writer.SetInputData(htg)
   writer.Write()


if __name__ == "__main__":
   __main__()

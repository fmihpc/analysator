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
import os
import pickle
from operator import itemgetter
import numbers

import analysator as pt
try:
   import vtk
   import vtk.numpy_interface
   import vtk.numpy_interface.dataset_adapter
   import vtk.vtkCommonColor
   import vtkmodules.vtkCommonColor
   from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
except Exception as e:
   logging.warning("VTK import (VTK >= 9.2.0 required) did not succeed due to "+str(e))
   raise(e)
import time


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
      self.__cellIDtoIdx = {}

   def SetFileName(self, filename):
      if filename != self.__FileName:
         self.Modified()
         self.__FileName = filename
         if self.__FileName is not None:
            self.__reader = pt.vlsvfile.VlsvReader(self.__FileName)
            fn = os.path.basename(self.__FileName)
            self.__metafile = os.path.join(self.__reader.get_cache_folder(),"vlsvvtkcache.pkl")

   def GetFileName(self):
      return self.__FileName

   def buildDescriptor(self):
      f = self.__reader
      fileindex_for_cellid = f.get_cellid_locations()
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

      def children(cid, level):
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

      if hasattr(self, "__descriptor") and hasattr(self, "__idxToFileIndex"):
         return self.__descriptor, self.__idxToFileIndex
      try:
         if reinit:
            raise Exception("reinit = True")
         with open(self.__metafile,'rb') as mfile:
            data = pickle.load(mfile)
            self.__idxToFileIndex = data['idxToFileIndexMap']
            self.__descriptor = data['descr']
      except Exception as e:
         print("Re-initializing HTG, no metadata accessible because of ", e)
         self.__descriptor, self.__idxToFileIndex = self.buildDescriptor()
         if not os.path.isdir(os.path.dirname(self.__metafile)):
            os.makedirs(os.path.dirname(self.__metafile), exist_ok=True)
         with open(self.__metafile,'wb') as mfile:
            pickle.dump({"descr":self.__descriptor, "idxToFileIndexMap":self.__idxToFileIndex}, mfile)

      return self.__descriptor, self.__idxToFileIndex
   
   def getCellIDtoIdxMap(self):
      if len(self.__cellIDtoIdx) == 0:
         FileIndexToID = {v:k for k,v in self.getDescriptor()[1].items()}
         for c,fi in self.__reader.get_cellid_locations().items():
            self.__cellIDtoIdx[c] = FileIndexToID[fi]

      return self.__cellIDtoIdx

   def getHTG(self):
      if self.__htg is not None:
         return self.__htg
      
      f = self.__reader

      descr, idxToFileIndex = self.getDescriptor()
      basegridsize = f.get_spatial_mesh_size()

      src = vtk.vtkHyperTreeGridSource()
      src.SetMaxDepth(f.get_max_refinement_level()+1)
      src.SetDimensions(int(basegridsize[0]+1),int(basegridsize[1]+1),int(basegridsize[2]+1))
      src.SetGridScale(f._VlsvReader__dx,f._VlsvReader__dy,f._VlsvReader__dz)
      src.SetBranchFactor(2)
      

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

      if self.__htg is None:
         print("Creating htg via RequestInformation")
         self.__htg = self.getHTG()
         vars = self.findVariablesFromVlsv(getReducers=True)
         for name in vars:
            if ("vg" in name.lower()) or (name.lower() == "cellid"):
               # print(name)
               self._arrayselection.AddArray(name)
               self._arrayselection.DisableArray(name)
               self.__cellarrays.append(name)
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

      htg = self.getHTG()
      # Do not re-add an already existing array
      if htg.GetCellData().HasArray(varname):
         print("skipped existing array")
         return True
      
      print("varname ", varname)

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

      print("VlsvVtkReader RequestData:")
      print(outInfo)

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

# TODO
   # def GetDataArraySelection(self):
   #    return self._arrayselection
   
   # def GetIdx(self, cellid):
   #    cidToIdx = self.getCellIDtoIdxMap()
   #    if isinstance(cellid,numbers.Number):
   #       return cidToIdx[cellid]
   #    else:
   #       return itemgetter(*cellid)(cidToIdx)


def __main__():
   import analysator as pt
   # This initializes a hypertreegrid from the given reader.
   reader = VlsvVtkReader()
   reader.SetFileName("/home/mjalho/Downloads/bulk.0002189.vlsv")

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


   htg = reader.GetOutputDataObject(0)


   writer = vtk.vtkXMLHyperTreeGridWriter()
   writer.SetFileName("output_FID.htg")
   
   writer.SetInputData(htg)
   writer.Write()


if __name__ == "__main__":
   __main__()

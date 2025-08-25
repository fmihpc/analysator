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

''' Utilities for caching VLSV data and metadata.
'''

import logging
import os
import sys
import warnings
import numbers
import numpy as np
from operator import itemgetter
import re
import rtree 

class VariableCache:
    ''' Class for handling in-memory variable/reducer caching.
    '''
    def __init__(self, reader):
        self.__varcache = {} # {(varname, operator):data}
        self.__reader = reader

    def keys(self):
        return self.__varcache.keys()
    
    def __getitem__(self, key):
       return self.__varcache[key]
       
    def __setitem__(self, key, value):
       self.__varcache[key] = value
       

    def read_variable_from_cache(self, name, cellids, operator):
      ''' Read variable from cache instead of the vlsv file.
         :param name: Name of the variable
         :param cellids: a value of -1 reads all data
         :param operator: Datareduction operator. "pass" does no operation on data
         :returns: numpy array with the data, same format as read_variable

         .. seealso:: :func:`read_variable`
      '''

      var_data = self.__varcache[(name,operator)]
      if var_data.ndim == 2:
         value_len = var_data.shape[1]
      else:
         value_len = 1
         
      if isinstance(cellids, numbers.Number):
         if cellids == -1:
            return var_data
         else:
            return var_data[self.__reader.get_cellid_locations()[cellids]]
      else:
         if(len(cellids) > 0):
            indices = np.array(itemgetter(*cellids)(self.__reader.get_cellid_locations()),dtype=np.int64)
         else:
            indices = np.array([],dtype=np.int64)
         if value_len == 1:
            return var_data[indices]
         else:
            return var_data[indices,:]

class FileCache:
   ''' Top-level class for caching to file.
   '''

   def __init__(self, reader) -> None:
      self.__reader = reader

      self.__rtree_index_files = []
      self.__rtree_index = None
      self.__rtree_idxfile = os.path.join(self.get_cache_folder(),"rtree.idx")
      self.__rtree_datfile = os.path.join(self.get_cache_folder(),"rtree.dat")
      self.__rtree_properties = rtree.index.Property()
      self.__rtree_properties.dimension = 3
      self.__rtree_properties.overwrite=True


   def get_cache_folder(self):
      fn = self.__reader.file_name

      head,tail = os.path.split(fn)
      path = head
      numslist = re.findall(r'\d+(?=\.vlsv)', tail)

      if(len(numslist) == 0):
         path = os.path.join(path,"vlsvcache",tail[:-5])
      else:
         nums = numslist[-1]
         head, tail = tail.split(nums)

         leading_zero = True
         path = os.path.join(path,"vlsvcache",head[:-1])
         for i,n in enumerate(nums):
            if n == '0' and leading_zero: continue
            if leading_zero:
               fmt = "{:07d}"
            else:
                fmt = "{:0"+str(7-i)+"d}"
            path = os.path.join(path, fmt.format(int(n)*10**(len(nums)-i-1)))
            leading_zero = False

      return path

   def clear_cache_folder(self):
      path = self.get_cache_folder()
      import shutil
      shutil.rmtree(path)

   def set_cellid_spatial_index(self, force = False):
      if not os.path.exists(self.get_cache_folder()):
         os.makedirs(self.get_cache_folder())

      if(force or (not os.path.isfile(self.__rtree_idxfile) or not os.path.isfile(self.__rtree_datfile))):
         
         
         
         bboxes = self.__reader.get_mesh_domain_extents("SpatialGrid")
         bboxes = bboxes.reshape((-1,6), order='C')
         print(bboxes.shape)

         self.__rtree_index = rtree.index.Index(self.__rtree_idxfile[:-4],properties=rtree_properties, interleaved=False)
         for rank, bbox in enumerate(bboxes):
            print(rank, bbox)
            self.__rtree_index.insert(rank, bbox)

         print("index set")
      else:
         print("index exists")

   def get_cellid_spatial_index(self, force = False):
      if self.__rtree_index == None:
         if(force or (not os.path.isfile(self.__rtree_idxfile) or not os.path.isfile(self.__rtree_datfile))):
            self.set_cellid_spatial_index(force)
         else:
            self.__rtree_index = rtree.index.Index(self.__rtree_idxfile[:-4], properties=self.__rtree_properties, interleaved=False)

      return self.__rtree_index


class MetadataFileCache(FileCache):
   ''' File caching class for storing "lightweight" metadata.
   '''
   pass
   # superclass constructor called instead if no __init__ here
   # def __init__(self, reader) -> None:
   #    super(MetadataFileCache, self).__init__(reader)



class VariableFileCache(FileCache):
   ''' File caching class for storing intermediate data, such as 
   gradient terms that are more expensive to compute, over whole grids with
   more HDD footprint.
   '''
   pass
   # superclass constructor called instead if no __init__ here
   # def __init__(self, reader) -> None:
   #    super(MetadataFileCache, self).__init__(reader)

class PicklableFile(object):
   ''' Picklable file pointer object.
   '''
   def __init__(self, fileobj):
      self.fileobj = fileobj

   def __getattr__(self, key):
      return getattr(self.fileobj, key)

   def __getstate__(self):
      ret = self.__dict__.copy()
      ret['_file_name'] = self.fileobj.name
      ret['_file_mode'] = self.fileobj.mode
      if self.fileobj.closed:
         ret['_file_pos'] = 0
      else:
         ret['_file_pos'] = self.fileobj.tell()
      del ret['fileobj']
      return ret

   def __setstate__(self, dict):
      self.fileobj = open(dict['_file_name'], dict['_file_mode'])
      self.fileobj.seek(dict['_file_pos'])
      del dict['_file_name']
      del dict['_file_mode']
      del dict['_file_pos']
      self.__dict__.update(dict)

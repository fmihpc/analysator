# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
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

import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
import os
import numbers
import vlsvvariables
from reduction import datareducers,multipopdatareducers,data_operators,v5reducers,multipopv5reducers
from collections import Iterable,OrderedDict
from vlsvwriter import VlsvWriter
from variable import get_data

class VlsvReader(object):
   ''' Class for reading VLSV files
   ''' 


   ''' Meshinfo is an information container for multiple meshes.
       Implemented as an empty class.
   '''
   class MeshInfo:
      pass

   file_name=""
   def __init__(self, file_name):
      ''' Initializes the vlsv file (opens the file, reads the file footer and reads in some parameters)

          :param file_name:     Name of the vlsv file
      '''
      # Make sure the path is set in file name: 
      file_name = os.path.abspath(file_name)

      self.file_name = file_name
      self.__fptr = open(self.file_name,"rb")
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}
      self.__max_spatial_amr_level = -1

      self.use_dict_for_blocks = False
      self.__fileindex_for_cellid_blocks={} # [0] is index, [1] is blockcount
      self.__cells_with_blocks = {} # per-pop
      self.__blocks_per_cell = {} # per-pop
      self.__blocks_per_cell_offsets = {} # per-pop
      self.__order_for_cellid_blocks = {} # per-pop
      
      self.__read_xml_footer()
      # Check if the file is using new or old vlsv format
      # Read parameters (Note: Reading the spatial cell locations and
      # storing them will anyway take the most time and memory):

      meshName="SpatialGrid"
      bbox = self.read(tag="MESH_BBOX", mesh=meshName)
      if bbox is None:
         #read in older vlsv files where the mesh is defined with parameters
         self.__xcells = (int)(self.read_parameter("xcells_ini"))
         self.__ycells = (int)(self.read_parameter("ycells_ini"))
         self.__zcells = (int)(self.read_parameter("zcells_ini"))
         self.__xblock_size = 1
         self.__yblock_size = 1
         self.__zblock_size = 1
         self.__xmin = self.read_parameter("xmin")
         self.__ymin = self.read_parameter("ymin")
         self.__zmin = self.read_parameter("zmin")
         self.__xmax = self.read_parameter("xmax")
         self.__ymax = self.read_parameter("ymax")
         self.__zmax = self.read_parameter("zmax")
      else:
         #new style vlsv file with 
         nodeCoordinatesX = self.read(tag="MESH_NODE_CRDS_X", mesh=meshName)   
         nodeCoordinatesY = self.read(tag="MESH_NODE_CRDS_Y", mesh=meshName)   
         nodeCoordinatesZ = self.read(tag="MESH_NODE_CRDS_Z", mesh=meshName)   
         self.__xcells = bbox[0]
         self.__ycells = bbox[1]
         self.__zcells = bbox[2]
         self.__xblock_size = bbox[3]
         self.__yblock_size = bbox[4]
         self.__zblock_size = bbox[5]
         self.__xmin = nodeCoordinatesX[0]
         self.__ymin = nodeCoordinatesY[0]
         self.__zmin = nodeCoordinatesZ[0]
         self.__xmax = nodeCoordinatesX[-1]
         self.__ymax = nodeCoordinatesY[-1]
         self.__zmax = nodeCoordinatesZ[-1]

      self.__dx = (self.__xmax - self.__xmin) / (float)(self.__xcells)
      self.__dy = (self.__ymax - self.__ymin) / (float)(self.__ycells)
      self.__dz = (self.__zmax - self.__zmin) / (float)(self.__zcells)

      self.__meshes = {}

      # Iterate through the XML tree, find all populations
      # (identified by their BLOCKIDS tag)
      self.active_populations=[]
      for child in self.__xml_root:
          if child.tag == "BLOCKIDS":
              if "name" in child.attrib:
                  popname = child.attrib["name"] 
              else:
                  popname = "avgs"

              # Create a new (empty) MeshInfo-object for this population
              pop = self.MeshInfo()
              
              # Update list of active populations
              if not popname in self.active_populations: self.active_populations.append(popname)

              bbox = self.read(tag="MESH_BBOX", mesh=popname)
              if bbox is None:
                 if self.read_parameter("vxblocks_ini") is not None:
                    #read in older vlsv files where the mesh is defined with
                    #parameters (only one possible)
                    pop.__vxblocks = (int)(self.read_parameter("vxblocks_ini"))
                    pop.__vyblocks = (int)(self.read_parameter("vyblocks_ini"))
                    pop.__vzblocks = (int)(self.read_parameter("vzblocks_ini"))
                    pop.__vxblock_size = 4
                    pop.__vyblock_size = 4
                    pop.__vzblock_size = 4
                    pop.__vxmin = self.read_parameter("vxmin")
                    pop.__vymin = self.read_parameter("vymin")
                    pop.__vzmin = self.read_parameter("vzmin")
                    pop.__vxmax = self.read_parameter("vxmax")
                    pop.__vymax = self.read_parameter("vymax")
                    pop.__vzmax = self.read_parameter("vzmax")
                    # Velocity cell lengths
                    pop.__dvx = ((pop.__vxmax - pop.__vxmin) / (float)(pop.__vxblocks)) / (float)(pop.__vxblock_size)
                    pop.__dvy = ((pop.__vymax - pop.__vymin) / (float)(pop.__vyblocks)) / (float)(pop.__vyblock_size)
                    pop.__dvz = ((pop.__vzmax - pop.__vzmin) / (float)(pop.__vzblocks)) / (float)(pop.__vzblock_size)

                 else:
                    #no velocity space in this file, e.g., file n ot written by Vlasiator 
                    pop.__vxblocks = 0
                    pop.__vyblocks = 0
                    pop.__vzblocks = 0
                    pop.__vxblock_size = 4
                    pop.__vyblock_size = 4
                    pop.__vzblock_size = 4
                    pop.__vxmin = 0
                    pop.__vymin = 0
                    pop.__vzmin = 0
                    pop.__vxmax = 0
                    pop.__vymax = 0
                    pop.__vzmax = 0
                    # Velocity cell lengths
                    pop.__dvx = 1
                    pop.__dvy = 1
                    pop.__dvz = 1

              else:
                 #new style vlsv file with bounding box
                 nodeCoordinatesX = self.read(tag="MESH_NODE_CRDS_X", mesh=popname)   
                 nodeCoordinatesY = self.read(tag="MESH_NODE_CRDS_Y", mesh=popname)   
                 nodeCoordinatesZ = self.read(tag="MESH_NODE_CRDS_Z", mesh=popname)   
                 pop.__vxblocks = bbox[0]
                 pop.__vyblocks = bbox[1]
                 pop.__vzblocks = bbox[2]
                 pop.__vxblock_size = bbox[3]
                 pop.__vyblock_size = bbox[4]
                 pop.__vzblock_size = bbox[5]
                 pop.__vxmin = nodeCoordinatesX[0]
                 pop.__vymin = nodeCoordinatesY[0]
                 pop.__vzmin = nodeCoordinatesZ[0]
                 pop.__vxmax = nodeCoordinatesX[-1]
                 pop.__vymax = nodeCoordinatesY[-1]
                 pop.__vzmax = nodeCoordinatesZ[-1]
                 # Velocity cell lengths
                 pop.__dvx = ((pop.__vxmax - pop.__vxmin) / (float)(pop.__vxblocks)) / (float)(pop.__vxblock_size)
                 pop.__dvy = ((pop.__vymax - pop.__vymin) / (float)(pop.__vyblocks)) / (float)(pop.__vyblock_size)
                 pop.__dvz = ((pop.__vzmax - pop.__vzmin) / (float)(pop.__vzblocks)) / (float)(pop.__vzblock_size)

              self.__meshes[popname]=pop
              if os.getenv('PTNONINTERACTIVE') == None:
                 print("Found population " + popname)
              
              # Precipitation energy bins
              i = 0
              energybins = []
              binexists = True
              while binexists:
                 binexists = self.check_parameter("{}_PrecipitationCentreEnergy{}".format(popname,i))
                 if binexists:
                    binvalue = self.read_parameter("{}_PrecipitationCentreEnergy{}".format(popname,i))
                    energybins.append(binvalue)
                 i = i + 1
              if i > 1:
                 pop.__precipitation_centre_energy = np.asarray(energybins)
                 vlsvvariables.speciesprecipitationenergybins[popname] = energybins

      self.__fptr.close()


   def __read_xml_footer(self):
      ''' Reads in the XML footer of the VLSV file and store all the content
      ''' 
      max_xml_size = 1000000
      #(endianness,) = struct.unpack("c", fptr.read(1))
      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr
      # Eight first bytes indicate whether the system is big_endianness or something else
      endianness_offset = 8
      fptr.seek(endianness_offset)
      # Read 8 bytes as unsigned long long (uint64_t in this case) after endianness, this tells the offset of the XML file.
      uint64_byte_amount = 8
      (offset,) = struct.unpack("Q", fptr.read(uint64_byte_amount))
      # Move to the xml offset
      fptr.seek(offset)
      # Read the xml data
      xml_data = fptr.read(max_xml_size)
      # Read the xml as string
      (xml_string,) = struct.unpack("%ds" % len(xml_data), xml_data)
      # Input the xml data into xml_root
      self.__xml_root = ET.fromstring(xml_string)
      if self.__fptr.closed:
         fptr.close()

   def __read_fileindex_for_cellid(self):
      """ Read in the cell ids and create an internal dictionary to give the index of an arbitrary cellID
      """
      cellids=self.read(mesh="SpatialGrid",name="CellID", tag="VARIABLE")

      #Check if it is not iterable. If it is a scale then make it a list
      if(not isinstance(cellids, Iterable)):
         cellids=[ cellids ]
      for index,cellid in enumerate(cellids):
         self.__fileindex_for_cellid[cellid]=index
         

   def __read_blocks(self, cellid, pop="proton"):
      ''' Read raw velocity block data from the open file.
      
      :param cellid: Cell ID of the cell whose velocity blocks are read
      :returns: A numpy array with block ids and their data
      '''
      if self.use_dict_for_blocks: #deprecated version
         if not pop in self.__fileindex_for_cellid_blocks:
            self.__set_cell_offset_and_blocks(pop)

         if( (cellid in self.__fileindex_for_cellid_blocks[pop]) == False ):
            # Cell id has no blocks
            return []
         offset = self.__fileindex_for_cellid_blocks[pop][cellid][0]
         num_of_blocks = self.__fileindex_for_cellid_blocks[pop][cellid][1]
      else:
         # Uses arrays (much faster to initialize)
         if not pop in self.__cells_with_blocks:
            self.__set_cell_offset_and_blocks_nodict(pop) 
         # Check that cells has vspace
         try:
            cells_with_blocks_index = self.__order_for_cellid_blocks[pop][cellid]
         except:
            print("Cell does not have velocity distribution")
            return []
         offset = self.__blocks_per_cell_offsets[pop][cells_with_blocks_index]
         num_of_blocks = self.__blocks_per_cell[pop][cells_with_blocks_index]


      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in block values
         if ("name" in child.attrib) and (child.attrib["name"] == pop) and (child.tag == "BLOCKVARIABLE"):
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = int(offset * vector_size * element_size + ast.literal_eval(child.text))
#            for i in range(0, cells_with_blocks_index[0]):
#               offset_avgs += blocks_per_cell[i]*vector_size*element_size

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)

         # Read in block coordinates:
         # (note the special treatment in case the population is named 'avgs'
         if (pop == 'avgs' or ("name" in child.attrib) and (child.attrib["name"] == pop)) and (child.tag == "BLOCKIDS"):
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_ids = int(offset * vector_size * element_size + ast.literal_eval(child.text))

            fptr.seek(offset_block_ids)
            if datatype == "uint" and element_size == 4:
               data_block_ids = np.fromfile(fptr, dtype = np.uint32, count = vector_size*num_of_blocks)
            elif datatype == "uint" and element_size == 8:
               data_block_ids = np.fromfile(fptr, dtype = np.uint64, count = vector_size*num_of_blocks)
            else:
               print("Error! Bad block id data!")
               print("Data type: " + datatype + ", element size: " + str(element_size))
               return

            data_block_ids = np.reshape(data_block_ids, (len(data_block_ids),) )

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      print("data_avgs = " + str(data_avgs) + ", data_block_ids = " + str(data_block_ids))
      if len(data_avgs) != len(data_block_ids):
         print("BAD DATA SIZES")

      return [data_block_ids, data_avgs]

   def __set_cell_offset_and_blocks(self, pop="proton"):
      ''' Read blocks per cell and the offset in the velocity space arrays for
          every cell with blocks into a private dictionary.
          Deprecated in favor of below version.
      '''
      print("Getting offsets for population " + pop)
      if pop in self.__fileindex_for_cellid_blocks:
         # There's stuff already saved into the dictionary, don't save it again
         return
      #these two arrays are in the same order: 
      #list of cells for which dist function is saved
      cells_with_blocks = self.read(mesh="SpatialGrid",tag="CELLSWITHBLOCKS", name=pop)
      #number of blocks in each cell for which data is stored
      blocks_per_cell = self.read(mesh="SpatialGrid",tag="BLOCKSPERCELL", name=pop)

      # Navigate to the correct position:
      from copy import copy
      offset = 0
      #self.__fileindex_for_cellid_blocks[pop] = {}
      self.__fileindex_for_cellid_blocks[pop] = dict.fromkeys(cells_with_blocks) # should be faster but negligible difference
      for i in range(0, len(cells_with_blocks)):
         self.__fileindex_for_cellid_blocks[pop][cells_with_blocks[i]] = [copy(offset), copy(blocks_per_cell[i])]
         offset += blocks_per_cell[i]

   def __set_cell_offset_and_blocks_nodict(self, pop="proton"):
      ''' Read blocks per cell and the offset in the velocity space arrays for every cell with blocks.
          Stores them in arrays. Creates a private dictionary with addressing to the array.
          This method should be faster than the above function.
      '''
      if pop in self.__cells_with_blocks:
         # There's stuff already saved into the dictionary, don't save it again
         return

      print("Getting offsets for population " + pop)

      self.__cells_with_blocks[pop] = np.atleast_1d(self.read(mesh="SpatialGrid",tag="CELLSWITHBLOCKS", name=pop))
      self.__blocks_per_cell[pop] = np.atleast_1d(self.read(mesh="SpatialGrid",tag="BLOCKSPERCELL", name=pop))

      self.__blocks_per_cell_offsets[pop] = np.empty(len(self.__cells_with_blocks[pop]))
      self.__blocks_per_cell_offsets[pop][0] = 0.0
      self.__blocks_per_cell_offsets[pop][1:] = np.cumsum(self.__blocks_per_cell[pop][:-1])
      self.__order_for_cellid_blocks[pop] = {}
      for index,cellid in enumerate(self.__cells_with_blocks[pop]):
         self.__order_for_cellid_blocks[pop][cellid]=index

   def list(self, parameter=True, variable=True, mesh=False, datareducer=False, operator=False, other=False):
      ''' Print out a description of the content of the file. Useful
         for interactive usage. Default is to list parameters and variables, query selection can be adjusted with keywords:

         Default and supported keywords:

         parameter=True 
         variable=True
         mesh=False
         datareducer=False 
         operator=False 
         other=False
      '''
      if parameter:
         print("tag = PARAMETER")
         for child in self.__xml_root:
            if child.tag == "PARAMETER" and "name" in child.attrib:
               print("   ", child.attrib["name"])
      if variable:
         print("tag = VARIABLE")
         for child in self.__xml_root:
            if child.tag == "VARIABLE" and "name" in child.attrib:
               print("   ", child.attrib["name"])
      if mesh:
         print("tag = MESH")
         for child in self.__xml_root:
            if child.tag == "MESH" and "name" in child.attrib:
               print("   ", child.attrib["name"])
      if datareducer:
         print("Datareducers:")
         for name in datareducers:
            print("   ",name, " based on ", datareducers[name].variables)
      if operator:
         print("Data operators:")
         for name in data_operators:
            if type(name) is str:
               if not name.isdigit():
                  print("   ",name)
      if other:
         print("Other:")
         for child in self.__xml_root:
            if child.tag != "PARAMETER" and child.tag != "VARIABLE" and child.tag != "MESH":
               print("    tag = ", child.tag, " mesh = ", child.attrib["mesh"])

   def check_parameter( self, name ):
      ''' Checks if a given parameter is in the vlsv reader

          :param name:             Name of the parameter
          :returns:                True if the parameter is in the vlsv file, false if not

          .. note:: This should be used for checking if a parameter exists, e.g. for different Vlasiator versions and time output

          .. code-block:: python

             # Example usage:
             vlsvReader = pt.vlsvfile.VlsvReader("test.vlsv")
             if vlsvReader.check_parameter( "time" ):
                time = vlsvReader.read_parameter("time")
             elif vlsvReader.check_parameter( "t" ):
                time = vlsvReader.read_parameter("t")
             else:
                time = None
      '''
      for child in self.__xml_root:
         if child.tag == "PARAMETER" and "name" in child.attrib:
            if child.attrib["name"].lower() == name.lower():
               return True
      return False

   def check_variable( self, name ):
      ''' Checks if a given variable is in the vlsv reader

          :param name:             Name of the variable
          :returns:                True if the variable is in the vlsv file, false if not

          .. note:: This should be used for checking if a variable exists in case a function behaves differently for ex. if B vector is in the vlsv and if not

          .. code-block:: python

             # Example usage:
             vlsvReader = pt.vlsvfile.VlsvReader("test.vlsv")
             if vlsvReader.check_variable( "B" ):
                # Variable can be plotted
                plot_B()
             else:
                # Variable not in the vlsv file
                plot_B_vol()
      '''
      for child in self.__xml_root:
         if child.tag == "VARIABLE" and "name" in child.attrib:
            if child.attrib["name"].lower() == name.lower():
               return True
      return False

   def check_population( self, popname ):
      ''' Checks if a given population is in the vlsv file

          :param name:             Name of the population
          :returns:                True if the population is in the vlsv file, false if not

          .. code-block:: python

             # Example usage:
             vlsvReader = pt.vlsvfile.VlsvReader("test.vlsv")
             if vlsvReader.check_population( "avgs" ):
                plot_population('avgs')
             else:
                if vlsvReader.check_population( "proton" ):
                   # File is newer with proton population
                   plot_population('proton')
      '''
      blockidsexist = False
      foundpop = False
      for child in self.__xml_root:
         if child.tag == "BLOCKIDS":
            if "name" in child.attrib:
               if popname.lower() == child.attrib["name"].lower():
                  foundpop = True
            else:
               blockidsexist = True
      if blockidsexist:
         for child in self.__xml_root:
            if child.tag == "BLOCKVARIABLE":
               if "name" in child.attrib:
                  if popname.lower() == child.attrib["name"].lower(): # avgs
                     foundpop = True
      return foundpop

   def get_all_variables( self ):
      ''' Returns all variables in the vlsv reader and the data reducer
          :returns:                List of variable is in the vlsv file
          .. code-block:: python
             # Example usage:
             vlsvReader = pt.vlsvfile.VlsvReader("test.vlsv")
             vars = vlsvReader.get_variables()
      '''
      varlist = [];
      for child in self.__xml_root:
         if child.tag == "VARIABLE" and "name" in child.attrib:
            varlist.append(child.attrib["name"])
      return varlist

   def get_cellid_locations(self):
      ''' Returns a dictionary with cell id as the key and the index of the cell id as the value. The index is used to locate the cell id's values in the arrays that this reader returns
      '''
      if len( self.__fileindex_for_cellid ) == 0:
         self.__read_fileindex_for_cellid()
      return self.__fileindex_for_cellid

   def read(self, name="", tag="", mesh="", operator="pass", cellids=-1):
      ''' Read data from the open vlsv file. 
      
      :param name: Name of the data array
      :param tag:  Tag of the data array.
      :param mesh: Mesh for the data array
      :param operator: Datareduction operator. "pass" does no operation on data.
      :param cellids:  If -1 then all data is read. If nonzero then only the vector
                       for the specified cell id or cellids is read
      :returns: numpy array with the data

      .. seealso:: :func:`read_variable` :func:`read_variable_info`
      '''

      if (len( self.__fileindex_for_cellid ) == 0):
         # Do we need to construct the cellid index?
         if isinstance(cellids, numbers.Number): # single or all cells
            if cellids >= 0: # single cell
               self.__read_fileindex_for_cellid()
         else: # list of cellids
            self.__read_fileindex_for_cellid()
               
      if tag == "" and name == "":
         print("Bad arguments at read")

      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      # Force lowercase name for internal checks
      name = name.lower()
         
      # Get population and variable names from data array name 
      if '/' in name:
         popname = name.split('/')[0]
         varname = name.split('/')[1]
      else:
         popname = 'pop'
         varname = name

      # Seek for requested data in VLSV file
      for child in self.__xml_root:
         if tag != "":
            if child.tag != tag:
               continue
         if name != "":
            if "name" in child.attrib and child.attrib["name"].lower() != name:
               continue
         if mesh != "":
            if "mesh" in child.attrib and child.attrib["mesh"] != mesh:
               continue
         if child.tag == tag:
            # Found the requested data entry in the file
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]
            variable_offset = ast.literal_eval(child.text)

            # Define efficient method to read data in
            try: # try-except to see how many cellids were given
               lencellids=len(cellids) 
               # Read multiple specified cells
               # If we're reading a large amount of single cells, it'll be faster to just read all
               # data from the file system and sort through it. For the CSC disk system, this
               # becomes more efficient for over ca. 5000 cellids.
               arraydata = []
               if lencellids>5000: 
                  result_size = len(cellids)
                  read_size = array_size
                  read_offsets = [0]
               else: # Read multiple cell ids one-by-one
                  result_size = len(cellids)
                  read_size = 1
                  read_offsets = [self.__fileindex_for_cellid[cid]*element_size*vector_size for cid in cellids]
            except: # single cell or all cells
               if cellids < 0: # -1, read all cells
                  result_size = array_size
                  read_size = array_size
                  read_offsets = [0]
               else: # single cell id
                  result_size = 1
                  read_size = 1
                  read_offsets = [self.__fileindex_for_cellid[cellids]*element_size*vector_size]
                  
            for r_offset in read_offsets:
               use_offset = int(variable_offset + r_offset)
               fptr.seek(use_offset)
               if datatype == "float" and element_size == 4:
                  data = np.fromfile(fptr, dtype = np.float32, count=vector_size*read_size)
               if datatype == "float" and element_size == 8:
                  data = np.fromfile(fptr, dtype=np.float64, count=vector_size*read_size)
               if datatype == "int" and element_size == 4:
                  data = np.fromfile(fptr, dtype=np.int32, count=vector_size*read_size)
               if datatype == "int" and element_size == 8:
                  data = np.fromfile(fptr, dtype=np.int64, count=vector_size*read_size)
               if datatype == "uint" and element_size == 4:
                  data = np.fromfile(fptr, dtype=np.uint32, count=vector_size*read_size)
               if datatype == "uint" and element_size == 8:
                  data = np.fromfile(fptr, dtype=np.uint64, count=vector_size*read_size)
               if len(read_offsets)!=1:
                  arraydata.append(data)
            
            if len(read_offsets)==1 and result_size<read_size:
               # Many single cell id's requested
               # Pick the elements corresponding to the requested cells
               for cid in cellids:
                  append_offset = self.__fileindex_for_cellid[cid]*vector_size
                  arraydata.append(data[append_offset:append_offset+vector_size])
               data = np.squeeze(np.array(arraydata))

            if len(read_offsets)!=1:
               # Not-so-many single cell id's requested
               data = np.squeeze(np.array(arraydata))

            if self.__fptr.closed:
               fptr.close()

            if vector_size > 1:
               data=data.reshape(result_size, vector_size)
            
            # If variable vector size is 1, and requested magnitude, change it to "absolute"
            if vector_size == 1 and operator=="magnitude":
               print("Data variable with vector size 1: Changed magnitude operation to absolute")
               operator="absolute"

            if result_size == 1:
               return data_operators[operator](data[0])
            else:
               return data_operators[operator](data)

      # Check which set of datareducers to use
      if varname[0:3]=="vg_":
         reducer_reg = v5reducers
         reducer_multipop = multipopv5reducers
      else:
         reducer_reg = datareducers
         reducer_multipop = multipopdatareducers
            
      # If this is a variable that can be summed over the populations (Ex. rho, PTensorDiagonal, ...)
      if self.check_variable(self.active_populations[0]+'/'+name): 
         tmp_vars = []
         for pname in self.active_populations:
            vlsvvariables.activepopulation = pname
            tmp_vars.append( self.read( pname+'/'+name, tag, mesh, "pass", cellids ) )
         return data_operators[operator](data_operators["sum"](tmp_vars))

      # Check if the name is in datareducers
      if name in reducer_reg:
         reducer = reducer_reg[name]
         # Read the necessary variables:

         # If variable vector size is 1, and requested magnitude, change it to "absolute"
         if reducer.vector_size == 1 and operator=="magnitude":
            print("Data reducer with vector size 1: Changed magnitude operation to absolute")
            operator="absolute"

         # Return the output of the datareducer
         if reducer.useVspace:
            actualcellids = self.read(mesh="SpatialGrid", name="CellID", tag="VARIABLE", operator=operator, cellids=cellids)
            output = np.zeros(len(actualcellids))
            index = 0
            for singlecellid in actualcellids:
               velocity_cell_data = self.read_velocity_cells(singlecellid)
               # Get cells:
               vcellids = list(velocity_cell_data.keys())
               # Get coordinates:
               velocity_coordinates = self.get_velocity_cell_coordinates(vcellids)
               tmp_vars = []
               for i in np.atleast_1d(reducer.variables):
                  tmp_vars.append( self.read( i, tag, mesh, "pass", singlecellid ) )
               output[index] = reducer.operation( tmp_vars , velocity_cell_data, velocity_coordinates )
               index+=1
               print(index,"/",len(actualcellids))
            return data_operators[operator](output)
         else:
            tmp_vars = []
            for i in np.atleast_1d(reducer.variables):
               tmp_vars.append( self.read( i, tag, mesh, "pass", cellids ) )
            return data_operators[operator](reducer.operation( tmp_vars ))

      # Check if the name is in multipop datareducers
      if 'pop/'+varname in reducer_multipop:

         reducer = reducer_multipop['pop/'+varname]
         # If variable vector size is 1, and requested magnitude, change it to "absolute"
         if reducer.vector_size == 1 and operator=="magnitude":
            print("Data reducer with vector size 1: Changed magnitude operation to absolute")
            operator="absolute"

         if reducer.useVspace:
            print("Error: useVspace flag is not implemented for multipop datareducers!") 
            return

         # sum over populations
         if popname=='pop':
            # Read the necessary variables:
            tmp_vars = []
            for pname in self.active_populations:
               vlsvvariables.activepopulation = pname
               tmp_vars.append( self.read( pname+'/'+varname, tag, mesh, "pass", cellids ) )
            return data_operators[operator](data_operators["sum"](tmp_vars))
         else:
            vlsvvariables.activepopulation = popname

         # Read the necessary variables:
         tmp_vars = []
         for i in np.atleast_1d(reducer.variables):
            if '/' not in i:
               tmp_vars.append( self.read( i, tag, mesh, "pass", cellids ) )
            else:
               tvar = i.split('/')[1]
               tmp_vars.append( self.read( popname+'/'+tvar, tag, mesh, "pass", cellids ) )
         return data_operators[operator](reducer.operation( tmp_vars ))

      if name!="":
         print("Error: variable "+name+"/"+tag+"/"+mesh+"/"+operator+" not found in .vlsv file or in data reducers!") 
      if self.__fptr.closed:
         fptr.close()


   def read_metadata(self, name="", tag="", mesh=""):
      ''' Read variable metadata from the open vlsv file. 
      
      :param name: Name of the data array
      :param tag:  Tag of the data array.
      :param mesh: Mesh for the data array
      :returns: four strings:
                the unit of the variable as a regular string
                the unit of the variable as a LaTeX-formatted string
                the description of the variable as a LaTeX-formatted string
                the conversion factor to SI units as a string                  
      '''

      if tag == "" and name == "":
         print("Bad arguments at read")

      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      # Force lowercase name for internal checks
      name = name.lower()
      
      # Seek for requested data in VLSV file
      for child in self.__xml_root:         
         if tag != "":
            if child.tag != tag:
               continue
         if name != "":
            if "name" in child.attrib and child.attrib["name"].lower() != name:
               continue
         if mesh != "":
            if "mesh" in child.attrib and child.attrib["mesh"] != mesh:
               continue
         if not "name" in child.attrib:
            continue
         # Found the requested data entry in the file
         unit = child.attrib["unit"]
         unitLaTeX = child.attrib["unitLaTeX"]
         variableLaTeX = child.attrib["variableLaTeX"]
         unitConversion = child.attrib["unitConversion"] 
         return unit, unitLaTeX, variableLaTeX, unitConversion
            
      if name!="":
         print("Error: variable "+name+"/"+tag+"/"+mesh+" not found in .vlsv file!" )
      if self.__fptr.closed:
         fptr.close()
      return -1
         

   def read_interpolated_fsgrid_variable(self, name, coordinates, operator="pass",periodic=["True", "True", "True"]):
      ''' Read a linearly interpolated FSgrid variable value from the open vlsv file.
      Arguments:
      :param name: Name of the (FSgrid) variable
      :param coords: Coordinates from which to read data 
      :param periodic: Periodicity of the system. Default is periodic in all dimension
      :param operator: Datareduction operator. "pass" does no operation on data
      :returns: numpy array with the data

      .. seealso:: :func:`read` :func:`read_variable_info`
      '''

      # At this stage, this function has not yet been implemented -- print a warning and exit
      print('Interpolation of FSgrid variables has not yet been implemented; exiting.')
      return -1


   def read_interpolated_variable(self, name, coordinates, operator="pass",periodic=["True", "True", "True"]):
      ''' Read a linearly interpolated variable value from the open vlsv file.
      Arguments:
      :param name: Name of the variable
      :param coords: Coordinates from which to read data 
      :param periodic: Periodicity of the system. Default is periodic in all dimension
      :param operator: Datareduction operator. "pass" does no operation on data
      :returns: numpy array with the data

      .. seealso:: :func:`read` :func:`read_variable_info`
      '''

      # First test whether the requested variable is on the FSgrid, and redirect to the dedicated function if needed
      if name[0:3] == 'fg_':
         return self.read_interpolated_fsgrid_variable(name, coordinates, operator, periodic)

      coordinates = get_data(coordinates)
      
      if len(np.shape(coordinates)) == 1:
         # Get closest id
         closest_cell_id=self.get_cellid(coordinates)
         if closest_cell_id == 0:
            return None
         closest_cell_coordinates=self.get_cell_coordinates(closest_cell_id)

         # Now identify the lower one of the 8 neighbor cells
         offset = [0 if coordinates[0] > closest_cell_coordinates[0] else -1,\
                   0 if coordinates[1] > closest_cell_coordinates[1] else -1,\
                   0 if coordinates[2] > closest_cell_coordinates[2] else -1]
         lower_cell_id = self.get_cell_neighbor(closest_cell_id, offset, periodic)
         lower_cell_coordinates=self.get_cell_coordinates(lower_cell_id)
         offset = [1,1,1]
         upper_cell_id = self.get_cell_neighbor(lower_cell_id, offset, periodic)
         upper_cell_coordinates=self.get_cell_coordinates(upper_cell_id)
         if (lower_cell_id<1 or upper_cell_id<1):
            print("Error: requested cell id for interpolation outside simulation domain")
            return -1
         # If the interpolation is done across two different refinement levels, issue a warning
         if self.get_amr_level(upper_cell_id) != self.get_amr_level(lower_cell_id):
            print("Warning: Interpolation across different AMR levels; this might lead to suboptimal results.")

         scaled_coordinates=np.zeros(3)
         for i in range(3):
            if lower_cell_coordinates[i] != upper_cell_coordinates[i]:
               scaled_coordinates[i]=(coordinates[i] - lower_cell_coordinates[i])/(upper_cell_coordinates[i] - lower_cell_coordinates[i])
            else:
               scaled_coordinates[i] = 0.0 # Special case for periodic systems with one cell in a dimension

         test_val=self.read_variable(name,lower_cell_id,operator)
         if isinstance(test_val, Iterable):
            try:
               value_length=len(test_val)
            except Exception as e:
               # Can't determine size, maybe some division by zero?
               value_length=1
         else:
            value_length=1
         
         # Now identify 8 cells, starting from the lower one
         ngbrvalues=np.zeros((2,2,2,value_length))
         for x in [0,1]:
            for y in [0,1]:
               for z  in [0,1]:
                  ngbrvalues[x,y,z,:] = self.read_variable(name, \
                                                           self.get_cell_neighbor(lower_cell_id, [x,y,z] , periodic), \
                                                           operator)

         c2d=np.zeros((2,2,value_length))
         for y in  [0,1]:
            for z in  [0,1]:
               c2d[y,z,:]=ngbrvalues[0,y,z,:]* (1- scaled_coordinates[0]) +  ngbrvalues[1,y,z,:]*scaled_coordinates[0]

         c1d=np.zeros((2,value_length))
         for z in [0,1]:
            c1d[z,:]=c2d[0,z,:]*(1 - scaled_coordinates[1]) + c2d[1,z,:] * scaled_coordinates[1]
            
         final_value=c1d[0,:] * (1 - scaled_coordinates[2]) + c1d[1,:] * scaled_coordinates[2]
         if len(final_value)==1:
            return final_value[0]
         else:
            return final_value

      else:
         # Multiple coordinates
         
         # Read all cells as we're anyway going to need this
         whole_cellids  = self.read_variable("CellID")
         whole_variable = self.read_variable(name,operator=operator)

         # Check one value for the length
         if isinstance(whole_variable[0], Iterable):
            value_length=len(whole_variable[0])
         else:
            value_length=1
         
         # Start iteration
         if value_length == 1:
            ret_array = np.zeros(len(coordinates))
         else:
            ret_array = np.zeros((len(coordinates), value_length))
         for i,cell_coords in enumerate(coordinates):
            closest_cell_id=self.get_cellid(cell_coords)
            if closest_cell_id == 0:
               if value_length==1:
                  ret_array[i]=None
               else:
                  ret_array[i,:] = None
               continue
            closest_cell_coordinates=self.get_cell_coordinates(closest_cell_id)
            
            # Now identify the lower one of the 8 neighbor cells
            offset = [0 if cell_coords[0] > closest_cell_coordinates[0] else -1,\
                      0 if cell_coords[1] > closest_cell_coordinates[1] else -1,\
                      0 if cell_coords[2] > closest_cell_coordinates[2] else -1]
            lower_cell_id = self.get_cell_neighbor(closest_cell_id, offset, periodic)
            lower_cell_coordinates=self.get_cell_coordinates(lower_cell_id)
            offset = [1,1,1]
            upper_cell_id = self.get_cell_neighbor(lower_cell_id, offset, periodic)
            upper_cell_coordinates=self.get_cell_coordinates(upper_cell_id)
           
            if self.get_amr_level(upper_cell_id) != self.get_amr_level(lower_cell_id):
               print("Warning: Interpolation across different AMR levels; this might lead to suboptimal results.")
 
            scaled_coordinates=np.zeros(3)
            for j in range(3):
               if lower_cell_coordinates[j] != upper_cell_coordinates[j]:
                  scaled_coordinates[j]=(cell_coords[j] - lower_cell_coordinates[j])/(upper_cell_coordinates[j] - lower_cell_coordinates[j])
               else:
                  scaled_coordinates[j] = 0.0 # Special case for periodic systems with one cell in a dimension
            
            ngbrvalues=np.zeros((2,2,2,value_length))
            for x in [0,1]:
               for y in [0,1]:
                  for z  in [0,1]:
                     cellid_neighbor = int(self.get_cell_neighbor(lower_cell_id, [x,y,z] , periodic))
                     ngbrvalues[x,y,z,:] = whole_variable[np.nonzero(whole_cellids==cellid_neighbor)[0][0]]
            
            c2d=np.zeros((2,2,value_length))
            for y in  [0,1]:
               for z in  [0,1]:
                  c2d[y,z,:]=ngbrvalues[0,y,z,:]* (1- scaled_coordinates[0]) +  ngbrvalues[1,y,z,:]*scaled_coordinates[0]
            
            c1d=np.zeros((2,value_length))
            for z in [0,1]:
               c1d[z,:]=c2d[0,z,:]*(1 - scaled_coordinates[1]) + c2d[1,z,:] * scaled_coordinates[1]
            
            final_value=c1d[0,:] * (1 - scaled_coordinates[2]) + c1d[1,:] * scaled_coordinates[2]
            if len(final_value)==1:
               ret_array[i] = final_value[0]
            else:
               ret_array[i, :] = final_value
         
         # Done.
         return ret_array

   def read_fsgrid_variable(self, name, operator="pass"):
       ''' Reads fsgrid variables from the open vlsv file.
       Arguments:
       :param name: Name of the variable
       :param operator: Datareduction operator. "pass" does no operation on data
       :returns: *ordered* numpy array with the data

       ... seealso:: :func:`read_variable`
       '''

       # Get fsgrid domain size (this can differ from vlasov grid size if refined)
       bbox = self.read(tag="MESH_BBOX", mesh="fsgrid")

       # Read the raw array data
       rawData = self.read(mesh='fsgrid', name=name, tag="VARIABLE", operator=operator)

       # Determine fsgrid domain decomposition
       numWritingRanks = self.read_parameter("numWritingRanks")
       if len(rawData.shape) > 1:
         orderedData = np.zeros([bbox[0],bbox[1],bbox[2],rawData.shape[1]])
       else:
         orderedData = np.zeros([bbox[0],bbox[1],bbox[2]])

       # Helper functions ported from c++ (fsgrid.hpp)
       def computeDomainDecomposition(globalsize, ntasks):
           processDomainDecomposition = [1,1,1]
           processBox = [0,0,0]
           optimValue = 999999999999999.
           for i in range(1,min(ntasks,globalsize[0]+1)):
               processBox[0] = max(1.*globalsize[0]/i,1)
               for j in range(1,min(ntasks,globalsize[1]+1)):
                   if(i * j > ntasks):
                       break
                   processBox[1] = max(1.*globalsize[1]/j,1)
                   for k in range(1,min(ntasks,globalsize[2]+1)):
                       if(i * j * k > ntasks):
                           continue
                       processBox[2] = max(1.*globalsize[2]/k,1)
                       value = 10 * processBox[0] * processBox[1] * processBox[2] + \
                        ((processBox[1] * processBox[2]) if i>1 else 0) + \
                        ((processBox[0] * processBox[2]) if j>1 else 0) + \
                        ((processBox[0] * processBox[1]) if k>1 else 0)
                       if i*j*k == ntasks:
                           if value < optimValue:
                              optimValue = value
                              processDomainDecomposition=[i,j,k]
           return processDomainDecomposition

       def calcLocalStart(globalCells, ntasks, my_n):
           n_per_task = globalCells//ntasks
           remainder = globalCells%ntasks
           if my_n < remainder:
               return my_n * (n_per_task+1)
           else:
               return my_n * n_per_task + remainder
       def calcLocalSize(globalCells, ntasks, my_n):
           n_per_task = globalCells//ntasks
           remainder = globalCells%ntasks
           if my_n < remainder:
               return n_per_task+1;
           else:
               return n_per_task

       currentOffset = 0;
       fsgridDecomposition = computeDomainDecomposition([bbox[0],bbox[1],bbox[2]],numWritingRanks)
       for i in range(0,numWritingRanks):
           x = (i // fsgridDecomposition[2]) // fsgridDecomposition[1]
           y = (i // fsgridDecomposition[2]) % fsgridDecomposition[1]
           z = i % fsgridDecomposition[2]
 	   
           thatTasksSize = [calcLocalSize(bbox[0], fsgridDecomposition[0], x), \
                            calcLocalSize(bbox[1], fsgridDecomposition[1], y), \
                            calcLocalSize(bbox[2], fsgridDecomposition[2], z)]
           thatTasksStart = [calcLocalStart(bbox[0], fsgridDecomposition[0], x), \
                             calcLocalStart(bbox[1], fsgridDecomposition[1], y), \
                             calcLocalStart(bbox[2], fsgridDecomposition[2], z)]
           thatTasksEnd = np.array(thatTasksStart) + np.array(thatTasksSize)

           totalSize = int(thatTasksSize[0]*thatTasksSize[1]*thatTasksSize[2])
           # Extract datacube of that task... 
           if len(rawData.shape) > 1:
               thatTasksData = rawData[currentOffset:currentOffset+totalSize,:]
               thatTasksData = thatTasksData.reshape([thatTasksSize[0],thatTasksSize[1],thatTasksSize[2],rawData.shape[1]],order='F')

               # ... and put it into place 
               orderedData[thatTasksStart[0]:thatTasksEnd[0],thatTasksStart[1]:thatTasksEnd[1],thatTasksStart[2]:thatTasksEnd[2],:] = thatTasksData
           else:
               # Special case for scalar data
               thatTasksData = rawData[currentOffset:currentOffset+totalSize]
               thatTasksData = thatTasksData.reshape([thatTasksSize[0],thatTasksSize[1],thatTasksSize[2]], order='F')

               # ... and put it into place 
               orderedData[thatTasksStart[0]:thatTasksEnd[0],thatTasksStart[1]:thatTasksEnd[1],thatTasksStart[2]:thatTasksEnd[2]] = thatTasksData

           currentOffset += totalSize

       return np.squeeze(orderedData)

   def read_variable(self, name, cellids=-1,operator="pass"):
      ''' Read variables from the open vlsv file. 
      Arguments:
      :param name: Name of the variable
      :param cellids: a value of -1 reads all data
      :param operator: Datareduction operator. "pass" does no operation on data
      :returns: numpy array with the data

      .. seealso:: :func:`read` :func:`read_variable_info`
      '''
      cellids = get_data(cellids)

      # Wrapper, check if requesting an fsgrid variable
      if (self.check_variable(name) and (name.lower()[0:3]=="fg_")):
         if not cellids == -1:
            print("Warning, CellID requests not supported for FSgrid variables! Aborting.")
            return False
         return self.read_fsgrid_variable(name=name, operator=operator)

      # Passes the list of cell id's onwards - optimization for reading is done in the lower level read() method
      return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", operator=operator, cellids=cellids)

   def read_variable_info(self, name, cellids=-1, operator="pass"):
      ''' Read variables from the open vlsv file and input the data into VariableInfo

      :param name: Name of the variable
      :param cellids: a value of -1 reads all data
      :param operator: Datareduction operator. "pass" does no operation on data
      :returns: numpy array with the data

      .. seealso:: :func:`read_variable`
      '''
      from variable import VariableInfo

      # Force lowercase
      name = name.lower()

      # Get population and variable names from data array name 
      if '/' in name:
         popname = name.split('/')[0]
         varname = name.split('/')[1]
      else:
         popname = "pop"
         varname = name

      # Check which set of datareducers to use
      if varname[0:3]=="vg_":
         reducer_reg = v5reducers
         reducer_multipop = multipopv5reducers
      else:
         reducer_reg = datareducers
         reducer_multipop = multipopdatareducers

      if (self.check_variable(name) and (varname[0:3]=="vg_" or varname[0:3]=="fg_")):
         # For Vlasiator 5 vlsv files, metadata is included
         units, latexunits, latex, conversion = self.read_metadata(name=name)
         # Correction for early version incorrect number density (extra backslash)
         if latex[0:3]==r"$\n":
            latex = r"$n"+latex[3:]
      elif (self.check_variable(name) and (name in vlsvvariables.unitsdict)):
         units = vlsvvariables.unitsdict[name]
         latex = vlsvvariables.latexdict[name]
         latexunits = vlsvvariables.latexunitsdict[name]            
      elif name in reducer_reg:
         units = reducer_reg[name].units
         latex = reducer_reg[name].latex
         latexunits = reducer_reg[name].latexunits
      elif 'pop/'+varname in reducer_multipop:
         poplatex='i'
         if popname in vlsvvariables.speciesdict:
            poplatex = vlsvvariables.speciesdict[popname]
         units = reducer_multipop['pop/'+varname].units
         latex = (reducer_multipop['pop/'+varname].latex).replace('REPLACEPOP',poplatex)
         latexunits = reducer_multipop['pop/'+varname].latexunits
      elif varname in vlsvvariables.unitsdict:
         poplatex='i'
         if popname in vlsvvariables.speciesdict:
            poplatex = vlsvvariables.speciesdict[popname]
         units = vlsvvariables.unitsdict[varname]
         latex = vlsvvariables.latexdictmultipop[varname].replace('REPLACEPOP',poplatex)
         latexunits = vlsvvariables.latexunitsdict[varname]
      else:
         units = ""
         latex = r""+name.replace("_","\_")
         latexunits = ""

      if name.startswith('fg_'):
          data = self.read_fsgrid_variable(name=name, operator=operator)
      else:
          data = self.read_variable(name=name, operator=operator, cellids=cellids)

      if operator != "pass":
         if operator=="magnitude":
            latex = r"$|$"+latex+r"$|$"
         else:
            latex = latex+r"${_{"+operator+r"}}$"
         return VariableInfo(data_array=data, name=name + "_" + operator, units=units, latex=latex, latexunits=latexunits)
      else:
         return VariableInfo(data_array=data, name=name, units=units, latex=latex, latexunits=latexunits)


   def get_max_refinement_level(self):
      ''' Returns the maximum refinement level of the AMR
      '''
      if self.__max_spatial_amr_level < 0:
         # Read the file index for cellid
         cellids=self.read(mesh="SpatialGrid",name="CellID", tag="VARIABLE")
         maxcellid = np.amax([cellids])

         AMR_count = 0
         while (maxcellid > 0):
            maxcellid -= 2**(3*(AMR_count))*(self.__xcells*self.__ycells*self.__zcells)
            AMR_count += 1
            
         self.__max_spatial_amr_level = AMR_count - 1
      return self.__max_spatial_amr_level
      
   def get_amr_level(self,cellid):
      '''Returns the AMR level of a given cell defined by its cellid
      
      :param cellid:        The cell's cellid
      :returns:             The cell's refinement level in the AMR
      '''
      AMR_count = 0
      while (cellid > 0):
         cellid -= 2**(3*(AMR_count))*(self.__xcells*self.__ycells*self.__zcells)
         AMR_count += 1
      return AMR_count - 1 

   def get_unique_cellids(self, coords):
      ''' Returns a list of cellids containing all the coordinates in coords,
          with no duplicate cellids. Relative order of elements is conserved.
      :param coords:         A list of coordinates
      :returns: a list of unique cell ids
      '''
      cids = [int(self.get_cellid(coord)) for coord in coords]

      #choose unique cids, keep ordering. This requires a bit of OrderedDict magic (python 2.7+)
      cidsout = list(OrderedDict.fromkeys(cids))
      return cidsout

   def get_cellid(self, coordinates):
      ''' Returns the cell id at given coordinates

      :param coordinates:        The cell's coordinates
      :returns: the cell id

      .. note:: Returns 0 if the cellid is out of bounds!
      '''
      # If needed, read the file index for cellid
      if len(self.__fileindex_for_cellid) == 0:
         self.__read_fileindex_for_cellid()

      # Check that the coordinates are not out of bounds:
      if (self.__xmax < coordinates[0]) or (self.__xmin >= coordinates[0]):
         return 0
      if (self.__ymax < coordinates[1]) or (self.__ymin >= coordinates[1]):
         return 0
      if (self.__zmax < coordinates[2]) or (self.__zmin >= coordinates[2]):
         return 0
      # Get cell lengths:
      cell_lengths = np.array([self.__dx, self.__dy, self.__dz])

      # Get cell indices:
      cellindices = np.array([(int)((coordinates[0] - self.__xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - self.__ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - self.__zmin)/(float)(cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * self.__xcells + cellindices[2] * self.__xcells * self.__ycells + 1

      # Going through AMR levels as needed
      AMR_count = 0
      ncells_lowerlevel = 0
      refmax = self.get_max_refinement_level()

      while AMR_count < refmax + 1:
          try:
              self.__fileindex_for_cellid[cellid]
              return cellid
          except:
              ncells_lowerlevel += 2**(3*AMR_count)*(self.__xcells*self.__ycells*self.__zcells) # Increment of cellID from lower lvl             
              AMR_count += 1
              # Get cell lengths:
              cell_lengths = np.array([self.__dx, self.__dy, self.__dz]) / 2**AMR_count # Check next AMR level

              # Get cell indices:
              cellindices = np.array([(int)((coordinates[0] - self.__xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - self.__ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - self.__zmin)/(float)(cell_lengths[2]))])
              # Get the cell id:
              cellid = ncells_lowerlevel + cellindices[0] + 2**(AMR_count)*self.__xcells*cellindices[1] + 4**(AMR_count)*self.__xcells*self.__ycells*cellindices[2] + 1

          if AMR_count == refmax + 1:
              raise Exception('CellID does not exist in any AMR level')

   def get_cell_coordinates(self, cellid):
      ''' Returns a given cell's coordinates as a numpy array

      :param cellid:            The cell's ID
      :returns: a numpy array with the coordinates

      .. seealso:: :func:`get_cellid`

      .. note:: The cell ids go from 1 .. max not from 0
      '''
      # Get cell lengths:
      xcells = self.__xcells
      ycells = self.__ycells
      zcells = self.__zcells
      cellid = (int)(cellid - 1)
      # Handle AMR
      cellscount = self.__xcells*self.__ycells*self.__zcells
      reflevel=0
      subtraction = (int)(cellscount * (2**reflevel)**3)
      while (cellid >= subtraction):
         cellid -= subtraction
         reflevel += 1
         subtraction = (int)(cellscount * (2**reflevel)**3)
         xcells *= 2
         ycells *= 2
         zcells *= 2
      # Get cell indices:
      cellindices = np.zeros(3)
      cellindices[0] = (int)(cellid)%(int)(xcells)
      cellindices[1] = ((int)(cellid)//(int)(xcells))%(int)(ycells)
      cellindices[2] = (int)(cellid)//(int)(xcells*ycells)
      # cellindices[0] = (int)(cellid)%(int)(self.__xcells)
      # cellindices[1] = ((int)(cellid)//(int)(self.__xcells))%(int)(self.__ycells)
      # cellindices[2] = (int)(cellid)//(int)(self.__xcells*self.__ycells)
   
      # Get cell coordinates:
      cell_lengths = np.array([(self.__xmax - self.__xmin)/(float)(xcells), (self.__ymax - self.__ymin)/(float)(ycells), (self.__zmax - self.__zmin)/(float)(zcells)])
      cellcoordinates = np.zeros(3)
      cellcoordinates[0] = self.__xmin + (cellindices[0] + 0.5) * cell_lengths[0]
      cellcoordinates[1] = self.__ymin + (cellindices[1] + 0.5) * cell_lengths[1]
      cellcoordinates[2] = self.__zmin + (cellindices[2] + 0.5) * cell_lengths[2]
      # Return the coordinates:
      return np.array(cellcoordinates)

   def get_cell_indices(self, cellid, reflevel):
      ''' Returns a given cell's indices as a numpy array

      :param cellid:            The cell's ID
      :param reflevel:          The cell's refinement level in the AMR
      :returns: a numpy array with the coordinates

      .. seealso:: :func:`get_cellid`

      .. note:: The cell ids go from 1 .. max not from 0
      '''
      # Calculating the index of the first cell at this reflevel
      index_at_reflevel = 0
      for i in range(0,reflevel):
         index_at_reflevel += 2**(3*i) * self.__xcells * self.__ycells * self.__zcells

      # Get cell indices:
      cellid = (int)(cellid - 1 - index_at_reflevel)
      cellindices = np.zeros(3)
      cellindices[0] = (int)(cellid)%(int)(2**reflevel*self.__xcells)
      cellindices[1] = ((int)(cellid)//(int)(2**reflevel*self.__xcells))%(int)(2**reflevel*self.__ycells)
      cellindices[2] = (int)(cellid)//(int)(4**reflevel*self.__xcells*self.__ycells)

      # Return the indices:
      return np.array(cellindices)

   def get_cell_neighbor(self, cellid, offset, periodic):
      ''' Returns a given cells neighbor at offset (in indices)

      :param cellid:            The cell's ID
      :param offset:            The offset to the neighbor in indices
      :param periodic:          For each dimension, is the system periodic
      :returns: the cellid of the neighbor

      .. note:: Returns 0 if the offset is out of bounds!

      '''
      reflevel = self.get_amr_level(cellid)
      indices = self.get_cell_indices(cellid, reflevel)

      # Special case if no offset      
      if ((offset[0]==0) * (offset[1]==0) * (offset[2]==0)):
         return cellid

      # Getting the neighbour at the same refinement level
      ngbr_indices = np.zeros(3)
      sys_size = [2**reflevel*self.__xcells, 2**reflevel*self.__ycells, 2**reflevel*self.__zcells]
      for i in range(3):
         ngbr_indices[i] = indices[i] + offset[i]
         if periodic[i]:
            for j in range(abs(offset[i])):
               #loop over offset abs as offset may be larger than the system size
               if ngbr_indices[i] < 0:
                  ngbr_indices[i] = ngbr_indices[i] + sys_size[i]
               elif ngbr_indices[i] >= sys_size[i]:
                  ngbr_indices[i] = ngbr_indices[i] - sys_size[i]
   
         elif ngbr_indices[i] < 0 or  ngbr_indices[i] >= sys_size[i]:
            print("Error in Vlsvreader get_cell_neighbor: out of bounds")
            return 0

      coord_neighbour = np.array([self.__xmin,self.__ymin,self.__zmin]) + (ngbr_indices + np.array((0.5,0.5,0.5))) * np.array([self.__dx,self.__dy,self.__dz])/2**reflevel
      cellid_neighbour = self.get_cellid(coord_neighbour)
      return cellid_neighbour

   def get_velocity_cell_ids(self, vcellcoord, pop="proton"):
      ''' Returns velocity cell ids of given coordinate

      Arguments:
      :param vcellcoords: One 3d coordinate
      :returns: Velocity cell id

      .. seealso:: :func:`get_velocity_cell_coordinates`
      '''
      vmin = np.array([self.__meshes[pop].__vxmin, self.__meshes[pop].__vymin, self.__meshes[pop].__vzmin])
      dv = np.array([self.__meshes[pop].__dvx, self.__meshes[pop].__dvy, self.__meshes[pop].__dvz])
      block_index = np.floor((vcellcoord - vmin) / (4 * dv))
      cell_index = np.floor(np.remainder(vcellcoord - vmin, 4 * dv) / dv)
      vcellid = int(block_index[0])
      vcellid += int(block_index[1] * self.__meshes[pop].__vxblocks)
      vcellid += int(block_index[2] * self.__meshes[pop].__vxblocks * self.__meshes[pop].__vyblocks)
      vcellid *= 64
      vcellid += int(cell_index[0])
      vcellid += int(cell_index[1] * 4)
      vcellid += int(cell_index[2] * 16)
      return vcellid

   def get_velocity_cell_coordinates(self, vcellids, pop="proton"):
      ''' Returns a given velocity cell's coordinates as a numpy array

      Arguments:
      :param vcellids:       The velocity cell's ID
      :returns: a numpy array with the coordinates

      .. seealso:: :func:`get_cell_coordinates` :func:`get_velocity_block_coordinates`
      '''
      vcellids = np.atleast_1d(vcellids)
      # Get block ids:
      blocks = vcellids.astype(int) // 64
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__meshes[pop].__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)//(int)(self.__meshes[pop].__vxblocks), (int)(self.__meshes[pop].__vyblocks))
      blockIndicesZ = blocks.astype(int)//(int)(self.__meshes[pop].__vxblocks*self.__meshes[pop].__vyblocks)
      blockCoordinatesX = blockIndicesX.astype(float) * self.__meshes[pop].__dvx * 4 + self.__meshes[pop].__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__meshes[pop].__dvy * 4 + self.__meshes[pop].__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__meshes[pop].__dvz * 4 + self.__meshes[pop].__vzmin
      # Get cell indices:
      cellids = np.remainder(vcellids.astype(int), (int)(64))
      cellIndicesX = np.remainder(cellids.astype(int), (int)(4))
      cellIndicesY = np.remainder((cellids.astype(int)//(int)(4)).astype(int), (int)(4))
      cellIndicesZ = cellids.astype(int)//(int)(16)
      # Get cell coordinates:
      cellCoordinates = np.array([blockCoordinatesX.astype(float) + (cellIndicesX.astype(float) + 0.5) * self.__meshes[pop].__dvx,
                                  blockCoordinatesY.astype(float) + (cellIndicesY.astype(float) + 0.5) * self.__meshes[pop].__dvy,
                                  blockCoordinatesZ.astype(float) + (cellIndicesZ.astype(float) + 0.5) * self.__meshes[pop].__dvz])

      return cellCoordinates.transpose()

   def get_velocity_block_coordinates( self, blocks, pop="proton"):
      ''' Returns the block coordinates of the given blocks in a numpy array

          :param blocks:         list of block ids
          :returns: a numpy array containing the block coordinates e.g. np.array([np.array([2,1,3]), np.array([5,6,6]), ..])

          .. seealso:: :func:`get_velocity_cell_coordinates`
      '''
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__meshes[pop].__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)//(int)(self.__meshes[pop].__vxblocks), (int)(self.__meshes[pop].__vyblocks))
      blockIndicesZ = blocks.astype(int)//(int)(self.__meshes[pop].__vxblocks*self.__meshes[pop].__vyblocks)
      blockCoordinatesX = blockIndicesX.astype(float) * self.__meshes[pop].__dvx * 4 + self.__meshes[pop].__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__meshes[pop].__dvy * 4 + self.__meshes[pop].__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__meshes[pop].__dvz * 4 + self.__meshes[pop].__vzmin
      # Return the coordinates:
      return np.array([blockCoordinatesX.astype(float),
                       blockCoordinatesY.astype(float),
                       blockCoordinatesZ.astype(float)]).transpose()

   def get_velocity_blocks( self, blockcoordinates, pop="proton" ):
      ''' Returns the block ids of the given block coordinates in a numpy array form

          :param blockcoordinates:         list of block coordinates e.g. np.array([np.array([2,1,3]), np.array([5,6,6]), ..])
          :returns: a numpy array containing the block ids e.g. np.array([4,2,56,44,2, ..])

          .. seealso:: :func:`get_velocity_block_coordinates`
      '''
      mins = np.array([self.__meshes[pop].__vxmin, self.__meshes[pop].__vymin, self.__meshes[pop].__vzmin]).astype(float)
      dvs = np.array([4*self.__meshes[pop].__dvx, 4*self.__meshes[pop].__dvy, 4*self.__meshes[pop].__dvz]).astype(float)
      multiplier = np.array([1, self.__meshes[pop].__vxblocks, self.__meshes[pop].__vxblocks * self.__meshes[pop].__vyblocks]).astype(float)
      velocity_block_ids = np.sum(np.floor(((blockCoordinates.astype(float) - mins) / dvs)) * multiplier, axis=-1)
      return velocity_block_ids

   def construct_velocity_cells( self, blocks ):
      ''' Returns velocity cells in given blocks

          :param blocks:         list of block ids
          :returns: a numpy array containing the velocity cell ids e.g. np.array([4,2,56,44,522, ..])
      '''
      return np.ravel(np.outer(np.array(blocks), np.ones(64)) + np.arange(64))

   def construct_velocity_cell_coordinates( self, blocks ):
      ''' Returns velocity cell coordinates in given blocks

          :param blocks:         list of block ids
          :returns: a numpy array containing the velocity cell ids e.g. np.array([4,2,56,44,522, ..])
      '''
      # Construct velocity cell coordinates from velocity cells and return them
      return self.get_velocity_cell_coordinates( self.construct_velocity_cells(blocks) )


   def construct_velocity_cell_nodes( self, blocks, pop="proton" ):
      ''' Returns velocity cell nodes in given blocks

          :param blocks:         list of block ids
          :returns: a numpy array containing velocity cell nodes and the keys for velocity cells

          .. note:: This is used for constructing velocity space inside the mayavi module

          .. seealso:: :mod:`grid`
      '''
      blocks = np.array(blocks)
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__meshes[pop].__vxblocks)).astype(np.uint16)
      blockIndicesY = np.remainder(blocks.astype(int)//(int)(self.__meshes[pop].__vxblocks), (int)(self.__meshes[pop].__vyblocks)).astype(np.uint16)
      blockIndicesZ = (blocks.astype(np.uint64)//(int)(self.__meshes[pop].__vxblocks*self.__meshes[pop].__vyblocks)).astype(np.uint16)

      cellsPerDirection = 4
      cellsPerBlock = 64

      # Get velocity cell min coordinates (per velocity block)
      vcellids = np.arange(cellsPerBlock).astype(np.uint32)
      cellIndicesX = np.remainder(vcellids.astype(int), (int)(cellsPerDirection)).astype(np.uint16)
      cellIndicesY = np.remainder((vcellids.astype(int)//(int)(cellsPerDirection)).astype(int), (int)(cellsPerDirection)).astype(np.uint16)
      cellIndicesZ = (vcellids.astype(int)//(int)(cellsPerDirection*cellsPerDirection)).astype(np.uint16)

      # Construct velocity cell node indices for every velocity cell per velocity block

      nodesPerCell = 8

      # NOTE: The ordering of the numpy array won't make sense to anyone who hasn't read VTK documentation. For further info check VTK_VOXEL. The numpy array is constructed according to VTK voxel's nodes
      cellNodeIndicesX = np.ravel(np.outer(cellIndicesX, np.ones(nodesPerCell)) + np.array([0, 1, 0, 1, 0, 1, 0, 1])).astype(np.uint16)
      cellNodeIndicesY = np.ravel(np.outer(cellIndicesY, np.ones(nodesPerCell)) + np.array([0, 0, 1, 1, 0, 0, 1, 1])).astype(np.uint16)
      cellNodeIndicesZ = np.ravel(np.outer(cellIndicesZ, np.ones(nodesPerCell)) + np.array([0, 0, 0, 0, 1, 1, 1, 1])).astype(np.uint16)

      nodeIndices_local = []
      nodesPerDirection = 5

      for i in range(nodesPerDirection):
         for j in range(nodesPerDirection):
            for k in range(nodesPerDirection):
               nodeIndices_local.append(np.array([i,j,k]))
      nodeIndices_local = np.array(nodeIndices_local).astype(np.uint16)

      nodesPerBlock = (int)(nodesPerDirection * nodesPerDirection * nodesPerDirection)


      def calculate_node_indices( self, blockIndicesX, blockIndicesY, blockIndicesZ, nodeIndices_local, nodesPerBlock, cellsPerDirection ):
         nodeIndicesX = np.ravel(np.outer(blockIndicesX, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,0])
         nodeIndicesY = np.ravel(np.outer(blockIndicesY, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,1])
         nodeIndicesZ = np.ravel(np.outer(blockIndicesZ, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,2])
   
         nodeIndices = np.transpose(np.array([nodeIndicesX, nodeIndicesY, nodeIndicesZ], copy=False))

         # Transform indices into unique keys
         nodeKeys = np.sum(nodeIndices * np.array([1, cellsPerDirection*self.__meshes[pop].__vxblocks+1, (cellsPerDirection*self.__meshes[pop].__vxblocks+1)*(cellsPerDirection*self.__meshes[pop].__vyblocks+1)]), axis=1)
         # Sort the keys and delete duplicates
         return np.unique(nodeKeys)
      #nodeIndices = calculate_node_indices( blockIndicesX, blockIndicesY, blockIndicesZ, nodeIndices_local, nodesPerBlock, cellsPerDirection )

      # Put the node indices into keys:
      nodeKeys = np.array([], dtype=np.uint64)
      N = 10
      for i in range(N):
         fromIndex = i*(len(blockIndicesX)//N)
         if i != N-1:
            toIndex = (i+1)*(len(blockIndicesX)//N)
         else:
            toIndex = len(blockIndicesX)
         nodeKeys = np.append(nodeKeys, calculate_node_indices( self, blockIndicesX[fromIndex:toIndex], blockIndicesY[fromIndex:toIndex], blockIndicesZ[fromIndex:toIndex], nodeIndices_local, nodesPerBlock, cellsPerDirection ) )


      # Delete duplicate nodes and sort the list:
      nodeKeys = np.unique(nodeKeys) #We now have all of the nodes in a list!




      def calc_global_cell_keys( self, blockIndicesX, blockIndicesY, blockIndicesZ, cellNodeIndicesX, cellNodeIndicesY, cellNodeIndicesZ, cellsPerBlock, nodesPerCell, cellsPerDirection, nodeKeys ):
         # reate node  indices for the cells
         globalCellIndicesX = np.ravel(np.outer(blockIndicesX, np.ones(cellsPerBlock * nodesPerCell).astype(np.uint16)) * cellsPerDirection + cellNodeIndicesX)
         globalCellIndicesY = np.ravel(np.outer(blockIndicesY, np.ones(cellsPerBlock * nodesPerCell).astype(np.uint16)) * cellsPerDirection + cellNodeIndicesY)
         globalCellIndicesZ = np.ravel(np.outer(blockIndicesZ, np.ones(cellsPerBlock * nodesPerCell).astype(np.uint16)) * cellsPerDirection + cellNodeIndicesZ)
   
         globalCellIndices = np.array([globalCellIndicesX, globalCellIndicesY, globalCellIndicesZ], copy=False)
         globalCellIndices = np.transpose(globalCellIndices)
         # Transform cell indices into unique keys
         globalCellIndices = np.sum(globalCellIndices * np.array([1, cellsPerDirection*self.__meshes[pop].__vxblocks+1, (cellsPerDirection*self.__meshes[pop].__vxblocks+1)*(cellsPerDirection*self.__meshes[pop].__vyblocks+1)]), axis=1)
         # Return cell nodes' indexes in the nodeKeys list
         return np.searchsorted(nodeKeys, globalCellIndices)


      # Create cellKeys
      cellKeys = np.zeros(len(blockIndicesX)*cellsPerBlock*nodesPerCell, dtype=np.uint32)
      N = 10
      # Append keys in cuts to save memory
      for i in range(N):
         fromIndex = i*(len(blockIndicesX)//N)
         if i != N-1:
            toIndex = (i+1)*(len(blockIndicesX)//N)
         else:
            toIndex = len(blockIndicesX)
         # Append cell keys
         cellKeys[fromIndex*cellsPerBlock*nodesPerCell:toIndex*cellsPerBlock*nodesPerCell] = calc_global_cell_keys( self, blockIndicesX[fromIndex:toIndex], blockIndicesY[fromIndex:toIndex], blockIndicesZ[fromIndex:toIndex], cellNodeIndicesX, cellNodeIndicesY, cellNodeIndicesZ, cellsPerBlock, nodesPerCell, cellsPerDirection, nodeKeys )

      cellKeys = np.reshape(cellKeys, (len(blocks)*64,8))

      # We now have all the cell keys and avgs values! (avgs is in the same order as cell keys)
      # Now transform node indices back into real indices
      nodeCoordinatesX = np.remainder(nodeKeys, (int)(cellsPerDirection*self.__meshes[pop].__vxblocks+1)).astype(np.float32) * self.__meshes[pop].__dvx + self.__meshes[pop].__vxmin
      nodeCoordinatesY = np.remainder(nodeKeys//(int)(cellsPerDirection*self.__meshes[pop].__vxblocks+1), cellsPerDirection*self.__meshes[pop].__vyblocks+1).astype(np.float32) * self.__meshes[pop].__dvy + self.__meshes[pop].__vymin
      nodeCoordinatesZ = ( nodeKeys // (int)((cellsPerDirection*self.__meshes[pop].__vxblocks+1) * (cellsPerDirection*self.__meshes[pop].__vyblocks+1)) ).astype(np.float32) * self.__meshes[pop].__dvz + self.__meshes[pop].__vzmin
      
      # Nodekeyss is no longer needed
      del nodeKeys

      nodes = np.array([nodeCoordinatesX, nodeCoordinatesY, nodeCoordinatesZ], copy=False)
      # Take a transpose
      nodes = np.transpose(nodes)

      return [nodes, cellKeys]





   def read_parameter(self, name):
      ''' Read a parameter from the vlsv file

      :param name:   Name of the parameter
      :returns: The parameter value

      .. seealso:: :func:`read_variable` :func:`read_variable_info`
      '''
      
      # Special handling for time
      if name=="time":
         if self.check_parameter(name="t"):
            return self.read(name="t", tag="PARAMETER")
      if name=="t":
         if self.check_parameter(name="time"):
            return self.read(name="time", tag="PARAMETER")

      return self.read(name=name, tag="PARAMETER")


   def read_velocity_cells(self, cellid, pop="proton"):
      ''' Read velocity cells from a spatial cell
      
      :param cellid: Cell ID of the cell whose velocity cells the function will read
      :returns: Map of velocity cell ids (unique for every velocity cell) and corresponding value

      #Example:

      example_cellid = 1111

      velocity_cell_map = vlsvReader.read_velocity_cells(example_cellid)
      velocity_cell_ids = velocity_cell_map.keys()
      velocity_cell_values = velocity_cell_map.values()

      random_index = 4 # Just some index
      random_velocity_cell_id = velocity_cell_ids[random_index]

      print "Velocity cell value at velocity cell id " + str(random_velocity_cell_id) + ": " + str(velocity_cell_map[random_velocity_cell_id])

      # Getting the corresponding coordinates might be more useful than having the velocity cell id so:
      velocity_cell_coordinates = vlsvReader.get_velocity_cell_coordinates(velocity_cell_ids) # Get velocity cell coordinates corresponding to each velocity cell id

      random_velocity_cell_coordinates = velocity_cell_ids[random_index]
      print "Velocity cell value at velocity cell id " + str(random_velocity_cell_id) + "and coordinates " + str(random_velocity_cell_coordinates) + ": " + str(velocity_cell_map[random_velocity_cell_id])

      .. seealso:: :func:`read_blocks`
      '''
      
      if self.use_dict_for_blocks: # old deprecated version, uses dict for blocks data
         if not pop in self.__fileindex_for_cellid_blocks:
            self.__set_cell_offset_and_blocks(pop) 
         # Check that cells has vspace
         if not cellid in self.__fileindex_for_cellid_blocks[pop]:
            print("Cell does not have velocity distribution")
            return []
         # Navigate to the correct position:
         offset = self.__fileindex_for_cellid_blocks[pop][cellid][0]
         num_of_blocks = self.__fileindex_for_cellid_blocks[pop][cellid][1]

      else:  # Uses arrays (much faster to initialize)
         if not pop in self.__cells_with_blocks:
            self.__set_cell_offset_and_blocks_nodict(pop) 
         # Check that cells has vspace
         try:
            cells_with_blocks_index = self.__order_for_cellid_blocks[pop][cellid]
         except:
            print("Cell does not have velocity distribution")
            return []
         # Navigate to the correct position:
         offset = self.__blocks_per_cell_offsets[pop][cells_with_blocks_index]
         num_of_blocks = self.__blocks_per_cell[pop][cells_with_blocks_index]


      if self.__fptr.closed:
         fptr = open(self.file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if "name" in child.attrib and (child.attrib["name"] == pop) and (child.tag == "BLOCKVARIABLE"):
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = int(offset * vector_size * element_size + ast.literal_eval(child.text))

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)
         # Read in block coordinates:
         if ("name" in child.attrib) and (child.attrib["name"] == pop) and (child.tag == "BLOCKIDS"):
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_ids = int(offset * vector_size * element_size + ast.literal_eval(child.text))

            fptr.seek(offset_block_ids)
            if datatype == "uint" and element_size == 4:
               data_block_ids = np.fromfile(fptr, dtype = np.uint32, count = vector_size*num_of_blocks)
            elif datatype == "uint" and element_size == 8:
               data_block_ids = np.fromfile(fptr, dtype = np.uint64, count = vector_size*num_of_blocks)
            else:
               print("Error! Bad data type in blocks!")
               return

         if (pop=="avgs") and (child.tag == "BLOCKIDS"): # Old avgs files did not have the name set for BLOCKIDS
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_ids = int(offset * vector_size * element_size + ast.literal_eval(child.text))

            fptr.seek(offset_block_ids)
            if datatype == "uint" and element_size == 4:
               data_block_ids = np.fromfile(fptr, dtype = np.uint32, count = vector_size*num_of_blocks)
            elif datatype == "uint" and element_size == 8:
               data_block_ids = np.fromfile(fptr, dtype = np.uint64, count = vector_size*num_of_blocks)
            else:
               print("Error! Bad data type in blocks!")
               return

            data_block_ids = data_block_ids.reshape(num_of_blocks, vector_size)

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_ids):
         print("BAD DATA SIZES")
      # Make a dictionary (hash map) out of velocity cell ids and avgs:
      velocity_cells = {}
      array_size = len(data_avgs)

      # Construct velocity cells:
      velocity_cell_ids = []
      for kv in range(4):
         for jv in range(4):
            for iv in range(4):
               velocity_cell_ids.append(kv*16 + jv*4 + iv)

      for i in range(array_size):
         velocity_block_id = data_block_ids[i]
         avgIndex = 0
         avgs = data_avgs[i]

         for j in velocity_cell_ids + 64*velocity_block_id:
            velocity_cells[(int)(j)] = avgs[avgIndex]
            avgIndex = avgIndex + 1
      return velocity_cells

   def get_spatial_mesh_size(self):
      ''' Read spatial mesh size
      
      :returns: Size of mesh in number of blocks, array with three elements
      '''
      return np.array([self.__xcells, self.__ycells, self.__zcells])

   def get_spatial_block_size(self):
      ''' Read spatial mesh block size
      
      :returns: Size of block in number of cells, array with three elements
      '''
      return np.array([self.__xblock_size, self.__yblock_size, self.__zblock_size])

   def get_spatial_mesh_extent(self):
      ''' Read spatial mesh extent
      
      :returns: Maximum and minimum coordinates of the mesh, [xmin, ymin, zmin, xmax, ymax, zmax]
      '''
      return np.array([self.__xmin, self.__ymin, self.__zmin, self.__xmax, self.__ymax, self.__zmax])

   def get_fsgrid_mesh_size(self):
      ''' Read fsgrid mesh size
      
      :returns: Size of mesh in number of cells, array with three elements
      '''
      # Get fsgrid domain size (this can differ from vlasov grid size if refined)
      bbox = self.read(tag="MESH_BBOX", mesh="fsgrid")
      return np.array(bbox[0:3])

   def get_fsgrid_mesh_extent(self):
      ''' Read fsgrid mesh extent
      
      :returns: Maximum and minimum coordinates of the mesh, [xmin, ymin, zmin, xmax, ymax, zmax]
      '''
      return np.array([self.__xmin, self.__ymin, self.__zmin, self.__xmax, self.__ymax, self.__zmax])

   def get_velocity_mesh_size(self, pop="proton"):
      ''' Read velocity mesh size
      
      :returns: Size of mesh in number of blocks, array with three elements
      '''
      return np.array([self.__meshes[pop].__vxblocks, self.__meshes[pop].__vyblocks, self.__meshes[pop].__vzblocks])

   def get_velocity_block_size(self, pop="proton"):
      ''' Read velocity mesh block size
      
      :returns: Size of block in number of cells, array with three elements
      '''
      return np.array([self.__meshes[pop].__vxblock_size, self.__meshes[pop].__vyblock_size, self.__meshes[pop].__vzblock_size])

   def get_velocity_mesh_extent(self, pop="proton"):
      ''' Read velocity mesh extent
      
      :returns: Maximum and minimum coordinates of the mesh, [vxmin, vymin, vzmin, vxmax, vymax, vzmax]
      '''
      return np.array([self.__meshes[pop].__vxmin, self.__meshes[pop].__vymin, self.__meshes[pop].__vzmin, self.__meshes[pop].__vxmax, self.__meshes[pop].__vymax, self.__meshes[pop].__vzmax])

   def read_blocks(self, cellid, pop="proton"):
      ''' Read raw block data from the open file and return the data along with block ids
      
      :param cellid: Cell ID of the cell whose velocity blocks are read
      :returns: A numpy array with block ids and data eg [array([2, 5, 6, 234, 21]), array([1.0e-8, 2.1e-8, 2.1e-8, 0, 4.0e-8])]

      .. seealso:: :func:`read_velocity_cells`
      '''
      # Uses new format
      return self.__read_blocks(cellid,pop)

      return []

   def get_precipitation_centre_energy(self, pop="proton"):
      ''' Read precipitation energy bins

      :returns: Array of centre energies
      '''
      return self.__meshes[pop].__precipitation_centre_energy

   def optimize_open_file(self):
      '''Opens the vlsv file for reading
         Files are opened and closed automatically upon reading and in the case of reading multiple times it will help to keep the file open with this command

         .. code-block: python

            #Example usage:
            variables = []
            vlsvReader.optimize_open_file()
            for i in range(1000):
               variables.append(vlsvReader.read_variable("rho", cellids=i))
            vlsvReader.optimize_close_file()

         .. note:: This should only be used for optimization purposes.
      '''
      self.__fptr = open(self.file_name,"rb")


   def optimize_close_file(self):
      '''Closes the vlsv file
         Files are opened and closed automatically upon reading and in the case of reading multiple times it will help to keep the file open with this command

         .. code-block: python

            # Example usage:
            variables = []
            vlsvReader.optimize_open_file()
            for i in range(1000):
               variables.append(vlsvReader.read_variable("rho", cellids=i))
            vlsvReader.optimize_close_file()

         .. note:: This should only be used for optimization purposes.
      '''
      if self.__fptr.closed:
         return
      else:
         self.__fptr.close()
         return

   def optimize_clear_fileindex_for_cellid_blocks(self):
      ''' Clears a private variable containing number of blocks and offsets for particular cell ids

         .. code-block: python

             # Example usage:
             vlsvReaders = []
             # Open a list of vlsv files
             for i in range(1000):
                vlsvReaders.append( VlsvReader("test" + str(i) + ".vlsv") )
             # Go through vlsv readers and print info:
             for vlsvReader in vlsvReaders:
                # Print something from the file on the screen
                print vlsvReader.read_blocks( cellid= 5021 ) # Stores info into a private variable
                # Upon reading from vlsvReader a private variable that contains info on cells that have blocks has been saved -- now clear it to save memory
                vlsvReader.optimize_clear_fileindex_for_cellid_blocks()

         .. note:: This should only be used for optimization purposes.
      '''
      self.__fileindex_for_cellid_blocks = {}
      self.__cells_with_blocks = {}
      self.__blocks_per_cell = {}
      self.__blocks_per_cell_offsets = {}
      self.__order_for_cellid_blocks = {}

   def optimize_clear_fileindex_for_cellid(self):
      ''' Clears a private variable containing cell ids and their locations

         .. code-block: python

             # Example usage:
             vlsvReaders = []
             # Open a list of vlsv files
             for i in range(1000):
                vlsvReaders.append( VlsvReader("test" + str(i) + ".vlsv") )
             # Go through vlsv readers and print info:
             for vlsvReader in vlsvReaders:
                # Print something from the file on the screen
                print vlsvReader.read_variable("B", cellids=2) # Stores info into a private variable
                # Upon reading from vlsvReader a private variable that contains info on cells that have blocks has been saved -- now clear it to save memory
                vlsvReader.optimize_clear_fileindex_for_cellid()

         .. note:: This should only be used for optimization purposes.
      '''
      self.__fileindex_for_cellid = {}



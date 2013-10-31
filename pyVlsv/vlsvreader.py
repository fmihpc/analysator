import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
from reduction import datareducers,data_operators
from collections import Iterable

class VlsvFile(object):
   ''' Class for reading VLSV files
   ''' 
   def __init__(self, file_name):
      self.__file_name = file_name
      self.__fptr = open(self.__file_name,"rb")
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}
      self.__fileindex_for_cellid_blocks={}
      self.__read_xml_footer()
      # Check if the file is using new or old vlsv format
      self.__uses_new_vlsv_format = False
      for child in self.__xml_root:
         if child.tag == "PARAMETER":
           if child.attrib["name"] == "version":
             self.__uses_new_vlsv_format = True
      #self.__read_fileindex_for_cellid()
      # Read parameters (Note: Reading the spatial cell locations and storing them will anyway take the most time and memory):
      self.__vxblocks = (int)(self.read_parameter("vxblocks_ini"))
      self.__vyblocks = (int)(self.read_parameter("vyblocks_ini"))
      self.__vzblocks = (int)(self.read_parameter("vzblocks_ini"))

      self.__xcells = (int)(self.read_parameter("xcells_ini"))
      self.__ycells = (int)(self.read_parameter("ycells_ini"))
      self.__zcells = (int)(self.read_parameter("zcells_ini"))

      self.__xmin = self.read_parameter("xmin")
      self.__ymin = self.read_parameter("ymin")
      self.__zmin = self.read_parameter("zmin")
      self.__xmax = self.read_parameter("xmax")
      self.__ymax = self.read_parameter("ymax")
      self.__zmax = self.read_parameter("zmax")

      self.__vxmin = self.read_parameter("vxmin")
      self.__vymin = self.read_parameter("vymin")
      self.__vzmin = self.read_parameter("vzmin")
      self.__vxmax = self.read_parameter("vxmax")
      self.__vymax = self.read_parameter("vymax")
      self.__vzmax = self.read_parameter("vzmax")

      self.__dx = (self.__xmax - self.__xmin) / (float)(self.__xcells)
      self.__dy = (self.__ymax - self.__ymin) / (float)(self.__ycells)
      self.__dz = (self.__zmax - self.__zmin) / (float)(self.__zcells)

      # Velocity cell lengths
      velocity_cells_per_direction = 4
      self.__dvx = ((self.__vxmax - self.__vxmin) / (float)(self.__vxblocks)) / (float)(velocity_cells_per_direction)
      self.__dvy = ((self.__vymax - self.__vymin) / (float)(self.__vyblocks)) / (float)(velocity_cells_per_direction)
      self.__dvz = ((self.__vzmax - self.__vzmin) / (float)(self.__vzblocks)) / (float)(velocity_cells_per_direction)

      self.__fptr.close()


   def __read_xml_footer(self):
      ''' Reads in the XML footer of the VLSV file and store all the content
      ''' 
      max_xml_size = 1000000
      #(endianness,) = struct.unpack("c", fptr.read(1))
      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr
      fptr.seek(8)
      (offset,) = struct.unpack("Q", fptr.read(8))
      fptr.seek(offset)
      xml_data = fptr.read(max_xml_size)
      (xml_string,) = struct.unpack("%ds" % len(xml_data), xml_data)
      self.__xml_root = ET.fromstring(xml_string)
      if self.__fptr.closed:
         fptr.close()

   def __read_fileindex_for_cellid(self):
      """ Read in the cell ids and create an internal dictionary to give the index of an arbitrary cellID
      """
      if self.__uses_new_vlsv_format == True:
         cellids=self.read(mesh="SpatialGrid",name="CellID", tag="VARIABLE")
      else:
         cellids=self.read(name="SpatialGrid",tag="MESH")
      #Check if it is not iterable. If it is a scale then make it a list
      if(not isinstance(cellids, Iterable)):
         cellids=[ cellids ]
      for index,cellid in enumerate(cellids):
         self.__fileindex_for_cellid[cellid]=index
         
   def __list_old_format(self):
      ''' Print out a description of the content of the file. Useful
         for interactive usage
      '''
      print "tag = PARAMETERS"
      for child in self.__xml_root:
         if child.tag == "PARAMETERS":
            print "   ", child.attrib["name"], " = ", child.attrib["value"]
      print "tag = VARIABLE"
      for child in self.__xml_root:
         if child.tag == "VARIABLE":
            print "   ", child.attrib["name"]
      print "Datareducers:"
      for name in datareducers:
         print "   ",name, " based on ", datareducers[name].variables
      print "Datareduction operators:"
      for name in data_operators:
         print "   ",name

      print "Other:"
      for child in self.__xml_root:
         if child.tag != "PARAMETERS" and child.tag != "VARIABLE":
            print "    tag = ", child.tag, " name = ", child.attrib["name"]
      return

   def __read_parameter_old_format(self, name):
      return self.read(name=name, tag="PARAMETERS")

   def __read_blocks_old_format(self, cellid):
      ''' Read raw block data from the open file.
      
      Arguments:
      :param cellid Cell ID of the cell whose velocity blocks are read
      :returns numpy array with block ids and their data
      '''
      if( len(self.__fileindex_for_cellid_blocks) == 0 ):
         self.__set_cell_offset_and_blocks()
      if( (cellid in self.__fileindex_for_cellid_blocks) == False ):
         # Cell id has no blocks
         return []
      offset = self.__fileindex_for_cellid_blocks[cellid][0]
      num_of_blocks = self.__fileindex_for_cellid_blocks[cellid][1]

      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = offset * vector_size * element_size + ast.literal_eval(child.text)
#            for i in range(0, cells_with_blocks_index[0]):
#               offset_avgs += blocks_per_cell[i]*vector_size*element_size

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)

         # Read in block coordinates:
         if child.attrib["name"] == "SpatialGrid" and child.tag == "BLOCKCOORDINATES":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_coordinates = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr.seek(offset_block_coordinates)
            if datatype == "float" and element_size == 4:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)

            data_block_coordinates = data_block_coordinates.reshape(num_of_blocks, vector_size)

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_coordinates):
         print "BAD DATA SIZES"
      # Make a dictionary (hash map) out of velocity cell ids and avgs:
      velocity_cells = {}
      array_size = len(data_avgs)

      mins = np.array([self.__vxmin, self.__vymin, self.__vzmin]).astype(float)
      dvs = np.array([4*self.__dvx, 4*self.__dvy, 4*self.__dvz]).astype(float)
      multiplier = np.array([1, self.__vxblocks, self.__vxblocks * self.__vyblocks]).astype(float)
      velocity_block_ids = np.sum(np.floor(((data_block_coordinates.astype(float)[:,0:3] - mins) / dvs)) * multiplier, axis=-1)

      return [velocity_block_ids.astype(int), data_avgs.astype(float)]

   def __read_blocks_new_format(self, cellid):
      ''' Read raw block data from the open file.
      
      Arguments:
      :param cellid Cell ID of the cell whose velocity blocks are read
      :returns numpy array with block ids and their data
      '''
      if( len(self.__fileindex_for_cellid_blocks) == 0 ):
         self.__set_cell_offset_and_blocks()

      if( (cellid in self.__fileindex_for_cellid_blocks) == False ):
         # Cell id has no blocks
         return []

      num_of_blocks = self.__fileindex_for_cellid_blocks[1]
      offset = self.__fileindex_for_cellid_blocks[0]

      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = offset * vector_size * element_size + ast.literal_eval(child.text)
#            for i in range(0, cells_with_blocks_index[0]):
#               offset_avgs += blocks_per_cell[i]*vector_size*element_size

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)

         # Read in block coordinates:
         if child.attrib["mesh"] == "SpatialGrid" and child.tag == "BLOCKIDS":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_ids = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr.seek(offset_block_ids)
            if datatype == "float" and element_size == 4:
               data_block_ids = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_block_ids = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)

            data_block_ids = data_block_ids.reshape(num_of_blocks, vector_size)

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_ids):
         print "BAD DATA SIZES"

      return [data_block_ids, data_avgs]



   def __read_velocity_cells_old_format( self, cellid, cells_with_blocks, blocks_per_cell, cells_with_blocks_index  ):
      # Read in the coordinates:
      #block_coordinates = self.read(name="",tag="BLOCKCOORDINATES")
      # Navigate to the correct position:
      offset = 0
      for i in xrange(0, cells_with_blocks_index[0]):
         offset += blocks_per_cell[i]

      num_of_blocks = np.atleast_1d(blocks_per_cell)[cells_with_blocks_index[0]]

      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = offset * vector_size * element_size + ast.literal_eval(child.text)
#            for i in range(0, cells_with_blocks_index[0]):
#               offset_avgs += blocks_per_cell[i]*vector_size*element_size

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)
         # Read in block coordinates:
         if child.attrib["name"] == "SpatialGrid" and child.tag == "BLOCKCOORDINATES":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_coordinates = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr.seek(offset_block_coordinates)
            if datatype == "float" and element_size == 4:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_block_coordinates = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)

            data_block_coordinates = data_block_coordinates.reshape(num_of_blocks, vector_size)

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_coordinates):
         print "BAD DATA SIZES"
      # Make a dictionary (hash map) out of velocity cell ids and avgs:
      velocity_cells = {}
      array_size = len(data_avgs)

      # Construct velocity cells:
      velocity_cell_ids = []
      for kv in xrange(4):
         for jv in xrange(4):
            for iv in xrange(4):
               velocity_cell_ids.append(kv*16 + jv*4 + iv)

      for i in xrange(array_size):
         block_coordinate = data_block_coordinates[i]
         # The minimum corner coordinates of the blocks
         vx = block_coordinate[0]
         vy = block_coordinate[1]
         vz = block_coordinate[2]
         # The diff in blocks
         avgs = data_avgs[i]
         # Get the velocity cell id (First transform coordinates to block indices, then block indices to block id and then block id to velocity cell ids):
         velocity_block_indices = np.array([np.floor((vx - self.__vxmin) / (4*self.__dvx)), np.floor((vy - self.__vymin) / (4*self.__dvy)), np.floor((vz - self.__vzmin) / (4*self.__dvz))])
         velocity_block_id = velocity_block_indices[0] + velocity_block_indices[1] * self.__vxblocks + velocity_block_indices[2] * self.__vxblocks * self.__vyblocks
         avgIndex = 0

         for j in velocity_cell_ids + 64*velocity_block_id:
            velocity_cells[(int)(j)] = avgs[avgIndex]
            avgIndex = avgIndex + 1
      return velocity_cells

   def __read_velocity_cells_new_format( self, cellid, cells_with_blocks, blocks_per_cell, cells_with_blocks_index  ):
      # Read in the coordinates:
      #block_ids = self.read(name="",tag="BLOCKCOORDINATES")
      # Navigate to the correct position:
      offset = 0
      for i in xrange(0, cells_with_blocks_index[0]):
         offset += blocks_per_cell[i]

      num_of_blocks = np.atleast_1d(blocks_per_cell)[cells_with_blocks_index[0]]

      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr

      # Read in avgs and velocity cell ids:
      for child in self.__xml_root:
         # Read in avgs
         if child.attrib["name"] == "avgs":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            # Navigate to the correct position
            offset_avgs = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr.seek(offset_avgs)
            if datatype == "float" and element_size == 4:
               data_avgs = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_avgs = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)
            data_avgs = data_avgs.reshape(num_of_blocks, vector_size)
         # Read in block coordinates:
         if child.attrib["mesh"] == "SpatialGrid" and child.tag == "BLOCKIDS":
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            #array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]

            offset_block_ids = offset * vector_size * element_size + ast.literal_eval(child.text)

            fptr.seek(offset_block_ids)
            if datatype == "float" and element_size == 4:
               data_block_ids = np.fromfile(fptr, dtype = np.float32, count = vector_size*num_of_blocks)
            if datatype == "float" and element_size == 8:
               data_block_ids = np.fromfile(fptr, dtype = np.float64, count = vector_size*num_of_blocks)

            data_block_ids = data_block_ids.reshape(num_of_blocks, vector_size)

      if self.__fptr.closed:
         fptr.close()

      # Check to make sure the sizes match (just some extra debugging)
      if len(data_avgs) != len(data_block_ids):
         print "BAD DATA SIZES"
      # Make a dictionary (hash map) out of velocity cell ids and avgs:
      velocity_cells = {}
      array_size = len(data_avgs)

      # Construct velocity cells:
      velocity_cell_ids = []
      for kv in xrange(4):
         for jv in xrange(4):
            for iv in xrange(4):
               velocity_cell_ids.append(kv*16 + jv*4 + iv)

      for i in xrange(array_size):
         velocity_block_id = data_block_ids[i]
         avgIndex = 0

         for j in velocity_cell_ids + 64*velocity_block_id:
            velocity_cells[(int)(j)] = avgs[avgIndex]
            avgIndex = avgIndex + 1
      return velocity_cells

   def __set_cell_offset_and_blocks(self):
      ''' Read blocks per cell and the offset in the velocity space arrays for every cell with blocks into a private dictionary
      '''
      if len(self.__fileindex_for_cellid_blocks) != 0:
         # There's stuff already saved into the dictionary, don't save it again
         return
      #these two arrays are in the same order: 
      if self.__uses_new_vlsv_format == False:
         #list of cells for which dist function is saved
         cells_with_blocks = self.read("SpatialGrid","CELLSWITHBLOCKS")
      else:
         #list of cells for which dist function is saved
         cells_with_blocks = self.read(mesh="SpatialGrid",tag="CELLSWITHBLOCKS")
      #number of blocks in each cell for which data is stored
      if self.__uses_new_vlsv_format == False:
         blocks_per_cell = self.read("SpatialGrid","BLOCKSPERCELL")
      else:
         blocks_per_cell = self.read(mesh="SpatialGrid",tag="BLOCKSPERCELL")

      # Navigate to the correct position:
      from copy import copy
      offset = 0
      self.__fileindex_for_cellid_blocks = {}
      for i in xrange(0, len(cells_with_blocks)):
         self.__fileindex_for_cellid_blocks[cells_with_blocks[i]] = [copy(offset), copy(blocks_per_cell[i])]
         offset += blocks_per_cell[i]

   def list(self):
      ''' Print out a description of the content of the file. Useful
         for interactive usage
      '''
      if self.__uses_new_vlsv_format == False:
         return self.__list_old_format()
      print "tag = PARAMETER"
      for child in self.__xml_root:
         if child.tag == "PARAMETER":
            print "   ", child.attrib["name"]
      print "tag = VARIABLE"
      for child in self.__xml_root:
         if child.tag == "VARIABLE":
            print "   ", child.attrib["name"]
      print "tag = MESH"
      for child in self.__xml_root:
         if child.tag == "MESH":
            print "   ", child.attrib["name"]
      print "Datareducers:"
      for name in datareducers:
         print "   ",name, " based on ", datareducers[name].variables
      print "Other:"
      for child in self.__xml_root:
         if child.tag != "PARAMETER" and child.tag != "VARIABLE" and child.tag != "MESH":
            print "    tag = ", child.tag, " mesh = ", child.attrib["mesh"]



   def get_cellid_locations(self):
      ''' Returns a dictionary with cell id as the key and the index of the cell id as the value. The index is used to locate the cell id's values in the arrays that this reader returns
      '''
      if len( self.__fileindex_for_cellid ) == 0:
         self.__read_fileindex_for_cellid()
      return self.__fileindex_for_cellid

   def read(self, name="", tag="", mesh="", operator="pass", read_single_cellid=-1):
      ''' Read data from the open vlsv file. 
      
      Arguments:
      :param name Name of the data array
      :param tag  Tag of the data array.
      :param mesh Mesh for the data array
      :param operator Datareduction operator. "pass" does no operation on data.
      :param read_single_cellid  If -1 then all data is read. If nonzero then only the vector for the specified cell id is read
      :returns numpy array with the data

      '''
      # Check first if the name is in datareducers
      if name in datareducers:
         reducer = datareducers[name]
         # Read the necessary variables:
         tmp_vars = []
         for i in np.atleast_1d(reducer.variables):
            tmp_vars.append( self.read( i, tag, mesh, "pass", read_single_cellid ) )
         # Return the output of the datareducer
         return data_operators[operator](reducer.operation( tmp_vars ))


      if (len( self.__fileindex_for_cellid ) == 0):
         if read_single_cellid >= 0:
            self.__read_fileindex_for_cellid()
      if tag == "" and name == "" and tag == "":
         print "Bad arguments at read"

      if self.__fptr.closed:
         fptr = open(self.__file_name,"rb")
      else:
         fptr = self.__fptr

      #TODO, read_single_cellid should perhaps be an list/numpy array with cellids that are read in. This could be more efficient to 
      #     study multiple cells, e.g., along a line
      for child in self.__xml_root:
         if tag != "":
            if child.tag != tag:
               continue
         if name != "":
            if child.attrib["name"] != name:
               continue
         if mesh != "":
            if child.attrib["mesh"] != mesh:
               continue
         if child.tag == tag and child.attrib["name"] == name:
            vector_size = ast.literal_eval(child.attrib["vectorsize"])
            array_size = ast.literal_eval(child.attrib["arraysize"])
            element_size = ast.literal_eval(child.attrib["datasize"])
            datatype = child.attrib["datatype"]            
            offset = ast.literal_eval(child.text)
            if read_single_cellid >= 0:
               offset=offset+self.__fileindex_for_cellid[read_single_cellid]*element_size*vector_size
               array_size=1

            fptr.seek(offset)

            if datatype == "float" and element_size == 4:
               data = np.fromfile(fptr, dtype = np.float32, count=vector_size*array_size)
            if datatype == "float" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.float64, count=vector_size*array_size)
            if datatype == "int" and element_size == 4:
               data = np.fromfile(fptr, dtype=np.int32, count=vector_size*array_size)
            if datatype == "int" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.int64, count=vector_size*array_size)
            if datatype == "uint" and element_size == 4:
               data = np.fromfile(fptr, dtype=np.uint32, count=vector_size*array_size)
            if datatype == "uint" and element_size == 8:
               data = np.fromfile(fptr, dtype=np.uint64, count=vector_size*array_size)

            if self.__fptr.closed:
               fptr.close()

            if vector_size > 1:
               data=data.reshape(array_size, vector_size)
            
            if array_size == 1:
               return data_operators[operator](data[0])
            else:
               return data_operators[operator](data)

      if self.__fptr.closed:
         fptr.close()

   def read_variable(self, name, operator="pass",cellids=-1):
      ''' Read variables from the open vlsv file. 
      Arguments:
      :param name Name of the variable
      :param operator Datareduction operator. "pass" does no operation on data      
      :param cellids, a value of -1 reads all data
      :returns numpy array with the data

      '''
      if isinstance(cellids, (int, float)):
         return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", operator=operator, read_single_cellid=cellids)
      else:
         # NOTE: Should the file read be optimized by opening the file here until all cellids have been read? It can be optimized by the user manually, as well
         variable = []
         for i in cellids:
            variable.append( self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", operator=operator, read_single_cellid=i) )
         return np.array(variable, copy=False)

   def read_variable_info(self, name, operator="pass",cellids=-1):
      ''' Read variables from the open vlsv file and input the data into VariableInfo
      Arguments:
      :param name Name of the variable
      :param operator Datareduction operator. "pass" does no operation on data
      :param cellids, a value of -1 reads all data
      :returns numpy array with the data
      '''
      data = self.read_variable(name=name, operator=operator, cellids=cellids)
      from variable import VariableInfo
      if name in datareducers:
         units = datareducers[name].units
      else:
         units = ""
      if operator != "pass":
         return VariableInfo(data_array=data, name=name + "_" + operator, units=units)
      else:
         return VariableInfo(data_array=data, name=name, units=units)


   def get_cellid(self, coordinates):
      ''' Returns the cell id at given coordinates

      Arguments:
      :param coordinates        The cell's coordinates
      :returns the cell id
      NOTE: Returns 0 if the cellid is out of bounds!
      '''
      # Check that the coordinates are not out of bounds:
      if (self.__xmax < coordinates[0]) or (self.__xmin > coordinates[0]):
         return 0
      if (self.__ymax < coordinates[1]) or (self.__ymin > coordinates[1]):
         return 0
      if (self.__zmax < coordinates[2]) or (self.__zmin > coordinates[2]):
         return 0
      # Get cell lengths:
      cell_lengths = np.array([self.__dx, self.__dy, self.__dz])
   
      # Get cell indices:
      cellindices = np.array([(int)((coordinates[0] - self.__xmin)/(float)(cell_lengths[0])), (int)((coordinates[1] - self.__ymin)/(float)(cell_lengths[1])), (int)((coordinates[2] - self.__zmin)/(float)(cell_lengths[2]))])
      # Get the cell id:
      cellid = cellindices[0] + cellindices[1] * self.__xcells + cellindices[2] * self.__xcells * self.__ycells + 1
      return cellid

   def get_cell_coordinates(self, cellid):
      ''' Returns a given cell's coordinates as a numpy array

      Arguments:
      :param cellid            The cell's ID
      :returns a numpy array with the coordinates
      '''
      # Get cell lengths:
      cell_lengths = np.array([(self.__xmax - self.__xmin)/(float)(self.__xcells), (self.__ymax - self.__ymin)/(float)(self.__ycells), (self.__zmax - self.__zmin)/(float)(self.__zcells)])
      # Get cell indices:
      cellid = (int)(cellid - 1)
      cellindices = np.zeros(3)
      cellindices[0] = (int)(cellid)%(int)(self.__xcells)
      cellindices[1] = ((int)(cellid)/(int)(self.__xcells))%(int)(self.__ycells)
      cellindices[2] = (int)(cellid)/(int)(self.__xcells*self.__ycells)
   
      # Get cell coordinates:
      cellcoordinates = np.zeros(3)
      cellcoordinates[0] = self.__xmin + (cellindices[0] + 0.5) * cell_lengths[0]
      cellcoordinates[1] = self.__ymin + (cellindices[1] + 0.5) * cell_lengths[1]
      cellcoordinates[2] = self.__zmin + (cellindices[2] + 0.5) * cell_lengths[2]
      # Return the coordinates:
      return np.array(cellcoordinates)

   def get_velocity_cell_coordinates(self, vcellids):
      ''' Returns a given velocity cell's coordinates as a numpy array

      Arguments:
      :param vcellid       The velocity cell's ID
      :return a numpy array with the coordinates
      '''
      vcellids = np.atleast_1d(vcellids)
      # Get block ids:
      blocks = vcellids.astype(int) / 64
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks))
      blockIndicesZ = blocks.astype(int)/(int)(self.__vxblocks*self.__vyblocks)
      blockCoordinatesX = blockIndicesX.astype(float) * self.__dvx * 4 + self.__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__dvy * 4 + self.__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__dvz * 4 + self.__vzmin
      # Get cell indices:
      cellids = np.remainder(vcellids.astype(int), (int)(64))
      cellIndicesX = np.remainder(cellids.astype(int), (int)(4))
      cellIndicesY = np.remainder((cellids.astype(int)/(int)(4)).astype(int), (int)(4))
      cellIndicesZ = cellids.astype(int)/(int)(16)
      # Get cell coordinates:
      cellCoordinates = np.array([blockCoordinatesX.astype(float) + (cellIndicesX.astype(float) + 0.5) * self.__dvx,
                                  blockCoordinatesY.astype(float) + (cellIndicesY.astype(float) + 0.5) * self.__dvy,
                                  blockCoordinatesZ.astype(float) + (cellIndicesZ.astype(float) + 0.5) * self.__dvz])

      return cellCoordinates.transpose()

   def get_velocity_block_coordinates( self, blocks ):
      ''' Returns the block coordinates of the given blocks in a numpy array

          :param blocks         list of block ids
          :returns a numpy array containing the block coordinates e.g. np.array([np.array([2,1,3]), np.array([5,6,6]), ..])
      '''
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks))
      blockIndicesZ = blocks.astype(int)/(int)(self.__vxblocks*self.__vyblocks)
      blockCoordinatesX = blockIndicesX.astype(float) * self.__dvx * 4 + self.__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__dvy * 4 + self.__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__dvz * 4 + self.__vzmin
      # Return the coordinates:
      return np.array([blockCoordinatesX.astype(float),
                       blockCoordinatesY.astype(float),
                       blockCoordinatesZ.astype(float)]).transpose()

   def get_velocity_blocks( self, blockcoordinates ):
      ''' Returns the block ids of the given block coordinates in a numpy array form

          :param blockCoordinates         list of block coordinates e.g. np.array([np.array([2,1,3]), np.array([5,6,6]), ..])
          :returns a numpy array containing the block ids e.g. np.array([4,2,56,44,2, ..])
      '''
      mins = np.array([self.__vxmin, self.__vymin, self.__vzmin]).astype(float)
      dvs = np.array([4*self.__dvx, 4*self.__dvy, 4*self.__dvz]).astype(float)
      multiplier = np.array([1, self.__vxblocks, self.__vxblocks * self.__vyblocks]).astype(float)
      velocity_block_ids = np.sum(np.floor(((blockCoordinates.astype(float) - mins) / dvs)) * multiplier, axis=-1)
      return velocity_block_ids

   def construct_velocity_cells( self, blocks ):
      ''' Returns velocity cells in given blocks

          :param blocks         list of block ids
          :returns a numpy array containing the velocity cell ids e.g. np.array([4,2,56,44,522, ..])
      '''
      return np.ravel(np.outer(np.array(blocks), np.ones(64)) + np.arange(64))

   def construct_velocity_cell_coordinates( self, blocks ):
      ''' Returns velocity cell coordinates in given blocks

          :param blocks         list of block ids
          :returns a numpy array containing the velocity cell ids e.g. np.array([4,2,56,44,522, ..])
      '''
      # Construct velocity cell coordinates from velocity cells and return them
      return self.get_velocity_cell_coordinates( self.construct_velocity_cells(blocks) )


   def construct_velocity_cell_nodes( self, blocks ):
      ''' Returns velocity cell nodes in given blocks

          :param blocks         list of block ids
          :param avgs           list of avgs values
          :returns a numpy array containing velocity cell nodes and the keys for velocity cells
          NOTE: THIS IS USED FOR CONSTRUCTING VELOCITY SPACES WITH VISUALIZATION PROGRAMS LIKE MAYAVI
      '''
      blocks = np.array(blocks)
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks)).astype(np.uint16)
      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks)).astype(np.uint16)
      blockIndicesZ = (blocks.astype(np.uint64)/(int)(self.__vxblocks*self.__vyblocks)).astype(np.uint16)

      cellsPerDirection = 4
      cellsPerBlock = 64

      # Get velocity cell min coordinates (per velocity block)
      vcellids = np.arange(cellsPerBlock).astype(np.uint32)
      cellIndicesX = np.remainder(vcellids.astype(int), (int)(cellsPerDirection)).astype(np.uint16)
      cellIndicesY = np.remainder((vcellids.astype(int)/(int)(cellsPerDirection)).astype(int), (int)(cellsPerDirection)).astype(np.uint16)
      cellIndicesZ = (vcellids.astype(int)/(int)(cellsPerDirection*cellsPerDirection)).astype(np.uint16)

      # Construct velocity cell node indices for every velocity cell per velocity block

      nodesPerCell = 8

      # NOTE: The ordering of the numpy array won't make sense to anyone who hasn't read VTK documentation. For further info check VTK_VOXEL. The numpy array is constructed according to VTK voxel's nodes
      cellNodeIndicesX = np.ravel(np.outer(cellIndicesX, np.ones(nodesPerCell)) + np.array([0, 1, 0, 1, 0, 1, 0, 1])).astype(np.uint16)
      cellNodeIndicesY = np.ravel(np.outer(cellIndicesY, np.ones(nodesPerCell)) + np.array([0, 0, 1, 1, 0, 0, 1, 1])).astype(np.uint16)
      cellNodeIndicesZ = np.ravel(np.outer(cellIndicesZ, np.ones(nodesPerCell)) + np.array([0, 0, 0, 0, 1, 1, 1, 1])).astype(np.uint16)

      nodeIndices_local = []
      nodesPerDirection = 5

      for i in xrange(nodesPerDirection):
         for j in xrange(nodesPerDirection):
            for k in xrange(nodesPerDirection):
               nodeIndices_local.append(np.array([i,j,k]))
      nodeIndices_local = np.array(nodeIndices_local).astype(np.uint16)

      nodesPerBlock = (int)(nodesPerDirection * nodesPerDirection * nodesPerDirection)


      def calculate_node_indices( self, blockIndicesX, blockIndicesY, blockIndicesZ, nodeIndices_local, nodesPerBlock, cellsPerDirection ):
         nodeIndicesX = np.ravel(np.outer(blockIndicesX, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,0])
         nodeIndicesY = np.ravel(np.outer(blockIndicesY, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,1])
         nodeIndicesZ = np.ravel(np.outer(blockIndicesZ, np.ones(nodesPerBlock).astype(np.uint16)) * cellsPerDirection + nodeIndices_local[:,2])
   
         nodeIndices = np.transpose(np.array([nodeIndicesX, nodeIndicesY, nodeIndicesZ], copy=False))

         # Transform indices into unique keys
         nodeKeys = np.sum(nodeIndices * np.array([1, cellsPerDirection*self.__vxblocks+1, (cellsPerDirection*self.__vxblocks+1)*(cellsPerDirection*self.__vyblocks+1)]), axis=1)
         # Sort the keys and delete duplicates
         return np.unique(nodeKeys)
      #nodeIndices = calculate_node_indices( blockIndicesX, blockIndicesY, blockIndicesZ, nodeIndices_local, nodesPerBlock, cellsPerDirection )

      # Put the node indices into keys:
      nodeKeys = np.array([], dtype=np.uint64)
      N = 10
      for i in xrange(N):
         fromIndex = i*(len(blockIndicesX)/N)
         if i != N-1:
            toIndex = (i+1)*(len(blockIndicesX)/N)
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
         globalCellIndices = np.sum(globalCellIndices * np.array([1, cellsPerDirection*self.__vxblocks+1, (cellsPerDirection*self.__vxblocks+1)*(cellsPerDirection*self.__vyblocks+1)]), axis=1)
         # Return cell nodes' indexes in the nodeKeys list
         return np.searchsorted(nodeKeys, globalCellIndices)


      # Create cellKeys
      cellKeys = np.zeros(len(blockIndicesX)*cellsPerBlock*nodesPerCell, dtype=np.uint32)
      N = 10
      # Append keys in cuts to save memory
      for i in xrange(N):
         fromIndex = i*(len(blockIndicesX)/N)
         if i != N-1:
            toIndex = (i+1)*(len(blockIndicesX)/N)
         else:
            toIndex = len(blockIndicesX)
         # Append cell keys
         cellKeys[fromIndex*cellsPerBlock*nodesPerCell:toIndex*cellsPerBlock*nodesPerCell] = calc_global_cell_keys( self, blockIndicesX[fromIndex:toIndex], blockIndicesY[fromIndex:toIndex], blockIndicesZ[fromIndex:toIndex], cellNodeIndicesX, cellNodeIndicesY, cellNodeIndicesZ, cellsPerBlock, nodesPerCell, cellsPerDirection, nodeKeys )

      cellKeys = np.reshape(cellKeys, (len(blocks)*64,8))

      # We now have all the cell keys and avgs values! (avgs is in the same order as cell keys)
      # Now transform node indices back into real indices
      nodeCoordinatesX = np.remainder(nodeKeys, (int)(cellsPerDirection*self.__vxblocks+1)).astype(np.float32) * self.__dvx + self.__vxmin
      nodeCoordinatesY = np.remainder(nodeKeys/(int)(cellsPerDirection*self.__vxblocks+1), cellsPerDirection*self.__vyblocks+1).astype(np.float32) * self.__dvy + self.__vymin
      nodeCoordinatesZ = ( nodeKeys / (int)((cellsPerDirection*self.__vxblocks+1) * (cellsPerDirection*self.__vyblocks+1)) ).astype(np.float32) * self.__dvz + self.__vzmin

      # Nodekeyss is no longer needed
      del nodeKeys

      nodes = np.array([nodeCoordinatesX, nodeCoordinatesY, nodeCoordinatesZ], copy=False)
      # Take a transpose
      nodes = np.transpose(nodes)

      return [nodes, cellKeys]





   def read_parameter(self, name):
      ''' Read a parameter from the vlsv file

      :param name   Name of the parameter
      '''
      if self.__uses_new_vlsv_format == True:
         return self.read(name=name, tag="PARAMETER")
      else:
         # Uses old format
         return self.__read_parameter_old_format(name=name)


   def read_velocity_cells(self, cellid):
      ''' Read velocity cells from a spatial cell
      
      Arguments:
      :param cellid Cell ID of the cell whose velocity cells are read
      :returns numpy array with blocks in the cell. Empty if cell has no stored blocks.
      '''
      #these two arrays are in the same order: 
      #list of cells for which dist function is saved
      if self.__uses_new_vlsv_format == False:
         cells_with_blocks = self.read("SpatialGrid","CELLSWITHBLOCKS")
      else:
         cells_with_blocks = self.read(mesh="SpatialGrid",tag="CELLSWITHBLOCKS")
      #number of blocks in each cell for which data is stored
      if self.__uses_new_vlsv_format == False:
         blocks_per_cell = self.read("SpatialGrid","BLOCKSPERCELL")
      else:
         blocks_per_cell = self.read(mesh="SpatialGrid",tag="BLOCKSPERCELL")
      (cells_with_blocks_index,) = np.where(cells_with_blocks == cellid)

      if len(cells_with_blocks_index) == 0:
         #block data did not exist
         print "Cell does not have velocity distribution"
         return []

      num_of_blocks = np.atleast_1d(blocks_per_cell)[cells_with_blocks_index[0]]

      # Check for the old library
      if self.__uses_new_vlsv_format == False:
         # Uses old format
         return self.__read_velocity_cells_old_format(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)
      else:
         # Uses new format:
         return self.__read_velocity_cells_new_format(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)


   def read_blocks(self, cellid):
      ''' Read raw block data from the open file and return the data along with block ids
      
      Arguments:
      :param cell_id Cell ID of the cell whose velocity blocks are read
      :returns numpy array with block ids and data eg [array([2, 5, 6, 234, 21]), array([1.0e-8, 2.1e-8, 2.1e-8, 0, 4.0e-8])]
      '''
      if( len(self.__fileindex_for_cellid_blocks) == 0 ):
         # Set the locations
         self.__set_cell_offset_and_blocks()

      if self.__uses_new_vlsv_format == False:
         # Uses old format
         return self.__read_blocks_old_format(cellid)
      else:
         # Uses new format
         return self.__read_blocks_new_format(cellid)

      return []

   def optimize_open_file(self):
      '''Opens the vlsv file for reading
         Files are opened and closed automatically upon reading and in the case of reading multiple times it will help to keep the file open with this command
         Example usage:
         variables = []
         vlsvReader.optimize_open_file()
         for i in xrange(1000):
            variables.append(vlsvReader.read_variable("rho", cellids=i))
         vlsvReader.optimize_close_file()
         NOTE: This should only be used for optimization purposes.
      '''
      self.__fptr = open(self.__file_name,"rb")


   def optimize_close_file(self):
      '''Closes the vlsv file
         Files are opened and closed automatically upon reading and in the case of reading multiple times it will help to keep the file open with this command
         Example usage:
         variables = []
         vlsvReader.optimize_open_file()
         for i in xrange(1000):
            variables.append(vlsvReader.read_variable("rho", cellids=i))
         vlsvReader.optimize_close_file()
         NOTE: This should only be used for optimization purposes.
      '''
      if self.__fptr.closed:
         return
      else:
         self.__fptr.close()
         return

   def optimize_clear_fileindex_for_cellid_blocks(self):
      ''' Clears a private variable containing number of blocks and offsets for particular cell ids
          NOTE: This should only be used for optimization purposes.
          Example usage:
          vlsvReaders = []
          # Open a list of vlsv files
          for i in xrange(1000):
             vlsvReaders.append( VlsvFile("test" + str(i) + ".vlsv") )
          # Go through vlsv readers and print info:
          for vlsvReader in vlsvReaders:
             # Print something from the file on the screen
             print vlsvReader.read_blocks( cellid= 5021 ) # Stores info into a private variable
             # Upon reading from vlsvReader a private variable that contains info on cells that have blocks has been saved -- now clear it to save memory
             vlsvReader.optimize_clear_fileindex_for_cellid_blocks()
      '''
      self.__fileindex_for_cellid_blocks = {}

   def optimize_clear_fileindex_for_cellid(self):
      ''' Clears a private variable containing cell ids and their locations
          NOTE: This should only be used for optimization purposes.
          Example usage:
          vlsvReaders = []
          # Open a list of vlsv files
          for i in xrange(1000):
             vlsvReaders.append( VlsvFile("test" + str(i) + ".vlsv") )
          # Go through vlsv readers and print info:
          for vlsvReader in vlsvReaders:
             # Print something from the file on the screen
             print vlsvReader.read_variable("B", cellids=2) # Stores info into a private variable
             # Upon reading from vlsvReader a private variable that contains info on cells that have blocks has been saved -- now clear it to save memory
             vlsvReader.optimize_clear_fileindex_for_cellid()
      '''
      self.__fileindex_for_cellid = {}



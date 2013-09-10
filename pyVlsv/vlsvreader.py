import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
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
         
   def __list_old(self):
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
      print "Other:"
      for child in self.__xml_root:
         if child.tag != "PARAMETERS" and child.tag != "VARIABLE":
            print "    tag = ", child.tag, " name = ", child.attrib["name"]
      return

   def __read_parameter_old(self, name):
      return self.read(name=name, tag="PARAMETERS")

   def __read_blocks_old(self, cellid):
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

   def __read_blocks_new(self, cellid):
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



   def __read_velocity_cells_old( self, cellid, cells_with_blocks, blocks_per_cell, cells_with_blocks_index  ):
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

   def __read_velocity_cells_new( self, cellid, cells_with_blocks, blocks_per_cell, cells_with_blocks_index  ):
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
         return self.__list_old()
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

   def read(self, name="", tag="", mesh="", read_single_cellid=-1):
      ''' Read data from the open vlsv file. 
      
      Arguments:
      :param name Name of the data array
      :param tag  Tag of the data array.
      :param read_single_cellid  If -1 then all data is read. If nonzero then only the vector for the specified cell id is read
      :returns numpy array with the data

      '''
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
               return data[0]
            else:
               return data

      if self.__fptr.closed:
         fptr.close()


   def read_variables(self, name):
      ''' Read variables from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :returns numpy array with the data

      '''
      return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", read_single_cellid=-1)

   def read_variables_for_cellids(self, name, cellids, index=3):
      ''' Read variables from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :param cellids List of cellids
      :returns numpy array with the data
      Note: Format of the numpy array:
      [variableForCellid1, variableForCellid2, variableForCellid3, ..]
      NOTE: THIS IS MAINLY USED FOR OPTIMIZATION PURPOSES
      '''
      if len( self.__fileindex_for_cellid ) == 0:
         self.__read_fileindex_for_cellid()
      # Read the variable:
      variablelist = self.read_variables(name)
      #make a list, if variablelist is a scalar
      if(not isinstance(variablelist, Iterable)):
         variablelist=[ variablelist ]
      #Pick the variables with the cell ids in the list:
      returnvariablelist = []
      for cellid in cellids:
         if index == 3:
            returnvariablelist.append(np.array(variablelist[self.__fileindex_for_cellid[cellid]]))
         else:
            returnvariablelist.append(np.array(variablelist[self.__fileindex_for_cellid[cellid]][index]))
      # Return the variables:
      return np.array(returnvariablelist)

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
      # Get xmax, xmin and xcells_ini
      xmax = self.read_parameter(name="xmax")
      xmin = self.read_parameter(name="xmin")
      xcells = (int)(self.read_parameter(name="xcells_ini"))
      # Do the same for y
      ymax = self.read_parameter(name="ymax")
      ymin = self.read_parameter(name="ymin")
      ycells = (int)(self.read_parameter(name="ycells_ini"))
      # And for z
      zmax = self.read_parameter(name="zmax")
      zmin = self.read_parameter(name="zmin")
      zcells = (int)(self.read_parameter(name="zcells_ini"))
      # Get cell lengths:
      cell_lengths = np.array([(xmax - xmin)/(float)(xcells), (ymax - ymin)/(float)(ycells), (zmax - zmin)/(float)(zcells)])
      # Get cell indices:
      cellid = (int)(cellid - 1)
      cellindices = np.zeros(3)
      cellindices[0] = (int)(cellid)%(int)(xcells)
      cellindices[1] = ((int)(cellid)/(int)(xcells))%(int)(ycells)
      cellindices[2] = (int)(cellid)/(int)(xcells*ycells)
   
      # Get cell coordinates:
      cellcoordinates = np.zeros(3)
      cellcoordinates[0] = xmin + (cellindices[0] + 0.5) * cell_lengths[0]
      cellcoordinates[1] = ymin + (cellindices[1] + 0.5) * cell_lengths[1]
      cellcoordinates[2] = zmin + (cellindices[2] + 0.5) * cell_lengths[2]
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

#      # Get block id:
#      block = (int)(vcellid) / 64
#      # Get block coordinates:
#      blockIndices = np.array( [(int)(block)%(int)(self.__vxblocks), (int)((int)(block)/(int)(self.__vxblocks))%(int)(self.__vyblocks), (int)(block)/(int)(self.__vxblocks*self.__vyblocks)] )
#      blockCoordinates = np.array([self.__vxmin + 4 * self.__dvx * blockIndices[0],
#                                   self.__vymin + 4 * self.__dvy * blockIndices[1],
#                                   self.__vzmin + 4 * self.__dvz * blockIndices[2]])
#      #vcellid = 64 * velocity_block_id + kv*4*4 + jv*4 + iv
#      # Get cell coordinates:
#      cellIndices = np.array([(int)((int)(vcellid)%64)%4, (int)(((int)((int)(vcellid)%64)/4))%4, (int)((int)(vcellid)%64)/(int)(4*4)])
#      cellCoordinates = np.array([(cellIndices[0] + 0.5) * self.__dvx, (cellIndices[1] + 0.5) * self.__dvy, (cellIndices[2] + 0.5) * self.__dvz])
#      # Get the coordinates:
#      vcellCoordinates = blockCoordinates + cellCoordinates
      # Return cell coordinates:
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

#   def construct_velocity_cell_nodes( blocks ):
#      ''' Returns velocity cell nodes in given blocks
#
#          :param blocks         list of block ids
#          :returns a numpy array containing velocity cell nodes
#          NOTE: THIS IS USED FOR CONSTRUCTING VELOCITY SPACES WITH VISUALIZATION PROGRAMS LIKE MAYAVI
#      '''
#      # Get block coordinates:
#      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks))
#      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks))
#      blockIndicesZ = blocks.astype(int)/(int)(self.__vxblocks*self.__vyblocks)
#      blockCoordinatesX = blockIndicesX.astype(float) * self.__dvx * 4 + self.__vxmin
#      blockCoordinatesY = blockIndicesY.astype(float) * self.__dvy * 4 + self.__vymin
#      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__dvz * 4 + self.__vzmin
#      # Get velocity cell min coordinates (per velocity block)
#      vcellids = np.arange(64)
#      minCellIndicesX = np.remainder(vcellids.astype(int), (int)(4))
#      minCellIndicesY = np.remainder((vcellids.astype(int)/(int)(4)).astype(int), (int)(4))
#      minCellIndicesZ = vcellids.astype(int)/(int)(16)
#      minCellCoordinatesX = minCellIndicesX.astype(float) * (float)(self.__dvx)
#      minCellCoordinatesY = minCellIndicesY.astype(float) * (float)(self.__dvy)
#      minCellCoordinatesZ = minCellIndicesZ.astype(float) * (float)(self.__dvz)
#      # Construct velocity cell nodes (per velocity block)
#      # Note: Every velocity cell has 8 nodes
#      cellNodesX = np.outer(minCellCoordinatesX, np.ones(8))
#      cellNodesY = np.outer(minCellCoordinatesY, np.ones(8))
#      cellNodesZ = np.outer(minCellCoordinatesZ, np.ones(8))
#      # Get the coordinates of the codes:
#      # Note: The arrays are a bit confusing but if you look up VTK_VOXEL in some vtk doc it will be clearer why the array looks like that
#      cellNodesX = cellNodesX + np.array([0, 1, 0, 1, 0, 1, 0, 1]) * self.__dvx
#      cellNodesY = cellNodesY + np.array([0, 0, 1, 1, 0, 0, 1, 1]) * self.__dvy
#      cellNodesZ = cellNodesZ + np.array([0, 0, 0, 0, 1, 1, 1, 1]) * self.__dvz
#      # The cellNodesX, y, z should now be a 64 * 8 matrix: Check to make sure:
#      if (len(cellNodesX) != 64) or (len(cellNodesX[0]) != 8):
#         print "BAD LENGTH CELL NODES X"
#         return []
#      if (len(cellNodesY) != 64) or (len(cellNodesY[0]) != 8):
#         print "BAD LENGTH CELL NODES Y"
#         return []
#      if (len(cellNodesZ) != 64) or (len(cellNodesZ[0]) != 8):
#         print "BAD LENGTH CELL NODES Z"
#         return []
#      # Ravel the cell nodes
#      cellNodesX = np.ravel(cellNodesX)
#      cellNodesY = np.ravel(cellNodesY)
#      cellNodesZ = np.ravel(cellNodesZ)
#      # Get all nodes:
#      # Transform blockCoordinatesX into len(blockCoordinatesX) * 64 * 8 matrix (for calculation) and then ravel:
#      nodeCoordinatesX = np.ravel(np.outer(blockCoordinatesX, np.ones(64*8)) + cellNodesX)
#      nodeCoordinatesY = np.ravel(np.outer(blockCoordinatesY, np.ones(64*8)) + cellNodesY)
#      nodeCoordinatesZ = np.ravel(np.outer(blockCoordinatesZ, np.ones(64*8)) + cellNodesZ)
#      # Return the cell coordinates:
#      return np.array([nodeCoordinatesX, nodeCoordinatesY, nodeCoordinatesZ]).transpose()

   def construct_velocity_cell_nodes( self, blocks ):
      ''' Returns velocity cell nodes in given blocks

          :param blocks         list of block ids
          :param avgs           list of avgs values
          :returns a numpy array containing velocity cell nodes and the keys for velocity cells
          NOTE: THIS IS USED FOR CONSTRUCTING VELOCITY SPACES WITH VISUALIZATION PROGRAMS LIKE MAYAVI
      '''
      blocks = np.array(blocks)
      # Get block coordinates:
      blockIndicesX = np.remainder(blocks.astype(int), (int)(self.__vxblocks))
      blockIndicesY = np.remainder(blocks.astype(int)/(int)(self.__vxblocks), (int)(self.__vyblocks))
      blockIndicesZ = blocks.astype(int)/(int)(self.__vxblocks*self.__vyblocks)

      cellsPerDirection = 4
      cellsPerBlock = 64

      blockCoordinatesX = blockIndicesX.astype(float) * self.__dvx * cellsPerDirection + self.__vxmin
      blockCoordinatesY = blockIndicesY.astype(float) * self.__dvy * cellsPerDirection + self.__vymin
      blockCoordinatesZ = blockIndicesZ.astype(float) * self.__dvz * cellsPerDirection + self.__vzmin
      # Get velocity cell min coordinates (per velocity block)
      vcellids = np.arange(cellsPerBlock)
      cellIndicesX = np.remainder(vcellids.astype(int), (int)(cellsPerDirection))
      cellIndicesY = np.remainder((vcellids.astype(int)/(int)(cellsPerDirection)).astype(int), (int)(cellsPerDirection))
      cellIndicesZ = vcellids.astype(int)/(int)(cellsPerDirection*cellsPerDirection)
      minCellCoordinatesX = cellIndicesX.astype(float) * (float)(self.__dvx)
      minCellCoordinatesY = cellIndicesY.astype(float) * (float)(self.__dvy)
      minCellCoordinatesZ = cellIndicesZ.astype(float) * (float)(self.__dvz)

      # Construct velocity cell node indices for every velocity cell per velocity block

      nodesPerCell = 8

      cellNodeIndicesX = np.ravel(np.outer(cellIndicesX, np.ones(nodesPerCell)) + np.array([0, 1, 0, 1, 0, 1, 0, 1]).astype(int))
      cellNodeIndicesY = np.ravel(np.outer(cellIndicesY, np.ones(nodesPerCell)) + np.array([0, 0, 1, 1, 0, 0, 1, 1]).astype(int))
      cellNodeIndicesZ = np.ravel(np.outer(cellIndicesZ, np.ones(nodesPerCell)) + np.array([0, 0, 0, 0, 1, 1, 1, 1]).astype(int))

      nodeIndices_local = []
      nodesPerDirection = 5

      for i in xrange(nodesPerDirection):
         for j in xrange(nodesPerDirection):
            for k in xrange(nodesPerDirection):
               nodeIndices_local.append(np.array([i,j,k]))
      nodeIndices_local = np.array(nodeIndices_local)

      nodesPerBlock = (int)(nodesPerDirection * nodesPerDirection * nodesPerDirection)

      nodeIndicesX = np.ravel(np.outer(blockIndicesX, np.ones(nodesPerBlock)) * cellsPerDirection + nodeIndices_local[:,0])
      nodeIndicesY = np.ravel(np.outer(blockIndicesY, np.ones(nodesPerBlock)) * cellsPerDirection + nodeIndices_local[:,1])
      nodeIndicesZ = np.ravel(np.outer(blockIndicesZ, np.ones(nodesPerBlock)) * cellsPerDirection + nodeIndices_local[:,2])


      nodeIndices = np.array([nodeIndicesX, nodeIndicesY, nodeIndicesZ]).transpose().astype(int)

      # Put the node indices into keys:
      
      nodeIndices = np.sum(nodeIndices * np.array([1, cellsPerDirection*self.__vxblocks+1, (cellsPerDirection*self.__vxblocks+1)*(cellsPerDirection*self.__vyblocks+1)]), axis=1)

      # Delete duplicate nodes and sort the list:
      nodeIndices = np.unique(nodeIndices) #We now have all of the nodes in a list!

      # Next create node  indices for the cells
      globalCellIndicesX = np.ravel(np.outer(blockIndicesX, np.ones(cellsPerBlock * nodesPerCell)) * cellsPerDirection + cellNodeIndicesX)
      globalCellIndicesY = np.ravel(np.outer(blockIndicesY, np.ones(cellsPerBlock * nodesPerCell)) * cellsPerDirection + cellNodeIndicesY)
      globalCellIndicesZ = np.ravel(np.outer(blockIndicesZ, np.ones(cellsPerBlock * nodesPerCell)) * cellsPerDirection + cellNodeIndicesZ)

      globalCellIndices = np.array([globalCellIndicesX, globalCellIndicesY, globalCellIndicesZ]).transpose().astype(int)

      # Put the indices into keys:
      cellKeys = globalCellIndices
      # Transform into keys
      cellKeys = np.sum(cellKeys * np.array([1, cellsPerDirection*self.__vxblocks+1, (cellsPerDirection*self.__vxblocks+1)*(cellsPerDirection*self.__vyblocks+1)]), axis=1).astype(int)
      # Get the proper indexes with the keys
      cellKeys = np.array(np.split(np.searchsorted(nodeIndices, cellKeys), 64*len(blocks))).astype(int)

      # We now have all the cell keys and avgs values! (avgs is in the same order as cell keys)
      # Now transform node indices back into real indices
      nodeCoordinatesX = np.remainder(nodeIndices, cellsPerDirection*self.__vxblocks+1) * self.__dvx + self.__vxmin
      nodeCoordinatesY = np.remainder(nodeIndices.astype(int)/(int)(cellsPerDirection*self.__vxblocks+1), cellsPerDirection*self.__vyblocks+1) * self.__dvy + self.__vymin
      nodeCoordinatesZ = ( nodeIndices.astype(int) / (int)((cellsPerDirection*self.__vxblocks+1) * (cellsPerDirection*self.__vyblocks+1)) ).astype(float) * self.__dvz + self.__vzmin

      nodes = np.array([nodeCoordinatesX, nodeCoordinatesY, nodeCoordinatesZ]).transpose()

      return [nodes, cellKeys]




   def read_variable(self, name, cellid):
      ''' Read a variable of a given cell from the open vlsv file. 
      
      Arguments:
      :param name Name of the variable
      :param cellid Cell's cell id
      :returns numpy array with the data

      '''
      return self.read(mesh="SpatialGrid", name=name, tag="VARIABLE", read_single_cellid=cellid)

   def read_parameter(self, name):
      ''' Read a parameter from the vlsv file

      :param name   Name of the parameter
      '''
      if self.__uses_new_vlsv_format == True:
         return self.read(name=name, tag="PARAMETER")
      else:
         # Uses old format
         return self.__read_parameter_old(name=name)


#getVelocityBlockCoordinates( cellStruct, blockId, blockCoordinates )
#void getVelocityBlockCoordinates(const CellStructure & cellStruct, const uint64_t block, array<Real, 3> & coordinates ) {
#   //First get indices:
#   array<uint64_t, 3> blockIndices;
#   blockIndices[0] = block % cellStruct.vcell_bounds[0];
#   blockIndices[1] = (block / cellStruct.vcell_bounds[0]) % cellStruct.vcell_bounds[1];
#   blockIndices[2] = block / (cellStruct.vcell_bounds[0] * cellStruct.vcell_bounds[1]);
#   //Store the coordinates:
#   for( int i = 0; i < 3; ++i ) {
#      coordinates[i] = cellStruct.min_vcoordinates[i] + cellStruct.vblock_length[i] * blockIndices[i];
#   }
#   return;
#}
#if( vlsvReader.readParameter( "vxblocks_ini", vcell_bounds[0] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#//y-direction
#if( vlsvReader.readParameter( "vyblocks_ini", vcell_bounds[1] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#//z-direction
#if( vlsvReader.readParameter( "vzblocks_ini", vcell_bounds[2] ) == false ) cerr << "FAILED TO READ PARAMETER AT " << __FILE__ << " " << __LINE__ << endl;
#for child in self.__xml_root:
#   if tag != "":
#      if child.tag != tag:
#         continue
#   if name != "":
#      if child.attrib["name"] != name:
#         continue
#   if mesh != "":
#      if child.attrib["mesh"] != mesh:
#         continue
#   if child.tag == tag and child.attrib["name"] == name:
#      vector_size = ast.literal_eval(child.attrib["vectorsize"])
#      array_size = ast.literal_eval(child.attrib["arraysize"])
#      element_size = ast.literal_eval(child.attrib["datasize"])
#      datatype = child.attrib["datatype"]            
#      offset = ast.literal_eval(child.text)
#      if read_single_cellid >= 0:
#         offset=offset+self.__fileindex_for_cellid[read_single_cellid]*element_size*vector_size
#         array_size=1
#
#      fptr = open(self.__file_name, "rb")
#      fptr.seek(offset)
#
#      if datatype == "float" and element_size == 4:
#         data = np.fromfile(fptr, dtype = np.float32, count=vector_size*array_size)
#      if datatype == "float" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.float64, count=vector_size*array_size)
#      if datatype == "int" and element_size == 4:
#         data = np.fromfile(fptr, dtype=np.int32, count=vector_size*array_size)
#      if datatype == "int" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.int64, count=vector_size*array_size)
#      if datatype == "uint" and element_size == 4:
#         data = np.fromfile(fptr, dtype=np.uint32, count=vector_size*array_size)
#      if datatype == "uint" and element_size == 8:
#         data = np.fromfile(fptr, dtype=np.uint64, count=vector_size*array_size)
#      fptr.close() 
#
#      if vector_size > 1:
#         data=data.reshape(array_size, vector_size)
#      
#      if array_size == 1:
#         return data[0]
#      else:
#         return data

      
#const velocity_block_indices_t indices = {{
#   (unsigned int) np.floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
#   (unsigned int) np.floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
#   (unsigned int) np.floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
#}};
#indices[0] = cell % block_vx_length;
#indices[1] = (cell / block_vx_length) % block_vy_length;
#indices[2] = cell / (block_vx_length * block_vy_length);
#static unsigned int get_velocity_block(
#   const Real vx,
#   const Real vy,
#   const Real vz
#) {
#   if (vx < SpatialCell::vx_min || vx >= SpatialCell::vx_max
#       || vy < SpatialCell::vy_min || vy >= SpatialCell::vy_max
#       || vz < SpatialCell::vz_min || vz >= SpatialCell::vz_max) {
#      return error_velocity_block;
#   }
#
#   const velocity_block_indices_t indices = {{
#      (unsigned int) np.floor((vx - SpatialCell::vx_min) / SpatialCell::block_dvx),
#      (unsigned int) np.floor((vy - SpatialCell::vy_min) / SpatialCell::block_dvy),
#      (unsigned int) np.floor((vz - SpatialCell::vz_min) / SpatialCell::block_dvz)
#   }};
#
#   return SpatialCell::get_velocity_block(indices);
#}

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
         return self.__read_velocity_cells_old(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)
      else:
         # Uses new format:
         return self.__read_velocity_cells_new(cellid=cellid, cells_with_blocks=cells_with_blocks, blocks_per_cell=blocks_per_cell, cells_with_blocks_index=cells_with_blocks_index)


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
         return self.__read_blocks_old(cellid)
      else:
         # Uses new format
         return self.__read_blocks_new(cellid)

      return []

   def optimize_open_file(self):
      '''Opens the vlsv file for reading
         NOTE: THIS SHOULD ONLY BE USED FOR OPTIMIZATION PURPOSES. FILES ARE OPENED AND CLOSED AUTOMATICALLY UPON READING BUT IN THE CASE OF HEAVY READING OPERATIONS OPEN_FILE WILL OPTIMIZE THE PROCESS
      '''
      self.__fptr = open(self.__file_name,"rb")


   def optimize_close_file(self):
      '''Closes the vlsv file
         NOTE: THIS SHOULD ONLY BE USED FOR OPTIMIZATION PURPOSES. FILES ARE OPENED AND CLOSED AUTOMATICALLY UPON READING BUT IN THE CASE OF HEAVY READING OPERATIONS OPEN_FILE WILL OPTIMIZE THE PROCESS
      '''
      if self.__fptr.closed:
         return
      else:
         self.__fptr.close()
         return

   def optimize_clear_fileindex_for_cellid_blocks(self):
      ''' Clears a private variable containing number of blocks and offsets for particular cell ids
          USED FOR OPTIMIZATION PURPOSES ONLY WHEN READING VELOCITY CELLS IS NO LONGER NEEDED
      '''
      self.__fileindex_for_cellid_blocks = {}

   def optimize_clear_fileindex_for_cellid(self):
      ''' Clears a private variable containing cell ids and their locations
          USED FOR OPTIMIZATION PURPOSES ONLY WHEN READING VARIABLES IS NO LONGER NEEDED
      '''
      self.__fileindex_for_cellid = {}



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
from reduction import datareducers,data_operators
from collections import Iterable


class VlsvWriter(object):
   ''' Class for reading VLSV files
   '''
   file_name = ""
   def __init__(self, vlsvReader, file_name ):
      ''' Initializes the vlsv file (opens the file, reads the file footer and reads in some parameters)

          :param vlsvReader:    Some open vlsv file for creating an XML footer as well as the grid
          :param file_name:     Name of the vlsv file where to input data
      '''
      self.file_name = os.path.abspath(file_name)
      self.__fptr = open(self.file_name,"wb")

      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}

      self.__offset = 0
      # Write endianness
      np.array(0, dtype=np.uint64).tofile(self.__fptr)
      # Write xml_offset, for now put this to zero:
      np.array(0, dtype=np.uint64).tofile(self.__fptr)

      self.__initialize( vlsvReader )

   def __initialize( self, vlsvReader ):
      ''' Writes the xml footer as well as the cell ids from the vlsvReader to the file and everything else needed for the grid
      '''
      # Get the xml sheet:
      xml_root = vlsvReader._VlsvReader__xml_root

      # Get list of tags to write:
      tags = {}
      tags['PARAMETER'] = ''
      tags['PARAMETERS'] = ''
      tags['MESH_NODE_CRDS_X'] = ''
      tags['MESH_NODE_CRDS_Y'] = ''
      tags['MESH_NODE_CRDS_Z'] = ''
      tags['MESH'] = ''
      tags['MESH_DOMAIN_SIZES'] = ''
      tags['CellID'] = ''
      tags['MESH_BBOX'] = ''
      tags['COORDS'] = ''

      # Copy the xml root
      for child in xml_root:
         if child.tag in tags:
            if 'name' in child.attrib: name = child.attrib['name']
            else: name = ''
            if 'mesh' in child.attrib: mesh = child.attrib['mesh']
            else: mesh = ''
            tag = child.tag
            extra_attribs = {}
            for i in child.attrib.iteritems():
               if i[0] != 'name' and i[0] != 'mesh':
                  extra_attribs[i[0]] = i[1]
            data = vlsvReader.read( name=name, tag=tag, mesh=mesh )
            # Write the data:
            self.write( data=data, name=name, tag=tag, mesh=mesh, extra_attribs=extra_attribs )


   def copy_variables( self, vlsvReader ):
      ''' Copies all variables from vlsv reader to the file

      '''
      # Get the xml sheet:
      xml_root = vlsvReader._VlsvReader__xml_root

      # Get list of tags to write:
      tags = {}
      tags['VARIABLE'] = ''

      # Copy the xml root and write variables
      for child in xml_root:
         if child.tag in tags:
            if 'name' in child.attrib:
                name = child.attrib['name']
            else:
                name = ''
            if 'mesh' in child.attrib:
                mesh = child.attrib['mesh']
            else:
                mesh = ''
            tag = child.tag
            # Copy extra attributes:
            extra_attribs = {}
            for i in child.attrib.iteritems():
               if i[0] != 'name' and i[0] != 'mesh':
                  extra_attribs[i[0]] = i[1]
            data = vlsvReader.read( name=name, tag=tag, mesh=mesh )
            # Write the data:
            self.write( data=data, name=name, tag=tag, mesh=mesh, extra_attribs=extra_attribs )
      return

   def write_velocity_space( self, vlsvReader, cellid, blocks_and_values ):
      ''' Writes given velocity space into vlsv file

          :param vlsvReader:        Some open vlsv reader file with velocity space in the given cell id
          :param cellid:            Given cellid
          :param blocks_and_values: Blocks and values in list format e.g. [[block1,block2,..], [block1_values, block2_values,..]] where block1_values are velocity block values (list length 64)
      '''

      # Get cells_with_blocks, blocks_per_cell etc
      cells_with_blocks = np.array([cellid])
      number_of_blocks  = len(blocks_and_values)
      blocks_per_cell    = np.array([number_of_blocks])

      # Write them out
      self.write( data=cells_with_blocks, name='', mesh="SpatialGrid", tag="CELLSWITHBLOCKS" )
      self.write( data=blocks_per_cell, name='', mesh="SpatialGrid", tag="BLOCKSPERCELL" )

      # Write blockids and values
      self.write( data=blocks_and_values[0], name='', mesh="SpatialGrid", tag="BLOCKIDS" )
      self.write( data=blocks_and_values[1], name='avgs', mesh="SpatialGrid", tag="BLOCKVARIABLE" )

   def write(self, data, name, tag, mesh, extra_attribs={}):
      ''' Writes an array into the vlsv file

      :param name: Name of the data array
      :param tag:  Tag of the data array.
      :param mesh: Mesh for the data array
      :param extra_attribs: Dictionary with whatever xml attributes that should be defined in the array that aren't name, tag, or mesh

      :returns: True if the data was written successfully

      '''
      # Make sure the data is in numpy array format:
      data = np.atleast_1d(data)
      fptr = self.__fptr

      datatype = ''

      # Add the data into the xml data:
      child = ET.SubElement(parent=self.__xml_root, tag=tag)
      child.attrib["name"] = name
      child.attrib["mesh"] = mesh
      child.attrib["arraysize"] = len(np.atleast_1d(data))
      if extra_attribs != '':
         for i in extra_attribs.iteritems():
            child.attrib[i[0]] = i[1]
      if len(np.shape(data)) == 2:
         child.attrib["vectorsize"] = np.shape(data)[1]
         datatype = str(type(data[0][0]))
      elif len(np.shape(data)) > 2:
         print "ERROR, np.shape returned len(np.shape(data)) > 2"
         return False
      else:
         child.attrib["vectorsize"] = 1
         datatype = str(type(data[0]))

      # Parse the data types:
      if 'uint' in datatype:
         child.attrib["datatype"] = "uint"
      elif 'int' in datatype:
         child.attrib["datatype"] = "int"
      elif 'float' in datatype:
         child.attrib["datatype"] = "float"
      else:
         print "BAD DATATYPE"
         return False

      if '64' in datatype:
         child.attrib["datasize"] = 8
      elif '32' in datatype:
         child.attrib["datasize"] = 4
      else:
         print "BAD DATASIZE"
         return False

      current_offset = fptr.tell()
      # Info the xml about the file offset for the data:
      child.text = str(current_offset)

      data.tofile(fptr)

      # write the xml footer:
      self.__write_xml_footer()

   def __write_xml_footer( self ):
      # Write the xml footer:
      max_xml_size = 1000000
      if self.__fptr.closed:
         fptr = open(self.file_name,"wb")
      else:
         fptr = self.__fptr
      current_offset = fptr.tell()
      #self.__xml_root.write( fptr )
      # Convert everything to string:
      for child in self.__xml_root:
         for i in child.attrib.iteritems():
            child.attrib[i[0]] = str(child.attrib[i[0]])
      tree = ET.ElementTree( self.__xml_root)
      tree.write(fptr)
      # Write the offset (first 8 bytes = endianness):
      offset_endianness = 8
      fptr.seek( offset_endianness )
      # Write the offset:
      np.array(current_offset, dtype=np.uint64).tofile(fptr)
      # Go back to the previous offset:
      fptr.seek(current_offset)

   def close( self ):
      self.__write_xml_footer()
      self.__fptr.close()


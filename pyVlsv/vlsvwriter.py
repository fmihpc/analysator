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

class VlsvWriter(object):
   ''' Class for reading VLSV files
   '''
   file_name = ""
   def __init__(self, vlsvReader, file_name, copy_meshes=None):
      ''' Initializes the vlsv file (opens the file, reads the file footer and reads in some parameters)

          :param vlsvReader:    Some open vlsv file for creating an XML footer as well as the grid
          :param file_name:     Name of the vlsv file where to input data
          :param copy_meshes:   list of mesh names to copy, default all
          
      '''
      self.file_name = os.path.abspath(file_name)
      try:
         self.__fptr = open(self.file_name,"wb")
      except FileNotFoundError as e:
         print("No such path: ", self.file_name)
         raise e
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}

      self.__offset = 0
      # Write endianness
      np.array(0, dtype=np.uint64).tofile(self.__fptr)
      # Write xml_offset, for now put this to zero:
      np.array(0, dtype=np.uint64).tofile(self.__fptr)

      self.__initialize( vlsvReader, copy_meshes )

   def __initialize( self, vlsvReader, copy_meshes=None ):
      ''' Writes the xml footer as well as the cell ids from the vlsvReader to the file and everything else needed for the grid
      '''
      # Get the xml sheet:
      xml_root = vlsvReader._VlsvReader__xml_root

      # Get list of tags to write:
      tags = {}
      tags['PARAMETER'] = ''
      tags['PARAMETERS'] = ''
      tags['MESH_NODE_CRDS'] = ''
      tags['MESH_NODE_CRDS_X'] = ''
      tags['MESH_NODE_CRDS_Y'] = ''
      tags['MESH_NODE_CRDS_Z'] = ''
      tags['MESH_OFFSETS'] = ''
      tags['MESH'] = ''
      tags['MESH_DOMAIN_SIZES'] = ''
      tags['MESH_GHOST_DOMAINS'] = ''
      tags['MESH_GHOST_LOCALIDS'] = ''
      tags['CellID'] = ''
      tags['MESH_BBOX'] = ''
      tags['COORDS'] = ''

      # Copy the xml root
      for child in xml_root:
         if child.tag in tags:
            if 'name' in child.attrib: name = child.attrib['name']
            else: name = ''
            if 'mesh' in child.attrib: mesh = child.attrib['mesh']
            else: mesh = None
            tag = child.tag

            if copy_meshes is not None:
               if mesh is not None and not mesh in copy_meshes:
                  continue
               if tag == "MESH" and not name in copy_meshes:
                  continue

            extra_attribs = {}
            for i in child.attrib.items():
               if i[0] != 'name' and i[0] != 'mesh':
                  extra_attribs[i[0]] = i[1]
            data = vlsvReader.read( name=name, tag=tag, mesh=mesh )
            # Write the data:
            #print("writing",name, tag)

            self.write( data=data, name=name, tag=tag, mesh=mesh, extra_attribs=extra_attribs )


   def copy_variables( self, vlsvReader, varlist=None ):
      ''' Copies variables from vlsv reader to the file.
           varlist = None: list of variables to copy; if no
           varlist is provided, copy all variables (default)
      '''

      # Delegate to the variable list handler
      if (varlist is not None):
         self.copy_variables_list(vlsvReader, varlist)
         return

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
                mesh = None
            tag = child.tag
            # Copy extra attributes:
            extra_attribs = {}
            for i in child.attrib.items():
               if i[0] != 'name' and i[0] != 'mesh':
                  extra_attribs[i[0]] = i[1]
            data = vlsvReader.read( name=name, tag=tag, mesh=mesh )
            # Write the data:
            self.write( data=data, name=name, tag=tag, mesh=mesh, extra_attribs=extra_attribs )
      return

   def copy_variables_list( self, vlsvReader, vars ):
      ''' Copies variables in the list vars from vlsv reader to the file

      '''
      # Get the xml sheet:
      xml_root = vlsvReader._VlsvReader__xml_root

      # Get list of tags to write:
      tags = {}
      tags['VARIABLE'] = ''
      found_vars = []

      # Copy the xml root and write variables
      for child in xml_root:
         if child.tag in tags:
            if 'name' in child.attrib:
               name = child.attrib['name']
               if not name in vars:
                  continue
               else:
                  found_vars.append(name)
            else:
                continue
            if 'mesh' in child.attrib:
                mesh = child.attrib['mesh']
            else:
               if tag in ['VARIABLE']:
                  print('MESH required')
                  return
               mesh = None
            tag = child.tag
            # Copy extra attributes:
            extra_attribs = {}
            for i in child.attrib.items():
               if i[0] != 'name' and i[0] != 'mesh':
                  extra_attribs[i[0]] = i[1]
            data = vlsvReader.read( name=name, tag=tag, mesh=mesh )
            # Write the data:
            self.write( data=data, name=name, tag=tag, mesh=mesh, extra_attribs=extra_attribs )

      for name in [varname for varname in vars if varname not in found_vars]:
         varinfo = vlsvReader.read_variable_info(name)
         self.write_variable_info(varinfo, 'SpatialGrid', 1)
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
      child = ET.SubElement(self.__xml_root, tag)
      child.attrib["name"] = name
      if mesh is not None:
         child.attrib["mesh"] = mesh
      child.attrib["arraysize"] = len(np.atleast_1d(data))

      if len(np.shape(data)) == 2:
         child.attrib["vectorsize"] = np.shape(data)[1]
         datatype = str(type(data[0][0]))
      elif len(np.shape(data)) > 2:
         print("ERROR, np.shape returned len(np.shape(data)) > 2")
         return False
      else:
         child.attrib["vectorsize"] = 1
         if(len(data) == 0):
            if tag=="MESH_GHOST_DOMAINS" or tag=="MESH_GHOST_LOCALIDS":
               datatype="int32"
            else:
               print("Trying to extract datatype from an empty array. I will fail as usual, since this is not the special case that is guarded against!")
               datatype = str(type(data[0]))
         else:
            datatype = str(type(data[0]))

      # Parse the data types:
      if 'uint' in datatype:
         child.attrib["datatype"] = "uint"
      elif 'int' in datatype:
         child.attrib["datatype"] = "int"
      elif 'float' in datatype:
         child.attrib["datatype"] = "float"
      else:
         print("BAD DATATYPE: " + datatype)
         return False

      if '64' in datatype:
         child.attrib["datasize"] = 8
      elif '32' in datatype:
         child.attrib["datasize"] = 4
      else:
         print("BAD DATASIZE")
         return False
      
      if (extra_attribs != '') and (extra_attribs is not None):
         for i in extra_attribs.items():
            child.attrib[i[0]] = i[1]

      current_offset = fptr.tell()
      # Info the xml about the file offset for the data:
      child.text = str(current_offset)

      try:
         data.tofile(fptr)
      except:
         np.ma.getdata(data).tofile(fptr) # numpy maskedarray tofile not implemented yet

      # write the xml footer:
      self.__write_xml_footer()

   def write_variable_info(self, varinfo, mesh, unitConversion, extra_attribs={}):
      ''' Writes an array into the vlsv file as a variable; requires input of metadata required by VlsvReader
      :param varinfo: VariableInfo object containing
         -data: The variable data (array)
         -name: Name of the data array
         -latex: LaTeX string representation of the variable name
         -units: plaintext string representation of the unit
         -latexunits: LaTeX string representation of the unit
      :param mesh: Mesh for the data array
      :param unitConversion: string representation of the unit conversion to get to SI
      :param extra_attribs: Dictionary with whatever xml attributes that should be defined in the array that aren't name, tag, or mesh,
        or contained in varinfo. Can be used to overwrite varinfo values besids name.

      :returns: True if the data was written successfully

      '''

      return self.write(varinfo.data, varinfo.name, 'VARIABLE', mesh, extra_attribs={'variableLaTeX':varinfo.latex, 'unit':varinfo.units, 'unitLaTeX':varinfo.latexunits, 'unitConversion':unitConversion}.update(extra_attribs))


   def write_fgarray_to_SpatialGrid(self, reader, data, name, extra_attribs={}):
      # get a reader for the target file
      #print(data.shape[0:3], reader.get_fsgrid_mesh_size(), (data.shape[0:3] == reader.get_fsgrid_mesh_size()))
      if not (data.shape[0:3] == reader.get_fsgrid_mesh_size()).all():
         print("Data shape does not match target fsgrid mesh")
         return
      vgdata = reader.fsgrid_array_to_vg(data)
      self.write(vgdata, name, "VARIABLE", "SpatialGrid",extra_attribs)

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
         for i in child.attrib.items():
            child.attrib[i[0]] = str(child.attrib[i[0]])
      tree = ET.ElementTree( self.__xml_root)
      xml_footer_indent( self.__xml_root)
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

def xml_footer_indent(elem, level=0):
   i = "\n" + level*"   "
   if len(elem):
      if not elem.text or not elem.text.strip():
            elem.text = i + "   "
      if not elem.tail or not elem.tail.strip():
            elem.tail = i
      for elem in elem:
            xml_footer_indent(elem, level+1)
      if not elem.tail or not elem.tail.strip():
            elem.tail = i
   else:
      if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

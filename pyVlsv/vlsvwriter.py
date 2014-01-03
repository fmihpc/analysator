import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
from reduction import datareducers,data_operators
from collections import Iterable


class VlsvWriter(object):
   ''' Class for reading VLSV files
   '''
   def __init__(self, file_name):
      ''' Initializes the vlsv file (opens the file, reads the file footer and reads in some parameters)

          :param file_name:     Name of the vlsv file
      '''
      self.__file_name = file_name
      self.__fptr = open(self.__file_name,"wb")
      self.__xml_root = ET.fromstring("<VLSV></VLSV>")
      self.__fileindex_for_cellid={}
      self.__fileindex_for_cellid_blocks={}
      

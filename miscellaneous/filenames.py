import glob
import os

'''@package docstring
This module imports glob and os and has some miscellaneous functions for basic file manipulation
'''

def get_sorted_file_names(name="*.vlsv"):
   '''Gets the file names in the current directory and sorts them.
      :param name            Name of the file(s), for example "*.vlsv"
      returns a list of file names in sorted order
   '''
   fileNames=glob.glob(name)
   fileNames.sort()
   return fileNames
   

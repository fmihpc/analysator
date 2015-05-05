''' Helper module for reading and locating files. Nothing interesting here..

'''

import os, sys

def sorted_filenames(name="*.vlsv"):
   '''Gets the file names in the current directory and sorts them.
      :param name:            Name of the file(s), for example "*.vlsv"
      :returns: a list of file names in sorted order
   '''
   import glob
   fileNames=glob.glob(name)
   fileNames.sort()
   return fileNames


import os, sys

def sorted_filenames(name="*.vlsv"):
   import glob
   fileNames=glob.glob(name)
   fileNames.sort()
   return fileNames


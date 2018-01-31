import filemanagement
import socket, re

# Input current folder's path
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)))
# Input folder paths
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "miscellaneous")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyCalculations")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyCellDataReduction")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyMayaVi")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyPlots")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyVisit")
filemanagement.sys.path.insert(0, filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__)) + "/" + "pyVlsv")

# Import modules
try:
   import calculations
except ImportError:
   print "Note: Did not import calculations module"

try:
   import vlsvfile
except ImportError:
   print "Note: Did not import vlsvfile module"

import os
import matplotlib.pyplot as plt
if os.getenv('PTINTERACTIVE') != None:
   try:
      import grid
   except ImportError:
      print "Note: Did not import grid module"
   try:
      plt.switch_backend('TkAgg')
   except:
      print "Note: Unable to switch to TkAgg backend"
else:
   try:
      plt.switch_backend('Agg')
   except:
      print "Note: Unable to switch to Agg backend"

try:
   import plot
except ImportError:
   print "Note: Did not import plot module"

try:
   import miscellaneous
except ImportError:
   print "Note: Did not import miscellaneous"


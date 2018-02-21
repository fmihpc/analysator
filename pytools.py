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
except ImportError as e:
   print "Note: Did not import calculations module: ", e

try:
   import vlsvfile
except ImportError as e:
   print "Note: Did not import vlsvfile module: ", e

import os
import matplotlib.pyplot as plt

if os.getenv('PTNONINTERACTIVE') != None:
   try:
      plt.switch_backend('Agg')
   except:
      print "Note: Unable to switch to Agg backend"
else:
   try:
      import grid
   except ImportError as e:
      print "Note: Did not import grid module: ", e
   try:
      plt.switch_backend('TkAgg')
   except:
      print "Note: Unable to switch to TkAgg backend"

try:
   import plot
except ImportError as e:
   print "Note: Did not import plot module: ", e

try:
   import miscellaneous
except ImportError as e:
   print "Note: Did not import miscellaneous: ", e


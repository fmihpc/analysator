import filemanagement
# Input paths:
fullPath = filemanagement.os.path.dirname(filemanagement.os.path.abspath(__file__))
# Input current folder's path
filemanagement.sys.path.insert(0, fullPath)
# Input folder paths
filemanagement.sys.path.insert(0, fullPath + "/" + "miscellaneous")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyCalculations")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyCellDataReduction")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyMayaVi")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyPlots")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyVisit")
filemanagement.sys.path.insert(0, fullPath + "/" + "pyVlsv")

# Import modules
try:
   import calculations
except ImportError:
   print "Note: Did not import calculations module"
try:
   import vlsvfile
except ImportError:
   print "Note: Did not import vlsvfile module"
try:
   import grid
except ImportError:
   print "Note: Did not import grid module"
try:
   import plot
except ImportError:
   print "Note: Did not import plot module"

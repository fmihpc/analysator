import filemanagement
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
try:
   import grid
except ImportError:
   print "Note: Did not import grid module"
try:
   import plot
except ImportError:
   print "Note: Did not import plot module"
try:
   import miscellaneous
except ImportError:
   print "Note: Did not import miscellaneous"


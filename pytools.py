import os, sys
# Input paths:
fullPath = os.path.dirname(os.path.abspath(__file__))
# Input current folder's path
sys.path.insert(0, fullPath)
# Input folder paths
sys.path.insert(0, fullPath + "/" + "miscellaneous")
sys.path.insert(0, fullPath + "/" + "pyCalculations")
sys.path.insert(0, fullPath + "/" + "pyCellDataReduction")
sys.path.insert(0, fullPath + "/" + "pyMayaVi")
sys.path.insert(0, fullPath + "/" + "pyPlots")
sys.path.insert(0, fullPath + "/" + "pyVisit")
sys.path.insert(0, fullPath + "/" + "pyVlsv")

# Import modules
import pycalculations
import vlsvreader
import creategrid

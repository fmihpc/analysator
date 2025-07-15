""" A rudimentary wrapper from the vlsvVtkInterface. Add and load from Paraview
plugins menu.

Built on Kitware example:
https://github.com/Kitware/ParaView/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py
This module demonstrates various ways of adding
VTKPythonAlgorithmBase subclasses as filters, sources, readers,
and writers in ParaView"""

# https://github.com/Kitware/ParaView/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py

# This is module to import. It provides VTKPythonAlgorithmBase, the base class
# for all python-based vtkAlgorithm subclasses in VTK and decorators used to
# 'register' the algorithm with ParaView along with information about UI.
from paraview.util.vtkAlgorithm import *
import pytools as pt

# Use the analysator reader, but include Paraview decorators via inheriting/chaining
# necessary methods etc

#------------------------------------------------------------------------------
# A reader example.
#------------------------------------------------------------------------------

# To add a reader, we can use the following decorators
#   @smproxy.source(name="PythonCSVReader", label="Python-based CSV Reader")
#   @smhint.xml("""<ReaderFactory extensions="csv" file_description="Numpy CSV files" />""")
# or directly use the "@reader" decorator.
@smproxy.reader(name="PythonVLSVReader", label="Python-based VLSV Reader, outputs htg",
                extensions="vlsv", file_description="VLSV files")
class PythonPvVLSVReader(pt.vlsvfile.VlsvVtkReader):
    """A reader that wraps an Analysator VLSV file reader.
    """    
    def __init__(self):
        super().__init__()


    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="vlsv", file_description="Vlasiator VLSV files")
    def SetFileName(self, filename):
        """Specify filename for the file to read."""
        if filename is None or filename == "None":
            return
        print(filename)
        super().SetFileName(filename)

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return None

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return super().GetDataArraySelection()


#------------------------------------------------------------------------------
# Todo: some tests
#------------------------------------------------------------------------------

def test_PythonVLSVReader(fname):
    pass

if __name__ == "__main__":
    pass

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
# A writer example.
#------------------------------------------------------------------------------
'''
@smproxy.writer(extensions="npz", file_description="NumPy Compressed Arrays", support_reload=False)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkTable"], composite_data_supported=False)
class NumpyWriter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=0, inputType='vtkTable')
        self._filename = None

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, fname):
        """Specify filename for the file to write."""
        if self._filename != fname:
            self._filename = fname
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkTable
        from vtkmodules.numpy_interface import dataset_adapter as dsa

        table = dsa.WrapDataObject(vtkTable.GetData(inInfoVec[0], 0))
        kwargs = {}
        for aname in table.RowData.keys():
            kwargs[aname] = table.RowData[aname]

        import numpy
        numpy.savez_compressed(self._filename, **kwargs)
        return 1

    def Write(self):
        self.Modified()
        self.Update()
'''
#------------------------------------------------------------------------------
# A filter example.
#------------------------------------------------------------------------------
'''
@smproxy.filter()
@smproperty.input(name="InputTable", port_index=1)
@smdomain.datatype(dataTypes=["vtkTable"], composite_data_supported=False)
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class ExampleTwoInputFilter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1, outputType="vtkPolyData")

    def FillInputPortInformation(self, port, info):
        if port == 0:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
        else:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkTable")
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkTable, vtkDataSet, vtkPolyData
        input0 = vtkDataSet.GetData(inInfoVec[0], 0)
        input1 = vtkDataSet.GetData(inInfoVec[1], 0)
        output = vtkPolyData.GetData(outInfoVec, 0)
        # do work
        print("Pretend work done!")
        return 1

@smproxy.filter()
@smproperty.input(name="Input")
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class PreserveInputTypeFilter(VTKPythonAlgorithmBase):
    """
    Example filter demonstrating how to write a filter that preserves the input
    dataset type.
    """
    def __init__(self):
        super().__init__(nInputPorts=1, nOutputPorts=1, outputType="vtkDataSet")

    def RequestDataObject(self, request, inInfo, outInfo):
        inData = self.GetInputData(inInfo, 0, 0)
        outData = self.GetOutputData(outInfo, 0)
        assert inData is not None
        if outData is None or (not outData.IsA(inData.GetClassName())):
            outData = inData.NewInstance()
            outInfo.GetInformationObject(0).Set(outData.DATA_OBJECT(), outData)
        return super().RequestDataObject(request, inInfo, outInfo)

    def RequestData(self, request, inInfo, outInfo):
        inData = self.GetInputData(inInfo, 0, 0)
        outData = self.GetOutputData(outInfo, 0)
        print("input type =", inData.GetClassName())
        print("output type =", outData.GetClassName())
        assert outData.IsA(inData.GetClassName())
        return 1
'''

def test_PythonVLSVReader(fname):
    reader = PythonPvVLSVReader()
    reader.SetFileName(fname)
    reader.Update()
    assert reader.GetOutputDataObject(0).GetNumberOfRows() > 0

if __name__ == "__main__":
    #test_PythonSuperquadricSource()
    test_PythonVLSVReader("/tmp/data.csv")

    from paraview.detail.pythonalgorithm import get_plugin_xmls
    from xml.dom.minidom import parseString
    for xml in get_plugin_xmls(globals()):
        dom = parseString(xml)
        print(dom.toprettyxml(" ","\n"))

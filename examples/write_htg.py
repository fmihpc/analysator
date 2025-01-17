import vtk
import pytools as pt

# This initializes a hypertreegrid from the given reader.
reader = pt.vlsvfile.VlsvVtkReader()
reader.SetFileName("../../bulk.0002189.vlsv")

# These functions grab one SpatialGrid variable and map that to 
# the hypertreegrid. Variable vector sizes of 1,2,3,4,9 supported.
reader.addArrayFromVlsv("vg_b_vol")
reader.addArrayFromVlsv("vg_beta_star")

reader.Update()

writer = vtk.vtkXMLHyperTreeGridWriter()
writer.SetFileName("output_EGE.htg")
writer.SetInputData(reader.GetOutputDataObject(0))
writer.Write()

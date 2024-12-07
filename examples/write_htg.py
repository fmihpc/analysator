import vtk
import pytools as pt

# This initializes a hypertreegrid from the given reader.
reader = pt.vlsvfile.VlsvReader("../../bulk.0002189.vlsv")
htg = pt.vlsvfile.vtkVlsvHyperTreeGrid(reader)

# These functions grab one SpatialGrid variable and map that to 
# the hypertreegrid. Variable vector sizes of 1,2,3,4,9 supported.
htg.addArrayFromVlsv("vg_b_vol")
htg.addArrayFromVlsv("vg_beta_star")

writer = vtk.vtkXMLHyperTreeGridWriter()
writer.SetFileName("output_EGE.htg")
writer.SetInputData(htg)
writer.Write()

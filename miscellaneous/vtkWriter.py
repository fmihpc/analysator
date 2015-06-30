#!/usr/bin/env python
 
import vtk
 
filename = "test.vtk"

# Create a line object
#######################################
# Create five points. 
origin = [0.0, 0.0, 0.0]
p0 = [1.0, 0.0, 0.0]
p1 = [0.0, 1.0, 0.0]
p2 = [0.0, 1.0, 2.0]
p3 = [1.0, 2.0, 3.0]

# Create a vtkPoints object and store the points in it
points = vtk.vtkPoints()
points.InsertNextPoint(origin)
points.InsertNextPoint(p0)
points.InsertNextPoint(p1)
points.InsertNextPoint(p2)
points.InsertNextPoint(p3)

# Create a cell array to store the lines in and add the lines to it
lines = vtk.vtkCellArray()

for i in range(3):
  line = vtk.vtkLine()
  line.GetPointIds().SetId(0,i)
  line.GetPointIds().SetId(1,i+1)
  lines.InsertNextCell(line)

# Create a polydata to store everything in
linesPolyData = vtk.vtkPolyData()

# Add the points to the dataset
linesPolyData.SetPoints(points)

# Add the lines to the dataset
linesPolyData.SetLines(lines)
#######################################


# Write the Line object
#######################################
polyDataWriter = vtk.vtkPolyDataWriter()
polyDataWriter.SetFileName(filename)
polyDataWriter.SetInput(linesPolyData)
polyDataWriter.Write()
#######################################

# Read and visualize the written line
#######################################
reader = vtk.vtkPolyDataReader()
reader.SetFileName(filename)
 
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(reader.GetOutput())
else:
    mapper.SetInputConnection(reader.GetOutputPort())
 
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
 
# Create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
# Assign actor to the renderer
ren.AddActor(actor)
 
# Enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()
#######################################

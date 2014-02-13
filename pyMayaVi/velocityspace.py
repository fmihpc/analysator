import mayavi.mlab
import mayavi
from tvtk.api import tvtk
import numpy as np

def generate_velocity_grid( vlsvReader, cellid, iso_surface=False ):
   '''Generates a velocity grid from a given spatial cell id
      :param cellid:           The spatial cell's ID
      :param iso_surface:      If true, plots the iso surface
   '''
   # Create nodes
   # Get velocity blocks and avgs:
   blocksAndAvgs = vlsvReader.read_blocks(cellid)
   if len(blocksAndAvgs) == 0:
      print "CELL " + str(cellid) + " HAS NO VELOCITY BLOCK"
      return False
   # Create a new scene
   #engine.new_scene()
   #mayavi.mlab.set_engine(engine)
   # Create a new figure
   #figure = mayavi.mlab.clf()
   #figure.scene.disable_render = True
   blocks = blocksAndAvgs[0]
   avgs = blocksAndAvgs[1]
   # Get nodes:
   nodesAndKeys = vlsvReader.construct_velocity_cell_nodes(blocks)
   # Create an unstructured grid:
   points = nodesAndKeys[0]
   tets = nodesAndKeys[1]
   tet_type=tvtk.Voxel().cell_type#VTK_VOXEL

   ug=tvtk.UnstructuredGrid(points=points)
   # Set up the cells
   ug.set_cells(tet_type,tets)
   # Input data
   values=np.ravel(avgs)
   ug.cell_data.scalars=values
   ug.cell_data.scalars.name='avgs'
   #figure.scene.disable_render = False
   d = mayavi.mlab.pipeline.add_dataset(ug)
   ptdata = mayavi.mlab.pipeline.cell_to_point_data(d)
   iso = mayavi.mlab.pipeline.iso_surface(ptdata, contours=[1e-15,1e-14,1e-12], opacity=0.3)
   mayavi.mlab.show()


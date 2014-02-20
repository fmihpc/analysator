from traits.api import HasTraits, Instance, Property, Button, Enum
from mayavi.core.ui.engine_view import EngineView
from traits.api import HasTraits, Range, Instance, \
                    on_trait_change
from traitsui.api import View, Item, HGroup, Group
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import \
                    MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
import vlsvfile
from numpy import mgrid, empty, sin, pi, ravel
import pylab as pl
from tvtk.api import tvtk
import traits.api
import mayavi.api
import mayavi.mlab
import numpy as np
import signal
import threading
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.modules.outline import Outline
from mayavi.modules.surface import Surface
from mayavi.modules.vectors import Vectors

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

class VelocitySpace(HasTraits):
   ''' A class for plotting velocity space with MayaVi

   '''
   picker = Enum('None',
                 "Cut_through")

   args = ""

   cut_through = []

   scene = Instance(MlabSceneModel, ())

   engine_view = Instance(EngineView)

   current_selection = Property

   # Define the view:
   view = View(
                HGroup(
                   Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                      height=250, width=300, show_label=False, resizable=True),
                   Group(
                      #'cell_pick',
                      'picker',
                      'args',
                      show_labels=True
                   ),
                ),
                resizable=True,
            )


   def __init__(self, vlsvReader, cellid, iso_surface=False, threaded=True, **traits):
      ''' Initializes the class and loads the mayavi grid

          :param vlsvReader:        Some vlsv reader with a file open
          :type vlsvReader:         :class:`vlsvfile.VlsvReader`
          :param threaded:          Boolean value for using threads or not using threads to draw the grid (threads enable interactive mode)
      '''
      HasTraits.__init__(self, **traits)
      self.__vlsvReader = vlsvReader
      self.engine_view = EngineView(engine=self.scene.engine)
      self.__engine = self.scene.engine
      self.__picker = []
      self.__mins = []
      self.__maxs = []
      self.__last_pick = []
      self.__structured_figures = []
      self.__unstructured_figures = []
      self.__thread = []
      self.__load_grid( cellid=cellid, iso_surface=iso_surface, threaded=threaded )

   def __load_grid( self, cellid, iso_surface, threaded=True ):
      ''' Creates a grid and inputs scalar variables from a vlsv file
          :param variable:        Name of the variable to plot
          :param operator:        Operator for the variable
          :param threaded:        Boolean value for using threads or not using threads to draw the grid (threads enable interactive mode)
      '''
      # Draw the grid:
      if threaded == True:
         thread = threading.Thread(target=self.__generate_grid, args=( cellid, iso_surface ))
         thread.start()
      else:
         self.__generate_grid( cellid, iso_surface )


   def __generate_grid( self, cellid, iso_surface=False ):
      '''Generates a velocity grid from a given spatial cell id
         :param cellid:           The spatial cell's ID
         :param iso_surface:      If true, plots the iso surface
      '''

      ''' Generates a grid from given data
          :param mins:           An array of minimum coordinates for the grid for ex. [-100, 0, 0]
          :param maxs:           An array of maximum coordinates for the grid for ex. [-100, 0, 0]
          :param cells:          An array of number of cells in x, y, z direction
          :param datas:          Scalar data for the grid e.g. array([ cell1Rho, cell2Rho, cell3Rho, cell4Rho, .., cellNRho ])
          :param names:          Name for the scalar data
      '''
      # Create nodes
      # Get velocity blocks and avgs:
      blocksAndAvgs = self.__vlsvReader.read_blocks(cellid)
      if len(blocksAndAvgs) == 0:
         print "CELL " + str(cellid) + " HAS NO VELOCITY BLOCK"
         return False

      blocks = blocksAndAvgs[0]
      avgs = blocksAndAvgs[1]
      # Get nodes:
      nodesAndKeys = self.__vlsvReader.construct_velocity_cell_nodes(blocks)
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

      # Plot B if possible:
      def plot_B( name ):
         ''' Helper function for plotting B vector (name can change from B_vol to B)
             :param name:         Name of the B vector ( "B_vol" or "B" )
         '''
         # Read B vector and plot it:
         B = self.__vlsvReader.read_variable(name=name,cellids=cellid)
         points2 = np.array([[0,0,0]])
         ug2 = tvtk.UnstructuredGrid(points=points2)
         ug2.point_data.vectors = [(B * 8000000000000) / np.linalg.norm( B )]
         ug2.point_data.vectors.name = 'B_vector'
         #src2 = VTKDataSource(data = ug2)
         d2 = self.scene.mlab.pipeline.add_dataset(ug2)
         #mayavi.mlab.add_module(Vectors())
         vec = self.scene.mlab.pipeline.vectors(d2)
         vec.glyph.mask_input_points = True
         vec.glyph.glyph.scale_factor = 100000

      if self.__vlsvReader.check_variable( "B" ) == True:
         plot_B( "B" )
      elif self.__vlsvReader.check_variable( "B_vol" ) == True:
         plot_B( "B_vol" )

      # Visualize
      d = self.scene.mlab.pipeline.add_dataset(ug)
      if iso_surface == False:
         iso = self.scene.mlab.pipeline.surface(d)
      else:
         ptdata = self.scene.mlab.pipeline.cell_to_point_data(d)
         iso = self.scene.mlab.pipeline.iso_surface(ptdata, contours=[1e-15,1e-14,1e-12], opacity=0.3)

      # Visualize the data
      d = self.scene.mlab.pipeline.add_dataset(ug)
      iso = self.scene.mlab.pipeline.surface(d)

      # Configure traits
      self.configure_traits()
      return


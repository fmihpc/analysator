from traits.api import HasTraits, Instance, Property, Button, Enum
from mayavi.core.ui.engine_view import EngineView
from traits.api import HasTraits, Range, Instance, \
                    on_trait_change
from traitsui.api import View, Item, HGroup, Group
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import \
                    MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene
import vlsvreader
from numpy import mgrid, empty, sin, pi, ravel
from tvtk.api import tvtk
import mayavi.api
import mayavi.mlab
import numpy as np
import signal

#Catch SIGINT as mayavi (VTK) has disabled the normal signal handler
def SigHandler(SIG, FRM):
    print "Ctrl+C"
    return
signal.signal(signal.SIGINT, SigHandler)


class MayaviPlots(HasTraits):
   '''Class for constructing plots with MayaVi
   '''
   cell_pick = Button('cell_pick')
   test_attribute2 = Button('test2')


   scene = Instance(MlabSceneModel, ())

   engine_view = Instance(EngineView)

   current_selection = Property

   # Define the view:
   view = View(
                HGroup(
                   Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                      height=250, width=300, show_label=False, resizable=True),
                   Group(
                      'cell_pick',
                      'test_attribute2',
                      show_labels=False
                   ),
                ),
                resizable=True,
            )


   def __init__(self, vlsvReader, **traits):
      HasTraits.__init__(self, **traits)
      print "Constructing mayavi plot"
      self.__vlsvReader = vlsvReader
      self.engine_view = EngineView(engine=self.scene.engine)
      self.__engine = self.scene.engine
      self.__pick_mode = True
      self.__picker = []

   def __picker_callback( self, picker ):
      """ This gets called when clicking on a cell
      """
      point_id = picker.cell_id
      print "CELL ID: " + str(point_id+1)
      # NOTE: In vlasiator cell ids start from 1, in mayavi they start from 0, hence the +1
      self.__generate_velocity_grid(point_id+1)
   
   def __generate_grid( self, mins, lengths, cells, datas, names ):
      ''' Generates a grid from given data
          :param mins           An array of minimum coordinates for the grid for ex. [-100, 0, 0]
          :param lengths        An array of cell lengths (the cell's lengths in x, y, z direction)
          :param cells          An array of number of cells in x, y, z direction
          :param datas          Scalar data for the grid e.g. array([ cell1Rho, cell2Rho, cell3Rho, cell4Rho, .., cellNRho ])
          :param names          Name for the scalar data
      '''
      self.figure = self.scene.mlab.gcf()
      # Create nodes
      x, y, z = mgrid[mins[0]:lengths[0]*(cells[0]+1):(cells[0]+1)*complex(0,1), mins[1]:lengths[1]*(cells[1]+1):(cells[1]+1)*complex(0,1), mins[2]:lengths[2]*(cells[2]+1):(cells[2]+1)*complex(0,1)]
      # Cell coordinates:
      x2 = 0.1*0.5 + np.arange(4)/4.0*0.1
      y2 = 0.1*0.5 + np.arange(4)/4.0*0.1
      z2 = 0.1*2.0/5.0*0.5 + np.arange(1)
      
      # Create points for the nodes:
      pts = empty(z.shape + (3,), dtype=float)
      pts[...,0] = x
      pts[...,1] = y
      pts[...,2] = z
      
      # Input scalars
      scalars = np.array(datas)
      
      # We reorder the points, scalars and vectors so this is as per VTK's
      # requirement of x first, y next and z last.
      pts = pts.transpose(2, 1, 0, 3).copy()
      pts.shape = pts.size/3, 3
      scalars = scalars.T.copy()
      
      # Create the dataset.
      sg = tvtk.StructuredGrid(dimensions=x.shape, points=pts)
      sg.cell_data.scalars = ravel(scalars.copy())
      sg.cell_data.scalars.name = names
      
      
      # Visualize the data
      d = self.scene.mlab.pipeline.add_dataset(sg)
      iso = self.scene.mlab.pipeline.surface(d)#CONTINUE

      # Configure traits
      self.configure_traits()
      

   def __generate_velocity_grid( self, cellid ):
      '''Generates a velocity grid from a given spatial cell id
         :param cellid           The spatial cell's ID
      '''
      # Create nodes
      # Get velocity blocks and avgs:
      blocksAndAvgs = self.__vlsvReader.read_blocks(cellid)
      if len(blocksAndAvgs) == 0:
         print "CELL " + str(cellid) + " HAS NO VELOCITY BLOCK"
         return False
      # Create a new scene
      self.__engine.new_scene()
      mayavi.mlab.set_engine(self.__engine)
      # Create a new figure
      figure = mayavi.mlab.gcf(engine=self.__engine)
      figure.scene.disable_render = True
      blocks = blocksAndAvgs[0]
      avgs = blocksAndAvgs[1]
      # Get nodes:
      nodesAndKeys = self.__vlsvReader.construct_velocity_cell_nodes(blocks)
      # Create an unstructured grid:
      points = nodesAndKeys[0]
      tets = nodesAndKeys[1]
      tet_type=tvtk.Voxel().cell_type#VTK_VOXEL

      ug=tvtk.UnstructuredGrid(points=points)
      #Thissetsupthecells.
      ug.set_cells(tet_type,tets)
      #Attributedata.
      values=np.ravel(avgs)
      ug.cell_data.scalars=values
      ug.cell_data.scalars.name='avgs'
      d = mayavi.mlab.pipeline.add_dataset(ug)
      iso = mayavi.mlab.pipeline.surface(d)
      figure.scene.disable_render = False
      return True

   def __do_nothing( self, picker ):
      return

   # Trait events:
   @on_trait_change('scene.activated')
   def set_mouse_click( self ):
      # Temporary bug fix (MayaVi needs a dummy pick to be able to remove cells callbacks from picker.. )
      self.figure.on_mouse_pick( self.__do_nothing, type='cell' )
      # Cell picker
      self.__picker = self.figure.on_mouse_pick( self.__picker_callback, type='cell' )

   @on_trait_change('cell_pick')
   def switch_cell_pick(self):
      if self.__pick_mode == True:
         self.figure.on_mouse_pick( self.__picker_callback, type='cell', remove=True )
         self.__pick_mode = False
         print "Turning cell pick off"
      else:
         self.figure.on_mouse_pick( self.__picker_callback, type='cell' )
         self.__pick_mode = True
         print "Turning cell pick on"

   def load_grid( self, variable ):
      ''' Creates a grid and inputs scalar variables from a vlsv file
          :param variable        Name of the variable to plot
      '''
      # Get the cell params:
      mins = np.array([self.__vlsvReader.read_parameter("xmin"), self.__vlsvReader.read_parameter("ymin"), self.__vlsvReader.read_parameter("zmin")])
      cells = np.array([self.__vlsvReader.read_parameter("xcells_ini"), self.__vlsvReader.read_parameter("ycells_ini"), self.__vlsvReader.read_parameter("zcells_ini")])
      maxs = np.array([self.__vlsvReader.read_parameter("xmax"), self.__vlsvReader.read_parameter("ymax"), self.__vlsvReader.read_parameter("zmax")])
      lengths = (maxs - mins) / cells.astype(float)
      # Get the variables:
      index_for_cellid_dict = self.__vlsvReader.get_cellid_locations()
      variable_array = self.__vlsvReader.read_variables( name=variable )
      # Sort the dictionary by cell id
      import operator
      sorted_index_for_cellid_dict = sorted(index_for_cellid_dict.iteritems(), key=operator.itemgetter(0))
      # Add the variable values:
      variable_array_sorted = []
      for i in sorted_index_for_cellid_dict:
         variable_array_sorted.append(variable_array[i[1]])
      # Draw the grid:
      self.__generate_grid( mins=mins, lengths=lengths, cells=cells, datas=variable_array_sorted, names=variable )






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
from variable import get_data, get_name, get_units
from mayavi.modules.labels import Labels
from mayavigrid import MayaviGrid

def kill_subprocess( pipe, wait ):
  import time
  time.sleep(wait)
  pipe.terminate()

def call_particle_pusher( particlepusherinterface, vlsvReader, coordinates, args ):
   ''' Calls the particle pusher with the given coordinates. See also the picker in the class MayaviGrid
   '''
   
   new_coordinates = []
   for coordinate in coordinates:
     new_coordinates.append( coordinate )
   
   # Read in the velocity coordinates:
   vlsvReader.read_variable
   for i in xrange(3):
     new_coordinates.append(0)
   
   particlepusherinterface.particle_coordinates.append(new_coordinates)
   
   if int(args) == len(particlepusherinterface.particle_coordinates):
     import subprocess
     parse_args = []
     
     parse_args.append("/home/otto/vlasiator/particle_post_pusher")
     parse_args.append("--run_config")
     parse_args.append("/home/otto/vlasiator/particles/particles.cfg")
     
     #for coordinate in coordinates:
     #  parse_args.append(str(coordinate))

     pipe = subprocess.Popen(parse_args,stdin=subprocess.PIPE, stderr=subprocess.PIPE)

     for coordinates in particlepusherinterface.particle_coordinates:
        str_coords = ""
        for coordinate in coordinates:
          str_coords = str_coords + str(coordinate)
          str_coords = str_coords + " "
        str_coords = str_coords+"\n"
        pipe.stdin.write(str_coords)
        print str_coords

     pipe.stdin.write(str_coords)
     pipe.stdin.close()
     import thread
     thread.start_new_thread( kill_subprocess, (pipe, 50*len(particlepusherinterface.particle_coordinates)) )
     particlepusherinterface.particle_coordinates = []

class Particlepusherinterface(MayaviGrid):
   ''' This class is used to plot the data in a vlsv file as a mayavi grid The following will bring up a new window and plot the grid in the vlsv file:

   .. code-block:: python

      grid = pt.grid.MayaviGrid(vlsvReader=f, variable="rho", operator='pass', threaded=False)

   Once you have the window open you can use the picker tool in the right-upper corner and use various point-click tools for analyzing data.

   Picker options:
   
   **None** Does nothing upon clicking somewhere in the grid
   
   **Velocity_space** Plots the velocity space at a specific position upon clicking somewhere in the grid Note: If the vlsv file does not have the velocity space at the position where you are clicking, this will not work
   
   **Velocity_space_iso_surface** Plots the velocity space at a specific position upon clicking somewhere in the grid in iso-surface plotting style Note: If the vlsv file does not have the velocity space at the position where you are clicking, this will not work
   
   **Velocity_space_nearest_cellid** Plots the velocity space of the closest cell id to the picking point Note: If the vlsv file does not have velocity space saved at all, this will not work
   
   **Velocity_space_nearest_cellid_iso_surface** Plots the velocity space of the closest cell id to the picking point in iso-surface plotting style Note: If the vlsv file does not have velocity space saved at all, this will not work
   
   **Pitch_angle** Plots the pitch angle distribution at the clicking position Note: If the vlsv file does not have the velocity space at the position where you are clicking, this will not work
   
   **Gyrophase_angle** Plots the gyrophase angle distribution at the clicking position Note: If the vlsv file does not have the velocity space at the position where you are clicking, this will not work
   
   **Cut_through** Is used to plot or save the cut-through between two clicking points. This option requires you to use the args section at top-left. To use the args section to plot variables you must write for example: **plot rho B,x E,y** Upon clicking at two points a new window would open with a cut-through plot of rho, x-component of B and y-component of E Alternatively, you can save the cut-through to a variable in the MayaviGrid class by typing instead: **rho B,x E,y** and then going to the terminal and typing

   **Particle_pusher**
   
   .. code-block:: python

      cut_through_data = grid.cut_through
      print cut_through_data

   '''
   picker = Enum('None',
                 'Velocity_space',
                 "Velocity_space_nearest_cellid",
                 'Velocity_space_iso_surface',
                 'Velocity_space_nearest_cellid_iso_surface',
                 "Pitch_angle",
                 "Gyrophase_angle",
                 "Cut_through",
                 "Particle_pusher")

   picker_functions = {"Particle_pusher": call_particle_pusher} # Picker functions (here we define what happens when the picker has the "Particle_pusher" option on and it clicks on something (In this case calls the __call_particle_pusher functions

   particle_coordinates = []

   pipe = []












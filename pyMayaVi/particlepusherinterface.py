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

def draw_streamlines( particlepusherinterface, coordinates_map ):
   ''' Draws streamlines into the particlepusherinterface
   '''
   # Loop through particles and their coordinates:
   for i in coordinates_map.iteritems():
      particle_number = i[0]
      coordinates_list = np.array(i[1])
      particlepusherinterface.draw_streamline_long( np.transpose(coordinates_list) )
   return

def update_coordinates( coordinates_map, output_parsed ):
   ''' Updated coordinates in the read_subprocess

   '''
   particle_number = int(output_parsed[0])
   coordinates = np.array([float(output_parsed[2]), float(output_parsed[3]), 0])
   if particle_number in coordinates_map:
      coordinates_map[particle_number].append(coordinates)
   else:
      coordinates_map[particle_number] = [coordinates]
   return coordinates_map

def read_subprocess( particlepusherinterface, pipe ):
   ''' Read in input from the particles pusher and processes it by drawing streamlines

       :param particlepusherinterface: The Particlepusherinterface class where we launched the particle pusher and which has a MayaVi grid open
       :param pipe: Popen process (subprocess) with stdout open
   '''
   # Read the output and draw streamlines:
   coordinates_map = {} # We want to keep track of coordinates for each particle in the form: {0: [coordinate1, coordinate2, ...], 2: [coordinate1, coordinate2, ..], ..}
   import re
   iterator = 1
   while True:
      output = pipe.stdout.readline()
      if output == "":
        break;
      output_parsed = re.split('\t| |\n', output)
      # Update coordinates:
      coordinates_map = update_coordinates( coordinates_map, output_parsed )
      print "Iteration " + str(iterator) + " done!"
      iterator = iterator + 1
   # Draw the streamlines:
   draw_streamlines(particlepusherinterface, coordinates_map)
   # Kill the proces
   pipe.terminate()

def call_particle_pusher( particlepusherinterface, vlsvReader, coordinates, args ):
   ''' Calls the particle pusher with the given coordinates. See also the picker in the class MayaviGrid

       :param particlepusherinterface: Particle pusher class
       :param vlsvReader: a vlsvReader file with a file open
       :param coordinates: Some coordinates on a mayavi grid
       :param args: Arguments (these are used to parse how many particles we want to launch currently)
   '''
   
   # Input new coordinates ( This is vx, vy, vz )
   new_coordinates = []
   for coordinate in coordinates:
     new_coordinates.append( coordinate )
   
   # Read in the velocity coordinates:
   bulk_V = vlsvReader.read_variable(name="v", cellids=vlsvReader.get_cellid(new_coordinates))
   for i in xrange(3):
     new_coordinates.append( bulk_V[i] )
   
   # Input the new coordinates into the particle pusher
   particlepusherinterface.particle_coordinates.append(new_coordinates)
   
   # Check if this is the amount of particles we want to input 
   user_defined_input = int(args)
   current_number_of_particles = len(particlepusherinterface.particle_coordinates)
   if user_defined_input <= current_number_of_particles:
     
     # Call the particle pusher
     import subprocess
     parse_args = []
     
     #Executable location
     parse_args.append("/home/otto/vlasiator/particle_post_pusher")
     # Options
     parse_args.append("--run_config")
     # CFG location
     parse_args.append("/home/otto/vlasiator/particles/particles.cfg")

     # Open a pipe for the process (get input, output and error output to the pipe)
     pipe = subprocess.Popen(parse_args,stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

     # Iterate through particles and write into the stdin (this is required by the particle pusher)
     for coordinates in particlepusherinterface.particle_coordinates:
        coordinates_to_output = ""
        for coordinate in coordinates:
          coordinates_to_output = coordinates_to_output + str(coordinate) #Append
          coordinates_to_output = coordinates_to_output + " "
        coordinates_to_output = coordinates_to_output+"\n"
        pipe.stdin.write(coordinates_to_output) # Write the coordinates to the std
     # Close the stdin
     pipe.stdin.close()
     # Read the output in a separate thread and kill the process at the end:
     import thread
     #thread.start_new_thread( read_subprocess, (particlepusherinterface, pipe) )
     read_subprocess(particlepusherinterface, pipe)
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












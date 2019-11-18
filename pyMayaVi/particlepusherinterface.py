# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

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
   for i in coordinates_map.items():
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
      print("Iteration " + str(iterator) + " done!")
      iterator = iterator + 1
   # Draw the streamlines:
   draw_streamlines(particlepusherinterface, coordinates_map)
   # Kill the proces
   pipe.terminate()

def get_file_name( name ):
   lower=-2; upper=-2;
   for i in range(len(name)):
      if name[i].isdigit():
         if upper < i-1:
            lower = i;
         upper = i
   return name[0:lower] + "%0" + str(upper+1-lower) + "i" + name[upper+1:]


def call_particle_pusher( particlepusherinterface, coordinates_list, args ):
   ''' Launches the particle pusher

   '''
   # Call the particle pusher
   import subprocess
   parse_args = []


   # Executable location
   parse_args.append(str(particlepusherinterface.particlepushercommand))
   # Options for the particle pusher:
   for i in range(1,len(args)):
      parse_args.append(args[i])
   # File location:
   parse_args.append("--particles.input_filename_pattern")
   print(particlepusherinterface.vlsvReader.file_name)
   parse_args.append(get_file_name(particlepusherinterface.vlsvReader.file_name))
   print(parse_args)

   # Open a pipe for the process (get input, output and error output to the pipe)
   pipe = subprocess.Popen(parse_args,stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

   # Iterate through particles and write into the stdin (this is required by the particle pusher)
   for coordinates in coordinates_list:
      coordinates_to_output = ""
      for coordinate in coordinates:
        coordinates_to_output = coordinates_to_output + str(coordinate) #Append
        coordinates_to_output = coordinates_to_output + " "
      coordinates_to_output = coordinates_to_output+"\n"
      pipe.stdin.write(coordinates_to_output) # Write the coordinates to the std
   # Close the stdin
   pipe.stdin.close()
   # Read the output in a separate thread and kill the process at the end:
   import _thread
   #thread.start_new_thread( read_subprocess, (particlepusherinterface, pipe) )
   read_subprocess(particlepusherinterface, pipe)
   particlepusherinterface.particle_coordinates = []




def call_particle_pusher_velocity_sampling( particlepusherinterface, vlsvReader, coordinates, args ):
   ''' Calls the particle pusher with the given coordinates and inputs particles in random places within the velocity distribution function in given coordinates

       :param particlepusherinterface: Particle pusher class
       :param vlsvReader: a vlsvReader file with a file open
       :param coordinates: Some coordinates on a mayavi grid
       :param args: Arguments (these are used to parse how many particles we want to launch currently)
   '''
   if len(args) <= 1:
      print("Bad args field, should be <number of velocity samples> <particle pusher options>")
      return
   print("COORDS:")
   print(coordinates)
   # Note: vlsvReader must be VlasiatorReader type
   vlasiatorReader = vlsvReader
   coordinates = vlasiatorReader.get_nearest_coordinates_with_distribution( coordinates )
   print("NEAREST COORDS:")
   print(coordinates)
   # Get the cellid:
   cellid = vlasiatorReader.get_cellid( coordinates )

   
   # Read in the velocity space
   velocity_cell_map = vlasiatorReader.read_velocity_cells( cellid )
   velocity_cell_coordinates = vlasiatorReader.get_velocity_cell_coordinates(list(velocity_cell_map.keys()))
   number_of_particles = int(args[0])
   step = int(float(len(velocity_cell_coordinates))/float(number_of_particles))
   # Input particles:
   new_coordinates = []
   for i in range(number_of_particles):
      index = i*step
      velocity_coordinate = velocity_cell_coordinates[index]
      new_coordinates.append([coordinates[0], coordinates[1], coordinates[2], velocity_coordinate[0], velocity_coordinate[1], velocity_coordinate[2]])

   call_particle_pusher( particlepusherinterface, new_coordinates, args )

def call_particle_pusher_bulk_v( particlepusherinterface, vlsvReader, coordinates, args ):
   ''' Calls the particle pusher with the given coordinates. See also the picker in the class MayaviGrid

       :param particlepusherinterface: Particle pusher class
       :param vlsvReader: a VlasiatorReader file with a file open
       :param coordinates: Some coordinates on a mayavi grid
       :param args: Arguments (these are used to parse how many particles we want to launch currently)
   '''
   if len(args) <= 1:
      print("Bad args field, should be <number of velocity samples> <particle pusher options>")
      return

   # Input new coordinates ( This is vx, vy, vz )
   new_coordinates = []
   for coordinate in coordinates:
     new_coordinates.append( coordinate )
   
   # Read in the velocity coordinates:
   bulk_V = vlsvReader.read_variable(name="v", cellids=vlsvReader.get_cellid(new_coordinates))
   for i in range(3):
     new_coordinates.append( bulk_V[i] )
   
   # Input the new coordinates into the particle pusher
   particlepusherinterface.particle_coordinates.append(new_coordinates)
   
   # Check if this is the amount of particles we want to input 
   user_defined_input = int(args[0])
   current_number_of_particles = len(particlepusherinterface.particle_coordinates)

   if user_defined_input <= current_number_of_particles:
     call_particle_pusher( particlepusherinterface, particlepusherinterface.particle_coordinates, args )

def draw_B_stream( particlepusherinterface, vlsvReader, coordinates, args ):
   ''' Draw B streamlines into the particlepusherinteraface
   '''
   # Draw the streamlines for vlsvreader:
   from fieldtracer import static_field_tracer
   stream=static_field_tracer( vlsvReader, x0=coordinates, max_iterations=10, dx=10000, direction='+' )
   print(stream)
   print(np.transpose(stream))
   particlepusherinterface.draw_streamline_long( np.transpose(stream) )


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
                 "Particle_pusher_bulk_v",
                 "Particle_pusher_velocity_sampling",
                 "Particle_pusher_B_stream")

   picker_functions = {"Particle_pusher_bulk_v": call_particle_pusher_bulk_v,
                       "Particle_pusher_velocity_sampling": call_particle_pusher_velocity_sampling, # Picker functions (here we define what happens when the picker has the "Particle_pusher_bulk_v" option on and it clicks on something (In this case calls the __call_particle_pusher_bulk_v functions
                       "Particle_pusher_B_stream": draw_B_stream}

   particle_coordinates = []

   pipe = []

   particlepushercommand = []

   def __init__(self, vlsvReader, variable, pushcommand, operator="pass", threaded=True, **traits):
      ''' Initializes the class and loads the mayavi grid

          :param vlsvReader:        Some vlsv reader with a file open
          :type vlsvReader:         :class:`vlsvfile.VlsvReader`
          :param variable:          Name of the variable
          :param pushcommand:       Command for calling the particle pusher
          :param operator:          Operator for the variable
          :param threaded:          Boolean value for using threads or not using threads to draw the grid (threads enable interactive mode)
      '''
      self.particlepushercommand=pushcommand
      super(Particlepusherinterface, self).__init__(vlsvReader, variable, operator, threaded, **traits)








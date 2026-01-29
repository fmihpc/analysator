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

# This file contains static spacecraft function for fetching variable data from multiple vlsv files and looking at the data in time

import numpy as np
import logging
import analysator as pt
try:
   from collections.abc import Iterable
except ImportError:
   from collections import Iterable

def cell_time_evolution( vlsvReader_list, variables, cellids, units="" ):
   ''' Returns variable data from a time evolution of some certain cell ids

       :param vlsvReader_list:         List containing VlsvReaders with a file open
       :type vlsvReader_list:          :class:`vlsvfile.VlsvReader`
       :param variables:               Name of the variables
       :param cellids:                 List of cell ids
       :param units:                   List of units for the variables (OPTIONAL)
       :returns: an array containing the data for the time evolution for every cell id

       .. code-block:: python

          import analysator as pt
          # Example of usage:
          time_data = pt.calculations.cell_time_evolution( vlsvReader_list=[VlsvReader("bulk.000.vlsv"), VlsvReader("bulk.001.vlsv"), VlsvReader("bulk.002.vlsv")], variables=["rho", "Pressure", "B"], cellids=[2,4], units=["N", "Pascal", "T"] )

          # Check output
          logging.info time_data

          # Now plot the results:
          time = time_data[0]
          rho = time_data[3]
          pt.plot.plot_variables(time, rho)

          # Do post processing:
          rho_data = rho.data
          non_existing_example_function(rho_data)
   

   '''
   vlsvReader_list = np.atleast_1d(vlsvReader_list)
   variables = np.atleast_1d(variables)
   cellids = np.atleast_1d(cellids)
   reader_0 = vlsvReader_list[0]

   # Check against legacy files with tstep instead of timestep:
   if reader_0.check_parameter("tstep"):
      parameters = ["t","tstep","fileIndex"]
      parameter_units=["s","",""]
   elif reader_0.check_parameter("timestep"):
      parameters = ["t","timestep","fileIndex"]
      parameter_units=["s","",""]
   else:
      logging.warning("Could not obtain tstep or timestep from readers. Returning only t and fileIndex.")
      parameters = ["t","fileIndex"]
      parameter_units=["s",""]
   #construct empty units, if none are given
   if (units == "") or (len(units) != len(variables)):
      units=[ "" for i in range(len(variables))]
   #construct data 
   data = [[] for i in range(len(parameters)+len(cellids)*len(variables))]
   for t in range(len(vlsvReader_list)):
      # Get the vlsv reader
      vlsvReader = vlsvReader_list[t]
      # Open the vlsv reader's file:
      vlsvReader.optimize_open_file()
      #go through parameters
      for j in range(len(parameters)):
         # Read the parameter
         # Save the data into the right slot in the data array:
         data[j].append(vlsvReader.read_parameter( parameters[j]))

      # Go through variables:
      for j in range(len(variables)):
         variable = variables[j]
         # Read the variable for all cell ids
         #variables_for_cellids = vlsvReader.read_variables_for_cellids( variable, cellids )
         # Save the data into the right slot in the data array:
         for i in range(len(cellids)):
            data[len(parameters)+i*len(variables)+j].append(vlsvReader.read_variable( variable, cellids[i] ))
            #TODO Use vectorization over cellids!

      # For optimization purposes we are now freeing vlsvReader's memory
      # Note: Upon reading data vlsvReader created an internal hash map that takes a lot of memory
      vlsvReader.optimize_clear_fileindex_for_cellid()
      # Close the vlsv reader's file:
      vlsvReader.optimize_close_file()
   from output import output_1d
   return output_1d( data, 
                     parameters +  [variables[(int)(i)%(int)(len(variables))] for i in range(len(data)-len(parameters))], 
                     parameter_units + [units[(int)(i)%(int)(len(units))] for i in range(len(data)-len(parameters))] )


def point_time_evolution( vlsvReader_list, variables, coordinates, units="", method='nearest'):
   ''' Returns variable data from a time evolution of given coordinates. Defaults to zeroth-order spatial interpolation.

       :param vlsvReader_list:         List containing VlsvReaders with a file open
       :type vlsvReader_list:          :class:`vlsvfile.VlsvReader`
       :param variables:               Name of the variables
       :param coordinates:             List of coordinates [n,3]
       :param units:                   List of units for the variables (OPTIONAL)
       :param method:                  name of interpolation method ['nearest','linear']
       :returns: an array containing the data for the time evolution for every coordinate

       .. code-block:: python

          import pytools as pt
          # Example of usage:
          time_data = pt.calculations.point_time_evolution( vlsvReader_list=[VlsvReader("bulk.000.vlsv"), VlsvReader("bulk.001.vlsv"), VlsvReader("bulk.002.vlsv")], variables=["rho", "Pressure", "B"], coordinates=[[1e8,0,0],[1.2e8,0,0]], units=["N", "Pascal", "T"] )

          # Check output
          logging.info time_data

          # Now plot the results:
          time = time_data[0]
          rho = time_data[3]
          pt.plot.plot_variables(time, rho)

          # Do post processing:
          rho_data = rho.data
          non_existing_example_function(rho_data)
   

   '''
   vlsvReader_list = np.atleast_1d(vlsvReader_list)
   variables = np.atleast_1d(variables)
   coordinates = np.array(coordinates)
   if coordinates.ndim == 1:
      coordinates = coordinates[np.newaxis,:]
   reader_0 = vlsvReader_list[0]

   # Check against legacy files with tstep instead of timestep:
   if reader_0.check_parameter("tstep"):
      parameters = ["t","tstep","fileIndex"]
      parameter_units=["s","",""]
   elif reader_0.check_parameter("timestep"):
      parameters = ["t","timestep","fileIndex"]
      parameter_units=["s","",""]
   else:
      logging.warning("Could not obtain tstep or timestep from readers. Returning only t and fileIndex.")
      parameters = ["t","fileIndex"]
      parameter_units=["s",""]
   #construct empty units, if none are given
   if (units == "") or (len(units) != len(variables)):
      units=[ "" for i in range(len(variables))]
   #construct data 
   data = [[] for i in range(len(parameters)+coordinates.shape[0]*len(variables))]
   for t in range(len(vlsvReader_list)):
      # Get the vlsv reader
      vlsvReader = vlsvReader_list[t]
      # Open the vlsv reader's file:
      vlsvReader.optimize_open_file()
      #go through parameters
      for j in range(len(parameters)):
         # Read the parameter
         # Save the data into the right slot in the data array:
         data[j].append(vlsvReader.read_parameter( parameters[j]))

      # Go through variables:
      for j in range(len(variables)):
         variable = variables[j]
         # Read the variable for all cell ids
         #variables_for_cellids = vlsvReader.read_variables_for_cellids( variable, cellids )
         # Save the data into the right slot in the data array:
         for i in range(coordinates.shape[0]):
            data[len(parameters)+i*len(variables)+j].append(vlsvReader.read_interpolated_variable( variable, coordinates[i,:], method=method ))
            #TODO Use vectorization over coordinates!

      # For optimization purposes we are now freeing vlsvReader's memory
      # Note: Upon reading data vlsvReader created an internal hash map that takes a lot of memory
      vlsvReader.optimize_clear_fileindex_for_cellid()
      vlsvReader.optimize_clear_fileindex_for_cellid_blocks()
      # Close the vlsv reader's file:
      vlsvReader.optimize_close_file()
   from output import output_1d
   return output_1d( data, 
                     parameters +  [variables[(int)(i)%(int)(len(variables))] for i in range(len(data)-len(parameters))], 
                     parameter_units + [units[(int)(i)%(int)(len(units))] for i in range(len(data)-len(parameters))] )

class VlsvTInterpolator:
   ''' Class for setting up a time-interpolation wrapper from a list of VLSV files.
   Requires either files_list or vlsvReaders_list as input.

       :param files_list:               List containing paths to VLSV files
       :type files_list:                :class:`str`
       :param vlsvReaders_list:         List containing VlsvReaders with a file open
       :type vlsvReaders_list:          :class:`vlsvfile.VlsvReader`

       :returns:  a callable object that can be used to interpolate values in time and space from the given files.
       :returns: an array containing the data for the time evolution for every coordinate

       .. code-block:: python

          import analysator as pt;
          # Example of usage:
          time_interpolator = pt.calculations.VlsvTInterpolator( vlsvReader_list=[VlsvReader("bulk.000.vlsv"), VlsvReader("bulk.001.vlsv"), VlsvReader("bulk.002.vlsv")])


          # Now plot the results:
          time = time_interpolator.ts
          rho = [time_interpolator(t, [1e8,0,0], "proton/vg_rho")]
          pt.plot.plot_variables(time, rho)
          pl.show()

   '''

   def __init__(self, files_list=[], vlsvReaders_list=[]):
      self.ts = []
      
      self.files = files_list
      self.readers = vlsvReaders_list
      if len(self.files) == 0 and len(self.readers) == 0:
         raise ValueError("Please provide either a list of vlsv files or a list of readers")
      if len(self.readers) == 0:
         for i,f in enumerate(self.files):
            reader = pt.vlsvfile.VlsvReader(f)
            self.readers.append(reader)
            self.ts.append(reader.read_parameter('time'))
      else:
         self.files = []
         for i,r in enumerate(self.readers):
            self.ts.append(r.read_parameter('time'))
            self.files.append(r.file_name)

         # self.readers.append(pt.vlsvfile.VlsvReader(f))
      sort_bunch = [(t, self.files[i], self.readers[i]) for i,t in enumerate(self.ts)]
      # self.readers = np.array(self.readers)
      sort_bunch.sort()
      self.ts = [b[0] for b in sort_bunch]
      self.files = [b[1] for b in sort_bunch]
      self.readers = [b[2] for b in sort_bunch]

      self.activeReaders = {'low':None, 'hi':None}

   def __call__(self, t, coordinates, variable, operator='pass', method="linear", extrapolate=False):
      ''' Spatiotemporal interpolation. 

       :param t:                 Time-coordinate of query
       :param coordinates:       List or numpy array of query coordinates, shape of [n,3] or [3]
       :param variable:          Name of the variable
       :param operator:          Name of operator (see VlsvReader/datareduction) ['pass']
       :param method:            Name of interpolation method ['nearest','linear'], only affects spatial interpolation
       :param extrapolate:       Allow time queries outside of reader range [False]. If enabled, queries beyond range return interpolants from the first/latest reader.
       
       '''
      # print("Calling for t ",t)
      crds = np.array(coordinates)
      stack = True
      if len(np.shape(crds)) == 1:
         stack = False
         crds = np.atleast_2d(crds)


      if not extrapolate:
         if t < np.min(self.ts) or t > np.max(self.ts):
            raise ValueError("Queried time "+str(t)+" is out of bounds for given files ["+str(np.min(self.ts))+", "+str(np.max(self.ts))+"]")
      # crds = np.reshape(coordinates,(len(coordinates)//3,3))
      # print(coordinates)
      rti = np.searchsorted(self.ts, t)
      lower_t = self.ts[max(rti-1,0)]
      upper_t = self.ts[min(rti,len(self.ts)-1)]
      if(upper_t==lower_t):
         alpha = 0
      else:
         alpha = (t-lower_t)/(upper_t - lower_t)
      # print(lower_t, upper_t, coordinates)
      if self.activeReaders['low']:
         if lower_t == self.activeReaders['low'].read_parameter('time'):
               pass
         else:
               self.activeReaders['low'] = self.readers[max(rti-1,0)]
      else:
         self.activeReaders['low'] = self.readers[max(rti-1,0)]

      if self.activeReaders['hi']:
         if lower_t == self.activeReaders['hi'].read_parameter('time'):
               pass
         else:
               self.activeReaders['hi'] = self.readers[min(rti,len(self.ts)-1)]
      else:
         self.activeReaders['hi'] = self.readers[min(rti,len(self.ts)-1)]
         
      lower_vals = self.activeReaders['low'].read_interpolated_variable(variable, crds, operator=operator, method=method)
      upper_vals = self.activeReaders['hi'].read_interpolated_variable(variable, crds, operator=operator, method=method)
      vals = (1-alpha)*lower_vals + alpha*upper_vals


      test_variable = self.activeReaders['low'].read_variable(variable,cellids=[self.activeReaders['low'].get_cellid(crds[0,:])],operator=operator)
      if isinstance(test_variable,np.ma.core.MaskedConstant):
         value_length=1
      elif isinstance(test_variable, Iterable):
         value_length=len(test_variable)
      else:
         value_length=1

      if stack:
         return vals.squeeze()
      else:
         if value_length == 1:
            return vals.squeeze()[()] # The only special case to return a scalar instead of an array
         else:
            return vals.squeeze()

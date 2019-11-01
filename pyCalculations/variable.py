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

# This file has a class "Variable" that holds all the important data for variables e.g. variable's name, the units and the data on the variable
import numpy as np

class VariableInfo:
   ''' A class/struct for holding variable info. This includes the variable data in array form, the name and the units.
       Variable info is in: data, name, units.

       .. code-block:: python

          LIST VARIABLES IN VARIBLE INFO:
          data              Data of the variable (array list)
          name              Name of the variable
          units             Units of the variable
          latex             Name of the variable in LaTeX
          latexunits        Units of the variable in LaTeX
   '''
   def __init__(self, data_array, name="", units="", latex="", latexunits=""):
      ''' Initializes variable info.

          :param data_array:         Raw data for the variable (an array for example)
          :param name:               Name of the variable
          :param units:              Name of the variable's units
          :param latex:              Name of the variable in LaTeX
          :param latexunits:         Name of the variable's units in LaTeX
      '''
      self.data = np.ma.asarray(data_array)
      self.name = name
      self.units = units
      self.latex = latex
      self.latexunits = latexunits

   def __repr__(self):
      return "VariableInfo(Name: \'" + str(self.name) + "\' Units: \'" + str(self.units) + "\')"

   def __str__(self):
      return self.__repr__()

   def get_variable(self, index ):
      ''' Creates a new variable with identical data but only certain index is included

          :param index:         Vector index
          :returns: an edited version of the variable

          .. note::

             E.g. if the data is an array of 3d vectors, get_variable(0) would return the variable with data[:,0] as the data
      '''
      import numpy as np
      if len(self.data) <= 0:
         print("BAD DATA LENGTH")
         return []
      if len(np.atleast_1d(self.data[0])) <= index:
         print("BAD INDEX, THE INDEX IS LARGER THAN VECTOR SIZE!")
         return []
      return VariableInfo( self.data[:,index], self.name, self.units, self.latex, self.latexunits )



def get_data( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance

       :param variable:           The variable as a VariableInfo instance or a list
       :returns: data of the variable
   '''
   if isinstance(variable, VariableInfo):
      return variable.data
   else:
      return variable

def get_name( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance

       :param variable:           The variable as a VariableInfo instance or a list
       :returns: the name of the variable or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.name
   else:
      return ""

def get_units( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance

       :param variable:           The variable as a VariableInfo instance or a list
       :returns: the units of the variable or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.units
   else:
      return ""

def get_latex( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance

       :param variable:           The variable as a VariableInfo instance or a list
       :returns: the variable in LaTeX format or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.latex
   else:
      return ""

def get_latexunits( variable ):
   ''' Function to use when not sure if variable is in raw form ( simply a list with data ), or a VariableInfo instance

       :param variable:           The variable as a VariableInfo instance or a list
       :returns: the units of the variable in LaTeX format or \"\" if not a VariableInfo instance
   '''
   if isinstance(variable, VariableInfo):
      return variable.latexunits
   else:
      return ""


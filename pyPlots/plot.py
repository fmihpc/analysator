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

''' The plot module has all the functions related to plotting variables.

   .. code-block:: python

      # Example:
      import pytools as pt
      pt.pt.
      #press [tab] -> get the functions

'''

from plot_variables import plot_variables, plot_multiple_variables
from plot2d_with_vspace import vlsv_plot2d_with_vspace

import colormaps
import plot_helpers
from plot_colormap import plot_colormap
from plot_vdf import plot_vdf
from plot_mxm import plot_colormapmxm
from plot_vdf import vSpaceReducer
from plot_colormap_with_VDF import plot_colormap_with_vdf
#from plot_colormap_with_spacecraft import plot_colormap_with_spacecraft
from plot_keogram import plot_keogram
from plot_precipitation import plot_prec_spectrum
from plot_precipitation import plot_prec_time_spectrum

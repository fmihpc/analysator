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

import colormaps
import plot_helpers
from plot_colormap import plot_colormap
from plot_vdf import plot_vdf
from plot_vdf_profiles import plot_vdf_profiles

from plot_colormap3dslice import plot_colormap3dslice

import numpy as np, os

decimalprecision_ax = 0
cb_linear = False

# Different style scientific format for colour bar ticks
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    # this should bring all colorbar ticks to the same horizontal position, but for
    # some reason it doesn't work. (signchar=r'\enspace')
    signchar=r'' 
    # Multiple braces for b take care of negative values in exponent
    # brackets around \times remove extra whitespace
    if os.getenv('PTNOLATEX') is None:
        # replaces minus sign with en-dash to fix big with latex descender value return
        if np.sign(x)<0: signchar=r'\mbox{\textbf{--}}'
        return r'$'+signchar+'{}'.format(abs(float(a)))+r'{\times}'+'10^{{{}}}$'.format(int(b))
    else:
        return r'$'+'{}'.format(float(a))+r'{\times}'+'10^{{{}}}$'.format(int(b))

# axisfmt replaces minus sign with en-dash to fix big with latex descender value return
def axisfmt(x, pos):
    f = r'{:.'+decimalprecision_ax+r'f}'
    if os.getenv('PTNOLATEX') is None:
        a = f.format(abs(x))
        if np.sign(x)<0: a = r'\mbox{\textbf{--}}'+a
        return r'$'+a+'$'
    else:
        return f.format(x)

# cbfmt replaces minus sign with en-dash to fix big with latex descender value return, used for colorbar
def cbfmt(x, pos):
    # Set required decimal precision
    a, b = '{:.1e}'.format(x).split('e')
    precision = '0'
    if (cb_linear is True):
        # for linear, use more precision
        if int(b)<1: precision = str(abs(-1+int(b)))
    else:
        if int(b)<1: precision = str(abs(int(b)))
    f = r'{:.'+precision+r'f}'
    if os.getenv('PTNOLATEX') is None:
        a = f.format(abs(x))
        if np.sign(x)<0: a = r'\mbox{\textbf{--}}'+a
        return r'$'+a+'$'
    else:
        return f.format(x)


# 
# This file is part of Analysator.
# Copyright 2010-2016 Finnish Meteorological Institute
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

# Default vmin and vmax values defined per variable and run
def loadrundefaults(run, var, op):
    vminuse=None
    vmaxuse=None

    if var == 'rho':
        vminuse = 1.e+6
        vmaxuse = 4.e+6
        if run == 'AEA' or run == 'ABA' or run == 'BCQ':
            vminuse = 8.e+5
            vmaxuse = 4.e+6
        elif run == 'AEC' or run == 'ABC':
            vminuse = 1.e+6
            vmaxuse = 4.e+6
        elif run == 'BEB':
            vminuse = 2.e+6
            vmaxuse = 2.e+7

    elif var == 'MA':
        vminuse = 1.0
        vmaxuse = 10.0

    elif var == 'Mms':
        vminuse = 1.0
        vmaxuse = 10.0

    elif var == 'B':        
        if op==None:
            vminuse = 3.e-9
            vmaxuse = 3.e-7
        else: #components
            vminuse = -0.15e-9
            vmaxuse = 0.15e-9

    elif var == 'E':
        if op==None:
            vminuse = 0.01
            vmaxuse = 0.1
        else: #components
            vminuse = -0.015
            vmaxuse = 0.015

    elif var == 'rhoBeam':
        vminuse = 1.e3
        vmaxuse = 1.e5

    elif var == 'V':
        if op==None:
            vminuse = 1.e5
            vmaxuse = 1.e6
        else: #components
            vminuse = -1.e6
            vmaxuse = 1.e6

    elif var == 'beta':
        vminuse = 0.8
        vmaxuse = 100

    elif var == 'temperature':
        vminuse = 0.5e6
        vmaxuse = 5.0e6

    elif var == 'vs':
        vminuse = 5e4
        vmaxuse = 5e5

    elif var == 'va':
        vminuse = 5e4
        vmaxuse = 5e5

    return(vminuse, vmaxuse)


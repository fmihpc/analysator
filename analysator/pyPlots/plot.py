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
      import analysator as pt
      pt.pt.
      #press [tab] -> get the functions

'''

from plot_variables import plot_variables, plot_multiple_variables

import analysator as pt

import logging
import matplotlib.pyplot as plt
import matplotlib
import colormaps as cmaps

import plot_helpers
from plot_colormap import plot_colormap
from plot_vdf import plot_vdf
from plot_vdfdiff import plot_vdfdiff
from plot_vdf_profiles import plot_vdf_profiles
from plot_colormap3dslice import plot_colormap3dslice
from plot_threeslice import plot_threeslice
from plot_ionosphere import plot_ionosphere
import platform
from packaging.version import Version
from numbers import Number

try:
    from plot_isosurface import plot_isosurface, plot_neutral_sheet
    if Version(platform.python_version()) < Version("3.7"):
       raise ImportError("Python>=3.7 required for RBFInterpolator.")
except:
    logging.warning("plot_isosurface not imported. To access it, use Python version >=3.7 and install scikit-image.")


import numpy as np, os

if Version(matplotlib.__version__) < Version("3.3.0"):
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    plt.register_cmap(name='viridis_r', cmap=cmaps.viridis_r)
    plt.register_cmap(name='plasma', cmap=cmaps.plasma)
    plt.register_cmap(name='plasma_r', cmap=cmaps.plasma_r)
    plt.register_cmap(name='inferno', cmap=cmaps.inferno)
    plt.register_cmap(name='inferno_r', cmap=cmaps.inferno_r)
    plt.register_cmap(name='magma', cmap=cmaps.magma)
    plt.register_cmap(name='magma_r', cmap=cmaps.magma_r)

# Register custom colourmaps
if Version(matplotlib.__version__) < Version("3.5.0"):
    plt.register_cmap(name='parula', cmap=cmaps.parula)
    plt.register_cmap(name='parula_r', cmap=cmaps.parula_r)
    plt.register_cmap(name='hot_desaturated', cmap=cmaps.hot_desaturated_colormap)
    plt.register_cmap(name='hot_desaturated_r', cmap=cmaps.hot_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step
    plt.register_cmap(name='pale_desaturated', cmap=cmaps.pale_desaturated_colormap)
    plt.register_cmap(name='pale_desaturated_r', cmap=cmaps.pale_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step
    plt.register_cmap(name='warhol', cmap=cmaps.warhol_colormap)
else:
    matplotlib.colormaps.register(name='parula', cmap=cmaps.parula)
    matplotlib.colormaps.register(name='parula_r', cmap=cmaps.parula_r)
    matplotlib.colormaps.register(name='hot_desaturated', cmap=cmaps.hot_desaturated_colormap)
    matplotlib.colormaps.register(name='hot_desaturated_r', cmap=cmaps.hot_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step
    matplotlib.colormaps.register(name='pale_desaturated', cmap=cmaps.pale_desaturated_colormap)
    matplotlib.colormaps.register(name='pale_desaturated_r', cmap=cmaps.pale_desaturated_colormap_r) # Listed colormap requires making reversed version at earlier step
    matplotlib.colormaps.register(name='warhol', cmap=cmaps.warhol_colormap)


decimalprecision_ax = 0
decimalprecision_cblin = 0
cb_linear = False


if matplotlib.__version__=="0.99.1.1" and np.__version__=="1.4.1":
   logging.info('Warning, according to loaded numpy and matplotlib versions, user appears to be')
   logging.info('either using csc.taito.fi without loading the mayavi2 module, or by invoking')
   logging.info('the system python interpeter by calling "./scriptname.py" instead of "python ./scriptname.py"')

# Run TeX typesetting through the full TeX engine instead of python's own mathtext. Allows
# for changing fonts, bold math symbols etc, but may cause trouble on some systems.
if not os.getenv('PTNOLATEX'):
   matplotlib.rc('text', usetex=True)
   matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
   # matplotlib.rcParams['mathtext.fontset'] = 'stix'
   # matplotlib.rcParams['font.family'] = 'STIXGeneral'
   # Matplotlib suppressed logging messages came out after enabling logging.INFO: font.family must be one of (serif, sans-serif, cursive, monospace) when text.usetex is True. serif will be used by default.
   matplotlib.rcParams['font.family'] = 'serif'
   logging.info("Using LaTeX formatting")
   # matplotlib.rcParams['text.dvipnghack'] = 'True' # This hack might fix it on some systems

# Set backends
if matplotlib.get_backend()[:9] == 'module://':
   logging.info("Using backend "+matplotlib.get_backend())
   backend_interactive = matplotlib.get_backend()
   backend_noninteractive = matplotlib.get_backend()
elif not os.getenv('PTBACKEND'):
   backend_interactive = 'TkAgg'
   backend_noninteractive = 'Agg'
else:
   backend_interactive = os.getenv('PTBACKEND')
   backend_noninteractive = os.getenv('PTBACKEND')






if os.getenv('PTNONINTERACTIVE') != None:
   # Non-interactive plotting mode
   try:
      plt.switch_backend(backend_noninteractive)
   except:
      logging.info("Note: Unable to switch to "+backend_noninteractive+" backend")
else:
   # Interactive plotting mode
   plt.ion()
   try:
      plt.switch_backend(backend_interactive)
   except:
      logging.info("Note: Unable to switch to "+backend_interactive+" backend")

pt.backend_interactive=backend_interactive
pt.backend_noninteractive=backend_noninteractive


# Output matplotlib version
logging.info("Using matplotlib version "+matplotlib.__version__)

def get_cmap(colormap):
    if Version(matplotlib.__version__) <= Version("3.6"):
        return matplotlib.cm.get_cmap(name=colormap)
    else:
        return matplotlib.colormaps.get_cmap(colormap)


# Default output directory for plots
defaultoutputdir=os.path.expandvars('$HOME/Plots/')
if os.getenv('PTOUTPUTDIR'):
    defaultoutputdir=os.getenv('PTOUTPUTDIR')

# axisfmt replaces minus sign with en-dash to fix bug with latex descender value return
# nb: axis ticks are never plotted with scientific format
def axisfmt(x, pos):
    f = r'{:.'+decimalprecision_ax+r'f}'
    if not os.getenv('PTNOLATEX'):
        a = f.format(abs(x))
        if np.sign(x)<0: a = r'\mbox{\textbf{--}}'+a
        return r'$'+a+'$'
    else:
        return f.format(x)

# cbfmtsci replaces minus sign with en-dash to fix bug with latex descender value return
# Scientific format for colour bar ticks
def cbfmtsci(x, pos):
    if (cb_linear is True):
        # for linear, use more precision
        a, b = ('{:.'+str(int(decimalprecision_cblin))+'e}').format(x).split('e')
        precisionvalue = int(decimalprecision_cblin)
        if int(b) < 0:
            precisionvalue += int(b)
        if abs(precisionvalue)<1:
            precisionvalue=1
        f = '{:.' + str(abs(precisionvalue)) + 'f}'
        number = f.format(abs(float(a)))+r'{\times}'+'10^{{{}}}'.format(int(b))
    else:
        a, b = '{:.1e}'.format(x).split('e')
        number = '{:.1f}'.format(abs(float(a)))+r'{\times}'+'10^{{{}}}'.format(int(b))
    signchar=r''
    # Multiple braces for b take care of negative values in exponent
    # brackets around \times remove extra whitespace
    if not os.getenv('PTNOLATEX'):
        # replaces minus sign with en-dash to fix big with latex descender value return
        if np.sign(x)<0: signchar=r'\mbox{\textbf{--}}'
    else:
        if np.sign(x)<0: signchar=r'-'
    # Final special treatment for zero value
    if x==0:
        number = r'0.0{\times}10^{{{0}}}'
    return r'$'+signchar+number+'$'

# cbfmt replaces minus sign with en-dash to fix bug with latex descender value return, used for colorbar
# nb: regular floating i.e. non-scientific format for colorbar ticks
def cbfmt(x, pos):
    # Special treatment for zero value
    if x==0:
        return r'$0.0$'
    # Set required decimal precision
    a, b = '{:.1e}'.format(x).split('e')
    # e.g. 9.0e-1 means we need precision 1
    if (cb_linear is True):
        # for linear, use more precision
        precision = str(int(decimalprecision_cblin))
        #if int(b)<1: precision = str(1+abs(int(b)))
    else:
        precision = '0'
        if int(b)<1: precision = str(abs(int(b)))
    f = r'{:.'+precision+'f}'
    if not os.getenv('PTNOLATEX'):
        a = f.format(abs(x))
        if np.sign(x)<0: a = r'\mbox{\textbf{--}}'+a
        return r'$'+a+'$'
    else:
        return f.format(x)


# Helper routines for latex output handling
def bfstring(string):
    if not os.getenv('PTNOLATEX'):
        if len(string)==0:
            return ''
        else:
            return r'\mathbf{'+string+'}'
    # LaTeX output off
    return string

def rmstring(string):
    if len(string)==0:
        return ''
    else:
        return r'\mathrm{'+string+'}'

def mathmode(string):
    if len(string)==0:
        return ''
    else:
        # First remove any internal possible dollar signs, then wrap whole string into math block
        result = string.replace('$','')
        if os.getenv('PTNOLATEX'):
            # Get rid of latex spaces
            result = result.replace(r'\,','~').replace(r'\qquad','~~~~~~')
        return r"$"+result+"$"

def textbfstring(string):
    if not os.getenv('PTNOLATEX'):
        if len(string)==0:
            return ''
        else:
            return r'\textbf{'+string+'}'
    # LaTex output off
    return string

#Draws sharp contour lines along cell edges (xmeshpass,ymeshpass), for some dpi+linewidth combinations antialiasing True causes issues
def cell_edgecontours(ax,XmeshPass,YmeshPass,heightmap,threshold=0,linewidth=0.5,linestyle='solid',colors='black',antialiased=False):
    x=np.array(XmeshPass[0])
    y=np.array([li[0] for li in YmeshPass])

    v = np.diff(heightmap > threshold, axis=1)
    h = np.diff(heightmap > threshold, axis=0)

    # From https://stackoverflow.com/questions/63458863/way-to-contour-outer-edge-of-selected-grid-region-in-python
    # Check that at least one spot exceeds threshold or the below will crash.
    if np.max(heightmap) > threshold:
        
        l = np.argwhere(v.T)
        vlines = np.array(list(zip(np.stack((x[l[:, 0]+1], y[l[:, 1]])).T,
                                    np.stack((x[l[:, 0] + 1], y[l[:, 1] + 1])).T)))
        l = np.argwhere(h.T)
        hlines = np.array(list(zip(np.stack((x[l[:, 0]], y[l[:, 1] + 1])).T,
                                    np.stack((x[l[:, 0] + 1], y[l[:, 1] + 1])).T)))
        lines = np.vstack((vlines, hlines))

        ax.add_collection(matplotlib.collections.LineCollection(lines, lw=linewidth, colors=colors, linestyle=linestyle,zorder=2,antialiased=antialiased))

    
    return 0
def output_path(outputfile,outputfile_default,outputdir,nooverwrite):
        check=False
        if not outputfile:
            if not outputfile_default:
                outputfile="plot.png"

            outputfile=outputfile_default
            if not outputdir:
                # Check later if outputfile contains path information and combine with defaultoutputdir
                check=True

        # Separate possible path information from outputfile
        outputprefixind = outputfile.rfind('/')
        if outputprefixind >= 0:            
            outputdir = outputfile[:outputprefixind+1]
            outputfile = outputfile[outputprefixind+1:]

        if not outputdir: # default initial path
            outputdir=defaultoutputdir
            check=False
    
        if check:
            #remove leading '/' on outputdir as this would cause issues with os.path.join
            outputprefixind = outputdir.find('/')
            if outputprefixind == 0:            
                outputdir = outputdir[1:]
            outputdir= os.path.join(defaultoutputdir,outputdir)
            outputfile=os.path.join(outputdir,outputfile)
        else:
            outputfile = os.path.join(outputdir,outputfile)

        # Ensure output directory exists
        if not os.path.exists(outputdir):
            try:
                os.makedirs(outputdir)
            except FileExistsError: 
                #Parallel jobs might try to create dir simultaneously.
                pass
            except Exception as e:
                raise IOError("Could not create output directory "+outputdir+": "+str(e))
        if not os.access(outputdir, os.W_OK):
            raise IOError("No write access for directory "+outputdir+" Exiting.")

        # Check if target file already exists and overwriting is disabled
        if (nooverwrite and os.path.exists(outputfile)):            
            if os.stat(outputfile).st_size > 0: # Also check that file is not empty
                raise IOError("Found existing file "+outputfile+" and nooverwrite=True.")
            else:
                logging.warning("Found existing file "+outputfile+" of size zero. Re-rendering.")

        return outputfile

def get_scaled_units(variable_info,vscale=None, env='EarthSpace', manualDict=None):
    ''' Return scaling metadata

        :param variable_info: VariableInfo object used for getting the units
        :param env:          A string to choose the scaling dictionary [default: EarthSpace]
        :param manualDict:   a dictionary of {units : {scalingparams}}; used to update the included dictionary
        :param vscale:       float, factor to scale the variable with
        
        :returns: (norming factor, scaledUnits, scaledLatexUnits)

    '''

    #
    if env=='EarthSpace':
        variable_info.scaleDict = {
                's'      : {'defaultScale':1,
                            1e6: {'scaledUnits':'us', 'scaledLatexUnit':r'$\mu\mathrm{s}$'},
                            1e3: {'scaledUnits':'ms', 'scaledLatexUnit':r'$\mathrm{ms}$'}
                        },
                'T'      : {'defaultScale':1e9,
                            1e9: {'scaledUnits':'nT', 'scaledLatexUnit':r'$\mathrm{nT}$'}
                        },
                'K'      : {'defaultScale':1e-6,
                            1e-6:{'scaledUnits':'MK', 'scaledLatexUnit':r'$\mathrm{MK}$'}
                        },
                'Pa'     : {'defaultScale':1e9,
                            1e9:{'scaledUnits':'nPa', 'scaledLatexUnit':r'$\mathrm{nPa}$'}
                        },
                '1/m'  :   {'defaultScale':1e3,
                            1e3:{'scaledUnits':'1/km', 'scaledLatexUnit':r'$\mathrm{km}^{-1}$'}
                        },
                '1/m^3'  : {'defaultScale':1e-6,
                            1e-6:{'scaledUnits':'1/cm^3', 'scaledLatexUnit':r'$\mathrm{cm}^{-3}$'}
                        },
                '1/m3'   : {'defaultScale':1e-6,
                            1e-6:{'scaledUnits':'1/cm^3', 'scaledLatexUnit':r'$\mathrm{cm}^{-3}$'}
                        },
                'm/s'    : {'defaultScale':1e-3,
                            1e-3:{'scaledUnits':'km/s', 'scaledLatexUnit':r'$\mathrm{km}\,\mathrm{s}^{-1}$'}
                        },
                'V/m'    : {'defaultScale':1e3,
                            1e3:{'scaledUnits':'mV/m', 'scaledLatexUnit':r'$\mathrm{mV}\,\mathrm{m}^{-1}$'}
                        },
                'eV/cm3' : {'defaultScale':1e-3,
                            1e-3:{'scaledUnits':'keV/cm^3', 'scaledLatexUnit':r'$\mathrm{keV}\,\mathrm{cm}^{-3}$'}
                        },
                'eV/cm^3': {'defaultScale':1e-3,
                            1e-3:{'scaledUnits':'keV/cm^3', 'scaledLatexUnit':r'$\mathrm{keV}\,\mathrm{cm}^{-3}$'}
                        },
                'T/m': {'defaultScale':1e-12,
                            1e-9:{'scaledUnits':'nT/m', 'scaledLatexUnit':r'$\mathrm{nT}\,\mathrm{m}^{-1}$'},
                            1e-12:{'scaledUnits':'nT/km', 'scaledLatexUnit':r'$\mathrm{nT}\,\mathrm{km}^{-1}$'}
                        },
                'kg/m3'  : {'defaultScale':5.97863741e26,
                            5.97863741e26:{'scaledUnits':'amu/m^3', 'scaledLatexUnit':r'$\mathrm{amu}\,\mathrm{m}^{-3}$'},
                            5.97863741e20:{'scaledUnits':'amu/cm^3', 'scaledLatexUnit':r'$\mathrm{amu}\,\mathrm{cm}^{-3}$'}
                        },
                'A/m^2' :   {'defaultScale':1,
                            1e-3:{'scaledUnits':'mA/m^2', 'scaledLatexUnit':r'$\mathrm{mA}\,\mathrm{m}^{-2}$'}
                        },
                'W/m^2' :   {'defaultScale':1,
                            1e-3:{'scaledUnits':'mW/m^2', 'scaledLatexUnit':r'$\mathrm{mW}\,\mathrm{m}^{-2}$'}
                        },
                'm^2' :   {'defaultScale':1,
                            1e-4:{'scaledUnits':'km^2', 'scaledLatexUnit':r'$\mathrm{km}^{2}$'}
                        },
                'V' :   {'defaultScale':1,
                            1e-3:{'scaledUnits':'mV', 'scaledLatexUnit':r'$\mathrm{mV}$'}
                        },
                'Degrees' :   {'defaultScale':1,
                        },
                'mho' :   {'defaultScale':1,
                        },
        }

    else:
        variable_info.scaleDict = {}
    if manualDict is not None:
        variable_info.scaleDict.update(manualDict)
    unitScale = 1.0
    scaledUnits = variable_info.units
    scaledLatexUnits = variable_info.latexunits

    if variable_info.units != '':
        dictKey = variable_info.units
        try:
            udict = variable_info.scaleDict[dictKey]
        except:
            if vscale is None:
                return 1.0, variable_info.units, variable_info.latexunits
            else:
                return vscale, variable_info.units, variable_info.latexunits
        if vscale is None:
            try:
                unitScale = udict['defaultScale']
            except:
                return 1.0, variable_info.units, variable_info.latexunits
        elif np.isclose(vscale, 1.0):
            return 1.0, variable_info.units, variable_info.latexunits
        else:
            unitScale = vscale

        if not any([np.isclose(unitScale, tryScale) for tryScale in udict.keys() if isinstance(tryScale, Number)]):
            #
            return vscale, variable_info.units+" x{vscale:e}".format(vscale=vscale), variable_info.latexunits+r"{\times}"+cbfmtsci(vscale,None)
        try:
        #above guarantees the list comprehension does not give an empty list
            unitScale = [scale for scale in udict.keys() if isinstance(scale, Number) and np.isclose(scale,unitScale)][0]
            scaledUnits = udict[unitScale]['scaledUnits']
        except KeyError:
            # logging.info('Missing scaledUnits in specialist dict for' + variable_info.units + ' for unitScale='+str(unitScale))
            return 1.0, variable_info.units, variable_info.latexunits
        try:
            scaledLatexUnits = udict[unitScale]['scaledLatexUnit']
        except:
            # logging.info('Missing scaledLatexUnits in specialist dict for ' + variable_info.units+ ' for unitScale='+str(unitScale))
            return 1.0, variable_info.units, variable_info.latexunits
    else:
        if vscale is None or np.isclose(vscale, 1.0):
            return 1.0, variable_info.units, variable_info.latexunits
        else:
            return vscale, variable_info.units+"x{vscale:e}".format(vscale=vscale), variable_info.latexunits+r"{\times}"+cbfmtsci(vscale,None)

    return unitScale, scaledUnits, scaledLatexUnits

# A utility to get variableinfo with corresponding units for simple plotting. Add "canonical" scalings as
# necessary, for default/other environments.
def get_scaled_var(variable_info,vscale=None, data=None, env='EarthSpace', manualDict=None):
    ''' Automatically scales the variableinfo data and adjusts the units correspondingly with the
        default dictionaries.

        :param variable_info: VariableInfo object used for getting the units
        :param data:         in case you wish to provide new data array (why, though?)
        :param env:          A string to choose the scaling dictionary [default: EarthSpace]
        :param manualDict:   a dictionary of {units : {scalingparams}}; used to update the included dictionary

        :returns: variable_info, with scaled units with pre-formatted units included in the varinfo.


    '''

    if data is None:
        data = variable_info.data
    else:
        variable_info.data = data

    unitScale, scaledUnits, scaledLatexUnits = get_scaled_units(vscale=vscale, env=env, manualDict=manualDict,variable_info=variable_info)
    if unitScale == 1: # no change, optimize out the calculation
        return variable_info

    variable_info.data = data*unitScale
    variable_info.units = scaledUnits
    variable_info.latexunits = scaledLatexUnits
    return variable_info


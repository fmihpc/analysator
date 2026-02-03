
import socket, re, os, tempfile, atexit, shutil, sys
import warnings
import logging
import importlib.util

#Initialize backend variables, pt.plot will assign these when imported
backend_noninteractive=""
backend_interactive=""


logging.basicConfig(format='%(levelname)s: %(message)s', level=os.environ.get('ANALYSATOR_LOG_LEVEL', 'INFO').upper())

# Input current folder's path
#sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Input folder paths

if os.getenv('PTMAYAVI2') != None:
   sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/" + "pyMayaVi")


# Make sure matplotlib has a unique temp directory
mpldir = tempfile.mkdtemp()
atexit.register(shutil.rmtree, mpldir)
os.environ['MPLCONFIGDIR']=mpldir

# Check if user is on taito.csc.fi without loading the mayavi2 module
import numpy as np
logging.getLogger('matplotlib').setLevel(os.environ.get('ANALYSATOR_MPL_LOG_LEVEL', 'WARNING').upper())


def lazyimport(module_name):
   try:
      return sys.modules[module_name]
   except KeyError:
      spec = importlib.util.find_spec(module_name)
      loader = importlib.util.LazyLoader(spec.loader)
      spec.loader = loader
      module = importlib.util.module_from_spec(spec)
      sys.modules[module_name]= module
      loader.exec_module(module)
      return module



# Import modules
try:
   calculations=lazyimport("analysator.calculations")
except ImportError as e:
   logging.info("Note: Did not import calculations module: " + str(e))

try:
   vlsvfile=lazyimport("analysator.vlsvfile")
except ImportError as e:
   logging.info("Note: Did not import vlsvfile module: " + str(e))

import os

if os.getenv('PTNONINTERACTIVE') == None: #was ineq

   #Only attempt loading MayaVi2 if requested

   if os.getenv('PTMAYAVI2') != None:
      try:
         import grid
      except ImportError as e:
         logging.info("Note: Did not import (outdated MayaVi2) grid module: " + str(e))

try:
   #import plot
   plot=lazyimport("analysator.plot")
except ImportError as e:
   logging.info("Note: Did not import plot module: " + str(e))

def register_colormaps():
    ''' Register included colormaps to matplotlib. Required if you are relying on colormaps included in analysator and you are not using analysator.plot functions.
    '''
    import analysator.plot.colormaps

try:
   #import miscellaneous
   miscellaneous=lazyimport("analysator.miscellaneous")
except ImportError as e:
   logging.info("Note: Did not import miscellaneous: " + str(e))

# from analysator.pytools import *
# import warnings
# warnings.filterwarnings("once", category=DeprecationWarning)
# warnings.filterwarnings("once", category=PendingDeprecationWarning)
# warnings.filterwarnings("once", category=FutureWarning)


# # from pytools import vslvfile

# from os import path as __path
# root = __path.dirname(__file__)
# with open(__path.join(root,'pytools.py'),'r') as f:
#     source = f.read()
#     exec(source)
#     f.close()
#     del f

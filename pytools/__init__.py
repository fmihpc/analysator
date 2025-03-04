
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

import os
import warnings
import sys
import logging
import traceback

logging.basicConfig(format='%(levelname)s: %(message)s', level=os.environ.get('ANALYSATOR_LOG_LEVEL', 'INFO').upper())
if logging.root.level <= logging.DEBUG:
   print("`import pytools` called from the script below - `import analysator` works as well!")
   traceback.print_stack()
else:
   logging.info("You can now also use `import analysator` instead of `import pytools`! Analysator has been renamed to... analysator(!) for consistency.\n\n`import pytools` and `import pytools as pt` will work via the dirty hack here in `pytools/__init__.py` until some time in the future (v.1 release?).\n Use `export ANALYSATOR_LOG_LEVEL=DEBUG` to find out where `import pytools` was called.")

# Slurp and exec the analysator package here to get all functionalities under the pytools module if needed
if 'analysator' not in sys.modules.keys():
   logging.info("Importing analysator to pytools")
   root = os.path.dirname(__file__)
   sys.path.append(os.path.join(root,"src/"))
   from analysator import *
else:
   logging.info("Analysator imported already. Pointing pytools to analysator.")
   sys.modules['pytools'] = sys.modules['analysator']

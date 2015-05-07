import struct
import xml.etree.ElementTree as ET
import ast
import numpy as np
import os
from reduction import datareducers,data_operators
from collections import Iterable
from vlsvwriter import VlsvWriter
from variable import get_data
from operator import itemgetter
from vlsvreader import VlsvReader

class VlasiatorReader(VlsvReader):
   ''' Class for reading VLSV files with support for Vlasiator velocity space and structures

   '''
   def get_nearest_cellid_with_distribution( self, cellid ):
      # Find the nearest cell id with distribution:
      # Read cell ids with velocity distribution in:
      cell_candidates = self.read("SpatialGrid","CELLSWITHBLOCKS")
      # Read in the coordinates of the cells:
      cell_candidate_coordinates = [self.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]
      # Read in the cell's coordinates:
      cell_coordinates = self.get_cell_coordinates(cellid)
      # Find the nearest:
      norms = np.sum((cell_candidate_coordinates - cell_coordinates)**2, axis=-1)**(1./2)
      norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))
      # Get the cell id:
      cellid = cell_candidates[i]
      return cellid

   def get_nearest_coordinates_with_distribution( self, coordinates ):
      # Find the nearest cell id with distribution:
      # Read cell ids with velocity distribution in:
      cell_candidates = self.read("SpatialGrid","CELLSWITHBLOCKS")
      # Read in the coordinates of the cells:
      cell_candidate_coordinates = [self.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]
      # Read in the cell's coordinates:
      cell_coordinates = coordinates
      # Find the nearest:
      norms = np.sum((cell_candidate_coordinates - cell_coordinates)**2, axis=-1)**(1./2)
      norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))
      # Get the cell id:
      cellid = cell_candidates[i]
      return self.get_cell_coordinates(cellid)


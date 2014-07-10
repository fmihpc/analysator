import numpy as np

def extract_velocity_cells_sphere( vlsvReader, cellid, origin, radius ):
   ''' Extracts the velocity cells inside some given sphere
       :param vlsvReader:         Some VlsvFile with a file open
       :param cellid:             The cellid whose rho to calculate
       :param origin:             Origin for the sphere
       :param radius:             Radius for the sphere
       :returns: Non backstream velocity cells and their avgs values as [vcellids, avgs]
   '''
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   # Get cells:
   vcellids = velocity_cell_data.keys()
   # Get avgs data:
   avgs = velocity_cell_data.values()
   # Get a list of velocity coordinates shifted by the solar wind bulk velocity:
   origin = np.array(origin)
   v = vlsvReader.get_velocity_cell_coordinates(vcellids) - origin
   # Get sum of radiuses:
   radiuses = np.sum(np.abs(v)**2,axis=-1)
   # Check radius condition
   radius2 = radius**2
   condition = (radiuses <= radius2)
   # Get the velocity cells of sphere
   vcellids_sphere = np.extract(condition, vcellids)
   # Get the avgs
   avgs_sphere = np.extract(condition, avgs)
   # Return
   return [vcellids_sphere, avgs_sphere]

def extract_velocity_cells_non_sphere( vlsvReader, cellid, origin, radius ):
   ''' Retrieves the velocity cells within a given sphere and returns the population outside the given sphere
       :param vlsvReader:         Some VlsvFile with a file open
       :param cellid:             The cellid whose rho to calculate
       :param origin:             Origin for the sphere
       :param radius:             Radius for the sphere
       :returns: Non backstream velocity cells and their avgs values as [vcellids, avgs]
   '''
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   # Get cells:
   vcellids = velocity_cell_data.keys()
   # Get avgs data:
   avgs = velocity_cell_data.values()
   # Get a list of velocity coordinates shifted by the solar wind bulk velocity:
   origin = np.array(origin)
   v = vlsvReader.get_velocity_cell_coordinates(vcellids) - origin
   # Get sum of radiuses:
   radiuses = np.sum(np.abs(v)**2,axis=-1)
   # Check radius condition
   radius2 = radius**2
   condition = (radiuses > radius2)
   # Get the velocity cells of nonsphere
   vcellids_nonsphere = np.extract(condition, vcellids)
   # Get the avgs
   avgs_nonsphere = np.extract(condition, avgs)
   # Return
   return [vcellids_nonsphere, avgs_nonsphere]

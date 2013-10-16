import numpy as np


def get_rho_nonbackstream( vlsvReader, cellid ):
   ''' Calculates the non backstream population's (solar wind's) contributions on rho
       :param vlsvReader          Some VlsvFile with a file open
       :param cellid              The cellid whose rho to calculate
       Returns the value of rho
   '''
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   # Get cells:
   vcellids = velocity_cell_data.keys()
   # Get avgs data:
   avgs = velocity_cell_data.values()
   # Get a list of velocity coordinates shifted by the solar wind bulk velocity:
   V_sw = np.array([-500000,0,0])
   v = vlsvReader.get_velocity_cell_coordinates(vcellids) - V_sw
   # Get sum of radiuses:
   radiuses = np.sum(np.abs(v)**2,axis=-1)
   # Go through velocity cells and calculate the data:
   radius2 = 468621**2
   rho_nonbackstream = 0
   for i in xrange(len(radiuses)):
      if radiuses[i] <= radius2:
         rho_nonbackstream = rho_nonbackstream + avgs[i]

   # Get the volume of a velocity cell:
   vxblocks = vlsvReader.read_parameter("vxblocks_ini")
   vxmin = vlsvReader.read_parameter("vxmin")
   vxmax = vlsvReader.read_parameter("vxmax")
   DV3 = ((vxmax-vxmin)/((float)(vxblocks)*4))**3
   print DV3
   return rho_nonbackstream*DV3

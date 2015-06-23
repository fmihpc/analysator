import numpy as np
import pylab as pl
from rotation import rotateVectorToVector

# Detector data obtained from the Themis ESA instrument paper
# http://dx.doi.org/10.1007/s11214-008-9440-2
detector_opening_angle = 6  # in degrees
detector_angle_bins = 32    # angular bins for 360 degress
detector_energy_bins = 27   # Number of energy bins
detector_min_speed = 125e3  # in m/s
detector_max_speed = 2000e3 # in m/s
detector_geometric_factor = 0.00000061 # m^2 sr E
detector_timestep = 0.003   # in s

proton_mass = 1.67e-27      # in kg

def get_dv(vlsvReader):
   # Get velocity grid sizes:
   vxcells = (int)(vlsvReader.read_parameter("vxblocks_ini"))*4
   vycells = (int)(vlsvReader.read_parameter("vyblocks_ini"))*4
   vzcells = (int)(vlsvReader.read_parameter("vzblocks_ini"))*4

   vxmin = vlsvReader.read_parameter("vxmin")
   vymin = vlsvReader.read_parameter("vymin")
   vzmin = vlsvReader.read_parameter("vzmin")
   vxmax = vlsvReader.read_parameter("vxmax")
   vymax = vlsvReader.read_parameter("vymax")
   vzmax = vlsvReader.read_parameter("vzmax")

   dvx = (vxmax - vxmin) / (float)(vxcells)
   dvy = (vymax - vymin) / (float)(vycells)
   dvz = (vzmax - vzmin) / (float)(vzcells)
   return [dvx,dvy,dvz]

def themis_observation_from_file( vlsvReader, cellid, spin_axis=np.array([0,1,0])):
   ''' Calculates artificial THEMIS EMS observation from the given cell
   :param vlsvReader:        Some VlsvReader class with a file open
   :type vlsvReader:         :class:`vlsvfile.VlsvReader`
   :param cellid:            The cell id where the distribution is supposet to be sampled NOTE: The cell id must have a velocity distribution!
   :param spin_axis: spacecraft spin axis direction
   :returns: detector bins angles, energies and counts [bin_thetas,energies,min,max,counts]
   
   .. code-block:: python
   
   # Example usage:
   vlsvReader = VlsvReader("fullf.0001.vlsv")
   angles, energies, vmin, vmax, values = themis_observation_from_file( vlsvReader=self.vlsvReader, cellid=cellid)
   # plot:
   grid_r, grid_theta = np.meshgrid(energies,angles)
   fig,ax=pl.subplots(subplot_kw=dict(projection="polar"),figsize=(12,10))
   ax.set_title("Detector view at cell " + str(cellid))
   cax = ax.pcolormesh(grid_theta,grid_r,values, norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
   pl.show()
   '''
   # Read the velocity cells:
   velocity_cell_data = vlsvReader.read_velocity_cells(cellid)
   if velocity_cell_data == []:
      # No velocity space data here, return empty result
      return [[],[],[]]
   # Calculate the detector histogram
   # Get cells:
   vcellids = velocity_cell_data.keys()
   # Get a list of velocity coordinates:
   velocity_coordinates = vlsvReader.get_velocity_cell_coordinates(vcellids)
   angles,energies,detector_values = themis_observation(velocity_cell_data, velocity_coordinates,spin_axis)

   # Scale values by velocity gridsize to calculate densities
   dvx,dvy,dvz = get_dv(vlsvReader)
   detector_values *= dvx*dvy*dvz

   # Calc min and max
   val_min = np.min(detector_values)
   if val_min == 0:
        val_min = 1e-18
   val_max = np.max(detector_values)
   if val_max <= val_min:
        val_max = val_min * 10.

   return [angles,energies,val_min,val_max,detector_values]



def themis_observation(velocity_cell_data, velocity_coordinates, spin_axis=np.array([0,0,1])):
   ''' Calculates artificial THEMIS EMS observation from the given velocity space data
   
   :param velocity_cell_data: velocity cell information as obtained from vlsvReader
   :param velocity_coordinates: coordinates associated with the cells
   :param spin_axis: spacecraft spin axis direction
   :returns: detector bins angles, energies and counts [bin_thetas,energies,counts]
   
   .. code-block:: python
   
   # Example usage:
   vlsvReader = VlsvReader("fullf.0001.vlsv")
   angles, energies, vmin, vmax, values = themis_observation_from_file( vlsvReader=self.vlsvReader, cellid=cellid)
   # plot:
   grid_r, grid_theta = np.meshgrid(energies,angles)
   fig,ax=pl.subplots(subplot_kw=dict(projection="polar"),figsize=(12,10))
   ax.set_title("Detector view at cell " + str(cellid))
   cax = ax.pcolormesh(grid_theta,grid_r,values, norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
   pl.show()
   '''
 
   # Get avgs data:
   avgs = velocity_cell_data.values()
   # Shift to plasma frame
   #if plasmaframe == True:
   #   velocity_coordinates = velocity_coordinates - bulk_velocity
   # Get velocity absolute values:
   v_abs = np.sum(np.abs(velocity_coordinates)**2,axis=-1)**(1./2)
   # Transform into detector's frame
   v_rotated = rotateVectorToVector(velocity_coordinates, spin_axis)
   v_rotated = np.asarray(v_rotated)

   # Calculate phi and theta angles
   phi_angles = (np.arctan2(v_rotated[:,0], v_rotated[:,1]) + np.pi) / (2*np.pi) * 360  # 0 .. 360
   theta_angles = np.arccos(v_rotated[:,0] / v_abs) / (2*np.pi) * 360                   # 0 .. 180
   
   # now map them all into detector bins
   detector_values = np.zeros([detector_angle_bins, detector_energy_bins])
   for i in range(0,v_abs.shape[0]):
       equator_angle = np.abs(90.-theta_angles[i])
       if equator_angle > detector_opening_angle:
           continue
       angle_bin = phi_angles[i] / 360 * detector_angle_bins - .5
       energy_bin = np.log(v_abs[i] / detector_min_speed)/np.log(detector_max_speed/detector_min_speed) * detector_energy_bins
       
       # Calculate actual countrate from phasespace density
       # (=rho*g/E*v*dt)
       # TODO: Missing sterads
       avgs[i] *= detector_geometric_factor * 2 / v_abs[i] / proton_mass * detector_timestep
       detector_values[angle_bin,energy_bin] += avgs[i]

   return_angles = np.linspace(0.,2*np.pi,detector_angle_bins)
   return_energies = np.logspace(np.log10(detector_min_speed/1e3),np.log10(detector_max_speed/1e3),detector_energy_bins)
   return [return_angles,return_energies,detector_values]

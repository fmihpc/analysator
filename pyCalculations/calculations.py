''' This is a module that includes most, if not all, calculations within the analysator. Please type "." and press tab to see which functions this module includes.

Example:
import pytools as pt

pt.calculations.
#press tab

#Output: list of functions

pt.calculations.fourier?
#Press Enter

'''

# List of functions and classes that should be imported into the interface
from intpol_file import vlsv_intpol_file
from intpol_points import vlsv_intpol_points
from cutthrough import cut_through
from fourier import fourier
from variable import VariableInfo
from timeevolution import cell_time_evolution
from pitchangle import pitch_angles
#from backstream import extract_velocity_cells_sphere, extract_velocity_cells_non_sphere
from gyrophaseangle import gyrophase_angles_from_file
from themis_observation import themis_observation_from_file
from themis_observation import themis_plot_detector
from themis_observation import themis_plot_phasespace_contour
from themis_observation import themis_plot_phasespace_helistyle
#from themis_observation import simulation_to_spacecraft_frame
from cut3d import cut3d
from lineout import lineout
import fit
from fieldtracer import static_field_tracer
from fieldtracer import dynamic_field_tracer
from temperature_anisotropy import expr_Tanisotropy

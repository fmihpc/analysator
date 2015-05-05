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
from cutthrough import cut_through
from fourier import fourier
from variable import VariableInfo
from timeevolution import cell_time_evolution
from pitchangle import pitch_angles
#from backstream import extract_velocity_cells_sphere, extract_velocity_cells_non_sphere
from gyrophaseangle import gyrophase_angles_from_file
from cut3d import cut3d
from lineout import lineout
import fit

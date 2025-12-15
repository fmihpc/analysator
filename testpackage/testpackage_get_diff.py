#from testpackage_helper import system_call
from create_env import system_call
import logging
branch='image_compare'
git_diff=system_call(f'git diff --name-only origin/{branch}...').split('\n')

#Dictionary that tells which testpackage runs to run (values) if changes were made to these files (keys).
#Checking uses 'in' operation, case insensitive
#If None, everything will be run

file_checks = {
"plot_threeslice.py":"plot_threeslice",
"plot_colormap.py":"plot_colormap",
"plot_colormap3dslice.py":"plot_colormap3dslice",
"plot_ionosphere.py":"plot_ionosphere",
"plot_isosurface.py":["plot_isosurface","plot_neutral_sheet"],
"plot_vdf.py":"plot_vdf",
"plot_vdf_profiles.py":"plot_vdf_profiles",
"plot_vdfdiff.py":"plot_vdfdiff",
"plot_variables.py":None,
"plot_helpers.py":None,
"plot.py":None,
"colormaps.py":None,
"calculations":None,
"vlsv":None,
"testpackage_":None,
"MayaVi":None,
"compare_images.yml":None,
"miscellaneous":None,
}
#Override if there are many changes as run all tests
logger = logging.getLogger(__name__)
logging.basicConfig(filename="diff_log.txt",level=logging.INFO)
run_all=False
testpackage_check=True
if len(git_diff)>6:
    logging.info(f"Multiple ({len(git_diff)}) changes, will run all tests")
    run_all=True

output=[]
for diff_line in git_diff:
    for key,val in file_checks.items():
        if key.lower() in diff_line.lower():
            if 'testpackage_' in key.lower() and testpackage_check:
                testpackage_check=False
                logging.warning(f'::warning:: Testpackage has changed in the current branch as compared to {branch}, make sure the test is still comparable with current verification_set!')
            if not val:
                run_all=True
            elif type(val)==list: 
                output.extend(val)
            else:
                output.append(val)



if run_all:
    quit()
elif output:
    print(" ".join(output))
else:
    print("pass")



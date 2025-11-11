#from testpackage_helper import system_call
from create_env import system_call

git_diff=system_call('git diff --name-only origin/master...').split('\n')
git_diff = ["plot_vdfdiff.py","plot_isosurface.py"]

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
"testpackage_commons.py":None,
"MayaVi":None,
"compare_images.yml":None
}

#Override if there are many changes as run all tests
if len(git_diff)>30:
    quit()

output=[]

for diff_line in git_diff:
    for key,val in file_checks.items():
        if key.lower() in diff_line.lower():
            if not val:
                #run all tests
                quit()
            elif type(val)==list: 
                output.extend(val)
            else:
                output.append(val)


if output:
    print(" ".join(output))
else:
    print("None")
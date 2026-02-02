from testpackage_helpers import system_call
import sys
branch=sys.argv[1]
if branch not in ["master","dev"]:
    raise SystemError("Pull request target not master or dev, this file should not even be running! Something is likely wrong with the github workflow.")

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
"testpackage_commons.py":None,
"testpackage_compare.py":None,
"testpackage_workflow.sh":None,
"testpackage_get_job_error.py":"plot_vdf",#Doubt this needs all, just some
"testpackage_get_diff.py":None,
"testpackage_helpers.py":None,
"testpackage_custom_expr.py":None,
"MayaVi":None,
"run_compare.sh":None,
"compare_images.yml":None,
"compare_images_full.yml":None,
"miscellaneous":None,
}

f=open("diff_log.txt","w")
#Override if there are many changes -> run all tests
run_all=False
testpackage_check=True
if len(git_diff)>30:
    f.write(f"Multiple ({len(git_diff)}) changes, will run all tests\n")
    run_all=True

output=[]
for diff_line in git_diff:
    if run_all:
        break
    for key,val in file_checks.items():
        if key.lower() in diff_line.lower():
            if 'testpackage_' in key.lower() and testpackage_check:
                testpackage_check=False
                f.write(f'::warning::Testpackage has changed in the current branch as compared to {branch}, make sure the test is still comparable with current verification_set!\n')
            if not val:
                run_all=True
            elif type(val) is list:
                output+=list(set(val)-set(output))
            elif val not in output:
                output.append(val)

f.close()

if run_all:
    quit()
elif output:
    print(" ".join(output))
else:
    print("pass")



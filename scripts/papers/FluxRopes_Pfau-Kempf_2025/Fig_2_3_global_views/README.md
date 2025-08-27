# Generating Figures 2 and 3

## Tracing the streamlines from regions of interest

Scripts to trace field lines from regions of interest and store them as VTK files for later plotting e.g. in VisIt.

Whereas `extract_vtk.py` contains both stages, their parallelization is different so it is recommended to run things in two stages.

In the current version, tracing is still using the [yt library](https://yt-project.org/) streamline tracing tool that is best parallelised with multiple threads per instance, e.g. using [[job_tracing.sh]]. That creates a pickledump file for each of forward and backward tracing, and for each time step. Note for future users: the newer Analysator field line tracing is more efficient and should be used instead of yt.

In a second stage, e.g. using [[job_vtk.sh]], the pickledumps are reimported and VTK files in `*.vtu` format are generated. 

Both stages are *very* intensive if run for multiple thresholds and/or hundreds of time steps of a large run, be mindful of other users when setting parameters of your array job, run a small benchmark to get an idea how long it takes per file/time step, and/or ask advice of current project managers/system administrators if unsure whether you will hog/crash the system on your first attempt.


## Plotting global views

Figs. 2 and 3 in the paper (and a lot of exploratory work) were done with [VisIt](https://visit-dav.github.io/visit-website/).

[[visit_generate.sh]] first runs VisIt twice to generate the raw frames at the three different cutoff values for both the dayside ([[visit_dayside_cutoff.py]]) and nightside ([[visit_nightside_cutoff.py]]) view. Then a lot of magick is used to composite them, overlay the annotations and arrows etc. This requires a working install of VisIt, as installed in the path indicated on Carrington.

Note the not-very-documented trick to fiddle with the legend in a VisIt script.


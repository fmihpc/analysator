# Generating Figure 4 and Animation 1

[[figure.sh]] and [[animation.sh]] are the one-step functions to generate the frames, [[makemovies.sh]] uses `ffmpeg` to generate the animation `mp4` file. 

[[plot_multi_with_xo.py]] uses `gridspec` to automatically layout the frames, in a parameterised way so that one can easily plot more or less slices at set interval and range.

[[plot_colormap3dslice_withcontours.py]] is a customised version of `plot_colormap3dslice.py` that adds contours and X/O-line markers as used in the Figure, and [[plot_ionosphere_with_footpoints.py]] is a customised version of `plot_ionosphere.py` including the footpoints and other tweaks. Notably, it contains the 2D colour scale feature used for the footpoints, courtesy of Lauri Pänkäläinen.
 

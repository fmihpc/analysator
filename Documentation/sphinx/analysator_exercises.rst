Analysator exercises
====================

Why we teach this lesson
------------------------

This page has a number of suggested exercises to help the
participants to learn of some of the possibilities of the Analysator
visualization toolkit. One should first familiarize oneself with the
Analysator introduction page and keep those lessons in mind. Remember
that most plotting routines accept **keywords** with certain pre-determined default values.

Outputs
-------

Running Analysator within a Jupyter interface will by default show the results on-screen, but will also store the image as a .png file in the default Analysator output directory. You can deactivate .png output by providing a plotting routine with the ``draw=True`` keyword.

VlsvReader
----------

It may be helpful to open a file for access by creating a VlsvReader object. Here are three sample files: one from a 5D run, one from a 6D run, and one from a coarse 6D test with an ionospheric inner boundary active.

.. code-block:: python

   f2d = pt.vlsvfile.VlsvReader("/scratch/project_465000693/example_data/AGD/bulk.0000904.vlsv")
   f3d = pt.vlsvfile.VlsvReader("/scratch/project_465000693/example_data/EGE/bulk.0001500.vlsv")
   fiono = pt.vlsvfile.VlsvReader("/scratch/project_465000693/example_data/IGtest/bulk.0000229.vlsv")

Listing available variables
---------------------------

Within python, you can list available variables for a given file as a list, or as an extended
list of all available data reducers and operators. The return value of this function will be grouped as PARAMETERs (one per file), VARIABLEs (one per simulation cell)

.. code-block:: python

   f.list()
   f.list(operator=True, datareducer=True)

It is important to note that the list of datareducers is *not* limited to those
supported by the available data, but rather lists all supported by the plotting library.
   
Plotting 2D/5D simulation data
------------------------------

The easiest manner to plot 2D data from a 5D simulation run is to use the pre-existing
plotting routines. These will sort, scale, and crop the data as requested and provide
suitable axes and color bars. You can access a given file via several approaches:

.. code-block:: python
                
   pt.plot.plot_colormap(filename="/scratch/project_465000693/example_data/AGD/bulk.0000904.vlsv")
   pt.plot.plot_colormap(vlsvobj=f2d)
   pt.plot.plot_colormap(filedir="/scratch/project_465000693/example_data/AGD/", step=904)
   
The first thing to trial is plotting different variables from the data. First, use ``f2d.list()`` to review which outputs were active during the simulation. PARAMETERs cannot be plotted as they are singular values, VARIABLEs can. For any given VLSV file you should see at least overall vg-variables and per-population vg-variables (e.g. proton values). Often, you can find also fieldsolver grid fg-variables and, only in ionospheric 6D runs, ionospheric grid ig-variables.  

Altering basic plot properties
******************************

Selecting a new output variable (replacing the default proton number density) is performed as, for example:

.. code-block:: python
                
   pt.plot.plot_colormap(vlsvobj=f2d, var='vg_blocks')

Select a handful of output variables from ``f2d.list()`` and try plotting them. Once you're familiar with how vg-variables work, you can also try proton-specific values and then fg-variables.

For scalar variables such as number densities, the default setting is to use a logarithmic colour scale. When selecting a vector variable such as ``proton/vg_v``, the default option is to reduce the value to the scalar magnitude of the vector. To view a singular component of a vector, you can select the component via the ``operator`` or ``op`` keyword:,  ``op='x'``, ``op='y'``, or  ``op='z'``. This will automatically swap the colour scale to be linear, facilitating both positive and negative values. If, on the other hand, you want to plot the magnitude of a component, you can fold the output into an absolute value with the keyword ``absolute=True``. As the default color map for a linear plot is a diverging map with white in the centre, it is possible to enforce white mapping to zero by forcing the range to be symmetric with ```symmetric=True``.

Next, try zooming around. You can adjust the plot domain by cropping it, selecting a plotting box either in metres or in units of Earth radius. You can select the crop box as either ``boxm=[xmin,xmax,ymin,ymax]`` in metres, or ``boxre=[xmin,xmax,ymin,ymax]`` in units of Earth radius. The default unit used for plot axes is Earth radii, but you can change them to metres with ``axisunit=0``, kilometres with ``axisunit=3``, or e.g. megametres with ``axisunit=6``. 

Altering the colour scale
*************************

Plotting is all about pretty pictures, and pretty pictures are all about colormaps. https://matplotlib.org/stable/users/explain/colors/colormaps.html contains a list of available matplotlib colormaps, supported by analysator. In addition to the standard colormaps, we have created some extra ones such as ``hot_desaturated`` and ``warhol``. Try plotting some of your variables with various colormaps e.g. with the keywords ``colormap='viridis'`` or ``colormap='plasma'``. 

Next, try plotting a variable which would make sense to plot on a linear scale, such as ``vg_rank``. To activate a linear scale, give the keyword ``lin=True``. To adjust the number of tick marks on the colourbar axis, you can provide them as e.g. ``lin=5``.

By default, the minimum and maximum values of the colourbar axis are selected from the visible data. If you which to adjust the limits, to over- or undersaturate regions, you can set the minimum and maximum values of the scale with e.g. ``vmin=1`` and ``vmax=100``. Now, as a next exercise, plot a magnetic field variable, such as ``fg_b``, and turn off the masking of the inner boundary region with ``nomask=True``. You will notice that the centre of the dipole field will dominate the plot. Next, adjust the minimum and maximum values along with the selected colour map to view the regions of interest of your plot. You may note that different regions of the magnetosphere (sheath, tail, foreshock) are best plotted with different colourbar ranges.

When plotting magnetospheric magnetic field values, you might note that having the field strength in Tesla is not always smart. By providing the keyword ``vscale=None``, you allow the routine to auto-scale to a suggested value, e.g. nanotesla. Setting the value manually as ``vscale=1e9`` provides the same result. Other scaling factors can also be used, and they may or may not offer suitable unit names.

The tick marks of the colourbar default to using scientific notation, which may be deactivated with the keyword ``usesci=False``.

Should you wish to evaluate a wide range of variability whist still allowing both positive and negative values, Analysator allows the use of the symmetric logarithmic colour scale. Try it by plotting the out-of-plane magnetic field component ``var=fg_b, op='z', symlog=0`` where giving symlog the value of zero allows it to self-determine the extent of the linear range at the centre between two mirrored logarithmic ranges.

Overlaying data on top of the plot
**********************************

In addition to plotting a single value as a colormap, it is possible to overlay other information. A Vlasiator watermark is available by setting the ``wmark`` (colour) or ``wmarkb`` (black) keyword, which takes a location string such as ``"NW","NE","SW"``, or ``"SW"``. Setting ``Earth=true`` plots an Earth symbol inside the inner boundary with a realistic radius.

Setting the keword ``fsaved=True`` will add contours delienating all spatial cells which contain saved-to-disk velocity distribution functions.

Vector quantities can be overplotted as streamlines or vector maps. Give some a go!
.. code-block:: python

   :kword vectors:     Set to a vector variable to overplot (unit length vectors, color displays variable magnitude)
   :kword vectordensity: Aim for how many vectors to show in plot window (default 100)
   :kword vectorcolormap: Colormap to use for overplotted vectors (default: gray)
   :kword vectorsize:  Scaling of vector sizes

   :kword streamlines: Set to a vector variable to overplot as streamlines
   :kword streamlinedensity: Set streamline density (default 1)
   :kword streamlinecolor: Set streamline color (default white)
   :kword streamlinethick: Set streamline thickness

If a pre-generated fluxfunction file is available (only for 5D simulations), magnetic field lines can be plotted
with greater accuracy than via streamlines. Flux functions are generated with the fluxfunction tool distributed as part of Vlasiator. A plotting example:

.. code-block:: python
                
   pt.plot.plot_colormap(filedir='/scratch/project_465000693/example_data/AGD/', var='proton/vg_v', boxre=[0,40,-20,20], fluxdir='/scratch/project_465000693/example_data/AGD/flux/',step=904,lin=True,vscale=None,vmin=400, flux_levels=50)

Fine-tuning of plot properties
******************************

Several keywords exist for fine-tuning Analysator plot properties. As usual, the description of these can be found by calling the help functionality with ``pt.plot.plot_colormap?``. Some examples are provided below

.. code-block:: python

   :kword noborder:    Plot figure edge-to-edge without borders (default off)
   :kword noxlabels:   Suppress x-axis labels and title
   :kword noylabels:   Suppress y-axis labels and title
   :kword scale:       Scale text size (default=1.0)
   :kword thick:       line and axis thickness, default=1.0
   :kword nocb:        Set to suppress drawing of colourbar
   :kword internalcb:  Set to draw colorbar inside plot instead of outside. If set to a text
                       string, tries to use that as the location, e.g. "NW","NE","SW","SW"
   :kword highres:     Creates the image in high resolution, scaled up by this value (suitable for print).
   :kword tickinterval: Interval at which to have ticks on axes (not colorbar)
   :kword title:       string to use as plot title instead of time.
                       Special case: Set to "msec" to plot time with millisecond accuracy or "musec"
                       for microsecond accuracy. "sec" is integer second accuracy.
   :kword cbtitle:     string to use as colorbar title instead of map name
   :kword nooverwrite: Set to only perform actions if the target output file does not yet exist                  
   
Plotting 3D/6D simulation data
------------------------------

For 3D/6D simulation data, interactive VisIt plotting is often the easiest and most powerful way to go. Analysator still provides several useful tools for plotting 3D data reduced into 3D images. These methods can be expensive, as the AMR MPIgrid data is resampled to the highest resolution. This allows handling of both vg-grid and fg-grid variables. The simplest method is to plot a Cartesian 2D slice out of the 3D domain:

.. code-block:: python
                
   pt.plot.plot_colormap3dslice(filename="/scratch/project_465000693/example_data/EGE/bulk.0001500.vlsv")
   pt.plot.plot_colormap3dslice(vlsvobj=f3d)
   pt.plot.plot_colormap3dslice(filedir="/scratch/project_465000693/example_data/EGE/", step=1500)

For these Cartesian slices, many of the options of regular 2D plots are accepted. In addition, you can select the direction normal to the plot slice with e.g. ``normal='y'`` (the default) or ``normal='z'``. The coordinate along this normal direction used for the slice is set with the keyword ``cutpoint`` in metres or ``cutpointre`` in Earth radii.

Another, much heavier option is the *threeslice* which intersects three Cartesian planes through the simulation domain and plots them all. Generating one of these images can take several tens of seconds and require significant memory from the Jupyter server.
   
.. code-block:: python

   pt.plot.plot_threeslice(filename="/scratch/project_465000693/example_data/EGE/bulk.0001500.vlsv",draw=1)

Another again heavy option is to plot an isosurface of a 3D plot, showing one variable at a surface where another variable is constant. Generating one of these images may also take several tens of seconds and require significant memory from the Jupyter server. A plotting example:

.. code-block:: python

   pt.plot.plot_isosurface(filename="/scratch/project_465000693/example_data/EGE/bulk.0001500.vlsv",surf_var='vg_b_vol',surf_op='x',surf_level=0,color_var='proton/vg_v',color_op='x',draw=1,usesci=0, vscale=None, boxre=[-40,0,-20,20,-20,20],angle=[20,210])

For isosurfaces, cropping the plotting region can help significantly with both image clarity and plot time. The ``angle`` keyword is used to define both the elevation angle and the azimuthal rotation around the z axis. The surface can be set to be opaque with ``transparent=False``.
   
Plotting ionospheric simulation data
------------------------------------

A separate routine exists for plotting ionospheric values flattened on a polar plot. For example:

.. code-block:: python
                
   pt.plot.plot_ionosphere(vlsvobj=fiono, var='ig_fac',viewdir=1, symlog=0, draw=1)

Plotting velocity distribution functions
----------------------------------------

For some spatial cells, Vlasiator stores VDFs in addition to reduced data such as moments. These velocity distribution functions can also be evaluated using Analysator. The VDF plotter routine accepts either direct CellID values as a list, or alternatively coordinates in units metres or Earth radii, from which it can search for the closest cell with a stored VDF. An example:

.. code-block:: python

   pt.plot.plot_vdf(filename="/scratch/project_465000693/example_data/AGD/bulk.0000904.vlsv",draw=1,
   coordre=[5,-10,0],center="peak",bpara1=1,slicethick=0)
                   
Since VDFs are a sparse blob of phase-space, which may or may not intersect the origin, it's usually a smart move to either plot a projection with ``slicethick=0`` and/or re-centre the plot on either the bulk velocity of the plasma or the peak position of phase-space density. Various rotations are available, providing directions evaluated from the bulk flow, the magnetic field, etc.

Advanced plotting methods
-------------------------

Advanced methods: axes, post-processing, overlaying
***************************************************

Several advanced methods exist beyond the scope of this tutorial.

An Analysator plot can be placed inside another image using the ``axes`` keyword. Examples can be found in ``examples/gridspec_plot.py`` and ``examples/multi_panel_plot.py``.

Variables read from the vlsv file can be passed to user-provided functions for post-processing. There are two types of functions supported:

Expression functions
--------------------

An expression takes arrays of variables, computes a new value from them (a scalar or a vector), and then returns it to the main plotting routine. The colormap variable is replaced with the result of the expression. Examples of more involved expressions can be found in ``examples/generate_panel.py``. A long list of simpler expressions can be found in ``pyPlots/plot_helpers.py``.

Expressions can also be instructed to be given several timesteps of data in order to facilitate time derivatives of variables.

External functions
------------------

External functions act in many ways like expressions, but they are also given the axes of the plot along with coordinate information. An external function does not replace the main variable of the plotting routine, but can be used to overlay information, variables, or other information on top of the main plot. 

Advanced methods 2: energy spectrograms, virtual spacecraft
***********************************************************

More involved data analysis is also outside the scope of this tutorial, but the interested reader should look at topics within the Analysator wiki.

Plotting a virtual spacecraft time profile involves deciding on the position for the spacecraft and choosing a time extent for the measurement. There is no quick tool for this, as for each case the required variables are likely different, but a simple example script can be found as ``examples/plot_VSC.py``

Plotting time-energy spectrograms is a highly requested feature, and rudimentary scripts do exist within the Analysator repository, e.g. ``scripts/plot_time_energy_spectrogram.py``. Future updates will most likely add a full-fledged plotting routine to the toolkit.

Interesting questions you might get
-----------------------------------


Q: How can I plot XXX or YYY?

A: look through the examples and the post-processing scripts in ``pyPlots/plot_helpers.py``. If you are not successful, you can look through those for insipration, code it yourself, and make a pull request!

A2: You can of course always ask the developers directly - perhaps the functionality exists in a non-released feature branch.

Q: Why does my plotting routine crash?

A: Analysator has formed organically over many years, with new tools and features coming in along with new data types. Sometimes mistakes slip through. Also, it's worth checking if you can find more up-to-date versions of scipy, matplotlib, or numpy - that may help.

Typical pitfalls
----------------

- Forgetting to pull the latest version of Analysator

- Directly editing the plotting scripts instead of using external or expression functionality, leading to conflicts when Analysator master gets updated

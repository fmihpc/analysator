Introduction to Analysator
==========================

Why we teach this lesson
------------------------

Here we introduce the use of Analysator python tools for plotting and accessing .vlsv data.

The Analysator tool package contains .vlsv file accessing routines, helper routines,
data post-processing libraries, and a handful of useful plotting scripts. The Analysator
files and routines are designed to act as a starting point for more involved python post-processing
of Vlasiator data, not as an end point, and as such, can and should be extended by the end-user.

Intended learning outcomes
--------------------------

You shall learn how to open a .vlsv file in Analysator, how to view file contents, and how to perform
plotting of Vlasiator spatial data in 2D and 3D. You will practice accessing Vlasov grid, FSgrid, and
ionospheric grid data. You will learn about some of the most useful Analysator plotting options.

Timing
------



Preparing exercises
-------------------

Clone the Analysator git repository onto the computer where you intend to use it.
Vlasiator .vlsv output files can be quite large, but at the same time you should have access to a
graphical interface on the machine.

For cloning, you can use the SSH or the HTTP method.

.. code-block:: bash

   cd $HOME
   git clone https://github.com/fmihpc/analysator.git
   git clone git@github.com:fmihpc/analysator.git

Next, you should set your ``PYTHONPATH`` variable to include your Analysator install directory. One way is to add this line to your startup script:
``export PYTHONPATH=$PYTHONPATH:$HOME/analysator``

Analysator on LUMI with jupyterhub
----------------------------------

For this introductory course, we will be using analysator in an interactive manner through the web interface at
`<https://www.lumi.csc.fi>`_

Log in and select a Jupyter session (*not* Jupyter for courses). Verify that you have the correct project selected (``project_465000693``) and the interactive partition. You may want to select more than just 1 CPU in order to have enough memory to open large files. In settings, select "advanced" and type the following into the window "Script to start":

.. code-block:: bash

   module use /appl/local/csc/modulefiles/
   module load pytorch/2.1
   export PTNOLATEX=1
   export PYTHONPATH=$PYTHONPATH:$HOME/analysator:

N.B. Remember to include the semicolon at the end of the last line! Next, launch the Jupyter session. Once your job has queued and launched, you can launch the Jupyter session and verify operation by importing Analysator. Type the command in the window and execute with shift+Enter.

.. code-block:: python

   [1]: import pytools as pt
   Using LaTeX formatting
   Using backend module://matplotlib_inline.backend_inline
   Using matplotlib version 3.8.1

   [2]:

Analysator required packages on other systems
---------------------------------------------

Analysator should work mostly equally well on both Python 2.7 and Python 3.x. However, versions prior to
Python 3.8 no longer receive security support and are thus not recommended. Use of iPython, jupyter, or
a similar interface is recommended for ease of use. 

Analysator only requires a small selection of python libraries, namely matplotlib, scipy>1.7 and numpy. On modern
systems, these should be pre-installed. Some legacy sections of Analysator also use MayaVi2 but those
are no longer updated or supported. 

To verify the availability of required libraries, it is suggested to try importing Analysator in python:

.. code-block:: python

   In [1]: import pytools as pt
   Using LaTeX formatting
   Using backend module://matplotlib_inline.backend_inline
   Using matplotlib version 3.8.1

   In [2]:

Other practical aspects
-----------------------
A TeX Live installation (or similar) is recommended for formatting of plotting text. If one is not available
on the target system, output can be forced to use TeX-like markup supported directly by matplotlib.
This is achieved by setting the system variable ``export PTNOLATEX=1``. This will negatively impact output
of e.g. bolded text, but is required on e.g. the LUMI web interface.

On systems without an x-windowing system such as compute nodes on a cluster (or if using it is
prohibitively slow due to e.g. network weather), Analysator can be set to ignore X-windowing and
use a non-interactive frontend by setting the system variable ``export PTNONINTERACTIVE=1``. In this
case, outputs are generated into .png files and should be transferred to another system for viewing. This is the suggested approach when using a batch job to generate several images/frames in order to e.g. build a movie.

If necessary, the matplotlib frontend can be declared manually with a system variable,
for example, ```export PTBACKEND=Qt5Agg```

The default directory for image file output for some Analysator plotting tools is ``$HOME/Plots``.
This setting can be altered with the system variable ``export PTOUTPUTDIR=/target/directory/``.

Analysator function options
---------------------------

The formalism of providing Analysator plotting functions with arguments is similar to matlab or IDL, utilizing keywords. Many keywords have a default value of e.g. None, which the code checks against. 

Interactive help
----------------

Most Analysator functions and classes contain up-to-date help, which is accessable in the python interpreter:

.. code-block::

   pt.plot.plot_colormap?

Interactive plots
-----------------

On some systems you can activate interactive backends in Jupyter notebooks by issuing the command ``%matplotlib ipympl`` or ``%matplotlib notebook`` before importing pytools. This is not supported on the LUMI web interface.

Reading data
------------

Access to Vlasiator output .vlsv files is handled through the Vlsvreader class. There are a number of
useful plotting routines which do not require editing the data directly, but for any in-depth scripting,
direct access routines are likely necessary.

VlsvReader
**********
Open a file for access by creating a VlsvReader object.

.. code-block:: python

   f=pt.vlsvfile.VlsvReader("/path/to/simulation/bulk.0001234.vlsv")

Listing available variables
***************************

Within python, you can list available variables as a concise list, or as a list of all available data reducers and operators:

.. code-block:: python

   f.list()
   f.list(datareducer=True,operator=True)

Reading in vlasov grid (MPIgrid) variables
******************************************

In older Vlasiator versions (before 5.0, simulation identifier second letter A through F) most
variables are saved on the MPIgrid and there is no identifying naming convention. Since version
5.0, with simulation version identifier letters starting from G, vlasov grid variables are
prepended with ``vg_``. Note that for per-population variables, this is placed after the population name.

Variables are read and returned as numpy arrays. MPIgrid (Vlasov grid) cell scalar variables are returned
as a simple 1-dimensional array. Vectors, tensors and so on have additional dimensions tacked on. Note that
the ordering of CellIDs (and thus, the corresponding order of proton number densities and all other MPIgrid
variables) will vary between files. The list of MPIgrid CellIDs and the corresponding proton number
densities can be found with

.. code-block:: python

   cellids = f.read_variable('cellid')
   rho = f.read_variable('proton/vg_rho')

In order to use the read data, it needs to be sorted and rearranged to correspond with the
spatial grid structure. If the grid is 2-D and AMR was not used, this is relatively straightforward.
Select the coordinate sizes to match the simulation domain.
   
.. code-block:: python

   [xsize, ysize, zsize] = f.get_spatial_mesh_size()
   rho_shaped = rho[cellids.argsort()].reshape([ysize,xsize])

For vector data, use

.. code-block:: python

   bvol = f.read_variable('vg_b_vol')
   bvol_shaped = bvol[cellids.argsort()].reshape([ysize,xsize,3])

Reading in vlasov grid (MPIgrid) AMR variables
**********************************************

Since the AMR mesh is not refined in blocks but rather as an octree-mesh, the cells
from which the refined mesh consists of does not directly translate to a 2D array.
Re-sampling the input data is a somewhat involved process, and the interested reader can
peruse the contents of e.g. the ``pyPlots/plot_colormap3dslice.py`` file for a working example.

Reading in field solver grid (FSgrid) variables
***********************************************

Since Vlasiator version 5.0, field solver grid (FSgrid) variables can be output and are
prepended with ``fg_``. FSgrid variables are returned as a numpy array, pre-sorted by the
reading routine, with dimensions matching the spatial dimensions and, if applicable, vector size.
For example, reading volumetric B-fields might yield an array of shape ``(1024, 736, 736, 3)``.
There is a separate routine for reading FSgrid variables, but the standard ``read_variable()``
routine will redirect to the FSgrid routine if an FSgrid variable is requested.

.. code-block:: python
                
   fg_b = f.read_fsgrid_variable('fg_b')

Please note that FSgrid variables do not support reading via CellID. Transforming CellIDs to coordinates
and to FSgrid file indices is possible via functions provided by ``pt.vlsvfile.VlsvReader`` but are outside
the scope of this introductory tutorial.

Reading variables with metadata
*******************************

Since Vlasiator 5.0, metadata is included for stored variables. The function ``read_variable_info`` returns
an object with the following fields: ``data`` (as per the ``read_variable`` or ``read_fsgrid_variable``
call), ``name``, ``units``, ``latex`` (LaTeX-formatted name), ``latexunits`` (LaTeX-formatted unit)

.. code-block:: python
                
   vg_b_vol_with_info = f.read_variable_info('vg_b_vol')

Reading spatial cut-throughs
****************************

Reading a spatial profile through the simulation can be achieved with the ``cut_through()`` method.
This supports only Vlasov grid data, not FSgrid data. AMR support is not yet included. Select the
starting and final positions and read the line profile with

.. code-block:: python

   cut=pt.calculations.cut_through(f,pos1,pos2);

here ``f`` is the .vlsv file used for reading, ``pos1`` and ``pos2`` are XYZ coordinates (in metres) and the
returned structure contains the relevant cellIDs (``cut[0]``) and position along the cut (``cut[1]``, in metres).
You can read the actual cut data with

.. code-block:: python
                
   variable=get_data(f.read_variable("vg_variablename",cut[0].data))

Plot the data with

.. code-block:: python
                
   ax.plot(cut[1].data/Re, variable)

Instead of reading all cells along a cut, there exists an alternative function which proceeds primarily along
the cut in the dominant cartesian direction and returns one cellID per row/column.

.. code-block:: python
                
   cut = pt.calculations.cut_through_step(f, pos1, pos2)

Writing Data with VlsvWriter
----------------------------

From time to time, one may wish to perform more involved operations on the grid, and re-use them later. ``VlsvWriter`` can be used to save derived data on SpatialGrid. It operates by copying the grid metadata and data layout from an existing file at initialization, and can thereafter be used to store the results of more involved processing.

.. code-block:: python

   f = pt.vlsvfile.VlsvReader(input_file)
   # Initialize, copy only the SpatialGrid mesh
   writer = pt.vlsvfile.VlsvWriter(f, output_file, copy_meshes=['SpatialGrid']) 
   # Copy some list of variables as a baseline. varlist accepts datareducer variables as well.
   writer.copy_variables(f,varlist=["CellID","proton/vg_rho","proton/vg_v","vg_b_vol","vg_e_vol","vg_beta","vg_beta_star"])

   # Do some heavy lifting that you don't want to repeat each time:
   orthogonality = lengthy_calculation_for_orthogonality(f)
   # Take care that this variable is compatible with the SpatialGrid variables,
   # and that is has the same memory layout as CellIDs!

   # Wrap the result with metadata
   varinfo = pt.calculations.VariableInfo(orthogonality, 
                                          name="vg_LMN_orthogonality",
                                          units="",
                                          latex=r"$|\hat{L}_\mathrm{MGA}\times\hat{L}_\mathrm{MDD}|$",latexunits=r"")

   # Write the result to SpatialGrid with the output_file writer
   writer.write_variable_info(varinfo, 'SpatialGrid', 1,)



Interesting questions you might get
-----------------------------------


Q: Why are the output formats so convoluted?

A: They are optimized for run-time performance, so that each MPI task can simply pour its data into
one contiguous region on-disk via MPI writes. 

A2: Evolution over time leads to interesting design choices.

Typical pitfalls
----------------

- Read Vlasov grid data and forget the order the cells based on CELLIDS

- Read FSGrid data and accidentally order that also according to CELLIDS

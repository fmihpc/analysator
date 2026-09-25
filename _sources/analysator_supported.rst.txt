Analysator supported data reducers and vlasiator variables
===========================================================

This documents is intended as a helpful but incomplete reference for vlasiator
outputs and post-processing data reducers supported by Analysator
   
Reading variables from .vlsv files
----------------------------------

Analysator supports reading multiple different types of variables:

- MPIgrid variables
- FSgrid variables
- variables directly available in the data files
- variables generated via receipes from the data file variables (named, somewhat incorrectly, datareducers)
                
Vlasiator variable naming scheme
--------------------------------

In Vlasiator versions before V5, variables were always saved on the DCCRG/MPI/Vlasov grid,
and as such, only a single naming scheme was required. With the introduction of Vlasiator V5,
variables could also be saved on the FSgrid. When the ionosphere grid was be implemented,
that can also be saved to on its own grid. To differentiate between the different grids, a
naming scheme was introduced where a prefix defines the grid type:

- No prefix: old pre-V5 data
- prefix ``vg_``: Vlasov/MPI/DCCRG grid
- prefix ``fg_``: fieldsolver/FSgrid
- prefix ``ig_``: ionosphere grid

Multipop naming (version 4 onwards)
-----------------------------------

Some variables are directly connected to a given particle species, with the variable name having a
particle species prefix. For example: ``proton/rho`` or ``proton/vg_rho`` for V4 and V5 outputs
respectively. Note that multipop outputs are always on the Vlasov grid, but still have the ``vg_``
prefix after the particle species name. Some variable such as ``rho`` and ``blocks`` which are
directly connected to particle species are availabe without the population identifier in datafiles
from version 3 and earlier, assuming the species to be protons.

Vlasiator variable correspondences
----------------------------------

The below tables list variable names used in different versions of Vlasiator proper (updated: 3.6.2021).
A blank entry in a corresponding column indicates this variable is not available.
Some debugging variables such as field derivatives have been omitted from the table.

.. csv-table:: FSgrid to pre-V5
   :file: ./supported_variables_1.csv
   :widths: 30, 70
   :header-rows: 1

.. csv-table:: V5 Vlasov grid to pre-V5
   :file: ./supported_variables_2.csv
   :widths: 30, 70
   :header-rows: 1

.. csv-table:: Multipop correspondence
   :file: ./supported_variables_3.csv
   :widths: 30, 60, 10
   :header-rows: 1

Analysator datareducers
-----------------------

Post-processed derived variables (called somewhat mistakenly datareducers) are available in analysator.
One limitation is that in current versions of analysator (as of March 2022) all datareducers require
Vlasov grid variables as they work with CellID lists. Also, datareducers do not support spatial or
temporal derivatives. However, both use of FSgrid variables and spatial and temporal derivatives in
plotting routines is possible via expressions and external functions as called by plot_colormap
and other such plotting routines (see also plot_helpers.py).

The up-to-date list of datareducers can always be found in pyVlsv/reduction.py

List of datareducers for Vlasiator versions 1...4
-------------------------------------------------

Note: If several populations exist for a v4 multipop run, temperature reducers are incorrect.
If the value is available directly instead of as a datareducer, that name is listed instead in
parathentheses. If it is available through multipop datareducers, it is shown with "mpop".

.. csv-table:: Generic datareducers
   :file: ./supported_datareducers_1.csv
   :widths: 30, 60, 10
   :header-rows: 1

List of per-population datareducers for Vlasiator versions 4 and 5.
-------------------------------------------------------------------

For these multipop datareducers, replace ``pop`` with required population name. If ``pop/`` is
omitted, a sum over per-population values is provided.

.. csv-table:: Multipop datareduces
   :file: ./supported_datareducers_2.csv
   :widths: 40, 30, 30
   :header-rows: 1


                 

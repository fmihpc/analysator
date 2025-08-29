Magnetosphere regions and bow shock: How to find
================================================


Bow shock
---------

Plasma properties for estimating bow shock position:

* plasma compression:
    * :math:`n_p > 2n_{p, sw}` [Battarbee_et_al_2020]_ (Vlasiator)
* solar wind core heating:
    * :math:`T_{core} > 4T_{sw}` [Battarbee_et_al_2020]_ (Vlasiator)
    * :math:`T_{core} = 3T_{sw}` [Suni_et_al_2021]_ (Vlasiator)
* magnetosonic Mach number:
    * :math:`M_{ms} < 1` [Battarbee_et_al_2020]_ (Vlasiator)


**In analysator:**

:mod:`regions` has an option to find the bow shock. Default method uses 1.5*solar wind density as limit.
Usage example:

.. code-block:: python

    datafile = "vlsvbulkfile.vlsv"
    outfilen = "bowshock.vlsv"
    RegionFlags(datafile, outfilen, regions=["bowshock"])


Magnetosheath
-------------

properties:

* density:
    * :math:`8 cm^{-3}` [Hudges_Introduction_to_space_physics_Ch_9]_
* temperature:
    * ion: :math:`150 eV` [Hudges_Introduction_to_space_physics_Ch_9]_
    * electron: :math:`25 eV` [Hudges_Introduction_to_space_physics_Ch_9]_
* magnetic field:
    * :math:`15 nT` [Hudges_Introduction_to_space_physics_Ch_9]_
* plasma :math:`\beta`:
    * 2.5  [Hudges_Introduction_to_space_physics_Ch_9]_

**In analysator:**

:mod:`regions` has an option to find the magnetosheath using bow shock and magnetopause:
Usage example:

.. code-block:: python

    datafile = "vlsvbulkfile.vlsv"
    outfilen = "magnetosheath.vlsv"
    RegionFlags(datafile, outfilen, regions=["magnetosheath"])


Polar cusps
-----------

*Properties:*

* plasma density: high density in comparison to solar wind
    * Ion density :math:`\geq` solar wind ion density [Pitout_et_al_2006]_ (Cluster spacecraft data)
* ion energy:
    * mean ion energy :math:`~2-3 keV` [Pitout_et_al_2006]_ /[Stenuit_et_al_2001]_
* energy flux:
    * energy flux



**In analysator:**

:mod:`regions` has an option to find cusps using convex hull of the magnetosphere.
Usage example:

.. code-block:: python

    datafile = "vlsvbulkfile.vlsv"
    outfilen = "cusps.vlsv"
    RegionFlags(datafile, outfilen, regions=["cusps"])


Tail lobes
----------

* plasma density: low
    * below :math:`0.03 cm^{-3}` [Grison_et_al_2025]_ (Cluster spacecraft data)
    * :math:`0.01 cm^{-3}` [Koskinen_Space_Storms]_ p.38
    * less than :math:`0.1 cm^{-3}` [Wolf_Introduction_to_space_physics_Ch_10]_ p.291
* plasma :math:`\beta`: low
    * typically around :math:`0.05` [Grison_et_al_2025]_ (Cluster spacecraft data)
    * :math:`3e-3` [Koskinen_Space_Storms]_ p.38
* temperature:
    * ion temperature :math:`300 eV` [Koskinen_Space_Storms]_ p.38
    * electron temperature :math:`50 eV` [Koskinen_Space_Storms]_ p.38
* magnetic field:
    * :math:`20 nT` [Koskinen_Space_Storms]_ p.38
* open magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10]_ p.291
* strong and stable magnetic field towards the Earth (northern lobe) and away from the Earth (southern lobe) [Coxon_et_al_2016]_

Separated from the plasma sheet by the plasma sheet boundary layer (PSBL)


**In analysator:**

:mod:`regions` has an option to find tail lobes.
Usage example:

.. code-block:: python

    datafile = "vlsvbulkfile.vlsv"
    outfilen = "lobes.vlsv"
    RegionFlags(datafile, outfilen, regions=["lobes"])




Low-latitude boundary layer (LLBL)
----------------------------------



Properties:

* density:
    * ion number densities between those of magnetosphere and magnetosheath [Hudges_Introduction_to_space_physics_Ch_9]_ p.267
* temperature
    * ion temperatures between those of magnetosphere and magnetosheath [Hudges_Introduction_to_space_physics_Ch_9]_ p.267
* unknown field line configuration, probably a mix of open and closed field lines [Hudges_Introduction_to_space_physics_Ch_9]_ p.262



High-latitude boundary layer (HLBL)
-----------------------------------

Includes the plasma mantle on the tail side and the entry layer on the dayside

Properties:

* open magnetic field lines [Hudges_Introduction_to_space_physics_Ch_9]_ p.261





Plasma sheet boundary layer (PSBL)
----------------------------------

The plasma sheet boundary layer is a very thin boundary layer separating the tail lobes from the tail plasma sheet [Koskinen_Johdatus]_

*Properties:*

* density:
    * :math:`0.1 cm^{-3}` [Koskinen_Space_Storms]_ p.38
* temperature:
    * ion temperature :math:`1000 eV` [Koskinen_Space_Storms]_ p.38
    * electron temperature :math:`150 eV` [Koskinen_Space_Storms]_ p.38
* magnetic field:
    * :math:`20 nT` [Koskinen_Space_Storms]_ p.38
* plasma :math:`\beta` :
    * :math:`0.1` [Koskinen_Space_Storms]_ p.38
* probably closed magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10]_ p.291




Central plasma sheet
--------------------


*Properties:*

* density:
    * :math:`0.3 cm^{-3}` [Koskinen_Space_Storms]_ p.38
    * :math:`0.1-1 cm^{-3}` [Wolf_Introduction_to_space_physics_Ch_10]_ p.291
* temperature: hot
    * ion temperature :math:`4200 eV` [Koskinen_Space_Storms]_ p.38
    * electron temperature :math:`600 eV` [Koskinen_Space_Storms]_ p.38
* magnetic field:
    * :math:`10 nT` [Koskinen_Space_Storms]_ p.38, [Hudges_Introduction_to_space_physics_Ch_9]_
* plasma :math:`\beta`: high
    * :math:`6` [Koskinen_Space_Storms]_ p.38
* Mostly closed magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10]_

Inner plasma sheet: unusually low plasma beta may exist (e.g., cold tenuous plasma near the neutral sheet after long periods of northward IMF) [Boakes_et_al_2014]_, (Cluster spacecraft data)


**In analysator:**

:mod:`regions` has an option to find the central plasma sheet.
Note that the default parameters may catch some dayside cells in addition to plasma sheet.
Usage example:

.. code-block:: python

    datafile = "vlsvbulkfile.vlsv"
    outfilen = "CPS.vlsv"
    RegionFlags(datafile, outfilen, regions=["central_plasma_sheet"])




.. automodule:: regions
    :members:


------------

References 

.. [Battarbee_et_al_2020] Battarbee, M., Ganse, U., Pfau-Kempf, Y., Turc, L., Brito, T., Grandin, M., Koskela, T., and Palmroth, M.: Non-locality of Earth's quasi-parallel bow shock: injection of thermal protons in a hybrid-Vlasov simulation, Ann. Geophys., 38, 625-643, https://doi.org/10.5194/angeo-38-625-2020, 2020
.. [Suni_et_al_2021] Suni, J., Palmroth, M., Turc, L., Battarbee, M., Johlander, A., Tarvus, V., et al. (2021). Connection between foreshock structures and the generation of magnetosheath jets: Vlasiator results. Geophysical Research Letters, 48, e2021GL095655. https://doi. org/10.1029/2021GL095655
.. [Grison_et_al_2025] Grison, B., Darrouzet, F., Maggiolo, R. et al. Localization of the Cluster satellites in the geospace environment. Sci Data 12, 327 (2025). https://doi.org/10.1038/s41597-025-04639-z
.. [Koskinen_Johdatus] Koskinen, H. E. J. (2011). Johdatus plasmafysiikkaan ja sen avaruussovellutuksiin. Limes ry.
.. [Koskinen_Space_Storms] Koskinen, H. E. J. (2011). Physics of Space Storms: From the Solar Surface to the Earth. Springer-Verlag. https://doi.org/10.1007/978-3-642-00319-6
.. [Pitout_et_al_2006] Pitout, F., Escoubet, C. P., Klecker, B., and Rème, H.: Cluster survey of the mid-altitude cusp: 1. size, location, and dynamics, Ann. Geophys., 24, 3011–3026, https://doi.org/10.5194/angeo-24-3011-2006, 2006.
.. [Coxon_et_al_2016] Coxon,J.C.,C.M.Jackman, M. P. Freeman, C. Forsyth, and I. J. Rae (2016), Identifying the magnetotail lobes with Cluster magnetometer data, J. Geophys. Res. Space Physics, 121, 1436–1446, doi:10.1002/2015JA022020.
.. [Hudges_Introduction_to_space_physics_Ch_9] Hudges, W. J. (1995) The magnetopause, magnetotail and magnetic reconnection. In Kivelson, M. G., & Russell, C. T. (Eds.), Introduction to space physics (pp.227-287). Cambridge University Press.
.. [Wolf_Introduction_to_space_physics_Ch_10] Wolf, R. A. (1995) Magnetospheric configuration. In Kivelson, M. G., & Russell, C. T. (Eds.), Introduction to space physics (pp.288-329). Cambridge University Press.
.. [Boakes_et_al_2014] Boakes, P. D., Nakamura, R., Volwerk, M., and Milan, S. E. (2014). ECLAT Cluster Spacecraft Magnetotail Plasma Region Identifications (2001–2009). Dataset Papers in Science, 2014(1):684305. eprint: https://onlinelibrary.wiley.com/doi/pdf/10.1155/2014/684305
.. [Stenuit_et_al_2001] Stenuit, H., Sauvaud, J.-A., Delcourt, D. C., Mukai, T., Kokubun, S., Fujimoto, M., Buzulukova, N. Y., Kovrazhkin, R. A., Lin, R. P., and Lepping, R. P. (2001). A study of ion injections at the dawn and dusk polar edges of the auroral oval. Journal of Geophysical Research: Space Physics, 106(A12):29619–29631. eprint: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2001JA900060.

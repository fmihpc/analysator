Bow shock: How to find
======================


Plasma properties for estimating bow shock position:
----------------------------------------------------

* plasma compression 
    * :math:`n_p > 2n_{p, sw}` [Battarbee_et_al_2020]_ (Vlasiator)
* solar wind core heating
    * :math:`T_{core} > 4T_{sw}` [Battarbee_et_al_2020]_ (Vlasiator)
    * :math:`T_{core} = 3T_{sw}` [Suni_et_al_2021]_ (Vlasiator)
* magnetosonic Mach number
    * :math:`M_{ms} < 1` [Battarbee_et_al_2020]_ (Vlasiator)




In analysator:
--------------

regions.py: 

*method*: 

Convex hull from cells where :math:`n_p > 1.5 n_{p, sw}`
criteria

also velocity?


Foreshock:
----------

*properties:*

* larger fluctuations in the magnetic field, plasma velocity and plasma density than unperturbed solar wind [Grison_et_al_2025]_


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


Regions inside the bow shock:
-----------------------------

* magnetosheath: area inside bow shock but outside the magnetopause, see ...
* magnetosphere: area inside magnetopause, see again ... and ... for magnetosphere regions


------------

References:

.. [Battarbee_et_al_2020] Battarbee, M., Ganse, U., Pfau-Kempf, Y., Turc, L., Brito, T., Grandin, M., Koskela, T., and Palmroth, M.: Non-locality of Earth's quasi-parallel bow shock: injection of thermal protons in a hybrid-Vlasov simulation, Ann. Geophys., 38, 625-643, https://doi.org/10.5194/angeo-38-625-2020, 2020
.. [Suni_et_al_2021] Suni, J., Palmroth, M., Turc, L., Battarbee, M., Johlander, A., Tarvus, V., et al. (2021). Connection between foreshock structures and the generation of magnetosheath jets: Vlasiator results. Geophysical Research Letters, 48, e2021GL095655. https://doi. org/10.1029/2021GL095655
.. [Grison_et_al_2025] Grison, B., Darrouzet, F., Maggiolo, R. et al. Localization of the Cluster satellites in the geospace environment. Sci Data 12, 327 (2025). https://doi.org/10.1038/s41597-025-04639-z

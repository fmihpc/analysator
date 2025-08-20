Magnetosphere regions: How to find
==================================




Cusps
-----


*Properties:*

* plasma density: high density in comparison to solar wind
    * Ion density :math:`\geq` solar wind ion density [Pitout_et_al_2006_] (Cluster spacecraft data)
* ion energy:
    * mean ion energy :math:`~2-3 keV` [Pitout_et_al_2006_] /[Stenuit_et_al_2001_]
* energy flux:
    * energy flux



**In analysator:**


Tail lobes
----------

* plasma density: low
    * below :math:`0.03 cm^{-3}` [Grison_et_al_2025_] (Cluster spacecraft data)
    * :math:`0.01 cm^{-3}` [Koskinen_Space_Storms_ p.38]
    * less than :math:`0.1 cm^{-3}` [Wolf_Introduction_to_space_physics_Ch_10_ p.291]
* plasma :math:`\beta`: low
    * typically around :math:`0.05` [Grison_et_al_2025_] (Cluster spacecraft data)
    * :math:`3e-3` [Koskinen_Space_Storms_ p.38]
* temperature:
    * ion temperature :math:`300 eV` [Koskinen_Space_Storms_ p.38]
    * electron temperature :math:`50 eV` [Koskinen_Space_Storms_ p.38]
* magnetic field:
    * :math:`20 nT` [Koskinen_Space_Storms_ p.38]
* open magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10_ p.291]
* strong and stable magnetic field towards the Earth (northern lobe) and away from the Earth (southern lobe) [Coxon_et_al_2016_]

Separated from the plasma sheet by the plasma sheet boundary layer (PSBL)


**in analysator:**

regions.py

conditions:

* inside the magnetosphere
* plasma :math:`\beta` ....




Low-latitude boundary layer (LLBL)
----------------------------------



Properties:

* density:
    * ion number densities between those of magnetosphere and magnetosheath [Hudges_Introduction_to_space_physics_Ch_9_ p.267] 
* temperature
    * ion temperatures between those of magnetosphere and magnetosheath [Hudges_Introduction_to_space_physics_Ch_9_ p.267] 
* unknown field line configuration, probably a mix of open and closed field lines [Hudges_Introduction_to_space_physics_Ch_9_ p.262]



High-latitude boundary layer (HLBL)
-----------------------------------

Includes the plasma mantle on the tail side and the entry layer on the dayside [... cit.]

Properties:

* open magnetic field lines [Hudges_Introduction_to_space_physics_Ch_9_ p.261]





Plasma sheet boundary layer (PSBL)
----------------------------------

The plasma sheet boundary layer is a very thin boundary layer separating the tail lobes from the tail plasma sheet [Koskinen_Johdatus_]

*Properties:*

* density:
    * :math:`0.1 cm^{-3}` [Koskinen_Space_Storms_ p.38]
* temperature:
    * ion temperature :math:`1000 eV` [Koskinen_Space_Storms_ p.38]
    * electron temperature :math:`150 eV` [Koskinen_Space_Storms_ p.38]
* magnetic field:
    * :math:`20 nT` [Koskinen_Space_Storms_ p.38]
* plasma :math:`\beta` :
    * :math:`0.1` [Koskinen_Space_Storms_ p.38]
* probably closed magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10_ p.291]




Central plasma sheet
--------------------


*Properties:*

* density:
    * :math:`0.3 cm^{-3}` [Koskinen_Space_Storms_ p.38]
    * :math:`0.1-1 cm^{-3}` [Wolf_Introduction_to_space_physics_Ch_10_ p.291]
* temperature: hot
    * ion temperature :math:`4200 eV` [Koskinen_Space_Storms_ p.38]
    * electron temperature :math:`600 eV` [Koskinen_Space_Storms_ p.38]
* magnetic field:
    * :math:`10 nT` [Koskinen_Space_Storms_ p.38], [Hudges_Introduction_to_space_physics_Ch_9_]
* plasma :math:`\beta`: high
    * :math:`6` [Koskinen_Space_Storms_ p.38]
* Mostly closed magnetic field lines [Wolf_Introduction_to_space_physics_Ch_10_]




------------

References 

.. [Grison_et_al_2025] Grison, B., Darrouzet, F., Maggiolo, R. et al. Localization of the Cluster satellites in the geospace environment. Sci Data 12, 327 (2025). https://doi.org/10.1038/s41597-025-04639-z
.. [Koskinen_Johdatus] Koskinen, H. E. J. (2011). Johdatus plasmafysiikkaan ja sen avaruussovellutuksiin. Limes ry.
.. [Koskinen_Space_Storms] Koskinen, H. E. J. (2011). Physics of Space Storms: From the Solar Surface to the Earth. Springer-Verlag. https://doi.org/10.1007/978-3-642-00319-6
.. [Pitout_et_al_2006] Pitout, F., Escoubet, C. P., Klecker, B., and Rème, H.: Cluster survey of the mid-altitude cusp: 1. size, location, and dynamics, Ann. Geophys., 24, 3011–3026, https://doi.org/10.5194/angeo-24-3011-2006, 2006.
.. [Coxon_et_al_2016] Coxon,J.C.,C.M.Jackman, M. P. Freeman, C. Forsyth, and I. J. Rae (2016), Identifying the magnetotail lobes with Cluster magnetometer data, J. Geophys. Res. Space Physics, 121, 1436–1446, doi:10.1002/2015JA022020.
.. [Hudges_Introduction_to_space_physics_Ch_9] Hudges, W. J. (1995) The magnetopause, magnetotail and magnetic reconnection. In Kivelson, M. G., & Russell, C. T. (Eds.), Introduction to space physics (pp.227-287). Cambridge University Press.
.. [Wolf_Introduction_to_space_physics_Ch_10] Wolf, R. A. (1995) Magnetospheric configuration. In Kivelson, M. G., & Russell, C. T. (Eds.), Introduction to space physics (pp.288-329). Cambridge University Press.
.. [Sckopke_et_al_1981] Sckopke, N., Paschmann, G., Haerendel, G., Sonnerup, B. U. , Bame, S. J., Forbes, T. G., Hones Jr., E. W., and Russell, C. T. (1981). Structure of the low-latitude boundary layer. Journal of Geophysical Research: Space Physics, 86(A4):2099–2110. eprint: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/JA086iA04p02099
.. [Boakes_et_al_2014] Boakes, P. D., Nakamura, R., Volwerk, M., and Milan, S. E. (2014). ECLAT Cluster Spacecraft Magnetotail Plasma Region Identifications (2001–2009). Dataset Papers in Science, 2014(1):684305. eprint: https://onlinelibrary.wiley.com/doi/pdf/10.1155/2014/684305
Magnetopause: how to find
=========================

Magnetopause 

Plasma beta, beta*
------------------


Modified plasma beta

.. math:: \beta * = \dfrac{P_{th}P_{dyn}}{B^2/2\mu_0}

[Xu_et_al_2016]_, [Brenner_et_al_2021]_ 


The beta* gets values below 1 inside the magnetosphere near magnetopause and can be used to create a magnetopause surface.

Caveats: magnetotail current sheet has beta* :math:`>` 1

**in analysator:**

datareducer: beta_star, vg_beta_star

regions.py: cells with beta* value within a certain treshold are chosen, and a convex hull is constructed with vtk to represent the magnetopause.
Ideal values of beta* for magnetopause construction might be run-dependent, and the surface construction works best with a small-ish range of beta* but still big enough to gather cells evenly from all sides.


.. [Xu_et_al_2016] Xu, S., M. W. Liemohn, C. Dong, D. L. Mitchell, S. W. Bougher, and Y. Ma (2016), Pressure and ion composition boundaries at Mars, J. Geophys. Res. Space Physics, 121,  6417–6429, doi:10.1002/2016JA022644.
.. [Brenner_et_al_2021] Brenner A, Pulkkinen TI, Al Shidi Q and Toth G (2021) Stormtime Energetics: Energy Transport Across the Magnetopause in a Global MHD Simulation. Front. Astron. Space Sci. 8:756732. doi: 10.3389/fspas.2021.756732



Field line connectivity
-----------------------




Solar wind flow
---------------

Method used in e.g. [Parlmroth_et_al_2003_] 

Streamlines of velocity field that are traced from outside the bow shock curve around the magnetopause.

Caveats: sometimes some streamlines can curve into the magnetotail 


**In analysator:**

magnetopause_sw_streamline_2d.py
magnetopause_sw_streamline_3d.py

Streamlines are traced from outside the bow shock towards Earth. A subsolar point for the magnetopause is chosen to be where streamlines get closest to Earth in x-axis [y/z~0]. 

From subsolar point towards Earth, ....


From the Earth towards negative x the space is divided into yz-planes. 
Each yz-plane is then divided into sectors and magnetopause is marked to be in the middle of the sector with the radius of the closest streamline to the x-axis. 

For subsolar point, radial dayside and -x planes the closest streamline point can be changed to be n:th closest by setting keyword "ignore", in which case n-1 closest streamline points are not taken into acocunt.


2d:
streamlines: [...]
magnetopause: [...]

example usage in [...]


3d:
streamlines: [...]
magnetopause: [...]
surface triangulation: [...]

example usage in [...]






.. [Parlmroth_et_al_2003] Palmroth, M., T. I. Pulkkinen, P. Janhunen, and C.-C. Wu (2003), Stormtime energy transfer in global MHD simulation, J. Geophys. Res., 108, 1048, doi:10.1029/2002JA009446, A1.




Shue et al. (1997)
------------------

The Shue et al. (1997) [Shue_et_al_1997]_ mangnetopause model:

.. math::

    r_0 &= (11.4 + 0.013 B_z)(D_p)^{-\frac{1}{6.6}}, for B_z \geq 0 \\
        &= (11.4 + 0.14 B_z)(D_p)^{-\frac{1}{6.6}}, for B_z \leq 0 

.. math::

    \alpha = (0.58-0.010B_z)(1+0.010 D_p)


where :math:`D_p` is the dynamic pressure of the solar wind and :math:`B_z` the magnetic field z-component magnitude.
:math:`r_0` is the magnetopause standoff distance and :math:`\alpha` is the level of tail flaring. 


The magnetopause radius as a function of angle :math:`\theta` is 

.. math::
    r =  r_0 (\frac{2}{1+\cos \theta})^\alpha

**In analysator:** 

*shue.py* in *scripts*:





.. [Shue_et_al_1997] Shue, J.-H., Chao, J. K., Fu, H. C., Russell, C. T., Song, P., Khurana, K. K., and Singer, H. J. (1997). A new functional form to study the solar wind control of the magnetopause size and shape. Journal of Geophysical Research: Space Physics, 102(A5):9497–9511. eprint: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/97JA00196.
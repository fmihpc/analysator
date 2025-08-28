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

see :func:`magnetopause.magnetopause` "method" keyword options "beta_star", "beta_star_with_connectivity"

datareducer: beta_star, vg_beta_star

.. [Xu_et_al_2016] Xu, S., M. W. Liemohn, C. Dong, D. L. Mitchell, S. W. Bougher, and Y. Ma (2016), Pressure and ion composition boundaries at Mars, J. Geophys. Res. Space Physics, 121,  6417–6429, doi:10.1002/2016JA022644.
.. [Brenner_et_al_2021] Brenner A, Pulkkinen TI, Al Shidi Q and Toth G (2021) Stormtime Energetics: Energy Transport Across the Magnetopause in a Global MHD Simulation. Front. Astron. Space Sci. 8:756732. doi: 10.3389/fspas.2021.756732



Field line connectivity
-----------------------





Solar wind flow
---------------

Method used in e.g. [Palmroth_et_al_2003]_ 

Streamlines of velocity field that are traced from outside the bow shock curve around the magnetopause.

Caveats: sometimes some streamlines can curve into the magnetotail or dayside magnetoshpere


**In analysator:**

see
:func:`calculations.find_magnetopause_sw_streamline_2d`
:func:`calculations.find_magnetopause_sw_streamline_3d`

Streamlines are traced from outside the bow shock towards Earth. A subsolar point for the magnetopause is chosen to be where streamlines get closest to Earth in x-axis [y/z~0]. 

From subsolar point towards Earth, space is divided radially by spherical coordiante (theta from x-axis) angles theta and phi, and magnetopause is located by looking at streamline point distances from origo and marked to be at middle of the sector

From the Earth towards negative x the space is divided into yz-planes. 
Each yz-plane is then divided into 2d sectors and magnetopause is marked to be in the middle of the sector with the radius of the n:th closest streamline to the x-axis. 


For subsolar point, radial dayside, and -x planes the closest streamline point can be changed to be n:th closest by setting keyword *ignore*, in which case *ignore* closest streamline points are not taken into account.


2d:
Note: no radial dayside

3d:
After the magnetopause points are chosen, they are made into a surface by setting connnection lines between vertices (magnetopause points) to make surface triangles. 
This surface is then made into a vtkPolyData and returned as vtkDataSetSurfaceFilter that can be written into e.g. .vtp file.



.. [Palmroth_et_al_2003] Palmroth, M., T. I. Pulkkinen, P. Janhunen, and C.-C. Wu (2003), Stormtime energy transfer in global MHD simulation, J. Geophys. Res., 108, 1048, doi:10.1029/2002JA009446, A1.




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

see :mod:`shue`


.. [Shue_et_al_1997] Shue, J.-H., Chao, J. K., Fu, H. C., Russell, C. T., Song, P., Khurana, K. K., and Singer, H. J. (1997). A new functional form to study the solar wind control of the magnetopause size and shape. Journal of Geophysical Research: Space Physics, 102(A5):9497–9511. eprint: https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/97JA00196.





**In analysator:**
------------------

*magnetopause.py* in scripts for 3d runs
Constructs the magentopause surface with vtk's vtkDelaunay3d triangulation with optional alpha to make the surface non-convex.
Uses regions.py functions.

Important: SDF of non-convex surface might not always work

options (magnetopause() method keyword) and some notes: 

* solar wind flow ("streamlines")
    * uses *magnetopause_sw_streamline_3d.py* from pyCalculations
    * if streamlines make a turn so that the velocity points sunwards (to pos. x), the streamline is ignored from that point onwards
        * this fixes some issues with funky and inwards-turning streamlines but not all
    * very dependent on how solar wind streamlines behave
    * streamline seeds and other options can greatly affect the resulting magnetopause

* beta* ("beta_star")
    * beta* treshold might need tweaking as sometimes there are small low beta* areas in the magnetosheath that get taken in distorting the magnetopause shape at nose
    * convex hull (Delaunay_alpha=None) usually makes a nice rough magnetopause but goes over any inward dips (like polar cusps)
    * alpha shape (Delaunay_alpha= e.g. 1*R_E) does a better job at cusps and delicate shapes like vortices but might fail at SDF inside the magnetosphere
    * Delaynay3d has an easier time if the treshold is something like [0.4, 0.5] and not [0.1, 0.5]

* beta* with magnetic field line connectivity ("beta_star_with_connectivity")
    * includes closed-closed magnetic field line areas if available, otherwise like "beta_star"
    * can help with nose shape as beta* can be set lower to exclude magnetosheath low beta* areas while still getting the full dayside from field lines

* Shue et al. 1997 ("shue")
    * uses *shue.py* from scripts
    * a rough theoretical magnetopause using Shue et al. 1997 method based on B_z, solar wind density, and solar wind velocity

* user-defined parameter tresholds ("dict")
    * creates a magnetopause (or other area) using the Delaunay3d triangulation of some area where user-defined tresholds given as dictionary
    * dictionary key is data name in datafile and value is treshold used, if dictionary has multiple conditions, they all need to be fulfilled
    * dictionary example: {"vg_rho": [None, 1e5]} makes a magnetopause using cells where density is less than 1e5


.. automodule:: magnetopause
    :members:
# Flux ropes

The tools in this folder were used to post-process data from the [FHA run](http://urn.fi/urn:nbn:fi:att:3ce0f038-2c69-4c7c-8f67-7a71e9e57b56) and produce the figures of

Global evolution of flux transfer events along the magnetopause from the dayside to the far tail
Yann Pfau-Kempf, Konstantinos Papadakis, Markku Alho, Markus Battarbee, Giulia Cozzani, Lauri Pänkäläinen, Urs Ganse, Fasil Kebede, Jonas Suni, Konstantinos Horaites, Maxime Grandin, and Minna Palmroth
https://doi.org/10.5194/angeo-43-469-2025
Annales Geophysicae, Volume 43, Issue 2, Pages 469–488, 2025

## Global views
[[Fig_2_3_global_views]] contains the python and job scripts to trace the field lines from the regions where `vg_fluxrope` is above the chosen threshold. The resulting files can then be used in VisIt to plot the global views of Figs. 2 and 3.

## Volume of detection
[[Fig_4_volume]] contains the scripts to produce Fig. 4.

## Slices and footprints
[[Fig_5_Animation_1_slices_and_footprints]] contains the scripts to produce Fig. 5 and Animation 1. In particular, the customised versions of `plot_colormap3dslice.py` and `plot_ionosphere.py` are here. Features should be ported to the main branch but yours truly did not have the time to complete that unfortunately.

## Slices and virtual spacecraft traces in minimum variance frame
[[Fig_6_7_Animation_2_slices_and_vsc_in_MVA_frame]] contains the scripts to produce Figs. 6 and 7 using a customised version of `plot_colormap3dslice.py`. Same comment about feature-porting unfortunately. 

## Occurrence map
[[Fig_8_occurrence_map]] contains the scripts to produce Fig. 8.


# Notes about xkcd
Figures 5, 6, 7, and 8 use some [xkcd named colours](https://xkcd.com/color/rgb/) based on a [large survey](https://xkcd.com/color/rgb/). The colours are [included in maplotlib](https://matplotlib.org/stable/gallery/color/named_colors.html#xkcd-colors) out of the box.

The plotting scripts for Fig.s 4 to 8 include a commented line that, if activated,
```
plt.xkcd(scale=1.5, length=100, randomness=3)
```
creates the xkcd-style versions of the plots. Follow the [documentation](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.xkcd.html) to go full bananas and install the fonts too.


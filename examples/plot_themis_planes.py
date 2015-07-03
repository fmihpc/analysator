#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Given input file and cellID, plot the velocity space as it would be seen by Themis ESA detector.
# Creates output pngs in multiple cut planes:
#    x Direction         |   y Direction             |     filename
# -----------------------+---------------------------+---------------------
#   parallel to B        | perp. to B, in dir of v   | themis_B_v.png
#   parallel to B        | perp. to B, in dir of B×v | themis_B_Bxv.png
#   perp. to B, dir of v | perp. to B, in dir of Bxv | themis_B_perpperp.png
#   simulation x         | simulation y              | themis_x_y.png
#   simulation x         | simulation z              | themis_x_z.png
#   simulation y         | simulation z              | themis_y_z.png
#
# The plots where B does not directly specify an axis also have the magnetic
# field vector drawn in as an arrow.

import pytools as pt
import numpy as np
import pylab as pl
import sys

if len(sys.argv) < 2:
    sys.stderr.write("Syntax: plot_themis_planes.py <vlsvfile> cellID\n")
    sys.exit()

filename = sys.argv[1]
f = pt.vlsvfile.VlasiatorReader(filename)

cellid = int(sys.argv[2])

B = f.read_variable("B",cellid)
v = f.read_variable("v",cellid)

Bxv = np.cross(B,v)
BxvxB = np.cross(Bxv,B)

pt.calculations.themis_plot_phasespace_helistyle(f, cellid, B,v,smooth=True,xlabel=u"v∥B",ylabel=u"v∥V",vmin_min=1e-16)
pl.savefig("themis_B_v.png")
pt.calculations.themis_plot_phasespace_helistyle(f, cellid, B,Bxv,smooth=True,xlabel=u"v∥B",ylabel=u"v∥B×V", vmin_min=1e-17)
pl.savefig("themis_B_Bxv.png")
pt.calculations.themis_plot_phasespace_helistyle(f, cellid, Bxv,BxvxB,smooth=True,ylabel=u"v⟂B (in Plane of V)",xlabel=u"v⟂B (in Plane of B×V)", vmin_min=1e-17)
pl.savefig("themis_B_perpperp.png")

Bn = B / np.linalg.norm(B)
ax = pt.calculations.themis_plot_phasespace_helistyle(f, cellid, np.array([1.,0,0]),np.array([0,1.,0.]),smooth=True,xlabel=u"vx",ylabel=u"vy", vmin_min=1e-17)
ax.arrow(0,0,2000*Bn[0],2000*Bn[1],head_width=100,head_length=200, edgecolor='black', facecolor='black')
pl.savefig("themis_x_y.png")
ax = pt.calculations.themis_plot_phasespace_helistyle(f, cellid, np.array([1.,0,0]),np.array([0,0.,1.]),smooth=True,xlabel=u"vx",ylabel=u"vz", vmin_min=1e-17)
ax.arrow(0,0,2000*Bn[0],2000*Bn[2],head_width=100,head_length=200, edgecolor='black', facecolor='black')
pl.savefig("themis_x_z.png")
ax = pt.calculations.themis_plot_phasespace_helistyle(f, cellid, np.array([0.,1.,0]),np.array([0,0.,1.]),smooth=True,xlabel=u"vy",ylabel=u"vz", vmin_min=1e-17)
ax.arrow(0,0,2000*Bn[1],2000*Bn[2],head_width=100,head_length=200, edgecolor='black', facecolor='black')
pl.savefig("themis_y_z.png")

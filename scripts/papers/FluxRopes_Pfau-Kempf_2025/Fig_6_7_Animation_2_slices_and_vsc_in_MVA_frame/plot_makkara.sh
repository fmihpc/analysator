#!/usr/bin/env python
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
import sys
import pytools as pt
from variable import get_data, get_name, get_units

#plt.xkcd(scale=1.5, length=100, randomness=3)

contourvariable = "vg_fluxrope"
contouroperator = "pass"
contourvalues = 7
contourcolors = "xkcd:swamp"
contourlinewidths = 5

thick=3
streamlinethick=2
streamlinearrowsize=2

Re = 6.371e+6 # Earth radius in m


numPlots = 4
panel_labels = ("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)")
panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)
time_label_bbox = dict(boxstyle='square,pad=0.1', fc='white', ec='black', alpha=0.85, lw=thick)
times = []
cutre_x = []
cutre_y = []
cutre_z = []
box_side_re = []

times.append(1300)
box_side_re.append(18.0)
cutre_x.append(-18.0)
cutre_y.append(-17.0)
cutre_z.append(-2.5)
times.append(1400)
box_side_re.append(15.0)
cutre_x.append(-28.0)
cutre_y.append(-17.5)
cutre_z.append(-3.5)
times.append(1500)
box_side_re.append(12.0)
cutre_x.append(-38.0)
cutre_y.append(-18.0)
cutre_z.append(-2.5)
times.append(1600)
box_side_re.append(9)
cutre_x.append(-48.0)
cutre_y.append(-18.2)
cutre_z.append(-4.5)

inputLocation="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"

gsunit=8
gshbuf=1
gsvbuf=2
inter=0
cbwid=2
dpi=150
matplotlib.rcParams.update({'font.size': 35})

# Init figure    
fig = pl.figure(figsize=(2*gsunit+cbwid+gshbuf+inter,numPlots*(gsunit+gsvbuf)), dpi=dpi)
gs = fig.add_gridspec(numPlots*(gsunit+gsvbuf), 2*gsunit+cbwid+gshbuf+inter)

for i in range(numPlots):
    fig.add_subplot(gs[i*(gsunit+gsvbuf):i*(gsunit+gsvbuf)+gsunit, :gsunit]) # left
    fig.add_subplot(gs[i*(gsunit+gsvbuf):i*(gsunit+gsvbuf)+gsunit, gsunit+gshbuf:2*gsunit+gshbuf]) # right
    fig.add_subplot(gs[i*(gsunit+gsvbuf):i*(gsunit+gsvbuf)+gsunit, 2*gsunit+gshbuf:2*gsunit+gshbuf+cbwid//2]) # colorbar

axes = fig.get_axes()
#    for label in ax.yaxis.get_ticklabels():
#        label.set_visible(False)


#fig.suptitle("Flight of the Makkara")


for i in range(numPlots):
    
    timefull=str(times[i]).rjust(7, '0')
    file_name=inputLocation+"bulk1."+timefull+".vlsv"
    f=pt.vlsvfile.VlsvReader(file_name=file_name)

    print("Slicing timestep "+str(times[i])+" s")

    # plot slice
    ax=axes[i*3]
    ax.text(cutre_y[i]-0.8*box_side_re[i], cutre_z[i]+0.57*box_side_re[i], "$t = "+str(times[i])+"\,\mathrm{s}$", fontsize=42, bbox=time_label_bbox)
    ax.text(cutre_y[i]-0.45*box_side_re[i], cutre_z[i]+0.45*box_side_re[i], panel_labels[i*2], fontsize=42, bbox=panel_label_bbox, ha="left", va="top")
    #cbax=None
    pt.plot.plot_colormap3dslice(
        vlsvobj=f,
        var="vg_b_vol",
        fluxrope=contourvalues,
        fluxropecolour=contourcolors,
        fluxropelinewidth=contourlinewidths,
        fluxropelinestyle="solid",
        operator="x",
        cbtitle="$B_x\,\mathrm{[nT]}$",
        boxre=[cutre_y[i]-0.5*box_side_re[i],cutre_y[i]+0.5*box_side_re[i],cutre_z[i]-0.5*box_side_re[i],cutre_z[i]+0.5*box_side_re[i]],
        cutpointre=cutre_x[i],
        normal='x',
        colormap='coolwarm',
        vmin=-100,vmax=100,
        symlog=True,
        usesci=None,
        vscale=1e9,
        step=times[i],
        streamlines = "vg_b_vol",
        streamlinedensity=2,
        streamlinethick=streamlinethick,
        streamlinecolor='k',
        streamlinearrowsize=streamlinearrowsize,
        axes=ax,
        nocb=True,
        scale=4,
        thick=thick)
    ax.yaxis.label.set_size(42)
    ax.hlines(-7, -20, -15, lw=7, color='k')
    ax=axes[i*3+1]
    cbax=axes[i*3+2]
    ax.text(cutre_x[i]-0.45*box_side_re[i], cutre_z[i]+0.45*box_side_re[i], panel_labels[i*2+1], fontsize=42, bbox=panel_label_bbox, ha="left", va="top")
    pt.plot.plot_colormap3dslice(
        vlsvobj=f,
        var="vg_b_vol",
        fluxrope=contourvalues,
        fluxropecolour=contourcolors,
        fluxropelinewidth=contourlinewidths,
        fluxropelinestyle="solid",
        operator="y",
        cbtitle="$B_\perp\,\mathrm{[nT]}$",
        boxre=[cutre_x[i]-0.5*box_side_re[i],cutre_x[i]+0.5*box_side_re[i],cutre_z[i]-0.5*box_side_re[i],cutre_z[i]+0.5*box_side_re[i]],
        cutpointre=cutre_y[i],
        normal='y',
        colormap='coolwarm',
        vmin=-100,vmax=100,
        symlog=True,
        usesci=None,
        vscale=1e9,
        step=times[i],
        streamlines = "vg_b_vol",
        streamlinedensity=2,
        streamlinethick=streamlinethick,
        streamlinecolor='k',
        streamlinearrowsize=streamlinearrowsize,
        axes=ax,
        cbaxes=cbax,
        scale=4,
        thick=thick,
        noylabels=True)
    ax.hlines(-7, -20-i*10, -15-i*10, lw=7, color='k')

for ax in axes[:-3]:
    ax.set_xlabel("")

axes[-3].xaxis.label.set_size(42)
axes[-2].xaxis.label.set_size(42)


fig.savefig("FHA_makkara_slices.png", bbox_inches='tight')

pl.close()

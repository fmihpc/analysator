#!/usr/bin/env python
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Circle, Wedge

import pylab as pl
import numpy as np
#import scipy
import sys
import pytools as pt
from variable import get_data, get_name, get_units

#plt.xkcd(scale=1.5, length=100, randomness=3)

Re = 6.371e+6 # Earth radius in m

if len(sys.argv)==9:
    time_start= int(sys.argv[1])
    time_stop = int(sys.argv[2])
    time_cadence = int(sys.argv[3])
    time_stride = int(sys.argv[4])
    cutoff_start= float(sys.argv[5])
    cutoff_stop = float(sys.argv[6])
    cutoff_cadence = float(sys.argv[7])
    cutoff_stride = int(sys.argv[8])
    
else:
    sys.stderr.write("Usage: volume.py <starting_index> <final_index+1> <cadence> <starting_cutoff> <final_cutoff+eps> <cadence> \n")
    sys.exit()


#cutoff_list = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
cutoff_list = np.arange(cutoff_start, cutoff_stop, cutoff_cadence)




inputLocation="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"
outputLocation="./"
outputfile = outputLocation + "FHA_fluxrope_volume.png"

thick = 3
thickseq = (1, 2, 4, 5, 6)
styleseq = ("dashed", "solid", "solid", "solid", "solid")
fontsize = 25

times=np.arange(time_start,time_stop,time_cadence)

time_cutoff_volume = np.zeros((len(times), len(cutoff_list)))


for t,time in enumerate(times):
    timefull=str(time).rjust(7, '0')
    file_name=inputLocation+"bulk1."+timefull+".vlsv"
    print("Extracting timestep "+str(time)+" s")
    f=pt.vlsvfile.VlsvReader(file_name=file_name)

    cid=f.read_variable("CellID")
    vg_fluxrope = f.read_variable("vg_fluxrope")[cid.argsort()]
    for c,cutoff in enumerate(cutoff_list):
        multiplier = np.copy(vg_fluxrope)
        multiplier[np.argwhere(multiplier > cutoff)] = 0
        multiplier[np.argwhere(multiplier != 0)] = 1
        time_cutoff_volume[t, c] = np.sum(np.prod(f.get_cell_dx(cid[cid.argsort()]), axis=-1) * multiplier)



fig = pl.figure(figsize=(24, 12), dpi=72)
gs = fig.add_gridspec(12, 24)

ax_time = fig.add_subplot(gs[:, :12])
ax_cutoff = fig.add_subplot(gs[:, 12:])

for c,cutoff in enumerate(cutoff_list):
    if c % cutoff_stride != 0:
        continue
    if abs(cutoff - 7) < 0.001:
        ax_time.plot(times, time_cutoff_volume[:, c] / Re**3, label="$R_\mathrm{cutoff}="+f"{cutoff:.0f}"+"\,R_\mathrm{C}$", lw=3, c='k')
    else:
        ax_time.plot(times, time_cutoff_volume[:, c] / Re**3, label="$R_\mathrm{cutoff}="+f"{cutoff:.0f}"+"\,R_\mathrm{C}$", lw=2)

for t,time in enumerate(times):
    if t % time_stride != 0:
        continue
    ax_cutoff.plot(cutoff_list, time_cutoff_volume[t, :] / Re**3, label="$t="+f"{time}"+"\,\mathrm{s}$")

ax_cutoff.axvline(7.0, c='k', lw=3)
ax_cutoff.set_xticks([2,4,6,7,8,10])

ax_time.legend(fontsize=fontsize, labelspacing=0.2)
ax_cutoff.legend(fontsize=fontsize, labelspacing=0.2)
#axes.set_xlim([-yminmax / Re,yminmax / Re])
#axes.set_ylim([-zminmax / Re,zminmax / Re])
#axes.set_aspect('equal')

for axis in ['top','bottom','left','right']:
    ax_time.spines[axis].set_linewidth(thick)
    ax_cutoff.spines[axis].set_linewidth(thick)
ax_time.xaxis.set_tick_params(width=thick,length=3*thick, which="major")
ax_time.yaxis.set_tick_params(width=thick,length=3*thick)
ax_cutoff.xaxis.set_tick_params(width=thick,length=3*thick, which="major")
ax_cutoff.yaxis.set_tick_params(width=thick,length=3*thick)

ax_cutoff.yaxis.tick_right()

for item in ax_time.get_xticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')
for item in ax_time.get_yticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')
for item in ax_cutoff.get_xticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')
for item in ax_cutoff.get_yticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')

ax_time.set_xlabel('$\mathrm{Time}\ [s]$',fontsize=1.2*fontsize,weight='black')
ax_cutoff.set_xlabel('$R_\mathrm{cutoff}\ [R_\mathrm{c}]$',fontsize=1.2*fontsize,weight='black')
ax_time.set_ylabel('$\mathrm{Flux\ rope\ detection\ volume}\ [R_\mathrm{E}^3]$',fontsize=1.2*fontsize,weight='black')

ax_time.xaxis.set_minor_locator(AutoMinorLocator(2))
ax_time.xaxis.set_tick_params(which="minor", width=0)
ax_cutoff.xaxis.set_minor_locator(AutoMinorLocator(2))
ax_cutoff.xaxis.set_tick_params(which="minor", width=0)


panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)
ax_time.text(1100, 2900, "(a)", fontsize=42, bbox=panel_label_bbox, ha="left", va="top")
ax_cutoff.text(10, 2900, "(b)", fontsize=42, bbox=panel_label_bbox, ha="right", va="top")

#ax_cutoff.set_xscale('log')
#ax_cutoff.set_yscale('log')


fig.savefig(outputfile, bbox_inches="tight")


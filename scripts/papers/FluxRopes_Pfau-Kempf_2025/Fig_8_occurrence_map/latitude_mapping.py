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

cutoff = 7.0

slice_list_re = np.asarray((-110,-70,-10,0,4,12))
tail_mpause_index = 2

yminmax = 35 * Re
zminmax = 15 * Re

if len(sys.argv)==4:
    start= int(sys.argv[1])
    stop = int(sys.argv[2])
    cadence = int(sys.argv[3])
else:
    sys.stderr.write("Usage: latitude_mapping.py <starting_index> <final_index+1> <cadence> \n")
    sys.exit()


inputLocation="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"
outputLocation="./"
outputfile = outputLocation + "latitude_mapping.png"

cmap = matplotlib.colormaps['plasma']
ncolors = len(slice_list_re) - 1
colorseq = cmap(np.linspace(0, 1, ncolors))
lobe_color = "xkcd:greenish grey"
lobe_alpha = 0.5
#colorseq = ("black", "xkcd:cool grey", "cyan", "xkcd:saffron", "xkcd:radioactive green", "xkcd:blood orange")#, "xkcd:barney")
thick = 3
thickseq = (1, 2, 4, 5, 6)
styleseq = ("dashed", "solid", "solid", "solid", "solid")
fontsize = 25

times=np.arange(start,stop,cadence)

f=pt.vlsvfile.VlsvReader(file_name=inputLocation+"bulk1."+str(start).rjust(7,'0')+".vlsv")
cid=f.read_variable("CellID")
cumulated_data = f.read_variable("vg_fluxrope")[cid.argsort()]
cumulated_data[np.argwhere(cumulated_data > cutoff)] = 0


[xsizefg, ysizefg, zsizefg] = f.get_fsgrid_mesh_size()
xsizefg = int(xsizefg)
ysizefg = int(ysizefg)
zsizefg = int(zsizefg)
[xminfg, yminfg, zminfg, xmaxfg, ymaxfg, zmaxfg] = f.get_fsgrid_mesh_extent()
cellsizefg = (xmaxfg-xminfg)/xsizefg

slice_list_indices = np.int64(np.floor(((slice_list_re * Re) - xminfg) / cellsizefg))
ymin_index = np.int64((-yminmax - yminfg) / cellsizefg)
ymax_index = np.int64(( yminmax - yminfg) / cellsizefg)+1
zmin_index = np.int64((-zminmax - zminfg) / cellsizefg)
zmax_index = np.int64(( zminmax - zminfg) / cellsizefg)+1

[XmeshXY,YmeshXY] = np.meshgrid(np.linspace(-yminmax / Re,yminmax / Re,num=ymax_index-ymin_index),np.linspace(-zminmax / Re,zminmax / Re,num=zmax_index-zmin_index))


for time in times[1:]:
    timefull=str(time).rjust(7, '0')
    file_name=inputLocation+"bulk1."+timefull+".vlsv"
    print("Extracting timestep "+str(time)+" s")
    f=pt.vlsvfile.VlsvReader(file_name=file_name)

    cid=f.read_variable("CellID")
    vg_fluxrope = f.read_variable("vg_fluxrope")[cid.argsort()]
    vg_fluxrope[np.argwhere(vg_fluxrope > cutoff)] = 0
    cumulated_data += vg_fluxrope

vg_to_fg_map = f.map_vg_onto_fg(sorted=True)
cumulated_data=cumulated_data[vg_to_fg_map]

vg_rho = f.read_variable("proton/vg_rho")[cid.argsort()]
vg_rho = vg_rho[vg_to_fg_map]


fig = pl.figure(figsize=(yminmax / Re / 2 + 4, zminmax / Re / 2), dpi=72)
axes = plt.gca()


for i,index in enumerate(slice_list_indices[:-1]):
    print(i,index)
    if i == tail_mpause_index:
        axes.contourf(XmeshXY,YmeshXY, vg_rho[slice_list_indices[2],ymin_index:ymax_index,zmin_index:zmax_index].T, [0,4e5], colors=[lobe_color], alpha=lobe_alpha)#, linewidths=[thick], linestyles=["--"])
    axes.contour(XmeshXY, YmeshXY, np.sum(cumulated_data[slice_list_indices[i]:slice_list_indices[i+1],ymin_index:ymax_index,zmin_index:zmax_index], axis=0).T, [1], colors=[colorseq[i]], linewidths=[thickseq[i]], linestyles=[styleseq[i]])

#circle = Circle((18, 0), 35.5, facecolor=None, fill=False, edgecolor='k', linewidth=6)
#axes.add_artist(circle)

color_lines = []
legend_labels = []
for i,color in enumerate(colorseq):
    if i == tail_mpause_index:
        color_lines.append(mpatches.Patch(color=lobe_color, alpha=lobe_alpha))
        legend_labels.append(r"$n_\mathrm{p}\leq0.4\,\mathrm{cm}^{-3}$"+"\n"+r"$ x="+str(slice_list_re[i])+"\,R_\mathrm{E}$"+"\n"+r"$ t="+str(times[-1])+"\,\mathrm{s}$")
    color_lines.append(mlines.Line2D([], [], color=color, linewidth=thickseq[i], linestyle=styleseq[i]))
    legend_labels.append("$x \in ["+str((slice_list_re[i]))+";"+str(int(slice_list_re[i+1]))+"]\,R_\mathrm{E}$")


axes.legend(color_lines, legend_labels, ncols=1, fontsize=fontsize, loc='upper right', bbox_to_anchor=(1.37, 1.03))
axes.set_xlim([-yminmax / Re,yminmax / Re])
axes.set_ylim([-zminmax / Re,zminmax / Re])
axes.set_aspect('equal')

for axis in ['top','bottom','left','right']:
    axes.spines[axis].set_linewidth(thick)
axes.xaxis.set_tick_params(width=thick,length=3*thick, which="major")
axes.yaxis.set_tick_params(width=thick,length=3*thick)

for item in axes.get_xticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')
for item in axes.get_yticklabels():
    item.set_fontsize(fontsize)
    item.set_fontweight('black')

axes.set_xlabel('$y [R_\mathrm{E}]$',fontsize=1.2*fontsize,weight='black')
axes.set_ylabel('$z [R_\mathrm{E}]$',fontsize=1.2*fontsize,weight='black')
#axes.set_title("$\mathrm{Occurrence~of~flux~ropes~over}~t=["+str(int(start))+";"+str(times[-1])+"]\,\mathrm{s}$", fontsize=1.5*fontsize)
axes.set_title("$t=["+str(int(start))+";"+str(times[-1])+"]\,\mathrm{s}$", fontsize=1.5*fontsize)

axes.xaxis.set_minor_locator(AutoMinorLocator(2))
axes.xaxis.set_tick_params(which="minor", width=0)

axes.grid(color="white", which="both", axis="both", lw=thick)
axes.set_axisbelow(True)
axes.set_facecolor("xkcd:light grey")

fig.savefig("FHA_fluxrope_latitude_mapping.png", bbox_inches="tight")


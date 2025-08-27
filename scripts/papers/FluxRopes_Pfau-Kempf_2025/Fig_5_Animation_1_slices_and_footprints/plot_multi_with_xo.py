import pytools as pt
import sys, os, socket
import numpy as np
import plot_colormap3dslice_withcontours
import plot_ionosphere_with_footpoints

import matplotlib.pyplot as plt
import matplotlib.path as mpath

#plt.xkcd(scale=1.5, length=100, randomness=3)

var="vg_v"
operator="z"

contourvariable="vg_beta_star"
contouroperator="pass"
contourvalues=[0.5]
contourcolors=["black"]
contourlinewidths=[8]

vmin = -400
vmax = 400
symlog = None
lin = 1
colormap="coolwarm"
#colormap="RdYlBu_r"
vscale = 1e-3
usesci = False

isymlog = None
ilin = 1
ivscale = 1
icontourlinewidth=15

boxre=(-2, 12, -14, 14)

nPlots = 4
step = 3
height = (boxre[3]-boxre[2])
width = (boxre[1]-boxre[0])
margin = 2
cbwidth = 14

scale = 10
thick = 6
dpi = 50

panel_labels = ("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)", "(m)")
panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)

star = mpath.Path.unit_regular_star(6)
circle = mpath.Path.unit_circle()
# concatenate the circle with an internal cutout of the star
cut_star = mpath.Path(
    vertices=np.concatenate([circle.vertices, star.vertices[::-1, ...]]),
    codes=np.concatenate([circle.codes, star.codes]))

fluxropelevel = float(sys.argv[4])
fluxropemarker = "o"
#fluxropemarkeredgecolor = "xkcd:dirty purple"
#fluxropemarkerfacecolor = "xkcd:saffron"
fluxropemarkeredgecolor = "xkcd:swamp"
fluxropemarkerfacecolor = "xkcd:swamp"
fluxropemarkersize = 250*scale
fluxropelinewidths = 0.6*scale
fluxropemaxcurvatureradius = 20e10

xmarker = "$\mathrm{X}$"
xmarkeredgecolor = "xkcd:velvet"
xmarkerfacecolor = "xkcd:velvet"
xmarkersize = 300*scale
xlinewidths = 1*scale

omarker = "s"
omarkeredgecolor = "None"
#omarkerfacecolor = "xkcd:dark blue grey"
omarkerfacecolor = "xkcd:saffron"
omarkersize = 70*scale
olinewidths = 0.8*scale

start = int(sys.argv[2])
end   = int(sys.argv[3])
if start == end:
  end=start+1

open_flux = np.load("open_flux.npy")

for id in np.arange(start, end):
  index = str(id).zfill(7)

  if id < 1000:
    ivmin = 0.49 # 0
    ivmax = 0.5  # 1
  else:
    ivmin = 0.9 # for contour # 1.49 # for colormap #1 open
    ivmax = 1.1 # for contour # 1.5 # for colormap #2 closed

  print(index)

  f=pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"+sys.argv[1]+"."+index+".vlsv")
  fxo=pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/pysidecar_sdf_"+sys.argv[1]+"."+index+".vlsv")
  outputdir="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/visualizations/movies/vg_v_cusp_fate/fluxrope_"+str(fluxropelevel)+"/"
  if start+1 == end:
    outputfile="FHA_multimap_y_v_"+str(fluxropelevel)+"_figure_"+index+".png"
  else:
    outputfile="FHA_multimap_y_v_"+str(fluxropelevel)+"_"+index+".png"

#  if(os.path.exists(outputdir+outputfile)):
#    continue

  title="t="+str(id)+" s"

  fig = plt.figure(figsize=(nPlots*(width+1*margin)+cbwidth, height + 2*width + 4*margin), dpi=dpi) # X, Y
  gs = fig.add_gridspec(height + 2*width + 4*margin, nPlots*(width+1*margin)+cbwidth) # Y, X
  cbaxes = []
  cbaxhint = []
  nocb = []
  noylabels = []
  cbaxes.append(fig.add_subplot(gs[margin:height+margin, -cbwidth+margin:-cbwidth+margin+2]))
#  cbaxes.append(fig.add_subplot(gs[:, -4]))
  cbaxhint.append(None)
#  cbaxhint.append("vectmap")
  nocb.append(None)
#  nocb.append(None)
  noylabels.append(None)
  for i in range(nPlots - 1):
    cbaxes.append(None)
    cbaxhint.append(None)
    nocb.append(1)
    noylabels.append(1)

  time_label_bbox = dict(boxstyle='square,pad=0.1', fc='white', ec='black', alpha=0.85, lw=thick)
  fig.suptitle(title, x=0.1, y=0.91, ha="left", va="top", fontsize=12*scale, bbox=time_label_bbox)

  for i in range(nPlots):
    # (left, bottom, width, height)
    ax1    = fig.add_subplot(gs[margin:height+margin,i*(width+1*margin)+margin:i*(width+1*margin)+width+margin])
    if start+1 == end:
      ax1.text(boxre[0]+1, boxre[3]-1, panel_labels[i], fontsize=12*scale, bbox=panel_label_bbox, ha="left", va="top")

    plot_colormap3dslice_withcontours.plot_colormap3dslice(
      vlsvobj=f,
      vlsvxo=fxo,
      outputdir=outputdir,
      nooverwrite=1,
      var=var,
      operator=operator,
      title="y="+str((2*i - nPlots + 1) * step * 0.5)+" RE",
      cbtitle="Vz (km/s)",
      draw=None,
      usesci=usesci,
      symlog=symlog,
      boxm=[],
      boxre=boxre,
      colormap=colormap,
      run="FHA",
      nocb=nocb[i],
      internalcb=None,
      wmark=None,
      wmarkb=True,
      axisunit=None,
      thick=thick,
      scale=scale,
      tickinterval=None,
      noborder=False,
      noxlabels=None,
      noylabels=noylabels[i],
      vmin=vmin,
      vmax=vmax,
      lin=lin,
      external=None,
      expression=None,
      vscale=vscale,
      pass_vars=None,
      pass_times=None,
      pass_full=None,
      fsaved=None,
      vectors=None, #"vg_poynting",
      vectordensity=10000000,
      vectorcolormap='gray',
      vectorsize=1.0,
      streamlines="vg_b_vol",
      streamlinedensityx=1,
      streamlinedensityy=1.5,
      streamlinecolor='black', 
      streamlinethick=8.0,
      axes=ax1,
      cbaxes=cbaxes[i],
      cbaxhint=cbaxhint[i],
      normal='y',
      cutpoint=None,
      cutpointre=(2*i - nPlots + 1) * step * 0.5,
      contourvariable=contourvariable,
      contouroperator=contouroperator,
      contourvalues=contourvalues,
      contourcolors=contourcolors,
      contourlinewidths=contourlinewidths,
      fluxropelevel=fluxropelevel,
      fluxropemarker=fluxropemarker, fluxropemarkeredgecolor=fluxropemarkeredgecolor, fluxropemarkerfacecolor=fluxropemarkerfacecolor,
      fluxropemarkersize=fluxropemarkersize,
      fluxropelinewidths=fluxropelinewidths,
      fluxropemaxcurvatureradius=fluxropemaxcurvatureradius,
      xmarker=xmarker, xmarkeredgecolor=xmarkeredgecolor, xmarkerfacecolor=xmarkerfacecolor,
      xmarkersize=xmarkersize,
      xlinewidths=xlinewidths,
      omarker=omarker, omarkeredgecolor=omarkeredgecolor, omarkerfacecolor=omarkerfacecolor,
      omarkersize=omarkersize,
      olinewidths=olinewidths)

  #directions=("North", "South")

  ax1 = fig.add_subplot(gs[height+3*margin:height+3*margin+2*width,margin:2*width+2*margin])

  if start+1 == end:
    f2=pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"+sys.argv[1]+".0001200.vlsv")
  else:
    f2=None

  plot_ionosphere_with_footpoints.plot_ionosphere(
    vlsvobj=f,
    outputdir=outputdir,
    nooverwrite=1,
    var=None,
    operator=None,
    contourvar="ig_openclosed",
    f2=f2,
    contourcmap="Greys",
    contourlevels=0,
    contourvmin=ivmin,
    contourvmax=ivmax,
    contourlinewidths=icontourlinewidth,
    minlatitude=66,
    cbtitle="Flux rope Z (RE)",
    title=" ",
    draw=None,
    usesci=False,
    symlog=isymlog,
    colormap="viridis",
    run="FHA",
    nocb=True,
    internalcb=None,
    wmark=None,
    wmarkb=True,
    thick=thick,
    scale=scale,
    lin=ilin,
    vmin=-10,
    vmax=10,
    vscale=ivscale,
    axes=ax1,
    cbaxes=None,
    viewdir=1)
  if start+1 == end:
    ax1.text(0.1, 0.48, panel_labels[nPlots], fontsize=12*scale, bbox=panel_label_bbox, ha="left", va="top", transform=fig.transFigure)

  ax1 = fig.add_subplot(gs[height+3*margin:height+3*margin+2*width,2*width+3*margin:4*width+4*margin])
  axc = fig.add_subplot(gs[height+np.int16(1.4*width)+3*margin:,-cbwidth:])

  plot_ionosphere_with_footpoints.plot_ionosphere(
    vlsvobj=f,
    outputdir=outputdir,
    nooverwrite=1,
    var=None,
    operator=None,
    contourvar="ig_openclosed",
    f2=f2,
    contourcmap="Greys",
    contourlevels=0,
    contourvmin=ivmin,
    contourvmax=ivmax,
    contourlinewidths=icontourlinewidth,
    minlatitude=66,
    cbtitle="flux rope Z (RE)",
    title=" ",
    draw=None,
    usesci=False,
    symlog=isymlog,
    colormap="viridis",
    run="FHA",
    nocb=None,
    internalcb=None,
    wmark=None,
    wmarkb=True,
    thick=thick,
    scale=scale,
    lin=ilin,
    vmin=-10,
    vmax=10,
    vscale=ivscale,
    axes=ax1,
    cbaxes=axc,
    viewdir=-1)
  if start+1 == end:
    ax1.text(0.42, 0.48, panel_labels[nPlots+1], fontsize=12*scale, bbox=panel_label_bbox, ha="left", va="top", zorder=-1000000000000, transform=fig.transFigure)

  plt.savefig(outputdir+outputfile, dpi=dpi, pad_inches=0.5, bbox_inches='tight')

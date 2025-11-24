                 
restartcalls = [
# Input and output methods, nooverwrite
"pt.plot.plot_vdf(figsize=[5,4],filename=fileLocation+bulkname, outputdir=outputLocation+'/'+REPLACEINDEX+'_', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1, coordre=REPLACECOORDRE)",

# cellids, coordinates
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACECELLID)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACECOORDINATES)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACEMULTIPLECELLID)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACEMULTIPLECOORDINATES)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACEMULTIPLECOORDRE)",


# Thickness, scale
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, thick=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, scale=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0, coordre=REPLACECOORDRE)",

# Tick interval
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1000, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=500, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=0.5,axisunit=6, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1,axisunit=6, coordre=REPLACECOORDRE)",

# msec musec titles
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='msec', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='musec', coordre=REPLACECOORDRE)",

# B vector
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",

# Zoom and units
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, box=[-2e6,2e6,-2e6,2e6], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, box=[-2e6,2e6,-2e6,2e6],axisunit=0, coordre=REPLACECOORDRE)",

# Watermarks
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmarkb=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NE', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NW', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SE', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SW', coordre=REPLACECOORDRE)",

# Biglabels
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='A', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='B', biglabloc=0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='C', biglabloc=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='D', biglabloc=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='E', biglabloc=3, coordre=REPLACECOORDRE)",

# title, axes, noborders
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title=r'$\mathcal{Title}$ and so forth $\odot$', cbtitle=r'$\mathcal{Color}$', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noxlabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,nocb=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1, coordre=REPLACECOORDRE)",

# slicethick
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=4, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1e3, coordre=REPLACECOORDRE)",

# cellsize
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=4, coordre=REPLACECOORDRE)",

# fmin, fmax, setThreshold
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14,fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-20, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-15, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=0, coordre=REPLACECOORDRE)",

# colormaps 
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='nipy_spectral', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='jet', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated_r', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='viridis', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='plasma', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='magma', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='bwr', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='PuOr', coordre=REPLACECOORDRE)",

# cbulk, center, bvector
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, cbulk=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[-7e5,0,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[2e5,2e5,2e5], coordre=REPLACECOORDRE)",

# wflux
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wflux=1, coordre=REPLACECOORDRE)",

# directions
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,0,5], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,1,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[-1,0,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[1,1,1], coordre=REPLACECOORDRE)",
]
nonrestartcalls = [
# Input and output methods, nooverwrite
"pt.plot.plot_vdf(figsize=[5,4],filename=fileLocation+bulkname, outputdir=outputLocation+'/'+REPLACEINDEX+'_', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1, coordre=REPLACECOORDRE)",

# cellids, coordinates
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACECELLID)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACECOORDINATES)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACEMULTIPLECELLID)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACEMULTIPLECOORDINATES)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACEMULTIPLECOORDRE)",


# Thickness, scale
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, thick=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, scale=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0, coordre=REPLACECOORDRE)",

# Tick interval
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1000, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=500, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=0.5,axisunit=6, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1,axisunit=6, coordre=REPLACECOORDRE)",

# msec musec titles
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='msec', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='musec', coordre=REPLACECOORDRE)",

# B vector
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",

# Zoom and units
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, box=[-2e6,2e6,-2e6,2e6], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, box=[-2e6,2e6,-2e6,2e6],axisunit=0, coordre=REPLACECOORDRE)",

# Watermarks
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmarkb=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NE', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NW', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SE', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SW', coordre=REPLACECOORDRE)",

# Biglabels
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='A', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='B', biglabloc=0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='C', biglabloc=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='D', biglabloc=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='E', biglabloc=3, coordre=REPLACECOORDRE)",

# title, axes, noborders
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title=r'$\mathcal{Title}$ and so forth $\odot$', cbtitle=r'$\mathcal{Color}$', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noxlabels=1,noborder=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,nocb=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1, coordre=REPLACECOORDRE)",

# slicethick
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=0, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=4, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1e3, coordre=REPLACECOORDRE)",

# cellsize
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=0.5, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=2, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=4, coordre=REPLACECOORDRE)",

# fmin, fmax, setThreshold
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14,fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-20, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-15, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=0, coordre=REPLACECOORDRE)",

# colormaps 
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='nipy_spectral', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='jet', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated_r', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='viridis', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='plasma', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='magma', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='bwr', coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='PuOr', coordre=REPLACECOORDRE)",

# cbulk, center, bvector
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, cbulk=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[-7e5,0,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[2e5,2e5,2e5], coordre=REPLACECOORDRE)",

# wflux
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, wflux=1, coordre=REPLACECOORDRE)",

# directions
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,0,5], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,1,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[-1,0,0], coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[1,1,1], coordre=REPLACECOORDRE)",
]

multipopcalls = [
"pt.plot.plot_vdf(figsize=[5,4],vlsvobj=f, run=verifydir+REPLACEINDEX, pop='REPLACEPOP', coordre=REPLACECOORDRE)"]

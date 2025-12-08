
v5restartcalls=[]

v5nonrestartcalls=[
#outputs, nooverwrite
"pt.plot.plot_ionosphere(filename=fileLocation+bulkname, outputdir=outputLocation,run=REPLACEINDEX)",
"pt.plot.plot_ionosphere(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",
"pt.plot.plot_ionosphere(vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1)",
#colormaps

"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='nipy_spectral')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='jet')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated_r')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='viridis')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='plasma')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='magma')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='bwr')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='PuOr')",

#colorbar
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, nocb=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, internalcb=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, symlog=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, symlog=0)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, lin=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, log=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, cbtitle='cbtitle')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, usesci=0)",

# Thickness, scale
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=0.5)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=5.0)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=0.5)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol',lin=1, vscale=1e-3)",


# Watermarks,title
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, wmarkb=1)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NE')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NW')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SE')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SW')",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, title='title_test')",

#vmin,vmax,Symmetric
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX, symmetric=True,vmin=2,vmax=40)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,vmin=2,vmax=40)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,vmin=0)",
"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,vscale=1e9)",
#operators
#"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,var='ig_electrontemp',op='x')",
#"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,var='ig_electrontemp',op='y')",
#"pt.plot.plot_ionosphere(vlsvobj=f, run=verifydir+REPLACEINDEX,var='ig_electrontemp',op='z')",

#variables
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_fac')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_fac',viewdir=-1)",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_fac',viewdir=10)",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_fac',absolute=1)",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_cellarea')",
#"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_deltaphi')", (Causes failure)
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_electrontemp')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_latitude')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_openclosed')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_potential')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_precipitation')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_rhon')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_sigmah')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_sigmap')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_sigmaparallel')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_source')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_residual')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_p')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_pp')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_z')",
"pt.plot.plot_ionosphere(vlsvobj=f,run=verifydir+REPLACEINDEX,var='ig_zz')",

]

v5multipopcalls=[
]


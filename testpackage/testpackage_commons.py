import analysator as pt
import sys, os
import numpy as np
import traceback
import inspect
import logging
import re

datalocation = "/wrk-vakka/group/spacephysics/vlasiator"
runs = []
runs.append( { 'name': 'FHA',
                 'verifydir': '/FHA/', 
                 'fileLocation': datalocation+'/3D/FHA/bulk1/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap3dslice','plot_ionosphere'],
                 'pops': ['avgs'],
                 'time': 1000,
                 'singletime': False,
                 'skipped_vars':{'plot_ionosphere':'vg_'},
                 'filename': None, #restart file
                 'manualcall':False,
                 'nosubpops': True, # backstreaming / non-backstreaming
                 'vlasiator5': True,
                 'cavitonparams': [6.6e6,2.64e6,4.e-9,10] } )

runs.append( { 'name': 'BCQ',
                'verifydir': '/BCQ/', 
                'fileLocation': '/wrk-vakka/group/spacephysics/vlasiator/2D/BCQ/bulk/',
                'fluxLocation': None,
                'singletime': False,
                'pops': ['avgs'],
                'funcs': ['plot_colormap','plot_vdf'],
                'manualcall':False,
                'time': 1600,
                'skipped_vars':None,
                'filename': None,
                'vlasiator5':False,
                'nosubpops':False,
                'cavitonparams': [6.6e6,2.64e6,4.e-9,10]
                  } )
                     
# Custom expression function                                                               
def exprMA_cust(exprmaps, requestvariables=False):
    if vlasiator5 is True:
        reqs = ['vg_va']
    else:
        reqs = ['va']
    if requestvariables==True:
        return reqs
    # exprmaps is a dictionary of numpy arrays
    # Each array has 2 dimensions [xsize, ysize]
    # or 3 dimensions [xsize, ysize, components]
    # here the function returns the M_A with a preset bulk velocity
    custombulkspeed=1500000. # m/s
    va = exprmaps[reqs[0]][:,:]
    MA = np.ma.divide(custombulkspeed,va)
    return MA

# Second example of a more involved custom expression function
def expr_cav_cust(pass_maps, requestvariables=False):
    if vlasiator5 is True:
        reqs= ['vg_rho','vg_b_vol','vg_beta']
    else:
        reqs= ['rho', 'B', 'beta']
    if requestvariables==True:
        return reqs
    # pass_maps is a dictionary of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of dictionaries
    # Each dictionary contains an entry called 'dstep'
    # Which contains the relative time step position, i.e.
    # a value of 0 indicates no time offset.

    # This custom expression returns a map with values of
    # either 0 (solar wind), 0.5 (caviton), or 1.0 (SHFA), calculated against
    # time-averaged background values. This version doesn't do intelligent checks for the
    # format of the incoming data.
    if type(pass_maps) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("expr_cav_cust expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()

    # Multiple time steps were found
    ntimes = len(pass_maps)
    dsteps = [x['dstep'] for x in pass_maps]
    curri = dsteps.index(0)
    thesemaps = pass_maps[curri]
    
    thisrho = np.ma.masked_less_equal(thesemaps[reqs[0]][:,:], 0)
    thisB = np.ma.masked_less_equal(thesemaps[reqs[1]],0)
    thisbeta = np.ma.masked_less_equal(thesemaps[reqs[2]],0)
    thisBmag = np.linalg.norm(thisB, axis=-1)
        
    avgrho = np.zeros(np.array(thisrho.shape))
    avgBmag = np.zeros(np.array(thisrho.shape))
    # avgbeta = np.zeros(np.array(thisrho.shape))
    
    for i in range(ntimes):
        if i==curri: # Exclude current frame from background value averaging
            continue
        nowmaps = pass_maps[i]
        print(nowmaps['dstep'])
        avgrho = np.add(avgrho, nowmaps[reqs[0]])
        avgBcalc = np.linalg.norm(nowmaps[reqs[1]], axis=-1)
        avgBmag = np.add(avgBmag, avgBcalc)
        # avgbeta = np.add(avgbeta, nowmaps[2])

    avgrho = np.divide(np.ma.masked_less_equal(avgrho,0), np.array([ntimes-1]))
    avgBmag = np.divide(np.ma.masked_less_equal(avgBmag,0), np.array([ntimes-1]))
    #avgbeta = np.divide(np.ma.masked_less_equal(avgbeta,0), np.array([ntimes-1]))

    rhoratio = np.ma.divide(thisrho, avgrho)
    Bmagratio = np.ma.divide(thisBmag, avgBmag)
    #betaratio = np.divide(thisbeta, avgbeta)

    rhoratioreq = 0.9
    bmagratioreq = 0.9
    betashfareq = level_beta_SHFA

    # Calculations using masked arrays proved problematic so a less-than elegant method is used here.
    empty = np.zeros(np.array(thisrho.shape))+0.0
    half = empty + 0.5
    one = empty + 1.0
    
    caviton=np.zeros(empty.shape,dtype='float64')
    np.add(empty, one,out=caviton, where=(rhoratio<rhoratioreq))
    print("sum of cavitons rho ",caviton.sum())
    cavitonbmag=np.zeros(empty.shape,dtype='float64')
    np.add(caviton, one, out=cavitonbmag,where=(Bmagratio<bmagratioreq))
    print("sum of cavitons Bmag ",cavitonbmag.sum())
    shfa=np.zeros(empty.shape,dtype='float64')
    np.add(cavitonbmag, one,out=shfa, where=(thisbeta>betashfareq))
    print("sum of SHFA ",shfa.sum())

    combo=np.zeros(empty.shape,dtype='float64')
    np.add(empty, half,out=combo,where=(cavitonbmag>1.5))
    print("sum of combo ",combo.sum())
    combo2=np.zeros(empty.shape,dtype='float64')
    combo2 = np.add(empty, half,out=combo2, where=(shfa>2.5))
    print("sum of combo2 ",combo2.sum())
    combo3 = combo+combo2
    print("sum of combo3 ",combo3.sum())

    # Mask out anything that is inside the bow shock
    combo3 = np.ma.masked_where(thisrho>level_bow_shock, combo3)
    print("sum of combo3 upstream ",combo3.sum())

    return combo3


# Second example of a more involved custom expression function
def timesmooth(pass_maps):
    # pass_maps is a dictionary of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of dictionaries
    # Each dictionary contains an entry called 'dstep'
    # Which contains the relative time step position, i.e.
    # a value of 0 indicates no time offset.

    # consists of a single [ysize,xsize] array to smooth over time

    if type(pass_maps) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("timesmooth expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()
    ntimes = len(pass_maps)
    print("this many time steps ",ntimes)
    # Select first valid variable
    listofkeys = iter(pass_maps[0])
    while True:
        var = next(listofkeys)
        if var!="dstep": break
    
    # Multiple time steps were found
    avg = np.zeros(pass_maps[0][var].shape)
    for i in range(ntimes):
        avg = np.add(avg, pass_maps[i][var])
    return np.divide(avg, np.array([ntimes]))

# Helper function for drawing on existing panel
def extcontour(ax, XmeshXY,YmeshXY, extmaps, requestvariables=False):
    if vlasiator5 is True:
        reqs= ['vg_rho','vg_b_vol','vg_beta']
    else:
        reqs= ['rho', 'B', 'beta']
    if requestvariables==True:
        return reqs

    rho = extmaps[reqs[0]]
    beta = extmaps[reqs[2]]
    # take magnitude of B
    B = np.linalg.norm(extmaps[reqs[1]],axis=-1)

    # Colours to use
    color_cavitons = '#924900'
    color_SHFAs    = '#B66DFF'
    color_BS       = '#FFFF6D'

    # mask cavitons
    cavitons = np.ma.masked_greater_equal(B,level_B_caviton)
    cavitons.mask[rho > level_n_caviton] = True
    cavitons.fill_value = 0.
    cavitons[cavitons.mask == False] = 1.

    # mask SHFAs
    SHFAs = np.ma.masked_greater_equal(B,level_B_caviton)
    SHFAs.mask[rho > level_n_caviton] = True
    SHFAs.mask[beta < level_beta_SHFA] = True
    SHFAs.fill_value = 0.
    SHFAs[SHFAs.mask == False] = 1.
 
    # draw contours
    contour_shock = ax.contour(XmeshXY,YmeshXY,rho,[level_bow_shock],
                               linewidths=1.2, colors=color_BS)
    contour_cavitons = ax.contour(XmeshXY,YmeshXY,cavitons.filled(),[0.5], linewidths=1.5, colors=color_cavitons)
    contour_SHFAs = ax.contour(XmeshXY,YmeshXY,SHFAs.filled(),[0.5], linewidths=1.5, colors=color_SHFAs)

regularcalls=[
    # cellids, coordinates
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACECELLID)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACECOORDINATES)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACEMULTIPLECELLID)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACEMULTIPLECOORDINATES)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACEMULTIPLECOORDRE)",
]

nonrestartcalls = []


multipopcalls = []

v5restartcalls = [

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v')"

] 




v5nonrestartcalls = [


# Input and output methods, nooverwrite
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, outputdir=outputLocation+'/'+REPLACEINDEX+'_')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX)",
"pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1)",

# Thickness, scale
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, highres=True)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0, highres=True)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0, highres=True)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, highres=3)",

# Tick interval
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1)",

# msec musec titles
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='msec')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='musec')",

# Watermarks
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmarkb=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NE')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NW')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SE')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SW')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, Earth=True)",

# title, axes, noborders
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title=r'$\mathcal{Title}$ and so forth $\odot$', cbtitle=r'$\mathcal{Color}$')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noborder=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1,noborder=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noborder=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noxlabels=1,noborder=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,nocb=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1,highres=True)",

# Variables, operators, colormaps, usesci, lin, vscale
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', colormap='nipy_spectral')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', colormap='nipy_spectral', vscale=1e9)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', op='x', colormap='bwr')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', op='z', colormap='bwr')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', colormap='warhol',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', colormap='warhol',lin=1, vscale=1e-3)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='x', colormap='PuOr')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='y', colormap='PuOr',symlog=0, usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='z', colormap='PuOr',symlog=0, usesci=0)",

# Zoom and units
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=3)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6],axisunit=0)",

# Externals and expressions
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, pass_vars=['vg_rho','vg_b_vol','vg_beta'])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, pass_vars=['vg_va'], vmin=1, vmax=20,lin=1,usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, boxre=[0,30,-15,15], vmin=1, vmax=20,lin=1,usesci=0)", #why error with plot3d? keyeerror vg_va
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['vg_rho','vg_b_vol','vg_beta'],lin=1,colormap='bwr',usesci=0)",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['vg_rho'], boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['vg_beta'])",


# Everything at once
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15], expression=exprMA_cust, vmin=1, vmax=20,lin=1,usesci=0, fsaved=1, fluxfile=fluxLocation+fluxname)",

# Streamlines, vectors
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='x',vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='y',vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='z',vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectordensity=400, boxre=[-10,10,5,50],vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectordensity=20, boxre=[-10,10,5,50],vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol', vectorcolormap='viridis')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectordensity=400)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectordensity=20, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='x',vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='y',vectordensity=200)", 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='z',vectordensity=200)", 

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinecolor='black')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinecolor='gray')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=0.5, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=2, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5,streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2,streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5,streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2,streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2, boxre=[-10,10,5,50])",


# More data reducers
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_rho')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_nonthermal',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_thermal',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel_nonthermal', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular_nonthermal',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel_thermal', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular_thermal',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pressure')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pdyn',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pdynx',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_temperature')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",

#colorbar
#not yet in in master see PR #359
#"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cb_horizontal=True)",

#AMR, fsaved
#Does not work currently, fix for the AMR contour is in PR #364
#"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, amr=0.1,amrlinestyles='dashed',amrcolours='red',amrlinewidths=1"),
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved='red')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=0.5)",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[14,0],pass_vars=['vg_rho'])"

]

v5multipopcalls= [
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_v_parallel', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_v_perpendicular',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pressure')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pdyn',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pdynx',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_temperature')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular_nonthermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular_thermal')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy_nonthermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy_thermal', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_parallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_perpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_thermalvelocity',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_blocks')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_gyrotropy')"
]

manualcalls=[
    "pt.plot.colormap3dslice(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_gyrotropy')"
]
#keys: v5bulk,v5restart,bulk,restart,v5multipop,multipop

# For handier debugging, uncomment these to overwrite call lists and include only relevant calls
# regularcalls = []
# nonrestartcalls = ["pt.plot.plot_colormap3dslice(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['rho','B','beta'],lin=1,colormap='bwr',usesci=0)","pt.plot.plot_colormap3dslice(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
# ]
# multipopcalls = []
# v5regularcalls = []
# v5nonrestartcalls = []
# v5multipopcalls = []

def call_replace(call,func,skipped_vars):
    call=call.replace('REPLACEFUNC',func)

    #Get the arguments of the call
    args = re.search(r'\((.+)\)',call).group(1)
    args = [x[0] or x[1] or x[2] for x in re.findall(r'(\w+=[^,(\[]+)|(\w+=\(.+\))|(\w+=\[.+?\])',args)]
    #Get parameters of the func
    function_pars=inspect.getfullargspec(eval("pt.plot."+func)).args
    #Remove args that are not present as parameters for the func
    args_out=[]
    for arg in args:
        if arg:
            if skipped_vars and func in skipped_vars.keys() and skipped_vars[func] in arg.split("=")[1]:
                continue
            elif arg.split("=")[0] in function_pars:
                args_out.append(arg)
            else:
                logging.warning(f"Argument {arg} removed from call {call}")
        
    call=call[:call.rfind("(")+1]+",".join(args_out)+")"

    return call



#construct calls
calls = []
callrunids = []
callrunindex = []
funcids=[]
offset=0
for i,run in enumerate(runs):
    # bulk and restart files
    vlasiator5 = run['vlasiator5']
    filename = run['filename']
    fileLocation = run['fileLocation']
    singletime = run['singletime']
    nosubpops = run['nosubpops']
    fluxLocation = run['fluxLocation']
    functions = run['funcs']
    skipped_vars=run['skipped_vars']


    callindex = 0
    if run['manualcall']:
        for call in manualcalls:
            funcids.append(-1)
            callrunids.append(i)
            calls.append(call)
            callrunindex.append(callindex)
            callindex += 1

    for j,func in enumerate(functions):
        calls_in=v5restartcalls if vlasiator5 else regularcalls
        for call in calls_in:
            funcids.append(j)
            call=call_replace(call,func,skipped_vars)
            if func == "plot_vdf" and ('cellids' not in call and 'coordinates' not in call and 'coordsre' not in call):
                continue
            if not filename is None: 
                if vlasiator5:
                    call = call.replace("var='vg_v'","var='vg_restart_v'")
                else:
                    call = call.replace("var='V'","var='restart_V'")

            if call not in calls:
                callrunids.append(i)
                calls.append(call)
                callrunindex.append(callindex)
                callindex += 1

        # non-restart files
        if filename is None:
            #non restart calls
            calls_in=v5nonrestartcalls if vlasiator5 else nonrestartcalls
            for call in calls_in:
                funcids.append(j)
                call=call_replace(call,func,skipped_vars)
                # Skip flux function calls if no flux files
                if "flux" in call and fluxLocation is None:
                    continue
                # skip time integration if only one file available
                if "pass_times" in call and singletime:
                    continue
                # thermal / non-thermal subpopulations
                if vlasiator5 and (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                    continue
                elif (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                    continue
                if call not in calls:
                    callrunids.append(i)
                    calls.append(call)
                    callrunindex.append(callindex)
                    callindex += 1
                
            #multipop calls
            calls_in=v5multipopcalls if vlasiator5 else multipopcalls
            for pop in run['pops']:
                if pop != 'avgs':
                    for call in calls_in:
                        funcids.append(j)
                        call=call_replace(call,func,skipped_vars)
                        # Skip flux function calls if no flux files
                        if "flux" in call and fluxLocation is None:
                            continue
                        # skip time integration if only one file available
                        if "pass_times" in call and singletime:
                            continue
                        # thermal / non-thermal subpopulations
                        if vlasiator5 and (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                            continue
                        elif (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                            continue
                        mpcall = call.replace('REPLACEPOP',pop)
                        if call not in calls:
                            callrunids.append(i)
                            calls.append(call)
                            callrunindex.append(callindex)
                            callindex += 1
            

nteststot = len(callrunids)

# How many jobs? 
jobcount=int(sys.argv[1])
jobcurr=int(sys.argv[2])
increment = int(nteststot/jobcount)
remainder = nteststot - jobcount * increment
start=jobcurr * increment
end=start + increment
# Remainder frames are divvied out evenly among tasks
if jobcurr < remainder:
    start = start + jobcurr
    end = end + jobcurr + 1
else:
    start = start + remainder
    end = end + remainder


# Perform call
for j in range(start,end):
    # Calculate which run
    jrun = callrunindex[j]
    runid = callrunids[j]
    call = calls[j]

    funcid=funcids[jrun] #offset due to manualcalls
    
    runname = runs[runid]['name']
    
    func = runs[runid]['funcs'][funcid]
    verifydir = func+runs[runid]['verifydir']

    fileLocation = runs[runid]['fileLocation']
    fluxLocation = runs[runid]['fluxLocation']
    pops = runs[runid]['pops']
    time = runs[runid]['time']
    filename = runs[runid]['filename']
    vlasiator5 = runs[runid]['vlasiator5']
    singletime = runs[runid]['singletime']

    level_bow_shock = runs[runid]['cavitonparams'][0]
    level_n_caviton = runs[runid]['cavitonparams'][1]
    level_B_caviton = runs[runid]['cavitonparams'][2]
    level_beta_SHFA = runs[runid]['cavitonparams'][3]

    verifydir=os.path.join("testpackage_run",verifydir)
    outputLocation=os.path.join(pt.plot.defaultoutputdir,verifydir)

    if func == "plot_ionosphere" and "var='vg_" in call:
        print(f"skipped call {j}")
        continue

    # Source data files
    if filename is None:
        if '2D' not in fileLocation:
            bulkname = "bulk1."+str(time).rjust(7,'0')+".vlsv"
        else:
            bulkname = "bulk."+str(time).rjust(7,'0')+".vlsv"
    else:
        bulkname = filename
    if 'fluxprefix' in runs[runid]:
        fluxname = runs[runid]['fluxprefix']+str(time).rjust(7,'0')+".bin"
    else:
        fluxname = "flux."+str(time).rjust(7,'0')+".bin"


    if run=="ABC":
        bulkname = "distributions."+str(time).rjust(7,'0')+".vlsv"

           
    call = call.replace('REPLACEPREVINDEX',"'"+str(jrun-1).rjust(4,'0')+"'")
    call = call.replace('REPLACEINDEX',"'"+str(jrun).rjust(4,'0')+"'")
    call = call.replace('REPLACETIME',"'"+str(time)+"'")

    call = call.replace('REPLACECELLID','1')
    call = call.replace('REPLACECOORDRE','[10,0,0]')
    call = call.replace('REPLACECOORDINATES','[6.371e7,0,0]')
    call = call.replace('REPLACEMULTIPLECELLID','[1,51,101]')
    call = call.replace('REPLACEMULTIPLECOORDRE','[[10,0,0],[15,0,0],[20,0,0]]')
    call = call.replace('REPLACEMULTIPLECOORDINATES','[[6.371e7,0,0],[9.5565e7,0,0],[12.742e7,0,0]]')


    # Many different plots
    print(j, runid, jrun, call,fileLocation+bulkname)
    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    try:
        exec(call)
    except Exception as e:
        print("----------------------------\nFAILURE DURING CALL ",j," \n```\n"+call+"```\n", repr(e))
        
        traceback.print_exc()
        print("END TRACE for call",j,"\n----------------------------")

#add way to specify which function to test 
#add a way to add expections to variables etc easily.
#currenlty does multiple calls

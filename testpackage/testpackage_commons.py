import analysator as pt
import sys, os
import numpy as np
import traceback
import inspect
import logging
import re
import argparse

#add manualcals filtering or something since rightnow it tries to do it for 3D data? does it work?

argp=argparse.ArgumentParser(
    prog='Analysator Testpackage',
    description='Outputs test plots'
)
argp.add_argument("jobcount",type=int)
argp.add_argument("jobindex",type=int)

argp.add_argument('funcs',type=str,help="function/list of functions to test, if none give does all.",nargs='*')
cmd_args=argp.parse_args()
funcs_to_use=cmd_args.funcs
if "pass" in funcs_to_use:
    quit()

datalocation = "/wrk-vakka/group/spacephysics/vlasiator"
runs = []


#required args for functions, lists are handled as OR statements, tuples within lists as AND
#add a way to add required args automatically

#list of tuples, first element is the list of required arguments and second is the defaults if argument is not found, leaving it as None skips defaults

#maybe add overide to this in runs append
required_args ={
    "plot_vdf":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"])],
    "plot_vdf_profiles":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"])],
    "plot_isosurface":[([("surf_step","surf_var")],["surf_step=10","surf_var='vg_rho'"])]

}

'''
runs.append( { 'name': 'FHA',
                 'verifydir': '/FHA/', 
                 'fileLocation': datalocation+'/3D/FHA/bulk1/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap3dslice','plot_ionosphere','plot_isosurface'],
                 'pops': ['avgs'],
                 'time': 1000,
                 'singletime': False,

                  #Uses 'in' operator for the values of the inner dict, so skipping arg completely can be done with {'var':''}
                 'skipped_args':{'plot_ionosphere':{'var':'vg_'},
                                'plot_colormap3dslice':{'expression':['expr_cav_cust','exprMA_cust']}
                                },
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
                'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles'],
                'manualcall':False,
                'time': 1600,
                'skipped_args':None,
                'filename': None,
                'vlasiator5':False,
                'nosubpops':False,
                'cavitonparams': [6.6e6,2.64e6,4.e-9,10]
                  } )
'''                     
runs.append( { 'name': 'BCQr',
                 'verifydir': '/BCQr/', 
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles'],
                 'fileLocation': '/wrk-vakka/group/spacephysics/vlasiator/2D/BCQ/restart/',
                 'pops': ['avgs'],
                'fluxLocation': None,
                 'singletime': True, # neighboring bulk files not available
                 'time': 0,
                 'skipped_args':None,
                 'manualcall':False,
                 'vlasiator5': False,
                 'nosubpops': False, # thermal / non-thermal
                 'filename': 'restart.0001361.vlsv',
                'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
'''
runs.append( { 'name': 'BGA',
                 'verifydir': '/BGA/', 
                 'fileLocation': datalocation+'/2D/BGA/zero_ehall_layers_23/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles'],
                 'pops': ['proton'],
                 'skipped_args':{'plot_vdf':{'normal':''}},
                 'time': 380,
                 'manualcall':False,
                 'singletime': True, # neighboring bulk files not available
                 'filename': None,
                 'vlasiator5': True,
                 'nosubpops': True, # thermal / non-thermal
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
                     
runs.append( { 'name': 'BFD',
                 'verifydir': '/BFD/', 
                 'fileLocation': datalocation+'/2D/BFD/bulk/',
                 'fluxLocation': datalocation+'/2D/BFD/fluxfunction/',
                 'fluxprefix': 'bulk.',
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles'],
                 'skipped_args':{'plot_vdf':{'normal':''}},
                 'pops': ['proton','helium'],
                 'time': 2000,
                 'singletime': False,
                 'filename': None,
                  'manualcall':False,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
'''
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

'''

    DISCLAIMER

    Add to the end of the list always so the output file names remain the same



'''




#This can be v4 or v5, these both get added into nonrestartcalls and v5nonrestartcalls, here for cleanliness sake
agnostic_call = [


# Input and output methods, nooverwrite
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, outputdir=outputLocation+'/'+REPLACEINDEX+'_')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX)",
"pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1)",
"pt.plot.REPLACEFUNC(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",

# cellids, coordinates
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACECELLID)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACECOORDINATES)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACEMULTIPLECELLID)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACEMULTIPLECOORDINATES)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACEMULTIPLECOORDRE)",



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
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=0.5,axisunit=6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1,axisunit=6)",

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



# Overplots and flux lines
"pt.plot.REPLACEFUNC(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxfile=fluxLocation+fluxname, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, step=REPLACETIME, fluxthick=0.5, fluxlines=10)",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, fluxthick=5, fluxlines=2)",

# Vscale
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e3, vmax=7.e6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e-3, vmax=7.e0, vscale=1e-6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e6, vmax=7.e9, vscale=1e3)",


# Zoom and units
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=3)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6],axisunit=0)",


# slicethick (Mostly for vdf,vdf_prof,vdfdiff)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=0, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=2, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=4, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1e3, coordre=REPLACECOORDRE)",
# cellsize (same)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=0.5, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=2, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=4, coordre=REPLACECOORDRE)",
# fmin, fmax, setThreshold (same)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14,fmax=1.e-12, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-20, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-15, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=0, coordre=REPLACECOORDRE)",
# Biglabels(vdf,vdfdiff)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='A', coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='B', biglabloc=0, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='C', biglabloc=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='D', biglabloc=2, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='E', biglabloc=3, coordre=REPLACECOORDRE)",
# cbulk, center, bvector(same)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, cbulk=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[-7e5,0,0], coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[2e5,2e5,2e5], coordre=REPLACECOORDRE)",

# wflux(same)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wflux=1, coordre=REPLACECOORDRE)",
# directions(about same)
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,0,5], coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,1,0], coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[-1,0,0], coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1,coordswap=1, coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[1,1,1], coordre=REPLACECOORDRE)",



# colormaps 
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='nipy_spectral')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='jet')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated_r')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='viridis')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='plasma')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='magma')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='bwr')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='PuOr')",

#colorbar
#not yet in in master see PR #359
#"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cb_horizontal=True)",

#AMR, fsaved
#Does not work currently, fix for the AMR contour is in PR #364
#"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, amr=0.1,amrlinestyles='dashed',amrcolours='red',amrlinewidths=1"),
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved='red')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=0.5)",
]


#This is here so the agnostic calls are first and the order doesnt change if you add something to nonrestartcalls
nonrestartcalls=[]
v5nonrestartcalls=[]
restartcalls=[]
v5restartcalls=[]
nonrestartcalls.extend(agnostic_call)
v5nonrestartcalls.extend(agnostic_call)
restartcalls.extend(agnostic_call)
v5restartcalls.extend(agnostic_call)

v5restartcalls.extend([

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v')"

])

restartcalls.extend([
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='V')"

])

nonrestartcalls.extend([

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', colormap='hot_desaturated')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', colormap='hot_desaturated', vscale=1e3)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='x', colormap='RdBu',symlog=0, usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='y', colormap='RdBu',symlog=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='z', colormap='RdBu',symlog=0, usesci=0)",

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='RhoBackstream', colormap='jet')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='beta',lin=1, usesci=0, colormap='viridis',vmax=50)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='MA',lin=1,usesci=0,vmin=2,vmax=40, colormap='inferno_r')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Mms',lin=1,usesci=0, vmin=0, colormap='magma')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='va',lin=1,usesci=0, colormap='magma_r')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vms',lin=1, colormap='plasma_r')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vs',lin=1, colormap='viridis_r')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='max_v_dt', vscale=1e6)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='max_v_dt', vscale=1e3)",

"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Temperature', colormap='plasma', vscale=1e-6,lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pressure', vscale=1e9)",

# Symlog and vscale
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-9)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-12)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=0,vscale=1e9)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1,vscale=1e9)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-3,vscale=1e9)",

# Externals and expressions
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, pass_vars=['rho','B','beta'])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, pass_vars=['va'], vmin=1, vmax=20,lin=1,usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, boxre=[0,30,-15,15], vmin=1, vmax=20,lin=1,usesci=0)",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['rho','B','beta'],lin=1,colormap='bwr',usesci=0)",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['rho'], boxre=[0,30,-15,15])",
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['beta'])",

# Everything at once
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15], expression=exprMA_cust, vmin=1, vmax=20,lin=1,usesci=0, fsaved=1, fluxfile=fluxLocation+fluxname)",

# Streamlines, vectors
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B', vectorcolormap='viridis')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B', vectorcolormap='magma')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=20)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=400)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=20, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=400, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20,vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400,vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50],vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50],vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50],vectorsize=1.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20,vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400,vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50],vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50],vectorsize=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50],vectorsize=0.5)",


"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinecolor='black')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinecolor='gray')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=0.5, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=2, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5,streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2,streamlinethick=2)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5,streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2,streamlinethick=0.5)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5, boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', boxre=[-10,10,5,50])",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2, boxre=[-10,10,5,50])",

# More data reducers
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VBackstream',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VNonBackstream',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallel', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicular',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallelBackstream', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicularBackstream',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallelNonBackstream', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicularNonBackstream',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pressure')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallelBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicularBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallelNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicularNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pdyn',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pdynx',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Temperature')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallelBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicularBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallelNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicularNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='beta')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vBeam',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vBeamRatio')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Thermalvelocity',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Blocks')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='gyrotropy')"])







v5nonrestartcalls.extend([
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
"pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[14,0],pass_vars=['vg_rho'])"

])


v5multipopcalls=[
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


multipopcalls=[
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, pop='REPLACEPOP', coordre=REPLACECOORDRE)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/VParallel', op='magnitude',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/VPerpendicular',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pressure')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallelBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicularBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallelNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicularNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pdyn',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pdynx',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Temperature')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallelBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicularBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallelNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicularNonBackstream')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverPar', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverParBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverParNonBackstream', vmin=0.1, vmax=10)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/beta')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaParallel')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpendicular')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Thermalvelocity',lin=1)",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Blocks')",
"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/gyrotropy')"]

manualcalls=[
    "pt.plot.colormap3dslice(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_gyrotropy')",
    #plot_vdf manual calls
    "pt.plot.plot_vdf(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)"

]
#keys: v5bulk,v5restart,bulk,restart,v5multipop,multipop

# For handier debugging, uncomment these to overwrite call lists and include only relevant calls
# restartcalls = []
# nonrestartcalls = ["pt.plot.plot_colormap3dslice(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['rho','B','beta'],lin=1,colormap='bwr',usesci=0)","pt.plot.plot_colormap3dslice(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
# ]
# multipopcalls = []
# v5restartcalls = []
# v5nonrestartcalls = []
# v5multipopcalls = []

def call_replace(call,func,skipped_args):
    #This is kind of scuffed maybe

    call=call.replace('REPLACEFUNC',func)

    #Get the arguments of the call
    args = re.search(r'\((.+)\)',call).group(1)
    args = [x[0] or x[1] or x[2] for x in re.findall(r'(\w+=[^,(\[]+)|(\w+=\(.+\))|(\w+=\[.+?\])',args)]
    named_parameters=[arg.split("=")[0] for arg in args]

    #Get parameters of the func
    function_pars=inspect.getfullargspec(eval("pt.plot."+func)).args
    #Remove args that are not present as parameters for the func
    args_out=[]
    #check that all required func args are set
    if required_args and func in required_args.keys():
        for required_tuple in required_args[func]:
            required_params=required_tuple[0]
            default_params=required_tuple[1]
            check=False
            for param in required_params:
                if any((all(r in named_parameters for r in param),(param in named_parameters))):
                    check=True
                    break
            if not check:
                #print("REQUIRED",call,named_parameters,required_params)

                #Add parameters if there are default_params
                if default_params:
                    for param in default_params:
                        if param not in named_parameters:
                            args_out.append(param)
  #                      print("ADDED",param,call)
                else:
                    #print("NOT ADDED",call)
                    return None


    #add execption for tuple?
    for arg in args:
        if arg:
            if skipped_args and func in skipped_args.keys():
                skipped_args_dict=skipped_args[func]
                call_args=arg.split("=")
                if call_args[0] in skipped_args_dict.keys():
                    if type(skipped_args_dict[call_args[0]])==str and skipped_args_dict[call_args[0]] in call_args[1]:
                        continue
                    elif type(skipped_args_dict[call_args[0]])==list and any(arg_skip in call_args[1] for arg_skip in skipped_args_dict[call_args[0]]):
                        continue
            if arg.split("=")[0] in function_pars:
                args_out.append(arg)
            #else:
            #    logging.warning(f"Argument {arg} removed from call {call}")
    
    if not args_out:
        return None
    call=call[:call.rfind("(")+1]+",".join(args_out)+")"
#    print("BEFORE",call)
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

    if not funcs_to_use:
        functions = run['funcs']
    else:
        functions = list(set(run['funcs']) & set(funcs_to_use))
        if not functions:
            continue

    skipped_args=run['skipped_args']


    callindex = 0
    if run['manualcall']:
        for call in manualcalls:
            funcids.append(-1)
            callrunids.append(i)
            calls.append(call)
            callrunindex.append(callindex)
            callindex += 1

    for j,func in enumerate(functions):
        if filename is not None:
            calls_in=v5restartcalls if vlasiator5 else restartcalls
            for call in calls_in:
                call=call_replace(call,func,skipped_args)
                if not call:
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
                    funcids.append(j)


        # non-restart files
        if filename is None:
            #non restart calls
            calls_in=v5nonrestartcalls if vlasiator5 else nonrestartcalls
            for call in calls_in:

                call=call_replace(call,func,skipped_args)
                if not call:
                    continue
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
                    funcids.append(j)
                
            #multipop calls
            calls_in=v5multipopcalls if vlasiator5 else multipopcalls
            for pop in run['pops']:
                if pop != 'avgs':
                    for call in calls_in:

                        call=call_replace(call,func,skipped_args)
                        if not call:
                            continue
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
                        call = call.replace('REPLACEPOP',pop)
                        if call not in calls:
                            callrunids.append(i)
                            calls.append(call)
                            callrunindex.append(callindex)
                            callindex += 1
                            funcids.append(j)            

nteststot = len(callrunids)

# How many jobs? 
jobcount=cmd_args.jobcount
jobcurr=cmd_args.jobindex

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

    funcid=funcids[j] 
    
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
        #print(call)
    except Exception as e:
        print("----------------------------\nFAILURE DURING CALL ",j," \n```\n"+call+"```\n", repr(e))
        
        traceback.print_exc()
        print("END TRACE for call",j,"\n----------------------------")

#add way to specify which function to test 
#add a way to add expections to variables etc easily. (DONE)
#currenlty does multiple calls (Fixed with list but still needs better implementation as we waste bit of time going through multiple things)
#add manual calls (DONE)
#why spend time going through all calls on all threads?
#v5 vdf?? (post 2019 are v5)


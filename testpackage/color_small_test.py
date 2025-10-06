import analysator as pt
import sys, os
import numpy as np
import traceback

datalocation = "/wrk/group/spacephysics/vlasator"
runs = []
runs.append( { 'name': 'ABC',
                 'verifydir': 'testpackage_colormap/ABC/', 
                 'fileLocation': datalocation+'/2D/ABC/bulk/',
                 'fluxLocation': datalocation+'/2D/ABC/flux/',
                 'pops': ['avgs'],
                 'time': 1000,
                 'singletime': False,
                 'filename': None,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [6.6e6,2.64e6,4.e-9,10] } )
runs.append( { 'name': 'BCQ',
                 'verifydir': 'testpackage_colormap/BCQ/', 
                 'fileLocation': datalocation+'/2D/BCQ/bulk/',
                 'fluxLocation': None, #datalocation+'/2D/BCQ/flux/', # missing on Vorna
                 'pops': ['avgs'],
                 'time': 2000,
                 'singletime': False,
                 'filename': None,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
runs.append( { 'name': 'BFD',
                 'verifydir': 'testpackage_colormap/BFD/', 
                 'fileLocation': datalocation+'/2D/BFD/bulk/',
                 'fluxLocation': datalocation+'/2D/BFD/fluxfunction/',
                 'fluxprefix': 'bulk.',
                 'pops': ['proton','helium'],
                 'time': 2000,
                 'singletime': False,
                 'filename': None,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
# runs.append( { 'name': 'BCQr',
#                  'verifydir': 'testpackage_colormap/BCQr/', 
#                  'fileLocation': datalocation+'/2D/BCQ/',
#                  'fluxLocation': None,
#                  'pops': ['avgs'],
#                  'time': 0,
#                  'singletime': False,
#                  'filename': 'restart.0001361.vlsv',
#                  'nosubpops': False, # backstreaming / non-backstreaming
#                  'vlasiator5': False,
#                  'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
# runs.append( { 'name': 'AFC',
#                  'verifydir': 'testpackage_colormap/AFC/', 
#                  'fileLocation': datalocation+'/2D/AFC/production_halfres/',
#                  'fluxLocation': datalocation+'/2D/AFC/production_halfres/flux/',
#                  'pops': ['proton'],
#                  'time': 1000,
#                  'singletime': False,
#                  'filename': None,
#                  'nosubpops': False, # backstreaming / non-backstreaming
#                  'vlasiator5': False,
#                  'cavitonparams': [24.0e6,9.6e6,11.27e-9,10] } )
# runs.append( { 'name': 'AFCr',
#                  'verifydir': 'testpackage_colormap/AFCr/', 
#                  'fileLocation': datalocation+'/2D/AFC/production_halfres/',
#                  'fluxLocation': None,
#                  'pops': ['proton'],
#                  'time': 0,
#                  'singletime': False,
#                  'filename': 'restart.0000592.2019-04-19_07-58-12.vlsv',
#                  'nosubpops': False, # backstreaming / non-backstreaming
#                  'vlasiator5': False,
#                  'cavitonparams': [24.0e6,9.6e6,11.27e-9,10] } )
# runs.append( { 'name': 'BFDr',
#                  'verifydir': 'testpackage_colormap/BFDr/', 
#                  'fileLocation': datalocation+'/2D/BFD/restart/',
#                  'fluxLocation': None,
#                  'pops': ['avgs'],
#                  'time': 0,
#                  'singletime': False,
#                  'filename': 'restart.0001126.2018-06-03_21-34-16.vlsv',
#                  'nosubpops': False, # backstreaming / non-backstreaming
#                  'vlasiator5': False,
#                  'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
runs.append( { 'name': 'BGA',
                 'verifydir': 'testpackage_colormap/BGA/', 
                 'fileLocation': datalocation+'/2D/BGA/zero_ehall_layers_23/',
                 'fluxLocation': None,
                 'pops': ['proton'],
                 'time': 380,
                 'singletime': True, # neighboring bulk files not available
                 'filename': None,
                 'vlasiator5': True,
                 'nosubpops': True, # thermal / non-thermal
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
                     
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
    caviton = np.add(empty, one, where=(rhoratio<rhoratioreq))
    print("sum of cavitons rho ",caviton.sum())
    caviton = np.add(caviton, one, where=(Bmagratio<bmagratioreq))
    print("sum of cavitons Bmag ",caviton.sum())
    shfa = np.add(caviton, one, where=(thisbeta>betashfareq))
    print("sum of SHFA ",shfa.sum())

    combo = np.add(empty, half, where=(caviton>1.5))
    print("sum of combo ",combo.sum())
    combo2 = np.add(empty, half, where=(shfa>2.5))
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

regularcalls = []

nonrestartcalls = [
# Overplots and flux lines
"pt.plot.plot_colormap(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",
"pt.plot.plot_colormap(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1)",
"pt.plot.plot_colormap(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1, boxre=[-10,10,5,50])",
"pt.plot.plot_colormap(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxfile=fluxLocation+fluxname, boxre=[-10,10,5,50])",
"pt.plot.plot_colormap(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, step=REPLACETIME, fluxthick=0.5, fluxlines=10)",
"pt.plot.plot_colormap(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, fluxthick=5, fluxlines=2)",
]


multipopcalls = []

v5regularcalls = []

v5nonrestartcalls = []

v5multipopcalls = []
# For handier debugging, uncomment these to overwrite call lists and include only relevant calls
# regularcalls = []
# nonrestartcalls = ["pt.plot.plot_colormap(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['rho','B','beta'],lin=1,colormap='bwr',usesci=0)","pt.plot.plot_colormap(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
# ]
# multipopcalls = []
# v5regularcalls = []
# v5nonrestartcalls = []
# v5multipopcalls = []

# Construct test list
calls = []
callrunids = []
callrunindex = []
for i,run in enumerate(runs):
    # bulk and restart files
    vlasiator5 = run['vlasiator5']
    filename = run['filename']
    fileLocation = run['fileLocation']
    singletime = run['singletime']
    nosubpops = run['nosubpops']
    fluxLocation = run['fluxLocation']

    callindex = 0
    if vlasiator5:
        for call in v5regularcalls:
            if not filename is None: 
                call = call.replace("var='V'","var='restart_V'")
            callrunids.append(i)
            calls.append(call)
            callrunindex.append(callindex)
            callindex += 1
    else:
        for call in regularcalls:
            if not filename is None: 
                call = call.replace("var='V'","var='restart_V'")
            callrunids.append(i)
            calls.append(call)                
            callrunindex.append(callindex)
            callindex += 1
    # non-restart files
    if filename is None:
        if vlasiator5:
            for call in v5nonrestartcalls:
                # Skip flux function calls if no flux files
                if "flux" in call and fluxLocation is None:
                    continue
                # skip time integration if only one file available
                if "pass_times" in call and singletime:
                    continue
                # thermal / non-thermal subpopulations
                if (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                    continue
                callrunids.append(i)
                calls.append(call)
                callrunindex.append(callindex)
                callindex += 1
            for pop in run['pops']:
                if pop != 'avgs':
                    for call in v5multipopcalls:
                        # Skip flux function calls if no flux files
                        if "flux" in call and fluxLocation is None:
                            continue
                        # skip time integration if only one file available
                        if "pass_times" in call and singletime:
                            continue
                        # thermal / non-thermal subpopulations
                        if (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                            continue
                        mpcall = call.replace('REPLACEPOP',pop)
                        callrunids.append(i)
                        calls.append(mpcall)
                        callrunindex.append(callindex)
                        callindex += 1
        else:
            for call in nonrestartcalls:
                # Skip flux function calls if no flux files
                if "flux" in call and fluxLocation is None:
                    continue
                # skip time integration if only one file available
                if "pass_times" in call and singletime:
                    continue
                # thermal / non-thermal subpopulations
                if (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                    continue
                callrunids.append(i)
                calls.append(call)
                callrunindex.append(callindex)
                callindex += 1
            for pop in run['pops']:
                if pop != 'avgs':
                    for call in multipopcalls:
                        # Skip flux function calls if no flux files
                        if "flux" in call and fluxLocation is None:
                            continue
                        # skip time integration if only one file available
                        if "pass_times" in call and singletime:
                            continue
                        # thermal / non-thermal subpopulations
                        if (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                            continue
                        mpcall = call.replace('REPLACEPOP',pop)
                        callrunids.append(i)
                        calls.append(mpcall)
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

    runname = runs[runid]['name']
    verifydir = runs[runid]['verifydir']
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

    outputLocation=pt.plot.defaultoutputdir+verifydir
    
    # Source data files
    if filename is None:
        bulkname = "bulk."+str(time).rjust(7,'0')+".vlsv"
    else:
        bulkname = filename
    if 'fluxprefix' in runs[runid]:
        fluxname = runs[runid]['fluxprefix']+str(time).rjust(7,'0')+".bin"
    else:
        fluxname = "flux."+str(time).rjust(7,'0')+".bin"
           
    call = call.replace('REPLACEPREVINDEX',"'"+str(jrun-1).rjust(4,'0')+"'")
    call = call.replace('REPLACEINDEX',"'"+str(jrun).rjust(4,'0')+"'")
    call = call.replace('REPLACETIME',"'"+str(time)+"'")

    # Many different plots
    print(j, runid, jrun, call)
    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    try:
        exec(call)
    except Exception as e:
        print("----------------------------\nFAILURE DURING CALL ",j," \n```\n"+call+"```\n", repr(e))
        
        traceback.print_exc()
        print("END TRACE for call",j,"\n----------------------------")

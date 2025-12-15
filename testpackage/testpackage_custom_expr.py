import numpy as np

vlasiator5=None
level_bow_shock = None
level_n_caviton = None
level_B_caviton = None
level_beta_SHFA = None



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

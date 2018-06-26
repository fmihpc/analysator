# Default vmin and vmax values defined per variable and run
def loadrundefaults(run, var, op):
    vminuse=None
    vmaxuse=None

    if var == 'rho':
        vminuse = 1.e+6
        vmaxuse = 4.e+6
        if run == 'AEA' or run == 'ABA' or run == 'BCQ':
            vminuse = 8.e+5
            vmaxuse = 4.e+6
        elif run == 'AEC' or run == 'ABC':
            vminuse = 1.e+6
            vmaxuse = 4.e+6
        elif run == 'BEB':
            vminuse = 2.e+6
            vmaxuse = 2.e+7

    elif var == 'MA':
        vminuse = 1.0
        vmaxuse = 10.0

    elif var == 'Mms':
        vminuse = 1.0
        vmaxuse = 10.0

    elif var == 'B':        
        if op==None:
            vminuse = 3.e-9
            vmaxuse = 3.e-7
        else: #components
            vminuse = -0.15e-9
            vmaxuse = 0.15e-9

    elif var == 'E':
        if op==None:
            vminuse = 0.01
            vmaxuse = 0.1
        else: #components
            vminuse = -0.015
            vmaxuse = 0.015

    elif var == 'rhoBeam':
        vminuse = 1.e3
        vmaxuse = 1.e5

    elif var == 'V':
        if op==None:
            vminuse = 1.e5
            vmaxuse = 1.e6
        else: #components
            vminuse = -1.e6
            vmaxuse = 1.e6

    elif var == 'beta':
        vminuse = 0.8
        vmaxuse = 100

    elif var == 'temperature':
        vminuse = 0.5e6
        vmaxuse = 5.0e6

    elif var == 'vs':
        vminuse = 5e4
        vmaxuse = 5e5

    elif var == 'va':
        vminuse = 5e4
        vmaxuse = 5e5

    return(vminuse, vmaxuse)


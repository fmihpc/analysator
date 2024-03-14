import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path
#from scipy.interpolate import interp2d, griddata

#x: 1D or 2D array of positions
#t: 1D or 2D array of times
#var: 2D array of data, dimensions [nx, nt] (if x and t multidimensional, then assumed to have same dimensions as x and t)

# kwargs are passed to pcolormesh
# e.g. cmap=plasma, vmin=0, vmax=1d9, shading='auto', etc.

def keogram(x, t, var, filename='myplot.png', xlabel ='', ylabel ='', title = '', cbar_label = '', log=False, xlim=None, ylim=None, **kwargs):

    max_abs = np.nanmax([np.abs(var)])
    min_abs = np.nanmin([np.abs(var)])

    if (np.nanmin(var) < 0):
        if 'vmin' in kwargs.keys():
            vmin = kwargs['vmin']
        else:
            vmin = -max_abs
        if 'vmax' in kwargs.keys():
            vmax = kwargs['vmax']
        else:
            vmax = max_abs
    else:
        if 'vmin' in kwargs.keys():
            vmin = kwargs['vmin']
        else:
            vmin = min_abs
        if 'vmax' in kwargs.keys():
            vmax = kwargs['vmax']
        else:
            vmax = max_abs

    #if 'cmap' in kwargs.keys():
        #cmap=kwargs['cmap']
    #else:

    if 'cmap' not in kwargs.keys():
        if (np.nanmin(var) < 0):
            kwargs['cmap'] = 'bwr'
            #cmap='bwr'
        else:
            kwargs['cmap'] = 'plasma'
            #cmap='plasma'

    if log:
        if vmin <= 0:
        #if vmin <= 0:
    # Create a 'logarithmic' data normalization.
            linthresh = np.nanmax([min_abs, vmax / 5 ])  # linthresh must be > 0
            min_log_level = np.log10(min_abs)
            max_log_level = np.log10(max_abs)
            norm = mpl.colors.SymLogNorm(linthresh=linthresh, linscale=1, vmin=vmin, vmax=vmax, base=10)
        else:
            #if vmin == 0:
            #    vmin = np.nanmin(var[(var > 0)])
            norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        #norm = None

    print('keogram vmin ' + str(vmin)) 
    fig, ax = plt.subplots()
    #im = ax.pcolormesh(x, t, var, norm=norm, cmap = cmap, **kwargs)
    im = ax.pcolormesh(x, t, var, norm=norm, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
       #color bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
    mycbar.set_label(cbar_label)
    plt.tight_layout()
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()




# Assorted subroutines that are called by python programs in scripts/

import os
import numpy as np
from functools import wraps
from time import time
import hashlib
import re

global mu_0
mu_0 = 4e-7 * np.pi

try:
    from pyspedas.utilities.time_string import *
except:
    from pyspedas import time_string, time_double
    

'''
 show() creates uses display to display a plot created by tell(),
 which is saved in a default location
 the input, plt, is assumed to be a matplotlib plot

    ex. (on Turso login node, won't work on a compute note)
 import matplotlib.pyplot as plt
 from utils import *
 plt.plot([1,2,3])
 tell(plt)
 show()

Note that show (and tell?) won't work in an interactive session.

Another (more useful) option is to just type 'show' into the command line on turso

''' 

def get_bulklocation(run):
    # load data
    if run.upper() == 'EGI':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"
    elif run.upper() == 'EGL':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/"
    elif run.upper() == 'EGM':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGM/bulk/"
    elif run.upper() == 'EGN':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGM/"
    elif run.upper() == 'EGO':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGO/"
    elif run.upper() == 'EGP':
        location = "/wrk-vakka/users/horakons/vlasiator_data_files/"
        #location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGP/bulk1/"
        #location = "/wrk-vakka/group/spacephysics/vlasiator/3D/{}/bulk5/".format(run.upper())
    elif run.upper() == 'EGILIKE':
        location = "/wrk-vakka/group/spacephysics/vlasiator/temp/EGI_like/"
        #location = "/wrk-vakka/group/spacephysics/vlasiator/3D/{}/bulk5/".format(run.upper())
    elif run.upper() == 'EGILIKE2':
        location = "/wrk-vakka/users/ykempf/ionosphere/EGI/FAC_fgb/"
    elif run.upper() == 'FHA':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"
    elif run.upper() == 'FHAFGB':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk_with_fg_10/"
    elif run.upper() == 'FIA':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/FIA/bulk/"
    return location


def get_filename(run, fileIndex):
    # load data
    if run.upper() == 'EGI':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGL':
        filename = "bulk1.{}.{}.vlsv".format(run.lower(), str(fileIndex).zfill(7) )
    elif run.upper() == 'EGM':
        filename = "bulk.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGN':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGO':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGP':
        filename = "bulk1.egp.{}.vlsv".format(str(fileIndex).zfill(7) )
        #filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
        #filename = "bulk5.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGILIKE':   # test run
        filename = "bulk.{}.vlsv_fg".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGILIKE2':   # test run
        filename = "bulk.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'FHA':   # test run
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'FHAFGB':   # test run
        filename = "bulk_with_fg_10.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'FIA':   # test run
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    return filename


def timeseries(run, var_list, start, stop, filestem = None):     #, step?
    '''
        run: name of the Vlasiator run. ex.  'EGL', 'FHA'
        var_list: a list of variable names (strings)
            specify operators using 'variable.operator'
        start: first timestep
        stop: last timestep
        filestem: set this keyword to the filename, omitting the end 'XXXXXXX.vlsv'
    '''
    import pytools as pt
    if type(var_list) == str:
        var_list = [var_list]
    dct = {}
    for var in var_list:
        var_timeseries = []
        for fileIndex in range(start, stop+1):
            if filestem is None:
                filename = get_vlsvfile_fullpath(run, fileIndex)
            else:
                filename = filestem + str(fileIndex).zfill(7) + '.vlsv'
            f = pt.vlsvfile.VlsvReader(filename)
            if '.' in var:
                l = var.split('.')
                vartemp = l[0]; operator = l[1]
                var_timeseries.append(f.read_variable(vartemp, operator = operator))
            else:
                var_timeseries.append(f.read_variable(var))
        var_timeseries = np.array(var_timeseries)
        dct[var] = var_timeseries
    return dct


def get_vlsvfile_fullpath(run, fileIndex):
    return get_bulklocation(run) + get_filename(run, fileIndex)


def tell(plt):
   tmp_dir = '/wrk-vakka/users/horakons/temp/'
   plt.savefig(temp_dir+'plot.png', dpi=300, bbox_inches='tight')


def show():
   tmp_dir = '/wrk-vakka/users/horakons/temp/'
   os.system('display ' + tmp_dir + 'plot.png')


def cartesian_to_spherical(x, y, z):
    '''
    r > 0
    0 < theta < pi
    -pi < phi < pi
    all are assumed to be numpy arrays of equal dimensions

    returns:  r, theta, phi  [tuple]
    '''
    r = (x**2 + y**2 + z**2)**0.5
    theta = np.arccos( z / r )
    phi = np.arctan2(y, x)
    return r, theta, phi


def spherical_to_cartesian(r, theta, phi):
    '''
    r > 0
    0 < theta < pi
    -pi < phi < pi
    all are assumed to be numpy arrays of equal dimensions

    returns:  x, y, z   [tuple]
    '''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def cartesian_to_spherical_vector(vx, vy, vz, x, y, z):
    '''
    Convert cartesian vector(s) with coordinates (vx, vy, vz)
    at the position(s) theta, phi (note: position r does not affect the vector transformation)
    to spherical coordinates (v_r, v_theta, v_phi)

    dimensions of vx, vy, vz, x, y, z arrays must either match or be a single number
    '''
    r, theta, phi = cartesian_to_spherical(x, y, z)
    v_r =     vx * np.sin(theta) * np.cos(phi) + vy * np.sin(theta) * np.sin(phi) + vz * np.cos(theta)
    v_theta = vx * np.cos(theta) * np.cos(phi) + vy * np.cos(theta) * np.sin(phi) - vz * np.sin(theta)
    v_phi =  -vx * np.sin(phi) + vy * np.cos(phi)
    return v_r, v_theta, v_phi


def mkdir_path(path):
    filedir_list = path.split('/')
    filedir = path[:-len(filedir_list[-1])]
    if not(os.path.exists(filedir)):
         os.system('mkdir -p {}'.format(filedir))


def numcurl3d(inputarray, CELLSIZE_XYZ):
    # Assumes input array is of format [nx,ny,nz,3]
    # CELLSIZE_XYZ is 3-element array or list of grid spacings
    jac = numjacobian3d(inputarray, CELLSIZE_XYZ)
    # Output array is of format [nx,ny,nz,3]
    curl = np.zeros(inputarray.shape)
    curl[:,:,:,0] = jac[:,:,:,2,1]-jac[:,:,:,1,2]
    curl[:,:,:,1] = jac[:,:,:,0,2]-jac[:,:,:,2,0]
    curl[:,:,:,2] = jac[:,:,:,1,0]-jac[:,:,:,0,1]
    return curl


def numjacobian3d(inputarray, CELLSIZE_XYZ):
    # Assumes input array is of format [nx,ny,nz,3]
    # CELLSIZE_XYZ is 3-element array or list of grid spacings
    nx,ny,nz = inputarray[:,:,:,0].shape
    jac = np.zeros([nx,ny,nz,3,3])
    jac[:,:,:,0,0], jac[:,:,:,0,1], jac[:,:,:,0,2] = np.gradient(inputarray[:,:,:,0], CELLSIZE_XYZ[0])
    jac[:,:,:,1,0], jac[:,:,:,1,1], jac[:,:,:,1,2] = np.gradient(inputarray[:,:,:,1], CELLSIZE_XYZ[1])
    jac[:,:,:,2,0], jac[:,:,:,2,1], jac[:,:,:,2,2] = np.gradient(inputarray[:,:,:,2], CELLSIZE_XYZ[2])
    # Output array is of format [nx,ny,nz,3,3]
    #  :,:,component, derivativedirection
    # so for example: dAx/dx = :,:,:,0,0
    #                 DAy/dz = :,:,:,1,2
    return jac


def save(file,**kwargs):
    """
    Save the value of some data in a file.
    Usage: 
    a = 3; bvar = 'string'; test = [True, '300', 3.14]
    save('misdatos.pickle',a=a,b=bvar,test=test)
    """
    import pickle
    f=open(file,"wb")
    pickle.dump(kwargs,f,protocol=2)
    f.close


def restore(file):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pickle')    #datos is a dict with the saved variables
    """
    import pickle
    f=open(file,"rb")
    result = pickle.load(f)
    f.close
    return result


def timer(my_func, *args, **kwargs):
    '''
     timer is meant to be used as a python decorator.
     Use @timer before the function definition.
     Every the function is called afterwards,
     the runtime of the function will be displayed

      #Example 1:
    @timer
    def sq(x, keyw='test'):
        print(keyw)
        return(x**2)
    
    sq(3)

      #Example 2:
    len = timer(len)
    len([1,2,'asdf',None]) 
    '''
    @wraps(my_func)
    def wrapTheFunction(*args, **kwargs):
        if (my_func is not time) and (my_func is not print):
            t1 = time()
            output = my_func(*args, **kwargs)
            diff = time() - t1
            #print("{}(args={},kwargs={}) ran in {} seconds".format(my_func.__name__,args,kwargs, diff) )
            print("{}() ran in {} seconds".format(my_func.__name__, diff) )
            return output
    return wrapTheFunction


def sidecar(my_func, *args, **kwargs):
    '''
     sidecar is meant to be used as a python decorator.
     Use @sidecar before the function definition.
     When the decorated function is called, the the local sidecar/
     directory will be searched for a .pck (Pickle) file that
     contains the results of the function call for the given inputs
     If no such file is found, a new .pck file is created upon completion
     Basically behaves like a cache in permanent memory.
     
      #Example 1:
    @sidecar
    def sq(x, keyw='test'):
        print(keyw)
        return(x**2)
    
    sq(3)

      #Example 2:
    len = sidecar(len)
    len([1,2,'asdf',None]) 
    '''
    @wraps(my_func)
    def wrapTheFunction(*args, **kwargs):
        # Encode the signature as a string, then save it as a hash
        signature_string = "{}_{}_{}".format(my_func.__name__,args, str(sorted(kwargs))+str([kwargs[x] for x in sorted(kwargs)]))
        signature_string = re.sub(r'<.*?at.*?>', '', signature_string)   #get rid of references to specific memory addresses, that might change btw. calls
        signature_hash = hashlib.md5('{}'.format(signature_string).encode()).hexdigest()
        filename = "sidecar/{}/{}.pck".format(my_func.__name__,signature_hash) 
        mkdir_path(filename)  # make a sidecar/directory if it doesn't exist
        if os.path.exists(filename):
            print('restoring sidecar file for {}(): {}'.format(my_func.__name__, filename))
            tmp = restore(filename)
            output = tmp['output']
        else:
            print('saving sidecar file for {}(): {}'.format(my_func.__name__, filename) )       # note: if you try to pickle an I/O buffer, an exception will be thrown
            output = my_func(*args, **kwargs)
            save(filename, output=output)
        return output
    return wrapTheFunction




def save(file,**kwargs):
    """
    Save the value of some data in a file.
    Usage: 
    a = 3; bvar = 'string'; test = [True, '300', 3.14]
    save('misdatos.pickle',a=a,b=bvar,test=test)
    """
    import pickle
    f=open(file,"wb")
    pickle.dump(kwargs,f,protocol=2)
    f.close


def restore(file):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pickle')    #datos is a dict with the saved variables
           print(datos['test'])                  #e.g.
    """
    import pickle
    f=open(file,"rb")
    result = pickle.load(f)
    f.close
    return result



def t_string(t, epoch = False):
    if epoch:
        return time_string(t/1000 -  719528. * 24.* 3600.)
    else:
        return time_string(t)


def t_double(date, epoch = False):
    if epoch:
        return (time_float(date) + 719528. * 24.* 3600.) * 1000
    else:
        return time_float(date)


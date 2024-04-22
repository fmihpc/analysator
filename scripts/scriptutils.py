# Assorted subroutines that are called by python programs in scripts/

import os
import numpy as np
from functools import wraps
from time import time

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


def get_vlsvfile_fullpath(run, fileIndex):
    return get_bulklocation(run) + get_filename(run, fileIndex)


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


def mkdir_path(path):
    filedir_list = path.split('/')
    filedir = path[:-len(filedir_list[-1])]
    if not(os.path.exists(filedir)):
         os.system('mkdir -p {}'.format(filedir))


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



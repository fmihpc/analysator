# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import pytools as pt
import numpy as np
from rotation import rotateTensorToVector

PLANE = 'XY'
# or alternatively, 'XZ'

CELLSIZE = 300000.0 # cell size
DT = 0.5 # time step

def inplane(inputarray):
    # Assumes input array is of format [nx,ny,3]
    if PLANE=='XY':
        inputarray[:,:,2] = np.zeros(inputarray[:,:,2].shape)
    elif PLANE=='XZ':
        inputarray[:,:,1] = np.zeros(inputarray[:,:,1].shape)
    else:
        print("Error defining plane!")
        return -1
    return inputarray

def inplanevec(inputarray):
    # Assumes input array is of format [nx,ny,3]
    if PLANE=='XY':
        return inputarray[:,:,0:1]
    elif PLANE=='XZ':
        return inputarray[:,:,0:3:2]
    else:
        print("Error defining plane!")
        return -1

def numjacobian(inputarray):
    # Assumes input array is of format [nx,ny,3]
    nx,ny = inputarray[:,:,0].shape
    jac = np.zeros([nx,ny,3,3])
    if PLANE=='XY':
        jac[:,:,0,0], jac[:,:,0,1] = np.gradient(inputarray[:,:,0], CELLSIZE)
        jac[:,:,1,0], jac[:,:,1,1] = np.gradient(inputarray[:,:,1], CELLSIZE)
        jac[:,:,2,0], jac[:,:,2,1] = np.gradient(inputarray[:,:,2], CELLSIZE)
    elif PLANE=='XZ':
        jac[:,:,0,0], jac[:,:,0,2] = np.gradient(inputarray[:,:,0], CELLSIZE)
        jac[:,:,1,0], jac[:,:,1,2] = np.gradient(inputarray[:,:,1], CELLSIZE)
        jac[:,:,2,0], jac[:,:,2,2] = np.gradient(inputarray[:,:,2], CELLSIZE)
    else:
        print("Error defining plane!")
        return -1
    # Output array is of format [nx,ny,3,3]
    return jac

def numgradscalar(inputarray):
    # Assumes input array is of format [nx,ny]
    nx,ny = inputarray.shape
    grad = np.zeros([nx,ny,3])
    if PLANE=='XY':
        grad[:,:,0],grad[:,:,1] = np.gradient(inputarray, CELLSIZE)
    elif PLANE=='XZ':
        grad[:,:,0],grad[:,:,2] = np.gradient(inputarray, CELLSIZE)
    else:
        print("Error defining plane!")
        return -1
    # Output array is of format [nx,ny,3]
    return grad

def numdiv(inputarray):
    # Assumes input array is of format [nx,ny,3]
    jac = numjacobian(inputarray)
    # Output array is of format [nx,ny]
    return jac[:,:,0,0] + jac[:,:,1,1] + jac[:,:,2,2]

def numcrossproduct(inputvector1, inputvector2):
    # assumes inputvectors are of shape [nx,ny,3]
    # in fact nothing special here
    # Output array is of format [nx,ny,3]
    return np.cross(inputvector1,inputvector2) 

def numcurl(inputarray):
    # Assumes input array is of format [nx,ny,3]
    jac = numjacobian(inputarray)
    # Output array is of format [nx,ny,3]
    curl = np.zeros(inputarray.shape)
    curl[:,:,0] = jac[:,:,2,1]-jac[:,:,1,2]
    curl[:,:,1] = jac[:,:,0,2]-jac[:,:,2,0]
    curl[:,:,2] = jac[:,:,1,0]-jac[:,:,0,1]
    return curl

def numvecdottensor(inputarray, inputtensor):
    # Assumesinput array is of format [nx,ny,3]
    # assumes inputtensor is of shape [nx,ny,3,3]
    result = np.zeros(inputarray.shape)
    result[:,:,0] = (inputarray[:,:,0]*inputtensor[:,:,0,0] # is this indexing correcT?
                  + inputarray[:,:,1]*inputtensor[:,:,0,1]
                  + inputarray[:,:,2]*inputtensor[:,:,0,2])
    result[:,:,1] = (inputarray[:,:,0]*inputtensor[:,:,1,0]
                  + inputarray[:,:,1]*inputtensor[:,:,1,1]
                  + inputarray[:,:,2]*inputtensor[:,:,1,2])
    result[:,:,2] = (inputarray[:,:,0]*inputtensor[:,:,2,0]
                  + inputarray[:,:,1]*inputtensor[:,:,2,1]
                  + inputarray[:,:,2]*inputtensor[:,:,2,2])
    # Output array is of format [nx,ny,3]
    return result

def TransposeVectorArray(inputarray):
    # Assumes input array is of format [nA,nB,nV]
    # Output array is of format [nB,nA,nV]
    return np.transpose(inputarray, (1,0,2))

def rotateTensorArrayToVectorArray(inputtensor, inputvector):
    # assumes inputtensor is of shape [nx,ny,3,3] and 
    # inputvector is of shape [nx,ny,3]
    rotated = np.zeros(inputtensor.shape)
    nx, ny = inputtensor[:,:,0,0].shape
    # rotates so that tensor z-axis is parallel with vector
    for i in np.arange(nx):
        for j in np.arange(ny):
            rotated[i,j,:,:]= rotateTensorToVector(inputtensor[i,j,:,:],inputvector[i,j,:])
    # Output array is of format [nx,ny,3,3]
    return rotated

def TensorArrayParallelComponent(inputtensor):
    # assumes inputtensor is of shape [nx,ny,3,3]
    # Output array is of format [nx,ny]
    return inputtensor[:,:,2,2]

def TensorArrayPerpendicularComponent(inputtensor):
    # assumes inputtensor is of shape [nx,ny,3,3]
    # Output array is of format [nx,ny]
    return 0.5*(inputtensor[:,:,0,0]+inputtensor[:,:,1,1])

def TensorArrayAnisotropy(inputtensor):
    # assumes inputtensor is of shape [nx,ny,3,3]
    # Output array is of format [nx,ny]
    return np.divide(np.abs(TensorArrayPerpendicularComponent(inputtensor)), np.ma.masked_less_equal(np.abs(TensorArrayParallelComponent(inputtensor)),0))

def VectorArrayParallelComponent(inputvector, directionvector):
    # assumes inputvector and directionvector are of shape [nx,ny,3]
    dirnorm = np.divide(directionvector, np.linalg.norm(directionvector, axis=-1)[:,:,np.newaxis])
    # Need to perform dot product in smaller steps due to memory constraints
    # result = np.zeros(inputvector[:,:,0].shape)
    # for i in np.arange(len(inputvector[:,0,0])):
    #     result[i,:] = np.inner(inputvector[i,:,:],dirnorm[i,:,:])
    # Output array is of format [nx,ny]
    result = (inputvector*dirnorm).sum(-1) #dot product, alternatively numpy.einsum("ijk,ijk->ij",a,b)
    return result

def VectorArrayPerpendicularComponent(inputvector, directionvector):
    # assumes inputvector and directionvector are of shape [nx,ny,3]
    dirnorm = np.divide(directionvector, np.linalg.norm(directionvector, axis=-1)[:,:,np.newaxis])
    # Need to perform dot product in smaller steps due to memory constraints
    # paravector = dirnorm
    # for i in np.arange(len(inputvector[:,0,0])):
    #     paravector[i,:] = paravector[i,:] * np.inner(inputvector[i,:,:],dirnorm[i,:,:])[:,:,np.newaxis]
    #paravector = dirnorm * np.inner(inputvector, dirnorm)[:,:,np.newaxis]
    paravector = dirnorm * (inputvector*dirnorm).sum(-1)[:,:,np.newaxis] #dot product, alternatively numpy.einsum("ijk,ijk->ij",a,b)
    perpcomp = np.linalg.norm(inputvector - paravector, axis=-1) 
    # Output array is of format [nx,ny]
    return perpcomp

def VectorArrayAnisotropy(inputvector, directionvector):
    # assumes inputvector and directionvector are of shape [nx,ny,3]
    paracomp = VectorArrayParallelComponent(inputvector, directionvector)
    perpcomp = VectorArrayPerpendicularComponent(inputvector, directionvector)
    # Output array is of format [nx,ny]
    return np.divide(np.abs(perpcomp), np.ma.masked_less_equal(np.abs(paracomp),0))








def vec_MagneticPressureForce(inputarray):
    # assumes Magnetic field of shape [nx,ny,3]
    Bmap = inputarray
    Bnorm = np.linalg.norm(Bmap, axis=-1)
    mu0 = 1.25663706144e-6
    MagPressure = Bnorm*Bnorm/(2.*mu0)
    # Output array is of format [nx,ny,3]
    return (-1.)*numgradscalar(MagPressure)

def vec_MagneticTensionForce(inputarray):
    # Assumes input array is of format [nx,ny,3] with content B
    mu0 = 1.25663706144e-6
    jac = numjacobian(inputarray)
    nx, ny = inputarray[:,:,0].shape
    tensforce = np.zeros([nx,ny,3])
    tensforce[:,:,0] = (1./mu0) * (inputarray[:,:,0]*jac[:,:,0,0] + inputarray[:,:,1]*jac[:,:,0,1] + inputarray[:,:,2]*jac[:,:,0,2])
    tensforce[:,:,1] = (1./mu0) * (inputarray[:,:,0]*jac[:,:,1,0] + inputarray[:,:,1]*jac[:,:,1,1] + inputarray[:,:,2]*jac[:,:,1,2])
    tensforce[:,:,2] = (1./mu0) * (inputarray[:,:,0]*jac[:,:,2,0] + inputarray[:,:,1]*jac[:,:,2,1] + inputarray[:,:,2]*jac[:,:,2,2])
    # Output array is of format [nx,ny,3]
    return tensforce

def vec_ThermalPressureForce(inputarray):
    # assumes Thermal Pressure of shape [nx,ny]
    return (-1.)*numgradscalar(inputarray)

def vec_currentdensity(inputarray):
    # Assumes input array is of format [nx,ny,3] with content B or perb
    mu0 = 1.25663706144e-6
    # Output array is of format [nx,ny,3]
    return numcurl(inputarray) / mu0

def vec_Hallterm(currentdensity, magneticfield, numberdensity):
    # assumes current density of shape [nx,ny,3]
    # assumes Magnetic field of shape [nx,ny,3]
    # assumes number density of shape [nx,ny]
    unitcharge = 1.602177e-19
    crossp = np.cross(currentdensity, magneticfield)
    chargedensity = unitcharge * np.ma.masked_less_equal(numberdensity, 0)
    # Output array is of format [nx,ny,3]
    return np.divide(crossp,chargedensity)

def vec_ElectricFieldForce(electricfield, numberdensity):
    unitcharge = 1.602177e-19
    force = (-1.) * numberdensity[:,:,np.newaxis] * electricfield * unitcharge
    return force
    # force = np.zeros(np.array(electricfield.shape))
    # force[:,:,0] = (-1.) * numberdensity * electricfield[:,:,0] * unitcharge
    # force[:,:,1] = (-1.) * numberdensity * electricfield[:,:,1] * unitcharge
    # force[:,:,2] = (-1.) * numberdensity * electricfield[:,:,2] * unitcharge
    # return force





# pass_maps is a list of numpy arrays
# Each array has 2 dimensions [ysize, xsize]
# or 3 dimensions [ysize, xsize, components]
def expr_Hall_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','rho']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Rhomap = pass_maps['rho'].T # number density
    Jmap = vec_currentdensity(Bmap)
    Hallterm = vec_Hallterm(Jmap,Bmap,Rhomap)
    return np.linalg.norm(Hallterm, axis=-1).T

def expr_Hall_aniso(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','rho']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Rhomap = pass_maps['rho'].T # number density
    Jmap = vec_currentdensity(Bmap)
    Hallterm = vec_Hallterm(Jmap,Bmap,Rhomap)
    return VectorArrayAnisotropy(Hallterm,Bmap).T

def expr_J_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Jmap = vec_currentdensity(Bmap)
    return np.linalg.norm(Jmap, axis=-1).T

def expr_J_aniso(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Jmap = vec_currentdensity(Bmap)
    return VectorArrayAnisotropy(Jmap,Bmap).T

def expr_MagneticPressureForce_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    return np.linalg.norm(MagneticPressureForce, axis=-1).T

def expr_MagneticPressureForce_aniso(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    MagneticPressureForceAniso = VectorArrayAnisotropy(MagneticPressureForce, Bmap)
    return MagneticPressureForceAniso.T

# def expr_MagneticPressureForce_inplane_mag(pass_maps):
#     Bmap = TransposeVectorArray(pass_maps[0]) # Magnetic field
#     MagneticPressureForce = vec_MagneticPressureForce(Bmap)
#     return np.linalg.norm(inplane(MagneticPressureForce), axis=-1).T

# def expr_MagneticPressureForce_inplane_aniso(pass_maps):
#     Bmap = TransposeVectorArray(pass_maps[0]) # Magnetic field
#     MagneticPressureForce = vec_MagneticPressureForce(Bmap)
#     MagneticPressureForceIPAniso = VectorArrayAnisotropy(inplane(MagneticPressureForce), inplane(Bmap))
#     return MagneticPressureForceIPAniso.T

def expr_ThermalPressureForce_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['Pressure']
    Pmap = pass_maps['Pressure'].T #Pressure (scalar)
    return np.linalg.norm(vec_ThermalPressureForce(Pmap), axis=-1).T

def expr_ThermalPressureForce_aniso(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','Pressure']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Pmap = pass_maps['Pressure'].T #Pressure (scalar)
    ThermalPressureForceAniso = VectorArrayAnisotropy(vec_ThermalPressureForce(Pmap), Bmap)
    return ThermalPressureForceAniso.T

# def expr_ThermalPressureForce_inplane_aniso(pass_maps):
#     Bmap = TransposeVectorArray(pass_maps[0]) # Magnetic field
#     Pmap = pass_maps[1].T #Pressure (scalar)
#     ThermalPressureForce = vec_ThermalPressureForce(Pmap)
#     ThermalPressureForceIPAniso = VectorArrayAnisotropy(inplane(ThermalPressureForce), inplane(Bmap))
#     return ThermalPressureForceIPAniso.T

def expr_Eforce(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['E','rho']
    Emap = TransposeVectorArray(pass_maps['E']) # Electric field
    Rhomap = pass_maps['rho'].T # number density
    Efieldforce = vec_ElectricFieldForce(Emap, Rhomap)
    return np.linalg.norm(Efieldforce, axis=-1).T

def expr_Btension_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    MagTensionForce = vec_MagneticTensionForce(Bmap)
    return np.linalg.norm(MagTensionForce, axis=-1).T

# No expr_Btension_aniso as it's all in the perpendicular component

def expr_Bforces_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    MagTensionForce = vec_MagneticTensionForce(Bmap)
    return np.linalg.norm(MagneticPressureForce+MagTensionForce, axis=-1).T

def expr_Totforces_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','Pressure']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Pmap = pass_maps['Pressure'].T #Pressure (scalar)
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    MagTensionForce = vec_MagneticTensionForce(Bmap)
    ThermalPressureForce = vec_ThermalPressureForce(Pmap)
    return np.linalg.norm(ThermalPressureForce+MagneticPressureForce+MagTensionForce, axis=-1).T

def expr_Totforces_aniso(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','Pressure']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Pmap = pass_maps['Pressure'].T #Pressure (scalar)
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    MagTensionForce = vec_MagneticTensionForce(Bmap)
    ThermalPressureForce = vec_ThermalPressureForce(Pmap)
    return VectorArrayAnisotropy(ThermalPressureForce+MagneticPressureForce+MagTensionForce, Bmap).T

def expr_ratio_thermal_mag(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','Pressure']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Pmap = pass_maps['Pressure'].T #Pressure (scalar)
    MagneticPressureForce = vec_MagneticPressureForce(Bmap)
    MagTensionForce = vec_MagneticTensionForce(Bmap)
    ThermalPressureForce = vec_ThermalPressureForce(Pmap)
    MagForcesTot = np.ma.masked_less_equal(np.linalg.norm(MagneticPressureForce+MagTensionForce, axis=-1), 0)
    ThermalPressureForceTot = np.linalg.norm(ThermalPressureForce, axis=-1)
    return np.divide(ThermalPressureForceTot,MagForcesTot).T

def expr_E_parallel(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E']
    # No need for transposing in this one
    Bmap = pass_maps['B'] # Magnetic field
    Emap = pass_maps['E'] # Electric field
    return VectorArrayParallelComponent(Emap,Bmap)

def expr_E_perpendicular(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E']
    # No need for transposing in this one
    Bmap = pass_maps['B'] # Magnetic field
    Emap = pass_maps['E'] # Electric field
    return VectorArrayPerpendicularComponent(Emap,Bmap)

def expr_flowcompression(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['V']
    Vmap = TransposeVectorArray(pass_maps['V']) # Bulk flow
    return numdiv(Vmap).T

# def expr_gradB_aniso(pass_maps):
#     Bmap = TransposeVectorArray(pass_maps[0]) # Magnetic field
#     gradB = numjacobian(Bmap)
#     rotatedgradB = rotateTensorArrayToVectorArray(gradB,Bmap)
#     gradBaniso = TensorArrayAnisotropy(rotatedgradB)
#     return gradBaniso.T

# def expr_gradPTD_aniso(pass_maps):
#     Bmap = TransposeVectorArray(pass_maps[0]) # Magnetic field
#     PTDmap = TransposeVectorArray(pass_maps[1]) # PressureTensorDiagonal
#     gradP = numjacobian(PTDmap)
#     rotatedgradP = rotateTensorArrayToVectorArray(gradP,Bmap)
#     gradPaniso = TensorArrayAnisotropy(rotatedgradP)
#     return gradPaniso.T

# Pressure anisotropy is PPerpOverPar

# betaParallel, betaPerpendicular, betaPerpOverPar, rMirror



def expr_betatron(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E','PPerpendicular']
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of lists, where the top level
    # is 2N+1 timesteps with the middle one the requested time step

    # This custom expression returns a proxy for betatron acceleration
    if type(pass_maps[0]) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("expr_betatron expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()

    # Multiple time steps were found. This should be 3, for a time derivative.
    dsteps = [x['dstep'] for x in pass_maps]
    curri = dsteps.index(0)
    previ = dsteps.index(-1)

    thesemaps = pass_maps[curri]
    pastmaps = pass_maps[previ]

    thisB = TransposeVectorArray(thesemaps['B'])
    pastB = TransposeVectorArray(pastmaps['B'])
    thisBmag = np.linalg.norm(thisB, axis=-1)
    pastBmag = np.linalg.norm(pastB, axis=-1)
    dBdt = (thisBmag-pastBmag)/DT

    thisE = TransposeVectorArray(thesemaps['E'])
    thisPperp = thesemaps['PPerpendicular'].T

    # E x B drift = (E X B)/B^2
    B2 = (thisBmag*thisBmag)
    UExB = numcrossproduct(thisE,thisB)/B2[:,:,np.newaxis]
    gradB = numgradscalar(thisBmag)
    
    # dot product is (A*B).sum(-1)
    result = (thisPperp/thisBmag)*(dBdt + (UExB*gradB).sum(-1)) 
    return result.T

def expr_Fermi(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E','PParallel']
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Emap = TransposeVectorArray(pass_maps['E']) # Electric field
    Pparallel = pass_maps['PParallel'].T #Pressure (scalar)

    Bmag = np.linalg.norm(Bmap, axis=-1)
    B2 = (Bmag*Bmag)
    UExB = numcrossproduct(Emap,Bmap)/B2[:,:,np.newaxis]

    Bnorm = (Bmap/Bmag[:,:,np.newaxis])
    Bjac = numjacobian(Bnorm)
    kappa = numvecdottensor(Bnorm,Bjac)
    # (A*B).sum(-1) is dot product
    result = Pparallel*(UExB*kappa).sum(-1)
    return result.T

def expr_EJ_parallel(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E']
    # Needs B,E
    Bmap = TransposeVectorArray(pass_maps['B']) # Magnetic field
    Emap = TransposeVectorArray(pass_maps['E']) # Electric field

    Jmap = vec_currentdensity(Bmap)
    Jpara = VectorArrayParallelComponent(Jmap, Bmap)
    Epara = VectorArrayParallelComponent(Emap, Bmap)
    EJpara = (Jpara*Epara)
    return EJpara.T

def expr_Eacc_parallel(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','E','V','rho']
    # No need for transposing in this one
    # Needs B,E,V,rho
    Bmap = pass_maps['B'] # Magnetic field
    Emap = pass_maps['E'] # Electric field
    Vmap = pass_maps['V'] # Bulk flow
    rhomap = pass_maps['rho'] # density
    protoncharge = .16021773e-18

    Vpara = VectorArrayParallelComponent(Vmap, Bmap)
    Epara = VectorArrayParallelComponent(Emap, Bmap)
    qnEVpara = (Vpara*Epara*rhomap*protoncharge)
    return qnEVpara


def expr_diamagnetic(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','rho','V','Pressure']
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    # for time averaging, it's a list of lists, where the top level
    # is 2N+1 timesteps with the middle one the requested time step

    # This custom expression returns a proxy for betatron acceleration
    if type(pass_maps) is not list:
        # Not a list of time steps, calculating this value does not make sense.
        print("expr_diamagnetic expected a list of timesteps to average from, but got a single timestep. Exiting.")
        quit()

    # Multiple time steps were found. This should be 3, for a time derivative.
    dsteps = [x['dstep'] for x in pass_maps]
    curri = dsteps.index(0)
    previ = dsteps.index(-1)
    thesemaps = pass_maps[curri]
    pastmaps = pass_maps[previ]

    thisV = TransposeVectorArray(thesemaps['V'])
    pastV = TransposeVectorArray(pastmaps['V'])

    B = TransposeVectorArray(thesemaps['B'])
    rhom = thesemaps['rho'].T * 1.6726e-27
    Pres = thesemaps['Pressure'].T

    dUdt = (thisV-pastV)/DT
    rhodUdt = rhom[:,:,np.newaxis] * dUdt
    BperB2 =  B / (B*B).sum(-1)[:,:,np.newaxis]

    gradP = numgradscalar(Pres)

    term1 = numcrossproduct(BperB2, rhodUdt)
    term2 = numcrossproduct(BperB2, gradP)

    result = np.linalg.norm(term1+term2, axis=-1)
    return result.T

def expr_diamagnetic_noinertial(pass_maps, requestvariables=False):
    if requestvariables==True:
        return ['B','Pressure']
    # pass_maps is a list of numpy arrays
    # Each array has 2 dimensions [ysize, xsize]
    # or 3 dimensions [ysize, xsize, components]

    B = TransposeVectorArray(pass_maps['B'])
    Pres = pass_maps['Pressure'].T

    BperB2 =  B / (B*B).sum(-1)[:,:,np.newaxis]
    gradP = numgradscalar(Pres)

    # Dropping the first term i.e. inertial current
    term2 = numcrossproduct(BperB2, gradP)

    result = np.linalg.norm(term2, axis=-1)
    return result.T

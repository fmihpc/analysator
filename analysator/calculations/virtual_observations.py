
import numpy as np
def simulation_to_spacecraft_frame(spinvector, detector_axis, phi=0):
    ''' Builds a matrix to transform coordinates from simulation frame into spaceraft frame
    :param spinvector  Spacecraft spin axis, in simulation coordinates
    :param detector_axis Detector plane normal axis
    :param phi         Rotate spacecraft around spin axis after setting up coordinate system
    '''
    # TODO: Normalize vectors?
    y = np.cross(detector_axis,spinvector)
    z = np.cross(spinvector, y)
    yr = np.cos(phi)*y - np.sin(phi)*z
    zr = np.sin(phi)*y + np.cos(phi)*z
    m = np.array([spinvector, yr, zr])

    return m

def spacecraft_to_simulation_frame(spinvector, detector_axis, phi=0):
    ''' Builds a matrix to transform coordinates from spaceraft frame back to simulation frame
    :param spinvector  Spacecraft spin axis, in simulation coordinates
    :param detector_axis Detector plane normal axis
    :param phi         Rotate spacecraft around spin axis after setting up coordinate system
    '''
    return simulation_to_spacecraft_frame(spinvector,detector_axis,phi).T

def simulation_to_observation_frame(x_axis,y_axis):
    ''' Builds a 3x3 matrix to transform velocities into an observation plane
    :param x_axis:  x-axis of the observation plane (in simulation coordinates)
    :param y_axis:  y-axis of the observation plane (gets orthonormalized)
    '''
    xn = np.linalg.norm(x_axis)
    x_axis /= xn
    p = x_axis.dot(y_axis)
    y_axis -= p*x_axis
    yn = np.linalg.norm(y_axis)
    y_axis /= yn
    z_axis = np.cross(x_axis,y_axis)
    return np.array([x_axis,y_axis,z_axis])



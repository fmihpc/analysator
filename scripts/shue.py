'''

    See Shue et al. (1997): "A new functional form to study the solar wind control of the magnetopause size and shape"

    They describe a model where the radial magnetopause position r is treated as a function of angle wrt x-GSE axis

        r(theta) = r_0 * (2 / (1+ cos(theta)))^alpha

    r_0 is the subsolar standoff distance [in earth radii], alpha is a dimensionless parameter

    The parameters r_0, alpha depend on the solar wind parameters:
        B_z: z-component (GSE) of the magnetic field [nT]
        n_p: number density [cm^-3]
        v_sw: solar wind speed [km/sec]

    Different Vlasiator runs had 

Example usage (plot magnetopause for run 'EGI'):

    import numpy as np
    import matplotlib.pyplot as plt

    theta = np.linspace(0,np.pi / 2, 1000)
    r, r0, alpha = f_shue(theta, run='EGI')
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    plt.plot(x,y)

'''

import numpy as np
import logging

def _f_shue_parametrized(theta_polar, r_0, alpha):
    '''
        returns r as a function of theta,
        parameters: r_0, alpha ()
    '''
    return r_0 * (2 / (1 + np.cos(theta_polar)))**alpha

def _f_shue_parameters(run):
    '''
        Lookup function for Vlasiator run parameters, to be input into Shue model

        Inputs:
            run [string]: the name of the Vlasiator run
        Outputs:
            Solar wind driving parameters: B_z [nT], n_p [cm^-3], v_sw [km/s]
    '''
    if (run == 'EGI'):
        B_z = -5                  # nT
        n_p = 1                   # cm^-3
        v_sw = 750                # km/sec
    elif (run == 'EGK'):
        B_z = -20                  # ""
        n_p = 1                   
        v_sw = 750                
    elif (run == 'EGL'):          # note: initialized with EGI conditions
        B_z = -10
        n_p = 4
        v_sw = 750
    elif (run == 'EGM'):
        B_z = -5
        n_p = 1
        v_sw = 750
    elif (run == 'EGP'):          # note: B_x = -0.5 for this run
        B_z = -20
        n_p = 7
        v_sw = 1000
    else:
        B_z = 0
        n_p = 1
        v_sw = 500
        logging.info('VLASIATOR RUN NOT SPECIFIED!!!')     # error message
    return B_z, n_p, v_sw


def _r_0_alpha_shue(B_z, n_p, v_sw):
    '''
        Calculate r_0, alpha, depending on the input solar wind parameters.

        Inputs:
            B_z [nT], n_p [cm^-3], v_sw [km/s], as described above
        Output:
            r_0 [R_E], alpha [dimensionless]
    '''
    m_p  = 1.67262158e-27              # proton mass [kg]
    n_p_SI = n_p * 1e6             # m^-3
    v_sw_SI = v_sw * 1e3           # m/sec
    rho = n_p_SI * m_p
    D_p = (rho / 2) * (v_sw_SI**2) * 1e9  # dynamic pressure, in nanoPAscals
    if B_z >= 0:
        r_0 = (11.4 + 0.013 * B_z) * D_p**(-1 / 6.6)              # Eq. 12, Shue 1997
    elif B_z < 0:
        r_0 = (11.4 + 0.14 * B_z) * D_p**(-1 / 6.6)
    alpha = ( 0.58 - (0.010 * B_z) ) * (1 + (0.010 * D_p))         # Eq. 13, Shue 1997
    return r_0, alpha


def f_shue(theta_polar, run = None, B_z = None, n_p = None, v_sw = None):
    '''
     Shue et al. (1997): 'A new functional form to study the solar wind control of the magnetopause size and shape'`

     Inputs:
        theta_polar: 1D numpy array [radians], angle wrt x-GSE axis
            0 < theta_polar < 2pi

        keyword run: string containing Vlasiator run name (e.g. 'EGI') 
            If set, look up the solar wind parameters B_z, n_p, v_sw for this run

     Output:
        magnetopause position r(theta) [R_E]

     TODO: instead of using _f_shue_parameters(), use VlsvReader object's get_config() method to look up parameters.
    '''

    if run is not None:
        B_z, n_p, v_sw = _f_shue_parameters(run)
    r_0, alpha = _r_0_alpha_shue( B_z, n_p, v_sw )
    r = _f_shue_parametrized(theta_polar, r_0, alpha)
    return r, r_0, alpha


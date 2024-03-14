
#Example:
# theta = np.linspace(0,np.pi, 1000)
# r = f_shue(theta, run='EGI')
#
# x = r*np.cos(theta)
# y = r*np.sin(theta)
# plt.plot(x,y)


def f_shue_parametrized(theta_polar, r_0, alpha):
    ''' Shue et al. (1997): A new functional form to study the solar wind control
    returns r as a function of theta,
    parameters: r_0, alpha ()
    '''
    return r_0 * (2 / (1 + np.cos(theta_polar)))**alpha

def f_shue_parameters(run):
    m_p  = 1.67262158e-27              # proton mass [kg]
    if (run == 'EGI'):
        B_z = -5                  # nT
        n_p = 1                   # cm^-3
        v_sw = 750                # km/sec
    elif (run == 'EGK'):
        B_z = -20                  # nT
        n_p = 1                   # cm^-3
        v_sw = 750                # km/sec
    elif (run == 'EGL'):
        B_z = -10
        n_p = 4
        v_sw = 750
    elif (run == 'EGM'):
        B_z = -5
        n_p = 1
        v_sw = 750
    elif (run == 'EGP'):
        B_z = -20    # B_x = -0.5??
        n_p = 7
        v_sw = 1000
    else:
        B_z = 0
        n_p = 1
        v_sw = 500
        print('VLASIATOR RUN NOT SPECIFIED!!!')     # error message
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


def f_shue(theta_polar, run = None):
    # Shue et al. (1997): 'A new functional form to study the solar wind control
    # of the magnetopause size and shape'
    # INPUTS 0 < theta_polar < 2pi
    # theta_polar is a numpy array
    # outputs magnetopause position r(theta) [in R_E]
    r_0, alpha = f_shue_parameters(run)
    r = f_shue_parametrized(theta_polar, r_0, alpha)
    return r, r_0, alpha


def f_shue_make_func(r_0, alpha):
    return lambda theta_polar: r_0 * (2 / (1 + np.cos(theta_polar)))**alpha


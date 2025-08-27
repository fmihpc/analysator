import numpy as np
from scipy import interpolate


def mpos_oscillations(t, n, v, c, dt = 1e-3, t_cross = 0, R_0 = None, dRdt_0 = 0, f=2.44, dip_mom = 8.22e22):
    '''
    Implements formula for magnetopause oscillations.
    See Horaites et al. (2023) https://doi.org/10.1029/2023JA031374

    The function requires the following time series inputs,
    which are 1d arrays of equal size:
    INPUTS:
    t: time [seconds]
    n: density [m^-3]
    v: velocity [m/s]
    c: the inertial term (see Horaites et al., 2023) --- of order one

    KEYWORDS:
    dt: step size of the integration [seconds]
    t_crossing: the (assumed constant) fast mode magnetosheath crossing time [seconds]
    R_0: The initial standoff distance [in meters]
    dRdt_0: the initial rate of change [m/s]
    f: the "dipole compression factor". Typically 2-3 (2.44 fiducially). Note if Bz<0, then f is may be <2
    dip_mom: Earth's dipole moment, [A*m^2] (note in Vlasiator simulations, this may be different)

    The program first interpolates n(t), v(t), and c(t) to the required resolution.
    Then the evolution equation for the subsolar magnetopause standoff, R(t), is integrated.

    The initial standoff distance R_0 and rate of change wrt time dRdt_0 are calculated
    by default from the data. Or their values can be specified as keywords
    '''
    mu_0 = 4*np.pi*1e-7 #A/m^2
    m_p = 1.67262192e-27 #kg (proton mass)
    
    t_intp = np.arange(t[0], t[-1], step = dt)

    n_F = n[-1]
    v_F = v[-1]

    if type(c) == float or type(c) == int:
        c = np.zeros(t.size) + c

    if R_0 is None:
        R_0 = (f**2*mu_0*dip_mom**2/(32*np.pi**2*n[0]*m_p*v[0]**2))**(1/6.) # initial equilibrium R [m]

    R_F = (f**2*mu_0*dip_mom**2/(32*np.pi**2*n_F*m_p*v_F**2))**(1/6.)   # final equilibrium, best guess

    # Interpolate
    f_n = interpolate.interp1d(t + t_cross, n, fill_value='extrapolate', bounds_error = False)
    f_v = interpolate.interp1d(t + t_cross, v, fill_value='extrapolate', bounds_error = False)
    f_c = interpolate.interp1d(t, c, fill_value='extrapolate', bounds_error = False)

    n_intp = f_n(t_intp)
    v_intp = f_v(t_intp)
    c_intp = f_c(t_intp)
     
    eta_intp = n_intp / n_F     # eta ~ n / n_F

    nt = t_intp.size + 1
    R = np.zeros(nt)
    dRdt = np.zeros(nt)

    # initial conditions
    R[0] = R_0
    dRdt[0] = dRdt_0

    # Integrate to find R(t)
    for i, t_i in enumerate(t_intp):
        # dx/dt = v(t)
        R[i+1] = R[i] + (dt * dRdt[i])
        dRdt[i+1] = dRdt[i] - dt/(c_intp[i] * R_F) * ( eta_intp[i] * (v_intp[i]+dRdt[i])**2 - v_F**2*R_F**6/R[i]**6)

    return R[0:-1], t_intp


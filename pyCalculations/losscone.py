import numpy as np

def loss_cone_angle(cellcoord=None,cellcoordre=None,deg=False):
    ''' Calculates the value of the loss cone angle at a given location
        :kword cellcoord:    The coordinates (X,Y,Z) of the cell whose loss cone angle value is calculated [in m]
        :kword cellcoordre:  The coordinates (X,Y,Z) of the cell whose loss cone angle value is calculated [in Re]
        :kword deg:          True if user wants the angle in degrees

        :returns:            The value of the loss cone angle

        .. code-blocks:: python
    '''

    # Earth radius [m]
    R_E = 6.371e6

    # Convert coordinates to Re if needed
    if cellcoord!=None:
        X = cellcoord[0]/R_E
        Y = cellcoord[1]/R_E
        Z = cellcoord[2]/R_E

    else:
        X = cellcoordre[0]
        Y = cellcoordre[1]
        Z = cellcoordre[2]

    # Calculation of R and L
    R = np.sqrt(X**2+Y**2+Z**2)                   # Radial distance to Earth centre
    if np.sqrt(X**2+Y**2)!=0:                     # Magnetic latitude
        lat_m = np.arctan(Z/np.sqrt(X**2+Y**2))
    else:
        lat_m = np.sign(Z)*np.pi/2.
    L = R/np.cos(lat_m)**2                        # L-shell


    # Analytical formula for loss cone angle (dipole approximation)
    alph0 = np.arcsin(R**-1.5 * (4*L-3*R)**.25/(4*L-3.)**.25)

    # Conversion to degrees if needed
    if deg:
        alph0 = alph0*180./np.pi

    return alph0

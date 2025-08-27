'''
Calculates shock crossing values from Rankine-Hugoniot relations
in the deHoffmann-Teller (dHT) frame.
'''

import numpy as np
import math
import scipy.optimize
import logging

mu0 = 4*math.pi*1.e-7
mp = 1.67e-27
c = 299792458
k = 1.3806506e-23
Gamma = 5.0/3.0

def polynome(X, theta, V1sq, beta1, vA1sq, Gamma, vs1sq):
    return (V1sq-X*vA1sq)**2 * (X*vs1sq + 0.5*V1sq*np.cos(theta)**2*(X*(Gamma-1)-(Gamma+1))) + 0.5*vA1sq*V1sq*X*np.sin(theta)**2 * ((Gamma+X*(2-Gamma))*V1sq - X*vA1sq*((Gamma+1)-X*(Gamma-1))) #from Koskinen's Physics of Space Storms -book

def newtonmethod(theta, V1sq, beta1, vA1sq, Gamma, vs1sq):
    calctemp1 = 1. + 0.5*Gamma*beta1
    cos12 = np.cos(theta)**2
    sin12 = np.sin(theta)**2
    MA2=V1sq/vA1sq
    Ztry = max( ((0.5/cos12)*(calctemp1 + np.sqrt(calctemp1**2 - 2.*Gamma*beta1*cos12)) -1.),
                ((0.5/cos12)*(calctemp1 - np.sqrt(calctemp1**2 - 2.*Gamma*beta1*cos12)) -1.), 0.)
    # First (root for M**2) -1
    # Second and third (root for M**2) -1
    fform = (1. +Ztry)*((Ztry**2)*8.*cos12 +(3. -5.*Ztry)*sin12) -(Ztry**2)*5.*beta1
    gform = (1. +Ztry)*((Ztry**2)*2.*cos12 +(3. +Ztry)*sin12)
    Rtry = fform/gform
    M2try = (1. +Ztry)*Rtry
    rstep = 1.0
    compr = 0.
    while ((rstep >= 0.0001) and (compr < 0.001)):
        Ztry = max( ((0.5/cos12)*(calctemp1 + np.sqrt(calctemp1**2 - 2.*Gamma*beta1*cos12)) -1.),
                    ((0.5/cos12)*(calctemp1 - np.sqrt(calctemp1**2 - 2.*Gamma*beta1*cos12)) -1.), 0.)
        fform = (1. +Ztry)*((Ztry**2)*8.*cos12 +(3. -5.*Ztry)*sin12) -(Ztry**2)*5.*beta1
        gform = (1. +Ztry)*((Ztry**2)*2.*cos12 +(3. +Ztry)*sin12)
        Rtry = fform/gform
        M2try = (1. +Ztry)*Rtry

        while (abs(M2try - MA2) > 0.0001):
            fderi = (Ztry**2)*8.*cos12 +(3. -8.*Ztry)*sin12 -10.*Ztry*beta1 +(1. +Ztry)*(16.*Ztry*cos12 -5.*sin12)
            gderi = (Ztry**2)*2.*cos12 +(3. +Ztry)*sin12 +(1. +Ztry)*(4.*Ztry*cos12 +sin12)
            rderi = (gform*fderi-fform*gderi)/(gform**2)
            m2deri = (1. +Ztry)*rderi + Rtry

            # Newton step forward
            Ztry = Ztry + (MA2 - M2try)/m2deri * 0.5*rstep

            # Calculate new Rtry and M2try
            fform = (1. +Ztry)*((Ztry**2)*8.*cos12 +(3. -5.*Ztry)*sin12) -(Ztry**2)*5.*beta1
            gform = (1. +Ztry)*((Ztry**2)*2.*cos12 +(3. +Ztry)*sin12)
            Rtry = fform/gform
            M2try = (1. +Ztry)*Rtry

        if (Rtry <= MA2):
            compr = Rtry
        rstep = rstep * 0.1
    return compr

def rankine(Tu, rhou, V, B, n, Vsh):
    '''
    Call obliqueshock.rankine(Tu, rhou, V, B, n, Vsh) to compute the shock crossing values

        Inputs:
            Tu: upstream proton temperature [K] rhou: upstream proton number density [1/m3] V: 3-element upstream proton inflow velocity vector [m/s]
            B: 3-element upstream magnetic field vector [T] n: 3-element shock normal vector Vsh: 3-element shock velocity vector [m/s]

        Returns:
            The shock compression ratio
            X: Compression ratio from scipy.optimize X2: Compression ratio from the Newton method

        Example:
            obliqueshock.rankine(5e5, 1.0e6, [-750e3,0,0], [3.5355e-9,0,-3.5355e-9], [1,0,0], 0)\n
            -> Computes the shock crossing for Tu = 500 kK, rhou = 1/cc, V = [-750,0,0]km/s (inflow plasma speed is 750 km/s in -X),
            B = [3.5355e,0,-3.5355e]nT (upstream magnetic field is 5 nT at 45 degree angle),
            n = [1,0,0], Vsh = 0 (shock front points in +X direction and is stationary in input frame)
    '''

    # V, B, n are vectors
    V = np.array(V)
    B = np.array(B)
    n = np.array(n)

    # normalise n
    n = n/np.linalg.norm(n)
    Vshvect = n*Vsh

    Vu = V - Vshvect
    Bu = B

    pu = rhou * k * Tu
    vHT = np.cross(n, np.cross(Vu, Bu)) / np.dot(n,Bu)
    logging.info("The de Hoffmann Teller transformation velocity is [km/s]: " + str(np.linalg.norm(vHT)/1e3))

    VuHT = Vu - vHT
    BuHT = Bu
    #BuHT2 = Bu + 1/c**2 * np.cross(vHT, np.cross(-Vu,Bu)) # Eu = -Vu X Bu

    Emotional = -np.cross(VuHT,BuHT)
    logging.info("Verify motional E-field in dHT frame: " + str(Emotional))

    theta = np.arccos(np.dot(BuHT,n)/np.linalg.norm(BuHT))
    logging.info("Theta [degree]: " + str(np.rad2deg(theta)))

    nPlane = np.cross(Vu, Bu) / np.linalg.norm(np.cross(Vu, Bu)) # normal to the plane containing B and V
    t = np.cross(n, nPlane) / np.linalg.norm(np.cross(n, nPlane)) # transversal direction, i.e. normal to n and within the (B,V) plane (i.e. in turn normal to nPlane)

    # compute the normal and tangential components of upstream V and B
    VnuHT = np.dot(VuHT, n)
    VtuHT = np.dot(VuHT, t)

    BnuHT = np.dot(BuHT, n)
    BtuHT = np.dot(BuHT, t)

    # upstream plasma state parameters
    beta1 = 2.0*mu0*pu / np.dot(Bu,Bu)
    vAusq = np.dot(BuHT,BuHT)/(mu0*mp*rhou)
    vsusq = Gamma*pu/(mp*rhou)

    logging.info("vAu [km/s]: " + str(np.sqrt(vAusq)/1e3))
    logging.info("vSu [km/s]: " + str(np.sqrt(vsusq)/1e3))
    logging.info("Shock MA: " + str(np.linalg.norm(Vu)/np.sqrt(vAusq)))

    # compare the upstream initial speed to dHT speed
    logging.info("|Vu| [km/s]: " + str(np.linalg.norm(Vu)/1e3))
    logging.info("|VuHT| [km/s]: " + str(np.linalg.norm(VuHT)/1e3))
    # compression ratio
    sol = scipy.optimize.root(polynome, 2., args=(theta, np.dot(VuHT,VuHT), beta1, vAusq, Gamma, vsusq))
    X = sol.x
    # Alternative own Newton method
    X2 = newtonmethod(theta, np.dot(VuHT,VuHT), beta1, vAusq, Gamma, vsusq)
    ##logging.info("X " +str(X) + " X2 " + str(X2))

    VuHTsq = np.dot(VuHT,VuHT)
    VndHT = VnuHT / X
    VtdHT = VtuHT * (VuHTsq - vAusq)/(VuHTsq-X*vAusq)
    VdHTsq = VndHT**2 + VtdHT**2
    pd = pu * (X+((Gamma-1.)*X*VuHTsq)/(2.0*vsusq)*(1-VdHTsq/VuHTsq))
    rhod = rhou * X
    BndHT = BnuHT
    BtdHT = BtuHT * (VuHTsq-vAusq)*X/(VuHTsq-X*vAusq)

    BdHT = (BndHT * n) + (BtdHT * t)
    VdHT = (VndHT * n) + (VtdHT * t)

    Bd = BdHT
    #Bd2 = BdHT - 1/c**2 * np.cross(vHT, np.cross(-Vu,Bu))
    #logging.info("Bd "+str(Bd)+" alt "+str(Bd2))

    Vd = VdHT + vHT + Vshvect
    XB = np.linalg.norm(Bd)/np.linalg.norm(Bu)

    #//Td = pd / (Gamma * rhod * k)
    Td = pd / (rhod * k)

    #print compression ratios and upstream/downstream state
    logging.info("Gas compression ratio: " + str(X[0]))
    logging.info("Magnetic compression ratio: " + str(XB))

    logging.info("")
    logging.info("This determines the dHT upstream state")
    logging.info("Density [1/cm3]: " + str(rhou/1e6))
    logging.info("Temperature [K]: " + str(Tu))
    logging.info(" V [km/s]: " + str(VuHT/1e3))
    logging.info(" B [nT]: " + str(BuHT*1e9))
    logging.info(" |V| [km/s]: " + str(np.linalg.norm(VuHT)/1e3))
    logging.info(" |B| [nT]: " + str(np.linalg.norm(BuHT)*1e9))

    logging.info("")
    logging.info("This determines the dHT downstream state")
    logging.info("Density [1/cm3]: " + str(rhod[0]/1e6))
    logging.info("Temperature [K]: " + str(Td[0]))
    logging.info(" V [km/s]: " + str(VdHT/1e3))
    logging.info(" B [nT]: " + str(BdHT*1e9))
    logging.info(" |V| [km/s]: " + str(np.linalg.norm(VdHT)/1e3))
    logging.info(" |B| [nT]: " + str(np.linalg.norm(BdHT)*1e9))

    return X[0], X2

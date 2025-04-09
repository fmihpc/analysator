'''
Script calculates shock crossing values from Rankine-Hugoniot relations
in the normal incidence frame (NIF).
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

def polynome(X, theta, beta1, MA1):
    return (np.cos(theta)**2)*(2*(MA1**2) + 5*beta1*(np.cos(theta)**2))*(X**3) + (MA1**2)*(MA1**2 - (np.cos(theta)**2)*(5*(MA1**2) + 8 + 10*beta1))*(X**2) + (MA1**4)*(11*(np.cos(theta)**2) + 2*(MA1**2) + 5 + 5*beta1)*X - 8*(MA1**6)

def rankine(Tu, rhou, V, B):
    '''
    Call obliqueshock_nif.rankine(Tu, rhou, V, B, n, Vsh) to compute the shock crossing values

        Inputs:
            Tu: upstream proton temperature [K] rhou: upstream proton number density [1/m3] V: 3-element upstream proton inflow velocity vector [m/s]
            B: 3-element upstream magnetic field vector [T]\n
            Note: Give the plasma velocity relative to the shock speed, which is set to zero in the script.

        Returns:
            The shock compression ratio
            X: Plasma compression ratio XB: Magnetic field compression ratio

        Example:
            obliqueshock_nif.rankine(5e5, 1.0e6, [-750e3,0,0], [3.5355e-9,0,-3.5355e-9])\n
            -> Computes the shock crossing for Tu = 500 kK, rhou = 1/cc, V = [-750,0,0]km/s (inflow plasma speed is 750 km/s in -X),
            B = [3.5355e,0,-3.5355e]nT (upstream magnetic field is 5 nT at 45 degree angle),
    '''

    # V, B, n are vectors
    V = np.array(V)
    B = np.array(B)

    # normal and tangential vectors
    n = -V/np.linalg.norm(V)
    nPlane = np.cross(V, B)/np.linalg.norm(np.cross(V, B)); # normal to the plane containing B and V
    t = np.cross(n, nPlane)/np.linalg.norm(np.cross(n, nPlane));

    logging.info("Normal vector n: " + str(n))
    logging.info("Tangential vector t: " + str(t))
    # upstream parameters
    Vu = V
    Bu = B
    pu = rhou * k * Tu

    Vnu = np.dot(Vu, n)
    Vtu = np.dot(Vu, t)

    Bnu = np.dot(Bu, n)
    Btu = np.dot(Bu, t)

    # compute relevant variables: Upstream beta, vA, MA and shock angle (theta)
    theta = np.arccos(np.dot(Bu,n)/np.linalg.norm(Bu))
    logging.info("Theta [degree]: " + str(np.rad2deg(theta)))

    betau = 2.0*mu0*pu / np.dot(Bu,Bu)
    vAu = np.sqrt(np.dot(Bu,Bu) / (mu0*mp*rhou))
    MAu = np.linalg.norm(Vu) / vAu

    logging.info("vAu [km/s]: " + str(vAu/1e3))
    logging.info("Shock MA: " + str(np.linalg.norm(Vu)/vAu))
    logging.info("Upstream plasma beta: " + str(betau))
    # compression ratio
    sol = scipy.optimize.root(polynome, 2., args=(theta, betau, MAu))
    X = sol.x

    # get downstrean parameters
    Z1 = (MAu**2-np.cos(theta)**2) / (MAu**2-X*np.cos(theta)**2)
    Z2 = np.sin(theta)*np.cos(theta)*(X-1) / (MAu**2-X*np.cos(theta)**2)
    Z3 = (betau/MAu**2 + 2*(X-1)/X + np.sin(theta)**2*(1 - (X*Z1)**2 )/MAu**2) / (2*X)

    rhod = rhou * X

    Vnd = (Vnu/X)
    Vnt = Z2*Vnu
    Vd = (Vnd * n) + (Vnt * t)

    Bnd = Bnu
    Btd = X*Z1*Btu
    Bd = (Bnu * n) + (Btd * t)

    vTHd_sq = Z3*Vu[0]**2;
    Td = (mp*vTHd_sq/k);

    XB = np.linalg.norm(Bd)/np.linalg.norm(Bu)

    # print compression ratios and upstream/downstream state
    logging.info("Gas compression ratio: " + str(X[0]))
    logging.info("Magnetic compression ratio: " + str(XB))

    logging.info("")
    logging.info("This determines the NIF upstream state")
    logging.info("Density [1/cm3]: " + str(rhou/1e6))
    logging.info("Temperature [K]: " + str(Tu))
    logging.info(" V [km/s]: " + str(Vu/1e3))
    logging.info(" B [nT]: " + str(Bu*1e9))
    logging.info(" |V| [km/s]: " + str(np.linalg.norm(Vu)/1e3))
    logging.info(" |B| [nT]: " + str(np.linalg.norm(Bu)*1e9))

    logging.info("")
    logging.info("This determines the NIF downstream state")
    logging.info("Density [1/cm3]: " + str(rhod[0]/1e6))
    logging.info("Temperature [K]: " + str(Td[0]))
    logging.info(" V [km/s]: " + str(Vd/1e3))
    logging.info(" B [nT]: " + str(Bd*1e9))
    logging.info(" |V| [km/s]: " + str(np.linalg.norm(Vd)/1e3))
    logging.info(" |B| [nT]: " + str(np.linalg.norm(Bd)*1e9))

    return X[0], XB


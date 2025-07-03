# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 15:50:26 2021

@author: hmfausey

Module used by create_data_set to estimate luminosity distances of GRBs at 
different redshifts.

See https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html for a full general 
overview
"""

###############################################################################
####################################IMPORTS####################################
###############################################################################

import numpy as np

###############################################################################
###################################CONSTANTS###################################
###############################################################################

H0 = 67.36  # (km/(sMpc))
Omega_M = 0.3166
Omega_Lam = 0.6847
c = 2.998 * 10**5  # km/s
d_H = c / H0  # Mpc -- hubble distance

###############################################################################
###################################FUNCTIONS###################################
###############################################################################


def d_L(z):
    ##Luminosity Distance (Weinberg 1972, pp. 420-424; Weedman 1986, pp. 60-62,
    # or https://ned.ipac.caltech.edu/level5/Hogg/Hogg7.html)
    return (1 + z) * d_M(z)


def d_M(z):
    ##Transverse comoving distance -- equal to comoving distance for flat
    # universe (see https://ned.ipac.caltech.edu/level5/Hogg/Hogg5.html)
    return d_C(z)


def d_C(z):
    ##Comoving distance (see Peebles 2013 (pgs 310-321) or
    # https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html)
    # Performs integration of E(z) over z, and then multiplies by d_H(the hubble
    # distance)
     #Inputs:
         #z - redshift
     #Returns:
         #dC - comoving distance in Mpc
    
    step = 0.001
    zrange = np.arange(0, z + step, step)
    # For the function below, we set Omega_k to 0 since the Universe is
    # essentially flat
    func = 1 / np.sqrt(Omega_M * (1 + zrange) ** 3 + Omega_Lam)

    integration = 0

    for i in range(len(func) - 1):
        # Using trapezoid rule
        integration += 0.5 * (func[i] + func[i + 1]) * step

    dC = d_H * integration
    return dC

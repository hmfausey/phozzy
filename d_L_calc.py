# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 15:50:26 2021

@author: hmfausey

See https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html for a full general overview
"""

####################
#IMPORTS
####################

import numpy as np

####################
#CONSTANTS
####################

H0 = 67.36 #(km/(sMpc))
Omega_M = 0.3166
Omega_Lam = 0.6847
Omega_k = -0.011
c = 2.998*10**5 #km/s
d_H = c/H0 #Mpc

####################
#FUNCTIONS
####################

def d_L(z):
    ##Luminosity Distance (see Peebles 2013 (pgs 310-321) or https://ned.ipac.caltech.edu/level5/Hogg/Hogg4.html)
    return (1+z)*d_M(z)

def d_M(z):
    return d_C(z)

def d_C(z):
    step = 0.001
    zrange = np.arange(0, z+step, step)
    func = 1/np.sqrt(Om*(1+zrange)**3 + Ol)
    
    integration = 0
    
    for i in range(len(func)-1):
        #Using trapezoid rule
        integration += 0.5*(func[i] + func[i+1])*step
    
    return d_H*integration

'''
peak = np.exp(6.35065)

ans = peak*(((d_L(5))**2)/((d_L(6))**2))

print(np.log(ans))
'''
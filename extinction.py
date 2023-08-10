# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 09:29:52 2022

@author: hmfausey

The extinction module includes equations for the SMC, LMC, and Milky Way 
extinction models based on work from Pei(1998).

Yichuan C. Pei, "Interstellar Dust from the Milky Way to the Magellenic 
Clouds", The Astrophysical Journal, Volume 395, August 1992, Page 130, 
http://doi.org/10.1086/171637
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def xsi(params, lam_emit):
    ##Function that creates base function for extinction profiles.
     #Inputs:
         #params -- 2D array, holds all relevant parameters for the desired
          #extinction model
         #lam_emit -- wavelengths when emitted at the source/in the host galaxy
     #Returns:
         #xsi -- array, extinction profile at each corresponding wavelength for 
          #the current extinction model and input parameters.
          
    xsi = np.zeros(len(lam_emit))   #initialize array
    
    for i in range(len(xsi)):
        #Pei extinction curves rely on a sum over 6 terms each with a different
         #set of parameters
        for j in range(6):
            #get parameters for current sum term
            ai= params[j,0]
            lami = params[j,1]
            bi = params[j,2]
            ni = params[j,3]
            
            #add current xsi term for current wavelength
            xsi[i] += (ai)/((lam_emit[i]/lami)**(ni) + (lami/lam_emit[i])**(ni) + bi)
    
    return xsi

def A_lambda(xsi_vals, E_bv, R_v):
    ##Function determines A_lambda using xsi function
     #Inputs:
         #xsi -- array, list of xsi values determined by input wavelengths and 
          #extinction law
         #E_bv -- scalar, E_b-v value
         #R_v -- scalar, R_v value according to the desired extinction model
     #Returns:
         #A_lam -- array, A_lambda values at corresponding input xsi values
    
    A_lam = xsi_vals*E_bv*(1+R_v)
    
    return A_lam
    
def smc(lam_emit):
    ##Function applies SMC host galaxy extinction model to emitted spectrum
     #Inputs:
         #lam_emit -- wavelengths when emitted at the source/in the host galaxy
     #Returns:
         #xsi -- array, extinction profile at each corresponding wavelength for 
          #the current extinction model and input parameters.
         #R_v -- float, specific R_v = A_v/E_b-v value for the SMC extinction 
          #model
    params = np.array([[185, 0.042, 90, 2.0],[27, 0.08, 5.50, 4.0],[0.005, 0.22, -1.95, 2.0],[0.010, 9.7, -1.95, 2.0],[0.012, 18,-1.8, 2.0],[0.030, 25, 0.00, 2.0]])
    R_v = 2.93 ##For SMC
    
    xsi_array = xsi(params, lam_emit)
    
    return xsi_array, R_v

def lmc(lam_emit):
    ##Function applies LMC host galaxy extinction model to emitted spectrum
     #Inputs:
         #lam_emit -- array, wavelengths when emitted at the source/in the host
          #galaxy
     #Returns:
         #xsi -- array, extinction profile at each corresponding wavelength for 
          #the current extinction model and input parameters.
         #R_v -- float, specific R_v = A_v/E_b-v value for the LMC extinction
          #model
    
    params = np.array([[175, 0.046, 90, 2.0],[19, 0.08, 5.50, 4.5],[0.023, 0.22, -1.95, 2.0],[0.005, 9.7, -1.95, 2.0], [0.006, 18, -1.80, 2.0], [0.020, 25, 0.00, 2.0]])
    R_v = 3.16 #For LMC
    
    xsi_array = xsi(params, lam_emit)
    
    return xsi_array, R_v

def milkyway(lam_emit):
    ##Function applies Milky Way host galaxy extinction model to emitted 
     #spectrum
     #Inputs:
         #lam_emit -- wavelengths when emitted at the source/in the host galaxy
     #Returns:
         #xsi -- array, extinction profile at each corresponding wavelength for 
          #the current extinction model and input parameters.
         #R_v -- float, specific R_v = A_v/E_b-v value for the Milky Way 
          #extinction model
    
    params = np.array([[165, 0.047, 90, 2.0], [14, 0.08, 4.0, 6.5], [0.045, 0.22, -1.95, 2.0],[0.002, 9.7, -1.95, 2.0], [0.002, 18, -1.80, 2.0],[0.012, 25, 0.00, 2.0]])
    R_v = 3.08 #For Milky Way
    
    xsi_array = xsi(params, lam_emit)
    
    return xsi_array, R_v

def transmission(lam_emit, E_bv, extinction_law):
    ##Function calculates the transmission for an array of wavelengths and an
     #extinction law.
     #Inputs:
         #lam_emit -- array, wavelengths when emitted at the source/in the host
          #galaxy
         #E_bv -- scalar, E_{B-V} value
         #extinction_law -- string, indicates desired extinction law:
             #smc -- Small Magellenic Cloud
             #lmc -- Large Magellenic Cloud
             #mw -- Milky Way
          #Default is smc as this extinction law is most prevalent for GRBs
     #Returns:
         #transmission -- transmission coefficient
    
    #get xsi array and R_v values based on extinction law
    if extinction_law == 'smc':
        xsi_array, R_v = smc(lam_emit)
    elif extinction_law == 'lmc':
        xsi_array, R_v = lmc(lam_emit)
    elif extinction_law == 'mw':
        xsi_array, R_v = milkyway(lam_emit)
    else:
        raise Exception("Invalid extinction law. Choices are 'smc', 'lmc', or 'mw' and must be a string")
    
    #Get A_lambda array
    a_lam = A_lambda(xsi_array, E_bv, R_v)
    
    #Calculate transmission based on A_lambda values
    
    transmission = 10**(-0.4 * a_lam)
    
    return transmission
    

def get_extincted_curve(lam_emit, flux, E_bv, extinction_law):
    ##Function returns extincted curve based on E_bv value and desired
     #extinction law
     #Inputs:
         #lam_emit -- array, wavelengths when emitted at the source/in the host
          #galaxy
         #flux -- array, flux at each coresponding emitted wavelengths when
          #emitted at the source
         #E_bv -- scalar, E_{B-V} value
         #extinction_law -- string, indicates desired extinction law:
             #smc -- Small Magellenic Cloud
             #lmc -- Large Magellenic Cloud
             #mw -- Milky Way
          #Default is smc as this extinction law is most prevalent for GRBs
     #Returns:
         #lam_emit -- array, wavelengths when emitted at the source/in the host
          #galaxy
         #flux_new -- array, extincted flux at each corresponding emitted
          #wavelength based on the input extinction law
     
    
    transmission_coeff = transmission(lam_emit, E_bv, extinction_law)
    
    flux_new = transmission_coeff * flux
    
    return lam_emit, flux_new
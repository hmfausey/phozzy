# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:44:58 2022

@author: hmfausey

The attenuation module applies intergalactic attenuation to the GRB spectrum
using a model from Meiksin(2006), which is based off of work by Madau(1995)

Avery Meiksin, "Colour corrections for high-redshift objects due to 
intergalactic attenuation", Monthly Notices of the Royal Astronomical Society, 
Volume 365, Issue 3, January 2006, Pages 807–812, 
https://doi.org/10.1111/j.1365-2966.2005.09756.x

Piero Madau, "Radiative Transfer in a Clumpy Universe: The Colors of 
High-Redshift Galaxies", The Astrophysical Journal, Volume 441, March 1995, 
Page 18, http://doi.org/10.1086/175332
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

from scipy import special

import math

##############################################################################
###################################CONSTANTS##################################
##############################################################################

RYDBERG_H = 1.09678 * 10**(7)   #Rydberg constant for Hydrogen in 1/m
PI = np.pi
LYMAN_L = 0.0912 #In MICROmeters
N = 31      #Highest lyman series transition included in meiksin model. All
             #higher order transition are considered to be negligible.

##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def lyman_series(x):
    ##Function find lyman series out to the n=x->n=1 transition
     #Inputs:
         #x -- scalar, final wavelength transition. This value should not be 
               #less than 2
     #Returns:
         #lym_series -- array, lyman series out to the n=x->n=1 transition
    
    if x<2:
        #Make sure input is 2 or greater
        print("ERROR: lyman_series input must be 2 or greater")
        return "ERROR"
    
    #Initialize lyman series array
    lym_series = np.zeros(x-1)
    
    for i in range(2, x+1):
        #Determine transition wavelength from n=i to n=1
        transition_wavelength = (RYDBERG_H*(1 - (1/(i**2))))**(-1)
        #Add transition wavelength to lyman series array
        lym_series[i-2] = transition_wavelength * 10**(6)
    
    return lym_series

def upper_incomplete_gamma_function(s, x):
    ##For this function I use the built in complete gamma function
     #(scipy.special.gamma) and the normalized upper incomplete gamma function
     #(scipy.special.gammaincc) to comute the value for the unnormalized 
     #upper incomplete gamma function. 
     #NOTE: gammaincc is normalized, meaning that it is the upper incomplete
     #gamma function DIVIDED BY the complete gamma function, so in must be
     #multiplied by the complete gamma function to get the actual value of 
     #the upper incomplete gamma function.
     #Inputs:
         #s -- scalar, gamma function input
         #x -- scalar, the upper incomplete gamma functions lower bound for 
         #integration
     #Returns:
         #upper_incomplete_gamma-- scalar, value of the upper 
          #incomplete gamma function at s and x
    
    complete_gamma = special.gamma(0.5)
    norm_upper_incomplete_gamma = special.gammaincc(0.5, 1) 
    #The above norm_upper_incomplete_gamma_function is the equivalent of
     #(upper incomplete gamma function)/(complete gamma function)
    
    upper_incomplete_gamma = complete_gamma*norm_upper_incomplete_gamma
    
    return upper_incomplete_gamma
    

def teff_below_LLS(lam_obs, z):
    ##Finds optical depth below lyman limit for Lyman Limit Systems (LLS) #Uses 
     #constants as defined in Meiksin 2018
     #Inputs:
         #lam_obs -- array, wavelengths observed at the source
         #z -- scaler, redshift
     #Returns:
         #teff_LLS_array -- array, LLS optical depths at each observed wavelength

    #Define tau_LLS constants from paragraph below eq. 7 in Meiksin (2006)
    n_0 = 0.25
    beta_meiksin = 1.5  #Called the variable beta_meiksin to distinguish
                         #between it and the spectral index
    gamma = 1.5
    
    #initiallize teff array
    teff_LLS_array = np.zeros(len(lam_obs))
    
    i = 0   #initialize index for while loop
    while lam_obs[i] < LYMAN_L*(1+z):
        lam = lam_obs[i]
        
        #Break equation 7 into smaller terms
        
        #Factor in front of the first bracket in eq. 7 of Meiksin(2006)
        front_term = n_0/(4 + gamma - (3*beta_meiksin))
        
        #Initialize the terms for the first and second sum in meiksin equation 7
        sum1 = 0
        sum2 = 0
        
        #loop to compute sums out to 10 terms. Meiksin(2006) notes that the first
         #10 terms of each sum provide a sufficiently high level of convergence
        for n in range(10):
            #Note that sum 1 starts at n=0, but sum 2 starts at n=1, so we need to
             #account for that in our indexing
            sum1 += ((beta_meiksin-1)/(n + 1 - beta_meiksin)) * (((-1)**(n))/(math.factorial(n)))
            
            #Since sum 2 is such a large term, I am going to break it in two pieces
            sum2_term1 = ((beta_meiksin - 1)/((3*(n+1) - gamma - 1) * ((n+1) + 1 - beta_meiksin))) * (((-1)**(n+1))/(math.factorial(n+1)))
            sum2_term2 = ((1+z)**(gamma + 1 - (3*(n+1))) * (lam/LYMAN_L)**(3*(n+1))) - (lam/LYMAN_L)**(gamma + 1)
            sum2 += sum2_term1 * sum2_term2
     
        #Arrange terms by brackets
        
        bracket_1 = upper_incomplete_gamma_function(beta_meiksin - 1, 1) - np.exp(-1) - sum1
        bracket_2 = (((1+z)**(-3*(beta_meiksin-1) + gamma + 1))*((lam/LYMAN_L)**(3*(beta_meiksin-1)))) - ((lam/LYMAN_L)**(gamma + 1))
        
        #Now combining into the first large term, and the second large term
         #(Here I am splitting the two terms by the subtraction sign)
        
        large_term_1 = front_term * bracket_1 * bracket_2
        large_term_2 = n_0 * sum2
       
        #now input into the corresponding index of the LLS optical depth array
        
        teff_LLS_array[i] = large_term_1 - large_term_2
        
        i+= 1
        
    return teff_LLS_array

def teff_below_ISM(lam_obs, z):
    ##Finds optical depth below lyman limit for ISM
     #Inputs:
         #lam_obs -- array, observed wavelengths at the source
         #z -- scalar, redshift
     #Returns:
         #teff_IGM_array -- array, optical depth contributions from optically
          #thin system as given in equation 5 of Meiksin(2006)
    
    teff_IGM_array = np.zeros(len(lam_obs))
    
    i=0 
    while(lam_obs[i] < LYMAN_L*(1+z)):
        lam = lam_obs[i]
        z_L = lam/(LYMAN_L) - 1 #As defined above equation 5 in Meiksin(2006)
        
        teff_IGM_array[i] = 0.805*(1+z_L)**(3) * (1/(1+z_L) - 1/(1+z))
        i+= 1
    return teff_IGM_array

def teff_above(lam_obs, z, n_tot):
    ##Finds optical depth as a function of observed wavelength above 
     #lyman limit. 
     #Inputs:
         #lam_obs -- array, observed wavelengths
         #z, -- scalar, redshift
         #n_tot -- scalar, highest lyman transition used to calculate optical
                   #depth avoce the lyman limit
     #Returns:
         #teff_array -- array, optical depth contributions due to resonant 
                        #scattering by lyman transitions
    
    #Initialize optical depth array
    teff_array = np.zeros(len(lam_obs))
    #Get array of lyman series transitions out to n=n_tot -> n=1
    lym_series = lyman_series(n_tot)
    #Initialize optical depth parameter
    teff = 0 
        
    #below is an array of coefficents for the ratios of n=3-9 -> n=1
     #transition optical depths to the n=2 -> n=1 transition optical depth.
     #These can be found in table 1 of Meiksin(2006)
    transition_ratio_coeff = np.array([0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373, 0.0283])
    
    #Loop will determine the optical depth at each corresponding wavelength
    for i in range(len(lam_obs)):
        lam = lam_obs[i]
        #if lam >= LYMAN_L*(1+z):
        for n in range(2, N+1):
            #the z_n determines whether a given transition affect the optical depth 
             #at a given observed wavelength. The first one (z_2) determines
             #whether the lyman-alpha transition has an effect on the optical depth
             #at the given observed wavelength
            #z_n given in Meiksin 2006 caption before table 1
            lam_n = lym_series[n-2]
            z_n = (lam/lam_n) - 1
            
            #First need to find the tau_alpha(aka tau_2) optical depth (corresponds 
             #to the n=2 -> n=1 transition). Uses equations 2 and 3 from 
             #Meiksin(2006). Equation depends only on source redshift. This value will
             #remain the same for every wavelength we try to find the optical depth for,
             #so we will keep it outside of the loop and use the variable tau_alpha to
             #keep track of the value and apply it where necessary
            
            if (0 <= z < 4):
                tau_alpha = 0.00211 * (1+z_n)**(3.7)
            elif (z > 4):
                tau_alpha = 0.00058 * (1+z_n)**(4.5)
            elif (z==4):
                #While the z=4 value using the equations for z<4 and z>4 are very
                 #similar, there is still a slight discontinuity going from one
                 #function to the other, so here we will take the average of these two
                 #values in the case of z=4
                tau_alpha_low = 0.00211 * (1+z_n)**(3.7)
                tau_alpha_high = 0.00058 * (1+z_n)**(4.5)
                tau_alpha = (tau_alpha_low + tau_alpha_high)/2
            else:
                print("ERROR: z must be greater than or equal to 0")
                return "ERROR"
            
            #If a given z_n is greater than or equal to z, all higher order z_n
             #terms will also be greater than or equal to z, so this condition
             #allows the loop to close once the first z_n is not less than z in
             #an attempt to save computing time
            if z_n < z:
                #First determine whether or not tau_alpha will have an effect
                if n == 2:
                    teff += tau_alpha
                #Now, use equations for tau_n ratios for n=x -> n=1 transitions for x=3 to 
                 #x=9. These equations can be found in table 1 of Meiksin(2006)
                elif n <= 5:
                    #In this case the two functions listed for the n=3, 4, and 
                     #5 give the exact same value at z=3, so we don't have to 
                     #add an extra condition here
                    if z_n <= 3:
                        tau_n = ((transition_ratio_coeff[n-3]*(0.25*(1+z_n)))**(1/3))*tau_alpha
                        teff += tau_n
                    else:
                        tau_n = ((transition_ratio_coeff[n-3]*(0.25*(1+z_n)))**(1/6))*tau_alpha
                        teff += tau_n
                elif n <= 9:
                    #equations up through 9 are listed in table 1 of Meiksin(2006)
                    tau_n = ((transition_ratio_coeff[n-3]*(0.25*(1+z_n)))**(1/3))*tau_alpha
                    teff += tau_n
                    #Need to save tau_9 (aka tau_theta) for higher order terms
                    if n==9:
                        tau_theta = tau_n
                elif n <= N:
                    tau_n = (720/(n*((n**2) - 1)))*tau_theta
                    teff += tau_n
            else:
                break
            
        teff_array[i] = teff    #save optical depth in corresponding index
        teff = 0    #set teff back down to 0
    
    return teff_array
        
def teff_total(lam_emit, z):
    ##Returns the total optical depth at each wavelength for all possible
     #contributing factors
     #Inputs:
         #lam_emit -- array, emitted wavelengths
         #z -- scalar, redshift
     #Returns:
         #lam_obs -- array, observed wavelengths
         #teff -- array, optical depths at each observed wavelength
         
    lam_obs = lam_emit*(1+z)
    
    tau_above = teff_above(lam_obs, z, N)            
    tau_lls = teff_below_LLS(lam_obs, z)    
    tau_ism = teff_below_ISM(lam_obs, z)
    
    teff = tau_above + tau_lls + tau_ism
    
    return lam_obs, teff

def attenuate(lam_emit, flux, z):
    ##Returns attenuated curve for input redshift
     #Inputs:
         #lam_emit -- array, emitted wavelengths
         #flux -- array, flux at corresopnding emitted wavelengths
         #z -- scalar, redshift
     #Returns:
         #lam_obs -- array, observed wavelengths
         #atten_flux -- array, attenuated flux at each observed wavelength
         
         lam_obs, teff = teff_total(lam_emit, z)
         
         transmission = np.exp(-teff)
         
         atten_flux = transmission*flux
         
         return lam_obs, atten_flux

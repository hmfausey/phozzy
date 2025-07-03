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

For better speed, a lot of the Meiksin equations have been reduced to their 
numerical values (e.g., in Eq. 7, the first normalization (N0/(4 + gamma - 3*beta)
is just 0.25. The first sum is only made up of constants, so its value will
never change regardless of wavelength. Additionally, the upper incomplete gamma
function inputs will never change, so we can directly insert that value rather 
than having to calculate it. These are just a few examples, but a lot of the 
meiksin equations have been reduced into the corresponding values rather than
symbols and functions).

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

RYDBERG_H = 1.09678 * 10 ** (7)  # Rydberg constant for Hydrogen in 1/m
PI = np.pi
LYMAN_L = 0.0912  # In MICROmeters


##############################################################################
###################################FUNCTIONS##################################
##############################################################################


def lyman_series(x):
    #Find the wavelength associated with a specific lyman series transition
    if x < 2:
        # Make sure input is 2 or greater
        raise Exception("ERROR: lyman_series input must be 2 or greater")
    else:
        transition_wavelength = (RYDBERG_H * (1 - (1 / (x**2)))) ** (-1)
        
        return transition_wavelength * 10 ** (6)
        
        

def lyman_series_array(x):
    ##Function find lyman series out to the n=x->n=1 transition
    # Inputs:
    # x -- scalar, final wavelength transition. This value should not be
    # less than 2
    # Returns:
    # lym_series -- array, lyman series out to the n=x->n=1 transition

    if x < 2:
        # Make sure input is 2 or greater
        raise Exception("ERROR: lyman_series input must be 2 or greater")

    # Initialize lyman series array
    lym_series_array = np.zeros(x - 1)

    for i in range(2, x + 1):
        # Determine transition wavelength from n=i to n=1
        transition_wavelength = (RYDBERG_H * (1 - (1 / (i**2)))) ** (-1)
        # Add transition wavelength to lyman series array
        lym_series_array[i - 2] = transition_wavelength * 10 ** (6)

    return lym_series_array


def upper_incomplete_gamma_function(s, x):
    ##For this function I use the built in complete gamma function
    # (scipy.special.gamma) and the normalized upper incomplete gamma function
    # (scipy.special.gammaincc) to comute the value for the unnormalized
    # upper incomplete gamma function.
    # NOTE: gammaincc is normalized, meaning that it is the upper incomplete
    # gamma function DIVIDED BY the complete gamma function, so in must be
    # multiplied by the complete gamma function to get the actual value of
    # the upper incomplete gamma function.
    # OTHER NOTE: This function is no longer used in the actual intergalactic
    # attenuation calculation, but is useful for testing purposes, and for those
    # who would like to verify that the value of the gamma function in teff_below_LLS
    # is correct
    # Inputs:
    # s -- scalar, gamma function input
    # x -- scalar, the upper incomplete gamma functions lower bound for
    # integration
    # Returns:
    # upper_incomplete_gamma-- scalar, value of the upper
    # incomplete gamma function at s and x

    complete_gamma = special.gamma(s)
    norm_upper_incomplete_gamma = special.gammaincc(s, x)
    # The above norm_upper_incomplete_gamma_function is the equivalent of
    # (upper incomplete gamma function)/(complete gamma function)

    upper_incomplete_gamma = complete_gamma * norm_upper_incomplete_gamma

    return upper_incomplete_gamma

def sum_alpha():
    #In the Mieksin paper it mentions that both sums in equation 7 are sufficiently
    #converged after the first 10 terms, so we calculated this sum  for n = 0 - 9
    #NOTE: this function is not actively used in the attenuation calculation, but is
    #useful for testing purposes, and for those who would like to verify the sum
    
    sumAlpha = 0
    for n in range(10):
        sumAlpha += (0.5/(n-0.5)) * ((-1)**((n))/math.factorial(n))
        
    return sumAlpha

def diff_lambda_1(lam, z):
    #Calculates a piece of Meiksin (2006) eq. 7. 
    # Inputs:
    # lam_obs -- array, wavelengths observed at the source
    # z -- scaler, redshift
    # Returns:
    # term -- array, calculated peice of Meiksin equation for each input wavelenth
    
    #if you calculated out the exponent of the 1+z term, its just 1
    #The exponents for the two lambda terms are 1.5 and 2.5 respectively
    
    term = ((1+z) * (lam/LYMAN_L)**(1.5)) - (lam/LYMAN_L)**(2.5)
    
    return term

def sum_beta(lam, z):
    
    #Calculates a piece of Meiksin (2006) eq. 7. 
    # Inputs:
    # lam_obs -- array, wavelengths observed at the source
    # z -- scaler, redshift
    # Returns:
    # beta_sum -- array, calculated peice of Meiksin equation for each input wavelenth
    
    #This sum is a bit more complicated
    
    beta_sum = 0
    
    for n in range(1, 11):
        #Note, this sum starts at n=1 and goes to 10 (since the range upper bound is not
         #inclusive, its 11 in the function)
        front_term = (1/((3*n - 2.5)*(2*n - 1))) * ((-1)**n)/math.factorial(n)
        lam_difference_term = ((1+z)**(2.5 - 3*n))*((lam/LYMAN_L)**(3*n)) - (lam/LYMAN_L)**(2.5)
        
        beta_sum += front_term*lam_difference_term
    
    return beta_sum

def teff_below_LLS(lam_obs, z):
    ##Finds optical depth below lyman limit for Lyman Limit Systems (LLS) #Uses
    # constants as defined in Meiksin 2018
    # Inputs:
    # lam_obs -- array, wavelengths observed at the source
    # z -- scaler, redshift
    # Returns:
    # teff_LLS_array -- array, LLS optical depths at each observed wavelength

    # Define tau_LLS constants from paragraph below eq. 7 in Meiksin (2006)
    n_0 = 0.25
    beta_meiksin = 1.5  # Called the variable beta_meiksin to distinguish
    # between it and the spectral index
    gamma = 1.5
    
    #Converting equation 7 to be of the form:
        #n_0*([UIGF_val - exp(-1) - sumAlpha]*[diffLambda1] - sumBeta)
        #There is only one normalization here because N0/(4 + gamma - 3beta) = N0 = 0.25, so
        #we can pull these terms out and just multiply the whole thing by n_0 = 0.25
        #Sum beta is large and does have a dependence on wavelength (see below)
    
    
    UIGF_val = 0.27880558528065474 #calculated using the upper_incomplete_gamma_function function above
    
    sumAlpha = -1.861527720191936 #calculated using the sum_alpha function
    
    #diff_lambda_1 is line 2 of equation 7 (of the form [(1+z)^(value)*(lambda/lambda_L)^(value) - (lambda/lambda_L)^(value)]
    
    diffLambda1 = diff_lambda_1(lam_obs, z)
    
    sumBeta = sum_beta(lam_obs, z)
    
    teff_LLS_array = n_0 * ( ((UIGF_val - np.exp(-1) - sumAlpha)*diffLambda1) - sumBeta )
    
    end_contr = np.argwhere(teff_LLS_array < 0)
    
    if len(end_contr)==0:
        return teff_LLS_array
    
    end_arg = end_contr.min()
    teff_LLS_array[end_arg:] = np.zeros(len(teff_LLS_array[end_arg:]))
    
    return teff_LLS_array


def teff_below_ISM(lam_obs, z):
    ##Finds optical depth below lyman limit for ISM
    # Inputs:
    # lam_obs -- array, observed wavelengths at the source
    # z -- scalar, redshift
    # Returns:
    # teff_IGM_array -- array, optical depth contributions from optically
    # thin system as given in equation 5 of Meiksin(2006)
    
    z_L = (lam_obs/LYMAN_L) - 1
    
    teff_IGM_array = 0.805 * ((1 + z_L)**(3)) * ((1/(1 + z_L)) - (1/(1 + z)))
    
    end_contr = np.argwhere(teff_IGM_array < 0)
    
    if len(end_contr)==0:
        return teff_IGM_array
    
    end_arg = end_contr.min()
    teff_IGM_array[end_arg:] = np.zeros(len(teff_IGM_array[end_arg:]))
    
    return teff_IGM_array


def get_alpha_optical_depth(lam_obs, z):
    #Gets the optical depth from Ly-alpha absorption as a function of z_alpha(wavelength)
    #Inputs:
        #lam_obs - array, observed wavelength
        #z - float, redshift
    #Returns:
        #tau_alpha - optical depth at each wavelength
    lam_alpha = lyman_series(2)
    
    z_alpha = (lam_obs/lam_alpha) - 1
    tau_alpha = np.zeros(len(z_alpha))

    
    if z_alpha[0] >= z:

        #if the input redshift falls below the lowest wavelength in the spectrum
         #then there is no absorption from ly-alpha(or any other) transition
        return tau_alpha
    
    elif z_alpha[-1] <= 4:
        #if we never look at wavelengths corresponding to a redshift above 4, 
         #we never need to worry about it

        low_contr = np.argwhere(z_alpha < z)
   
        low_contr_index = low_contr.max()

        
        tau_alpha[:low_contr_index + 1] = 0.00211 * (1 + z_alpha[:low_contr_index + 1]) ** (3.7)

        return tau_alpha
    else:
        #This should cover all other cases
        
        #below this z_alpha is < 4
        low_index = np.argwhere(4 < z_alpha).min()
        #above this is z_alpha > z
        high_index = np.argwhere(z_alpha < z).max()
        
        #add below 
        tau_alpha[:low_index] = 0.00211 * (1 + z_alpha[:low_index]) ** (3.7)
        
        #add above 4 but below z
        tau_alpha[low_index:high_index+1] = 0.00058 * ((1 + z_alpha[low_index:high_index+1])**(4.5))
        
        #make sure >=z contributions are 0
        tau_alpha[high_index + 1:] = np.zeros(len(tau_alpha[high_index + 1:]))
        
    return tau_alpha


def get_transition_optical_depth_terms(lam_obs, z, n, tau_alpha):
    ##Finds the optical depth associated with each transition at each wavelength
    #Inputs:
        #lam_obs - array, observed wavelength
        #z - float, redshift
        #n - int, transition number
        #tau_alpha - array, optical depth due to ly-alpha transition
    #Returns:
        #tau_n - optical depth at each wavelength
    transition_ratio_coeff = np.array(
        [0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373, 0.0283]
    )
    
    lam_n = lyman_series(n)
    z_n = (lam_obs/lam_n) - 1
    
    tau_n = np.zeros(len(z_n))
    
    if (z_n[0] > z):
        #if the input redshift falls below the lowest wavelength in the spectrum
         #then there is no absorption from this transition
        return tau_n
    
    n_contr_index = np.argwhere(z_n < z).max()
    
    tau_alpha_val = tau_alpha[n_contr_index]
    tau_alpha_temp = np.ones(len(lam_obs))*tau_alpha_val
    
    
    if (n <= 5) and (z<=3):
        tau_n[:n_contr_index+1] = tau_alpha_temp[:n_contr_index+1] * transition_ratio_coeff[n-3]*((1 + z_n[:n_contr_index+1])**(1/6))
        return tau_n
    
    
    tau_n[:n_contr_index+1] = tau_alpha_temp[:n_contr_index+1] * transition_ratio_coeff[n-3]*((1 + z_n[:n_contr_index+1])**(1/3))

    return tau_n
        
    
def get_higher_order_optical_depth(lam_obs, z, n, tau_theta):
    ##Finds the optical depth associated with higher order transitions (n > 9) at each wavelength
    #Inputs:
        #lam_obs - array, observed wavelength
        #z - float, redshift
        #n - int, transition number
        #tau_theta - array, optical depth due to ly-theta transition
    #Returns:
        #tau_n - optical depth at each wavelength
        
        
    # the z_n term determines whether a given transition affects the optical depth
    # at a given observed wavelength. The first one (z_2) determines
    # whether the lyman-alpha transition has an effect on the optical depth
    # at the given observed wavelength z_n given in Meiksin 2006 caption before table 1
    lam_n = lyman_series(n)
    z_n = (lam_obs/lam_n) - 1
    
    tau_n = np.zeros(len(z_n))
    
    if (z_n[0] > z):
        #if the input redshift falls below the lowest wavelength in the spectrum
         #then there is no absorption from this transition
        return tau_n
    
    n_contr_index = np.argwhere(z_n < z).max()
    
    tau_n[:n_contr_index + 1] = tau_theta[:n_contr_index + 1] * (720/(n*(n**2 - 1)))
   
    return tau_n
    

def teff_above(lam_obs, z):
    ##Finds optical depth as a function of observed wavelength above
    # lyman limit.
    # Inputs:
    # lam_obs -- array, observed wavelengths
    # z, -- scalar, redshift
    # n_tot -- scalar, highest lyman transition used to calculate optical
    # depth avoce the lyman limit
    # Returns:
    # teff_array -- array, optical depth contributions due to resonant
    # scattering by lyman transitions

    
    N = 31  # Highest lyman series transition included in meiksin model. All
    # higher order transition are considered to be negligible.
    
    # Initialize optical depth array
    teff = np.zeros(len(lam_obs))

    # Get array of lyman series transitions out to n=n_tot -> n=1
    lym_series = lyman_series(N)
    
    #get tau_alpha (optical depth due to Ly-alpha)
    tau_alpha = get_alpha_optical_depth(lam_obs, z)

    
    if np.all(tau_alpha == 0):
        #if the ly-alpha transition is not contributing to the dropoff within our 
         #wavelength range, then neither is anything else and we can just return
         #an array of zero
        return teff
    
    #otherwise, add tau-alpha contribution
    teff += tau_alpha
    if z<=4:
        tau_alpha_val = 0.00211 * (1 + z) ** (3.7)
    else:
        tau_alpha_val = 0.00058 * ((1 + z)**(4.5))  
        
    tau_alpha_base = tau_alpha_val*np.ones(len(lam_obs))
    
    for n in range(3, N+1):

        #add contributions for each transition. If any come back as all zeros, 
         #end the loop
        if n < 9:
            tau_n = get_transition_optical_depth_terms(lam_obs, z, n, (tau_alpha))# + tau_alpha_base)/2)

        elif n == 9:
            tau_n = get_transition_optical_depth_terms(lam_obs, z, n, (tau_alpha))# + (tau_alpha_base)/125))
            tau_theta = tau_n
        else:
            tau_n = get_higher_order_optical_depth(lam_obs, z, n, tau_theta)
        
        #otherwise we want to add the tau_n contribution and continue the loop
        
        teff += tau_n
        
  
    return teff


def teff_total(lam_emit, z):
    ##Returns the total optical depth at each wavelength for all possible
    # contributing factors
    # Inputs:
    # lam_emit -- array, emitted wavelengths
    # z -- scalar, redshift
    # Returns:
    # lam_obs -- array, observed wavelengths
    # teff -- array, optical depths at each observed wavelength

    lam_obs = lam_emit * (1 + z)
    
    #Split the array and pass them to the correct areas 
    split_array = np.argwhere(lam_obs > LYMAN_L*(1+z))

    
    if len(split_array)==0:
        split = len(lam_obs)  
    else:   
        split = split_array.min()

    lam_low_contr = lam_obs[:split]
    
    tau_above = teff_above(lam_obs, z)
    
    tau_lls = teff_below_LLS(lam_low_contr, z)

    tau_ism = teff_below_ISM(lam_low_contr, z)

    low_total = tau_lls + tau_ism
    
    teff = np.zeros(len(tau_above))
    
    teff += tau_above
    teff[:split] += low_total
      
    
    return lam_obs, teff


def attenuate(lam_emit, flux, z):
    ##Returns attenuated curve for input redshift
    # Inputs:
    # lam_emit -- array, emitted wavelengths
    # flux -- array, flux at corresopnding emitted wavelengths
    # z -- scalar, redshift
    # Returns:
    # lam_obs -- array, observed wavelengths
    # atten_flux -- array, attenuated flux at each observed wavelength

    lam_obs, teff = teff_total(lam_emit, z)

    transmission = np.exp(-teff)

    atten_flux = transmission * flux

    return lam_obs, atten_flux

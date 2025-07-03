#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 15:21:42 2025

@author: hmfausey

Some of the attenuation models are too complex for testing, but there are tests included 
for everything that can be tested with a moderate amount of effort
"""

import pytest
import numpy as np
from phozzy import attenuation
import math
from phozzy import power_law

def test_lyman_series():
    real_lyman_series = [0.12156701, 0.10257220, 0.0972536, 0.094974287, 0.093780331, 0.0930748142, 0.0926225605, 0.0923150275, 0.0920963006, 0.0919351334]
    
    for i in range(2, 12):
        np.testing.assert_allclose(attenuation.lyman_series(i), real_lyman_series[i-2], 0.0001)
    
    np.testing.assert_raises_regex(Exception, "ERROR: lyman_series input must be 2 or greater", attenuation.lyman_series, 1)

def test_lyman_series_array():
    #first 10 terms of lyman series -- according to wikipedia
    real_lyman_series = [0.12156701, 0.10257220, 0.0972536, 0.094974287, 0.093780331, 0.0930748142, 0.0926225605, 0.0923150275, 0.0920963006, 0.0919351334]
    
    output_lyman_series = attenuation.lyman_series_array(11)
    
    #This should be pretty close, so using a small tolerance
    np.testing.assert_allclose(output_lyman_series, real_lyman_series, 0.0001)
    
    #also check the exception catches issues
    
    
    np.testing.assert_raises_regex(Exception, "ERROR: lyman_series input must be 2 or greater", attenuation.lyman_series_array, -5)

    
def test_upper_incomplete_gamma_function():
    
    #for UIG with inputs (s+1, 1) where z is a positive integer, the output should be floor(e * (s!))/e
    
    #First, checking this for values of s ranging from 1 to 9 (so inputs will be 2-10)
    
    s = np.arange(1,10, 1)
    
    expected = np.zeros(len(s))
    output = np.zeros(len(s))
    
    for i in range(len(expected)):
        expected[i] = float(math.floor(math.factorial(s[i])*np.e))/np.e
        output[i] = attenuation.upper_incomplete_gamma_function(s[i]+1, 1)
        
    np.testing.assert_allclose(expected, output, 0.001)
    
    
    #if x not equal to 1, can use a different sum. for s = 3, for example:
        #expect -- 2! * e^(-x) * (1 + x + x^2/2) = 2exp(-x)*(1 + x + (x^2/2))
        #can test this with a few values. x does NOT have to be an integer in this case
        
    x = np.linspace(0, 10, 100)
    
    expected = 2*np.exp(-x) * (1 + x + (x**2)/2)
    
    output = attenuation.upper_incomplete_gamma_function(3, x)
    
    np.testing.assert_allclose(expected, output, 0.0001)
    
    #even simpler, if s = 2...
    
    expected = np.exp(-x) * (1 + x)
    
    output = attenuation.upper_incomplete_gamma_function(2, x)
    
    np.testing.assert_allclose(expected, output, 0.0001)
    
    #Finally, for s = 1, should just equal e^-x
    
    x = np.linspace(0, 500, 25000)
    
    expected = np.exp(-x)
    output = attenuation.upper_incomplete_gamma_function(1, x)
    
    np.testing.assert_allclose(expected, output, 0.0001)
    
def test_sum_alpha():
    ##Here is a less convoluted way of writing out the sum
    
    expected = 0
    for n in range(10):
        expected += 1/(2*n - 1) * ((-1)**n)/math.factorial(n)
    
    #Lets see if it matches up with what we get when we use the form in the
     #meiksin paper
     
    np.testing.assert_equal(attenuation.sum_alpha(), expected)
    
def test_diff_lambda_1():
    lam = 0.0912 #the Lyman limit wavlength in MICROmeters
    z = np.linspace(0, 10, 100)
    output = attenuation.diff_lambda_1(lam, z)
    
    #Because of how things cancel, the output should just equal z
    #set the tolerance very tight here
    np.testing.assert_allclose(output, z, 1e-8)
    
    
    lam = 4*lam
    #For lam = 4*(lyman limit wavelength) should get (1+z)*8 - 32, aka 8z - 24
    
    output = attenuation.diff_lambda_1(lam, z)
    
    np.testing.assert_allclose(output, 8*z - 24, 1e-8)
    
def test_sum_beta():
    #This is what the sum used to look like. This checks that rewriting the equation
     #did not change how it functions
     
    LYMAN_L = 0.0912  # In MICROmeters
    lam = np.linspace(0.5, 2.4, 1000)
    z=3
    
    beta_meiksin = 1.5  # Called the variable beta_meiksin to distinguish
    # between it and the spectral index
    gamma = 1.5
    
    sum2 = 0

    # loop to compute sums out to 10 terms. Meiksin(2006) notes that the first
    # 10 terms of each sum provide a sufficiently high level of convergence
    for n in range(10):
        #Note that sum2 starts at 1, so accounting for that by replacing n with n+1
        # Since sum 2 is such a large term, I am going to break it in two pieces
        sum2_term1 = (
            (beta_meiksin - 1)
            / ((3 * (n + 1) - gamma - 1) * ((n + 1) + 1 - beta_meiksin))
        ) * (((-1) ** (n + 1)) / (math.factorial(n + 1)))
        sum2_term2 = (
            (1 + z) ** (gamma + 1 - (3 * (n + 1)))
            * (lam / LYMAN_L) ** (3 * (n + 1))
        ) - (lam / LYMAN_L) ** (gamma + 1)
        sum2 += sum2_term1 * sum2_term2
         
    
    output = attenuation.sum_beta(lam, z)
    
    #using a very tight tolerance
    np.testing.assert_allclose(sum2, output, 1e-8)
        
def test_teff_below_LLS():
    #This is what the teff_below_LLS function originally looked like when calculating
     #the sums and terms as they appear in Meiksin (2006) Eq. 7. This test will verify
     #that rewriting the function has not altered the results
     
    LYMAN_L = 0.0912  # In MICROmeters
    lam_obs = np.linspace(0.5, 2.4, 1000)
    z=10
    
    
    output = attenuation.teff_below_LLS(lam_obs, z)
    
    # Define tau_LLS constants from paragraph below eq. 7 in Meiksin (2006)
    n_0 = 0.25
    beta_meiksin = 1.5  # Called the variable beta_meiksin to distinguish
    # between it and the spectral index
    gamma = 1.5
    
    # initiallize teff array
    teff_LLS_array = np.zeros(len(lam_obs))

    i = 0  # initialize index for while loop
    while (lam_obs[i] < LYMAN_L * (1 + z)) and (i < len(lam_obs)):
        lam = lam_obs[i]

        # Break equation 7 into smaller terms

        # Factor in front of the first bracket in eq. 7 of Meiksin(2006)
        front_term = n_0 / (4 + gamma - (3 * beta_meiksin))

        # Initialize the terms for the first and second sum in meiksin equation 7
        sum1 = 0
        sum2 = 0

        # loop to compute sums out to 10 terms. Meiksin(2006) notes that the first
        # 10 terms of each sum provide a sufficiently high level of convergence
        for n in range(10):
            # Note that sum 1 starts at n=0, but sum 2 starts at n=1, so we need to
            # account for that in our indexing
            sum1 += ((beta_meiksin - 1) / (n + 1 - beta_meiksin)) * (
                ((-1) ** (n)) / (math.factorial(n))
            )

            # Since sum 2 is such a large term, I am going to break it in two pieces
            sum2_term1 = (
                (beta_meiksin - 1)
                / ((3 * (n + 1) - gamma - 1) * ((n + 1) + 1 - beta_meiksin))
            ) * (((-1) ** (n + 1)) / (math.factorial(n + 1)))
            sum2_term2 = (
                (1 + z) ** (gamma + 1 - (3 * (n + 1)))
                * (lam / LYMAN_L) ** (3 * (n + 1))
            ) - (lam / LYMAN_L) ** (gamma + 1)
            sum2 += sum2_term1 * sum2_term2

        # Arrange terms by brackets

        bracket_1 = (
            attenuation.upper_incomplete_gamma_function(beta_meiksin - 1, 1) - np.exp(-1) - sum1
        )
        bracket_2 = (
            ((1 + z) ** (-3 * (beta_meiksin - 1) + gamma + 1))
            * ((lam / LYMAN_L) ** (3 * (beta_meiksin - 1)))
        ) - ((lam / LYMAN_L) ** (gamma + 1))

        # Now combining into the first large term, and the second large term
        # (Here I am splitting the two terms by the subtraction sign)

        large_term_1 = front_term * bracket_1 * bracket_2
        large_term_2 = n_0 * sum2

        # now input into the corresponding index of the LLS optical depth array

        teff_LLS_array[i] = large_term_1 - large_term_2
        i += 1
    
    expected = teff_LLS_array
    
    np.testing.assert_allclose(expected, output, 1e-8)
        
    
def test_teff_below_ISM():
    
    LYMAN_L = 0.0912  # In MICROmeters
    lam_obs = np.linspace(0.5, 2.4, 1000)
    z=10
    
    output = attenuation.teff_below_ISM(lam_obs, z)
    #What it used to look like:
    
    teff_IGM_array = np.zeros(len(lam_obs))

    i = 0
    while (lam_obs[i] < LYMAN_L * (1 + z)) and (i < len(lam_obs) - 1):
        lam = lam_obs[i]
        z_L = lam / (LYMAN_L) - 1  # As defined above equation 5 in Meiksin(2006)

        teff_IGM_array[i] = 0.805 * (1 + z_L) ** (3) * (1 / (1 + z_L) - 1 / (1 + z))
        i += 1
    
    np.testing.assert_allclose(output, teff_IGM_array, 1e-8)
 
        
def test_teff_total():
    #From Meiksin table
    
    lam_obs = np.array([0.173, 0.1731, 0.1732, 0.1733, 0.1734])
    
    #For z = 2
    z = 2
    
    lam_emit = lam_obs/(1+z)
    
    
    expected = np.array([0.150055, 0.150066, 0.150077, 0.150089, 0.150101])
    #Need to change the line below
    lam_obs, output = attenuation.teff_total(lam_emit, z)
    
    np.testing.assert_allclose(expected, np.exp(-output), 0.005)
    
    #For z = 3.5
    z = 3.5
    
    lam_emit = lam_obs/(1+z)
    
    
    expected = np.array([0.020407, 0.020369, 0.02033 , 0.020292, 0.020253])
    #Need to change the line below
    lam_obs, output = attenuation.teff_total(lam_emit, z)
    
    np.testing.assert_allclose(expected, np.exp(-output), 0.005)
    
    #For z = 5
    
    z = 5
    
    lam_emit = lam_obs/(1+z)
    
    
    expected = np.array([0.003277, 0.003266, 0.003254, 0.003243, 0.003231])
    #Need to change the line below
    lam_obs, output = attenuation.teff_total(lam_emit, z)
    
    #for very small values I am relaxing the relative tolerance a bit
    np.testing.assert_allclose(expected, np.exp(-output), 0.01)
        
    #for z = 7
    
    z = 7
    
    lam_emit = lam_obs/(1+z)
    
    
    expected = np.array([0.000310, 0.000308, 0.000306, 0.000304, 0.000303])
    #Need to change the line below
    lam_obs, output = attenuation.teff_total(lam_emit, z)
    
    np.testing.assert_allclose(expected, np.exp(-output), 0.01)
         
def test_attenuate():
    
    lam_obs = np.linspace(0.5, 2.4, 1000)
    
    z = 0
    
    lam_emit, flux = power_law.fv_obs(lam_obs, 2.05, 100, 0.7, z)
    
    lam_obs_new, flux_atten = attenuation.attenuate(lam_emit, flux, z)
    
    np.testing.assert_equal(lam_obs, lam_obs_new)
    np.testing.assert_equal(flux, flux_atten)
    
    z = 20
    lam_emit, flux = power_law.fv_obs(lam_obs, 2.05, 100, 0.7, z)
    
    lam_obs_new, flux_atten = attenuation.attenuate(lam_emit, flux, z)
    
    np.testing.assert_allclose(lam_obs, lam_obs_new, 1e-10)
    #This is set to absolute instead of relative tolerance because everthing is
     #100% different from 0
    np.testing.assert_allclose(np.zeros(len(flux_atten)), flux_atten, atol = 1e-10)
    
    z = 4.43 #This will catch higher order lyman terms on the edge of the wavelength
     #range
    
    lam_emit, flux = power_law.fv_obs(lam_obs, 2.05, 100, 0.7, z)
    
    lam_obs_new, flux_atten = attenuation.attenuate(lam_emit, flux, z)
    
    
    lam_obs, teff = attenuation.teff_total(lam_emit, z)
    flux_real = flux*np.exp(-teff)
    
    np.testing.assert_allclose(flux_atten, flux_real, 0.01)
    
 
def test_teff_above_specific():
    
    #Check that the function ends early if only Ly-alpha contributes to the 
     #optical depth within the wavelength range we are looking at
    
    #Also, getting right up the the edge of the wavelength range might iluminate
     #anything funky about how things are being calculated right around the redshift
     #or the edge of the wavelength range
    z = 3.5
    lam_obs = np.linspace(0.5, 2.4, 1000)
    
    teff = attenuation.teff_above(lam_obs, z)
    
    #Note that the get_alpha_optical_depth() function has already been covered 
     #in full by tests above which will run first, so if there is a problem 
     #with this function it shoud (hopefully) have been caught earlier
    teff_actual = attenuation.get_alpha_optical_depth(lam_obs, z)
    
    #For this wavelength range, it should have stopped calculating the optical
     #depth at the ly-beta contribution, so it should just be equal to teff_alpha
     
    np.testing.assert_equal(teff, teff_actual)
  
 
         
         
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:23:06 2024

@author: hmfausey
"""


import pytest
import numpy as np
import sys
sys.path.append("../")

from phozzy import power_law as pl



def test_fv_emit():
    #test output for generic values
    lam_emit = np.array([0.5, 1, 2])
    lam_0 = 1
    f_0 = 1
    beta = 1
    output = pl.fv_emit(lam_emit, lam_0, f_0, beta)
    expected = [0.5, 1, 2]
    np.testing.assert_equal(output, expected)
    
    #if f_0 is doubled, so should the flux
    f_0 = 2
    output = pl.fv_emit(lam_emit, lam_0, f_0, beta)
    expected = [1, 2, 4]
    
    np.testing.assert_equal(output, expected)
    
    #if normalized at 2, the ratio of the flux density at one point vs another should stay 
     #the same but the flux should be 1 where lambda = 2.
    
    f_0 = 1
    lam_0 = 2
    
    output = pl.fv_emit(lam_emit, lam_0, f_0, beta)
    expected = [0.25, 0.5, 1]
    
    np.testing.assert_equal(output, expected)
    
    lam_0 = 1
    beta = 2
    
    #For beta = 2, the original output fluxes from the first case should be squared
    
    output = pl.fv_emit(lam_emit, lam_0, f_0, beta)
    expected = [0.25, 1, 4]
    
    np.testing.assert_equal(output, expected)
       
    
    
def test_fv_obs():
    
    lam_obs = np.array([2, 4, 8])
    lam_0 = 1
    f_0 = 1
    beta = 1
    z = 1
    
    #For a redshift of 1, we would expect the following flux-density output:
    f_expected = [2, 4, 8]
    #^^This is because we are also adjusting the lam_0 wavelength so it should stay
        #the same
    
    #but lam emit should be different:
    lam_expected = [1, 2, 4]
    
    output_lam, output_f = pl.fv_obs(lam_obs, lam_0, f_0, beta, z)
    np.testing.assert_equal([output_lam, output_f], [lam_expected, f_expected])
    
    #For a redshift of 3 instead, we would expect
    z = 3
    f_expected = [2, 4, 8]
    lam_expected = [0.5, 1, 2]
    
    output_lam, output_f = pl.fv_obs(lam_obs, lam_0, f_0, beta, z)
    np.testing.assert_equal([output_lam, output_f], [lam_expected, f_expected])
    
    

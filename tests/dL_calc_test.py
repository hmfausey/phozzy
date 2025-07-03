#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 14:48:56 2025

@author: hfausey
"""

import numpy as np
import pytest
from phozzy import d_L_calc
def test_dc():
    #Based on Ned Wright Javascript calculator
     #all for H00 = 67.36 and Omega_matter = 0.3166, and a flat universe
     
    #for z = 3, expect comoving distance of 6496.2
    
    expected = 6496.2 #Mpc
    output = d_L_calc.d_C(3)
    
    np.testing.assert_allclose(output, expected, 0.01)
    
    #for z = 6
    
    expected = 8408.0 #Mpc
    output = d_L_calc.d_C(6)
    
    np.testing.assert_allclose(output, expected, 0.01)
    
    
def test_dm():
    #For a flat Universe, this should be the same as dc
    
    #for z = 3
    
    expected = 6496.2 #Mpc
    output = d_L_calc.d_M(3)
    
    np.testing.assert_allclose(output, expected, 0.01)
    
    #for z = 6
    
    expected = 8408.0 #Mpc
    output = d_L_calc.d_M(6)
    
    np.testing.assert_allclose(output, expected, 0.01)
    
def test_dL():
    #Luminosity distance expected values based on Ned Wright Javascript calculator
     #all for H00 = 67.36 and Omega_matter = 0.3166, and a flat universe
     
     #for z = 3
     
     expected = 25984.8 #Mpc
     output = d_L_calc.d_L(3)
     
     np.testing.assert_allclose(output, expected, 0.01)
     
     #for z = 6
     
     expected = 58856.0 #Mpc
     output = d_L_calc.d_L(6)
     
     np.testing.assert_allclose(output, expected, 0.01)
    
    
    
    
    
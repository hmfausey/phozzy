#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 13:43:49 2025

@author: hfausey
"""

import numpy as np
from phozzy import extinction

def test_extinction_smc():
    
    lam_minus = np.array([0.45, 0.61, 0.80])
    expected = np.array([-2.61, -2.47, -2.12])
    
    lam = 1/lam_minus #(in microns)
    
    Ebv = 1
    
    #E_lam-v = Alam - Av
    
    xsi_array, Rv = extinction.smc(lam)
    
    Alam = extinction.A_lambda(xsi_array, Ebv, Rv)
    
    Av = Rv*Ebv
    
    E_lam_v = Alam - Av
    
    print(E_lam_v)
    
    np.testing.assert_allclose(E_lam_v/Ebv, expected, 0.05)
    
    
    
def test_extinction_lmc():
    
    lam_minus = np.array([0.45, 0.61, 0.80])
    expected = np.array([-2.92, -2.59, -2.22])
    
    lam = 1/lam_minus #(in microns)
    
    Ebv = 1
    
    #E_lam-v = Alam - Av
    
    xsi_array, Rv = extinction.lmc(lam)
    
    Alam = extinction.A_lambda(xsi_array, Ebv, Rv)
    
    Av = Rv*Ebv
    
    E_lam_v = Alam - Av
    
    print(E_lam_v)
    
    np.testing.assert_allclose(E_lam_v/Ebv, expected, 0.05)
 
    
def test_extinction_mw():
    
    lam_minus = np.array([0.45, 0.61, 0.80])
    expected = np.array([-2.76, -2.58, -2.23])
    
    lam = 1/lam_minus #(in microns)
    
    Ebv = 1
    
    #E_lam-v = Alam - Av
    
    xsi_array, Rv = extinction.milkyway(lam)
    
    Alam = extinction.A_lambda(xsi_array, Ebv, Rv)
    
    Av = Rv*Ebv
    
    E_lam_v = Alam - Av
    
    print(E_lam_v)
    
    np.testing.assert_allclose(E_lam_v/Ebv, expected, 0.05)
    
    
    
    
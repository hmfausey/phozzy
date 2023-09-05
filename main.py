# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 15:10:00 2022

@author: hmfausey

main is the front end of phozzy where one can define values for all input
values for the build and fit of the simulated data
"""

import os

import phozzy

import numpy as np

if __name__ == "__main__":
    ##Define overworld parameters
    filter_edges = np.array([0.5, 0.64, 0.87, 1.2, 1.7, 2.4])


    save_string = 'special_2_clever_extinction_expected_z_input'
    
    phozzy.phozzy(500, filter_edges, save_string, z_input = 'expected', Ebv_input= 'clever', Ebv_prior='clever', Ebv_fitting = True, bottom = 5, top = 25, acc = 0.1)   







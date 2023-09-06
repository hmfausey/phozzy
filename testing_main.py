# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 15:10:00 2022

@author: hmfausey

main is the front end of phozzy where one can define values for all input
values for the build and fit of the simulated data
"""

import os

#import phozzy

import numpy as np

if __name__ == "__main__":
    ##Define overworld parameters
    '''
    filter_edges = np.array([0.5, 0.64, 0.87, 1.2, 1.7, 2.4])


    save_string = 'special_2_clever_extinction_expected_z_input'
    
    phozzy.phozzy(500, filter_edges, save_string, z_input = 'expected', Ebv_input= 'clever', Ebv_prior='clever', Ebv_fitting = True, bottom = 5, top = 25, acc = 0.1)   
    '''
    import mcmc
    
    x = np.array([0.6, 0.95, 1.5])
    y = np.array([0, 15, 95])
    yerr = np.sqrt(y) + np.ones(len(y))
    initial_guess = np.array([90, 0.8, 4, 0.1])
    GRB_params = initial_guess
    save_string = str("Test")
    filter_edges = np.array([[0.5,0.7],[0.7,1.2],[1.2,1.8]])
    
    mcmc.mcmc(x, y, yerr, initial_guess, GRB_params, save_string, filter_edges, parallel = False, nwalkers = 10, burnin = 20, produc = 40)






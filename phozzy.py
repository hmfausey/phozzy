# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 20:52:50 2022

@author: hmfausey

phozzy is the front end of the program, where I define values for all input
values for the build and fit of the simulated data
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

import create_data_set as make_data

import mcmc

import filters

import analysis


##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def phozzy(num, filter_edges, save_string, bottom = 5, top = 25, acc = 0.1, z_input = 'uniform', z_prior = 'uniform', Ebv_input = 'clever', Ebv_prior = 'clever', Ebv_fitting = True, uncertainty = 0.05, sig_noise = 3, extinction_law = 'smc'):
    ##Function serves as the main function for running an entire simulation
     #with desired number of runs, filter edges, different inputs, uncertainty,
     #instrument noise, etc.
     #Inputs:
         #num -- int, number of runs
         #filter_edges -- numpy array, edges of all filters
         #save_string -- desired string for all input and output data
         #z_input -- string, desired input redshift distribution. Options are 
          #'uniform' or 'expected' (default 'uniform')
         #z_prior -- string, desired redshift prior. Options are 'uniform' or
          #'expected'. 'uniform' is the default and is highly suggested
         #Ebv_input -- string, desired input E_{b-v} distribution. Options are
          #'uniform', 'basic', and 'clever'. (default 'clever')
         #Ebv_prior -- string, desired E_{b-v} prior. Options are 'uniform', 
          #'basic', and 'clever'. (default 'clever')
         #Ebv_fitting -- boolean, determines whether E_{b-v} is a free 
          #parameter or not. (default True)
         #uncertainty -- float, statistical uncertainty (default 0.05)
         #sig_noise -- float, 1-sigma instrument noise (default 3 micro-Jansky)
         #extinction_law -- string, extinction law model to be used. Choices
          #are 'smc' (for small magellenic cloud), 'lmc' (for large magellenic
          #cloud) or 'mw' (for milky way). 'smc' default
    

    make_data.build_set(num, filter_edges, save_string, z_input = z_input, z_prior = z_prior, Ebv_input = Ebv_input, Ebv_prior = Ebv_prior, Ebv_fitting = Ebv_fitting, uncertainty = uncertainty, sig_noise = sig_noise)
    
    GRB_params = np.loadtxt("GRB_parameters_"+save_string+".txt", delimiter = ' ')
    filter_obs = np.loadtxt("Filter_observations_"+save_string+".txt", delimiter = ' ')
    uncertainties = np.loadtxt("Uncertainties_"+save_string+".txt", delimiter = ' ')
    initial_guesses = np.loadtxt("Initial_guesses_"+save_string+".txt", delimiter = ' ')
    
    filter_centers = filters.get_filter_centers(filter_edges)    

    for i in range(num):

        print("For run number", str(i))
        result_string = save_string+'_results_'+str(i)
        
        y = filter_obs[i]
        yerr = uncertainties[i]
        initial_guess = initial_guesses[i]
        GRB_param_arr = GRB_params[i]
        
        mcmc.mcmc(filter_centers, y, yerr, initial_guess, GRB_param_arr, result_string, filter_edges, z_prior = z_prior, Ebv_prior = Ebv_prior, Ebv_fitting = Ebv_fitting, extinction_law = extinction_law)
    
    analysis_string = save_string+"_results"
    analysis.density_and_stats(num, analysis_string, bottom, top, acc)

    return None

# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 20:52:50 2022

@author: hmfausey

phozzy runs the enitre simulation from start to finish
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import os

import numpy as np

import create_data_set as make_data

import mcmc

import filters

import analysis


##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def phozzy(num, filter_edges, save_string, extinction_law = 'smc', uncertainty = 0.05, sig_noise = 3, nwalkers = 50, burnin = 250, produc=500, z_input = 'uniform', z_prior = 'uniform', Ebv_input = 'evolving', Ebv_prior = 'evolving', Ebv_fitting = True, Ebv_upper_limit = 0, flux_input = 'kann', parallel = False, cpu_num = int(3/4*os.cpu_count()), highz_threshold = 5, acc = 0.1):
    ##The main function for running the entire simulation with desired number 
     #of runs, filter edges, parameter inputs and priors, uncertainty,
     #instrument noise, etc.
     #Inputs:
         #num -- int, number of runs
         #filter_edges -- numpy array, edges of all filters
         #save_string -- desired string for all input and output data
         #extinction_law -- string, extinction law model to be used. Choices
          #are 'smc' (for small magellenic cloud), 'lmc' (for large magellenic
          #cloud) or 'mw' (for milky way). 'smc' default
         #uncertainty -- float, statistical uncertainty (default 0.05,
          #corresponding to 5% statistical uncertainty)
         #sig_noise -- float, 1-sigma instrument noise (default 3 micro-Jansky)
         #nwalkers -- int, the number of walkers used by the MCMC fitting method
          #(default 50)
         #burnin -- int, number of steps in the burn-in phase (default 250)
         #produc -- int, number of steps in the production phase (default 500)
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
         #Ebv_upper_limit -- int, specifies whether the user would like upper 
          #limits applied to the evolving extinction prior and which one to 
          #use. 0 corresponds to no upper limit while 1 and 2 correspond to 
          #upper limit 1 (less constraining) and upper limit 2 (more 
          #constraining) from the paper (see for more details)
         #flux_input -- string, desired flux input distribution. Options are
          #'kann' (based on Kann et al. 2024), and 'basic' (default 'kann')
         #parallel -- bool, indicates whether the user would like to
          #parellalize the code. (default False)
         #cpu_num -- int, number of CPUs to be used when parellalizing the 
          #fit (default floor(3/4 * cpu_count))
         #highz_threshold -- int, cutoff between low and high redshift
         #acc -- float, desired accuracy (eg. 0.10 for GRBs within 10% of the 
          #input reshift, 0.20 for GRBs within 20% of the input redshift)
         
    
    #Create a set of GRB spectra and determine their associated photometric 
     #band measurements, including perturbations due to statistical uncertainty
     #and instrument noise, while recording initial guesses, input parameter
     #values for each GRB, and uncertainties for each photometric band
     #measurement. The files holding this information are saved to the current
     #working directory
    make_data.build_set(num, filter_edges, save_string, extinction_law = extinction_law, uncertainty = uncertainty, sig_noise = sig_noise, z_input = z_input, z_prior = z_prior, Ebv_input = Ebv_input, Ebv_prior = Ebv_prior, Ebv_fitting = Ebv_fitting, flux_input=flux_input, upper_limit=Ebv_upper_limit)
    
    
    #Load the data back in
    GRB_params = np.loadtxt("GRB_parameters_"+save_string+".txt", delimiter = ' ')
    filter_obs = np.loadtxt("Filter_observations_"+save_string+".txt", delimiter = ' ')
    uncertainties = np.loadtxt("Uncertainties_"+save_string+".txt", delimiter = ' ')
    initial_guesses = np.loadtxt("Initial_guesses_"+save_string+".txt", delimiter = ' ')
    
    #Filter centers will be the x values for in the MCMC fitting method these 
     #will not change during the simulation
    x = filters.get_filter_centers(filter_edges)    
    
    #For each simulated GRB
    for i in range(num):
        
        print("For run number", str(i))
        #Create string for saving the results. Files will be numbered from 0 to
         #num-1
        result_string = save_string+'_results_'+str(i)
        
        #y values and y errors will be the filter observations and filter
         #observation uncertainties respectively
        y = filter_obs[i]
        yerr = uncertainties[i]
        #Initial guess will be starting point for MCMC method
        if num>1:
            #y values and y errors will be the filter observations and filter
             #observation uncertainties respectively
            y = filter_obs[i]
            yerr = uncertainties[i]
            #Initial guess will be starting point for MCMC method
            initial_guess = initial_guesses[i]
            #Keep input parameters for this run to compare to results later
            GRB_param_arr = GRB_params[i]
        else:
            #Workaround for if only running for one GRB (indexing error)
            y = filter_obs
            yerr = uncertainties
            initial_guess = initial_guesses
            GRB_param_arr = GRB_params

        #Run the MCMC fitting method
        mcmc.mcmc(x, y, yerr, initial_guess, GRB_param_arr, result_string, filter_edges, nwalkers=nwalkers, burnin=burnin, produc=produc, extinction_law = extinction_law, z_prior = z_prior, Ebv_prior = Ebv_prior, Ebv_fitting = Ebv_fitting, upper_limit=Ebv_upper_limit, parallel = parallel, cpu_num = cpu_num)

    #Setup string for saving the analysis results
    analysis_string = save_string+"_results"
    #Run the analysis
    analysis.input_output_density_plot(num, nwalkers, filter_edges, analysis_string, highz_threshold, acc, Ebv_fitting=Ebv_fitting)

    return None

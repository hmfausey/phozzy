# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 12:36:10 2022

@author: hmfausey

Module includes all necessary function for Markov-chain Monte Carlo (MCMC) 
simulations
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

from matplotlib import pyplot as plt

import emcee

import corner

import multiprocessing as mp

import os

import build_spectrum

import filters

##############################################################################
###################################CONSTANTS##################################
##############################################################################

PI = np.pi

##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def mcmc(x, y, yerr, initial_guess, GRB_params, save_string, filter_edges, z_prior = 'uniform', Ebv_prior = 'clever', Ebv_fitting = True, extinction_law = 'smc', parallel = True, cpu_num = int(3/4*os.cpu_count())):
   ##Performs a Markov-Chain Monte-Carlo fitting method for a set of simulated
    #GRB photometric band measurements and uncertainties, records the final
    #parameter results, and saves them to a file
    #Inputs:
        #x -- numpy array, contains the central wavelength for each 
         #photometric band
        #y -- numpy array, constains the flux measurement for each photometric
         #band
        #yerr -- numpy array, constains the uncertainty on the flux
         #measurement for each photometric band
        #initial_guess -- numpy array, initial guess for each parameter.
        #GRB_params -- numpy array, the true values for each GRB parameter. 
         #Used for plotting and comparison purposes
        #save_string -- desired string for all input and output data
        #filter_edges -- 2D numpy array, contains the upper and lower edges of 
         #each filter in order 
        #z_prior -- str, desired redshift prior. Options are 'uniform' or
         #'expected'. 'uniform' is the default and is highly suggested
        #Ebv_prior -- string, desired E_{b-v} prior. Options are 'uniform', 
         #'basic', and 'clever'. (default 'clever')
        #Ebv_fitting -- boolean, determines whether E_{b-v} is a free 
         #parameter or not. (default True)
        #extinction_law -- string, extinction law model to be used. Choices
         #are 'smc' (for small magellenic cloud), 'lmc' (for large magellenic
         #cloud) or 'mw' (for milky way). 'smc' default
        #parallel -- bool, indicates whether the user would like to
         #parellalize the code. True default
        #cpu_num -- int, number of CPUs to be used when parellalizing the 
         #fit
     #Returns:
         #None
        
    #set number of walkers, and lengths of burn-in and production chains
    nwalkers = 50
    burnin = 250
    produc = 500
    
    #dimensions set to number of parameters
    ndim = len(initial_guess)
    #perturb initial position of each walker
    p0 = [np.array(initial_guess) + 1e-5 * np.random.randn(ndim) for i in range(nwalkers)]
    
    #Determine whether Ebv is a fitting parameter
    if Ebv_fitting:
        #If Ebv is a fitting parameter, include it for parameter names and
         #labels
        param_names = ["F0", "beta", "z", "Eb-v"]
        labels = [r'$F_0$', r'$\beta$', r'z', r'$E_{b-v}$']
        f_real, beta_real, z_real, Ebv_real = GRB_params
    else:
        #If Ebv is not a fitting parameter, exclude it from parameter names and
         #labels
        param_names = ["F0", "beta", "z"]
        labels = [r'$F_0$', r'$\beta$', r'z']
        f_real, beta_real, z_real = GRB_params
        Ebv_real = 0
        
    #Calculate normalization wavelength with filter edges
    
    #Using multiprocessing tool with 3/4 of all cpus in use (rounded down)
    if parallel:
        with mp.Pool(processes = cpu_num) as pool:
            #set up sampler
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (y, yerr, filter_edges, z_prior, Ebv_prior, Ebv_fitting, extinction_law), pool = pool)
            
            print("Running burn-in...")
            p0, _, _ = sampler.run_mcmc(p0, burnin)
            
            sampler.reset()
        
            print("Running production...")
            sampler.run_mcmc(p0, produc)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (y, yerr, filter_edges, z_prior, Ebv_prior, Ebv_fitting, extinction_law))
            
        print("Running burn-in...")
        p0, _, _ = sampler.run_mcmc(p0, burnin)
        
        sampler.reset()
    
        print("Running production...")
        sampler.run_mcmc(p0, produc)
        
    samples = sampler.chain[:,:,:].reshape((-1, ndim))
    
    #plot production chains for each parameter
    xaxis = np.arange(produc)
    for i in range(ndim):
        for j in range(nwalkers):
            plt.plot(xaxis, sampler.chain[j, :, i])
        plt.title("Paramter "+param_names[i])
        plt.savefig(save_string+"_MCMC_produc_"+param_names[i]+".png")
        plt.close()
    
    #plot photometric band measurements with error bars along with the original 
     #GRB spectrum
    plt.errorbar(x, y, yerr = yerr, fmt = "ko", capsize=0)
    lam_obs, spectrum = build_spectrum.build(filter_edges, f_real, beta_real, z_real, Ebv_real, extinction_law=extinction_law)
    plt.plot(lam_obs, spectrum, "k-")
    
    #Initialize 2-D array for storing parameter results and true parameter 
     #values
    datastore = np.zeros((nwalkers, len(GRB_params)*2))
    
    samples = sampler.flatchain
    
    #Plot the spectrum associated with the final positions of each of the 
     #walkers
    for i in range(nwalkers):
        s = samples[-(i+1)]
        #Differentiating whether E_bv was used as a fitting parameter
        if Ebv_fitting:
            Fnew, betanew, znew, E_bvnew = s
        else:
            Fnew, betanew, znew = s
            E_bvnew = 0

        #Finding the spectrum associated with the final positions of each 
         #walkers the the (errorless) photometric band measurements associated
         #with it
        lam_obs, spectrum = build_spectrum.build(filter_edges, Fnew, betanew, znew, E_bvnew, extinction_law=extinction_law)
        filter_vals, _, _ = filters.filter_observations(lam_obs, spectrum, filter_edges)
        
        #Save results
        if Ebv_fitting:
            final_data = np.array([f_real, Fnew, beta_real, betanew, z_real, znew, Ebv_real, E_bvnew])
        else:
            final_data = np.array([f_real, Fnew, beta_real, betanew, z_real, znew])
        
        datastore[i][:] = final_data
        
        #Plot spectrum associated with the final position of each walker
        plt.plot(lam_obs, spectrum, "c-", alpha = 0.3)
    
    #Save
    plt.savefig(save_string+"_walker_plot.png")
    plt.show()
    plt.close()
    
    #Save 2-D array of parameter results and original values to text file (to
     #be used for analysis)
    np.savetxt(save_string+"_MCMC_testing_datastore.txt", datastore, delimiter = ' ')
    
    #Create corner plot of posteriors of each parameter and save it
    figure = corner.corner(samples, labels=labels)
    figure.savefig(save_string+"_corner_plot.png")
    plt.show(figure)
    plt.close(figure)    

def model(params, filter_edges, Ebv_fitting, extinction_law):
    if Ebv_fitting:
        f_0, beta, z, E_bv = params
        lam_obs, spectrum = build_spectrum.build(filter_edges, f_0, beta, z, E_bv, extinction_law=extinction_law)
        filter_vals, _, _ = filters.filter_observations(lam_obs, spectrum, filter_edges)
    else:
        f_0, beta, z = params
        lam_obs, spectrum = build_spectrum.build(filter_edges, f_0, beta, z, 0)
        filter_vals, _, _ = filters.filter_observations(lam_obs, spectrum, filter_edges)
    return filter_vals

def chi_squared(params, y, yerr, filter_edges, Ebv_fitting, extinction_law):
    model_filters = model(params, filter_edges, Ebv_fitting, extinction_law)
    
    chi_2 = np.sum(((y - model_filters)/(yerr))**2)
    
    return chi_2

def lnlikelihood(params, y, yerr, filter_edges, Ebv_fitting, extinction_law):
    chi2 = chi_squared(params, y, yerr, filter_edges, Ebv_fitting, extinction_law)
    
    return -chi2/2

def lnprior(params, z_prior, Ebv_prior, Ebv_fitting):
    if Ebv_fitting:
        f_0, beta, z, E_bv = params
        #figure out extinction prior if applicable

        if E_bv >= 0:

            if Ebv_prior == 'basic':
                Ebv_norm = 4.28032
                Ebv_prior = Ebv_norm * np.exp(-Ebv_norm * E_bv)
            elif Ebv_prior == 'clever':
                if z < 2:
                    constant = 6.9
                    if E_bv > 2.05:
                        return -np.inf
                elif 2 <= z < 4:
                    constant = 12.6
                    if E_bv > 1.02:
                        return -np.inf
                elif z > 4:
                    constant = 36.2
                    if E_bv > 0.17:
                        return -np.inf
                    
                Ebv_prior = constant * np.exp(-constant * E_bv)

            else:
                raise Exception("Invalid E_{b-v} prior. Ebv_prior must be either 'basic', 'clever', or 'uniform' and must be a string.")
        else:
            return -np.inf
    else:
        f_0, beta, z = params
        Ebv_prior = 1.0
    
    if (f_0 > 0) and (0 < beta <= 2.5) and (0 < z < 25):
        #spectral index prior
        mu = 0.7
        sig = 0.2
        beta_prior = 1/(sig*np.sqrt(2*PI)) * np.exp(-0.5*((beta - mu)/sig)**2)
    
        full_prior = beta_prior * Ebv_prior
        return np.log(full_prior)
    
    return -np.inf

def lnprob(params, y, yerr, filter_edges, z_prior, Ebv_prior, Ebv_fitting, extinction_law):
    prior = lnprior(params, z_prior, Ebv_prior, Ebv_fitting)
    
    if np.isfinite(prior):
        prob = prior + lnlikelihood(params, y, yerr, filter_edges, Ebv_fitting, extinction_law)
        
        return prob
    
    return -np.inf
    
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 20:46:24 2022

@author: hmfausey

create_data_set creates a set of GRB params, filter values, error bars, and
initial guesses to check that everything is alright across the board
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

from matplotlib import pyplot as plt

import random

import build_spectrum

import filters

import d_L_calc as dcal

##############################################################################
###################################FUNCTIONS##################################
##############################################################################


def build_set(
    num,
    filter_edges,
    save_string,
    extinction_law="smc",
    uncertainty=0.05,
    sig_noise=3,
    z_input="uniform",
    z_prior="uniform",
    Ebv_input="evolving",
    Ebv_prior="evolving",
    Ebv_fitting=True,
    flux_input="kann",
    upper_limit=0,
):
    ##Saves a set of GRB parameters, the corresponding photometric band
    # measurements (with error), uncertainties, and initial guesses to be used
    # later for fitting and analysis. Photometric band measurements are
    # determined by generating a GRB spectrum based on the input parameters,
    # determining the photometric band fluxes based on the filter edges, and
    # the perturbing the bands according to the statistical uncertainty and
    # instrumental noise. Initial guesses are generates using the parameter
    # priors(expected distributions) to be used as a starting point for the
    # MCMC fitting method.
    # Inputs:
    # num -- int, desired sample size
    # filter_edges -- 2D numpy array, contains the upper and lower edges of
    # each filter in order
    # save_string -- desired string for all input and output data
    # extinction_law -- string, extinction law model to be used. Choices
    # are 'smc' (for small magellenic cloud), 'lmc' (for large magellenic
    # cloud) or 'mw' (for milky way). 'smc' default
    # uncertainty -- float, statistical uncertainty (default 0.05)
    # sig_noise -- float, 1-sigma instrument noise (default 3 micro-Jansky)
    # z_input -- string, desired input redshift distribution. Options are
    #'uniform' or 'expected' (default 'uniform')
    # z_prior -- string, desired redshift prior. Options are 'uniform' or
    #'expected'. 'uniform' is the default and is highly suggested
    # Ebv_input -- string, desired input E_{b-v} distribution. Options are
    #'uniform', 'basic', and 'evolving'. (default 'evolving')
    # Ebv_prior -- string, desired E_{b-v} prior. Options are 'uniform',
    #'basic', and 'evolving'. (default 'evolving')
    # Ebv_fitting -- boolean, determines whether E_{b-v} is a free
    # parameter or not. (default True)
    # flux_input -- string, desired flux input distribution. Options are
    #'kann' (based on Kann et al. 2024), and 'basic' (default 'kann')
    # upper_limit -- int, specifies whether the user would like upper limits
    # applied to the evolving extinction prior and which one to use. 0
    # corresponds to no upper limit while 1 and 2 correspond to upper
    # limit 1 (less constraining) and upper limit 2 (more constraining)
    # from the paper (see url for more details)
    # Returns:
    # None

    # Initialize filter observations and uncertainties
    
    filter_obs = np.zeros((num, len(filter_edges)))
    uncertainties = np.zeros((num, len(filter_edges)))

    if Ebv_fitting:
        # Determine if there are 3 or 4 free parameters
        initial_guesses = np.zeros((num, 4))
        GRB_params = np.zeros((num, 4))
    elif Ebv_fitting == False:
        initial_guesses = np.zeros((num, 3))
        GRB_params = np.zeros((num, 3))
    else:
        raise Exception("Ebv_fitting must be a boolean (True or False)")

    # For flux distribution
    m = 6.18
    s = 2.65

    i = 0

    while i < num:
        # Create GRB parameters and initial guesses for the desired number of
        # samples

        # Get values for beta
        betaval = random.gauss(0.7, 0.2)
        betaguess = random.gauss(0.7, 0.2)

        # Get values for z
        if z_input == "uniform":
            zval = random.uniform(0, 20)
        elif z_input == "expected":
            zval = random.lognormvariate(0.93, 0.7)
        else:
            raise Exception(
                "Invalid redshift input distribution. Choices are 'uniform' or 'expected' and must be a string"
            )

        zguess = random.uniform(0, 20)

        if flux_input == "kann":
            # Get F_0 value(depends on z)
            # Pull from approximate flux distribution at z=10 and adjust for actual
            # redshift using ratio of luminosity distances
            F10 = np.exp(random.gauss(m, s))
            Fval = F10 * (((dcal.d_L(10)) ** 2) / ((dcal.d_L(zval)) ** 2))
        elif flux_input == "basic":
            Fval = random.gauss(sig_noise * 300, np.sqrt(sig_noise * 100))
        else:
            raise Exception(
                "Invalid flux input distribution. Choices are 'kann' or 'basic'."
            )

        # Get value for E_{b-v} if applicable(depends on z)
        if Ebv_fitting:
            # Different E_{b-v} values depending on the input distribution
            if Ebv_input == "none":
                Ebvval = 0
            elif Ebv_input == "uniform":
                Ebvval = random.uniform(0, 3)
            elif Ebv_input == "basic":
                lambd = 4.28032
                Ebvval = random.expovariate(lambd)
            elif Ebv_input == "evolving":
                if zval < 2:
                    constant = 6.9
                    Ebvval = random.expovariate(constant)
                    # Inculuding considerations for upper limits
                    if Ebvval > 2.05 and (upper_limit == 1 or upper_limit == 2):
                        while Ebvval > 2.05:
                            Ebvval = random.expovariate(constant)
                elif 2 <= zval < 4:
                    constant = 12.6
                    Ebvval = random.expovariate(constant)
                    if Ebvval > 1.02 and (upper_limit == 1 or upper_limit == 2):
                        while Ebvval > 1.02:
                            Ebvval = random.expovariate(constant)
                elif zval > 4:
                    constant = 36.2
                    Ebvval = random.expovariate(constant)
                    if Ebvval > 0.17 and upper_limit == 2:
                        while Ebvval > 0.17:
                            Ebvval = random.expovariate(constant)
                    elif Ebvval > 0.34 and upper_limit == 1:
                        while Ebvval > 0.34:
                            Ebvval = random.expovariate(constant)
            else:
                raise Exception(
                    "Invalid E_{b-v} input distribution. Choices are 'uniform', 'basic', or 'evolving'  or 'none' and must be a string"
                )
            # Different E_{B-V} values depending on the prior
            if Ebv_prior == "uniform":
                Ebvguess = random.uniform(0, 3)
            elif Ebv_prior == "basic":
                lambd = 4.28032
                Ebvguess = random.expovariate(lambd)
            elif Ebv_prior == "evolving":
                if zval < 2:
                    constant = 6.9
                    Ebvguess = random.expovariate(constant)
                    # Including considerations for upper limits
                    if Ebvguess > 2.05 and (upper_limit == 1 or upper_limit == 2):
                        while Ebvguess > 2.05:
                            Ebvguess = random.expovariate(constant)
                elif 2 <= zval < 4:
                    constant = 12.6
                    Ebvguess = random.expovariate(constant)
                    if Ebvguess > 1.02 and (upper_limit == 1 or upper_limit == 2):
                        while Ebvguess > 1.02:
                            Ebvguess = random.expovariate(constant)
                elif zval > 4:
                    constant = 36.2
                    Ebvguess = random.expovariate(constant)
                    if Ebvguess > 0.17 and upper_limit == 2:
                        while Ebvguess > 0.17:
                            Ebvguess = random.expovariate(constant)
                    elif Ebvguess > 0.34 and upper_limit == 1:
                        while Ebvguess > 0.34:
                            Ebvguess = random.expovariate(constant)
            else:
                raise Exception(
                    "Invalid E_{b-v} prior distribution. Choices are 'uniform', 'basic', or 'evolving' and must be a string"
                )

            # Save GRB original parameters
            OG_params = np.array([Fval, betaval, zval, Ebvval])

            filter_centers = filters.get_filter_centers(filter_edges)
            # Build the spectrum based on the GRB parameters
            lam_obs, spectrum = build_spectrum.build(
                filter_edges, Fval, betaval, zval, Ebvval, extinction_law=extinction_law
            )
            
            # determine filter values (perturbed according to uncertainty) and
            # their corresponding error
            _, filter_vals, quadrature = filters.filter_observations(
                lam_obs,
                spectrum,
                filter_edges,
                noisy=True,
                uncertainty=uncertainty,
                sig_noise=sig_noise,
            )
            
            # If the reddest filter is above the 5-sigma detection limit, keep
            # the current filters and uncertainties
            if filter_vals[-1] >= (5 * sig_noise):
                
                plt.plot(lam_obs, spectrum, "k-")
                plt.errorbar(filter_centers, filter_vals, yerr=quadrature, fmt="co")
                plt.show()
                plt.close()
                
                filter_obs[i, :] = filter_vals
                uncertainties[i, :] = quadrature

                # Initial guess for flux based on the flux of the reddest band
                Fguess = filter_vals[-1]
                # create initial guess array
                guess_array = np.array([Fguess, betaguess, zguess, Ebvguess])
                # add initial guesses to 2D array of all initial guesses
                initial_guesses[i, :] = guess_array
                # add GRB parameters to 2D array of all GRB parameter values
                GRB_params[i, :] = OG_params

                i += 1

        else:
            # Same as above, but excluding E_{b-v} as a free parameter
            OG_params = np.array([Fval, betaval, zval])

            filter_centers = filters.get_filter_centers(filter_edges)
            lam_obs, spectrum = build_spectrum.build(
                filter_edges, Fval, betaval, zval, 0
            )
            _, filter_vals, quadrature = filters.filter_observations(
                lam_obs,
                spectrum,
                filter_edges,
                noisy=True,
                uncertainty=uncertainty,
                sig_noise=sig_noise,
            )

            if filter_vals[-1] >= 5 * sig_noise:
                
                plt.plot(lam_obs, spectrum, "k-")
                plt.errorbar(filter_centers, filter_vals, yerr=quadrature, fmt="co")
                plt.show()
                plt.close()
                
                filter_obs[i, :] = filter_vals
                uncertainties[i, :] = quadrature

                Fguess = filter_vals[-1]

                guess_array = np.array([Fguess, betaguess, zguess])

                initial_guesses[i, :] = guess_array

                GRB_params[i, :] = OG_params

                i += 1

    # Save GRB parameters, filter observations, uncertainties and initial
    # guesses to file that can be pulled from by other modules for fitting and
    # analysis, or directly compared to results by the user in whatever way they
    # see fit.
    np.savetxt("GRB_parameters_" + save_string + ".txt", GRB_params, delimiter=" ")
    np.savetxt("Filter_observations_" + save_string + ".txt", filter_obs, delimiter=" ")
    np.savetxt("Uncertainties_" + save_string + ".txt", uncertainties, delimiter=" ")
    np.savetxt("Initial_guesses_" + save_string + ".txt", initial_guesses, delimiter=" ")

    return None

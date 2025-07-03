# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 11:39:49 2022

@author: hmfausey

Function takes a set of filter values and spectrum and determines average flux
through each wavelength band
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

##############################################################################
###################################FUNCTIONS##################################
##############################################################################


def get_filter_centers(filter_edges):
    ##Determines the centers of each filter given the filter edges
    # Inputs:
    # filter_edges -- 2D numpy array, contains the upper and lower edges of
    # each filter in order
    # (eg. filter_edges = [[lower_1, upper_1], [lower_2, upper_2], [lower_3, upper_3]])
    # Returns:
    # filter_centers -- array, centers of each filter

    # Initialize array for filter centers
    filter_centers = np.zeros(len(filter_edges))

    # Calculate filter centers
    for i in range(len(filter_centers)):
        filter_centers[i] = (filter_edges[i][0] + filter_edges[i][-1]) / 2

    return filter_centers


def filter_observations(
    lam_obs, spectrum, filter_edges, noisy=False, uncertainty=0, sig_noise=3
):
    ##Finds the average observed flux in each filter using trapezoid rule for
    # numerical integration
    # Inputs:
    # lam_obs -- numpy array, observed wavelengths
    # spectrum -- numpy array, observed fluxes at each wavelength
    # filter_edges -- 2D numpy array, contains the upper and lower edges of
    # each filter in order
    # (eg. filter_edges = [[lower_1, upper_1], [lower_2, upper_2], [lower_3, upper_3]])
    # noisy -- boolean, variable indicating whether testing noisy or quiet
    # data (default False)
    # uncertainty -- scalar, value of the uncertainty from transmission of
    # lthe spectrum (default 0)
    # sig_noise -- scalar, 1-sig noise from instrument in micro-Jansky
    # Returns:
    # filter_aves -- array, fluxes through each filter assuming no noise
    # perturbed_vals -- array, perturbed fluxes in each filter using
    # statistical and instrumental uncertainties added in quadrature
    # quadrature -- array, statistical and instrumental uncertainties
    # added in quadrature

    # Initialize array to hold filter fluxes
    filter_vals = np.zeros(len(filter_edges))

    filter_num = 0
    ##Using trapazoid rule for integration
    # The for loop below looks at each wavelength from shortest to longest,
    # determines if it is in the current filter or split between two filters.
    # The flux contribution for the wavelength range being examined is then
    # added to the correct filter(s)
    for i in range(len(lam_obs) - 1):
        if filter_num == (len(filter_edges)):
            # If reached the long edge of the longest filter, stop
            break

        elif lam_obs[i] < filter_edges[filter_num][0] < lam_obs[i + 1]:
            # If the wavelength range is partially in the short edge of the
            # filter, half the flux contribution to that filter (this
            # condition will only occur when looking at the first filter)
            filter_vals[filter_num] += np.abs(
                0.5
                * (
                    0.5
                    * (spectrum[i] + spectrum[i + 1])
                    * (lam_obs[i + 1] - lam_obs[i])
                )
            )

        elif (filter_edges[filter_num][0] <= lam_obs[i]) and (
            lam_obs[i + 1] <= filter_edges[filter_num][-1]
        ):
            # If the wavelength range being integrated is fully in the current
            # filter, full flux contribution is added to that filter
            filter_vals[filter_num] += np.abs(
                0.5 * (spectrum[i] + spectrum[i + 1]) * (lam_obs[i + 1] - lam_obs[i])
            )

        elif lam_obs[i] < filter_edges[filter_num][-1] < lam_obs[i + 1]:
            # If the wavelength range contains the upper edge of the current
            # filter, only half the flux contribution to that filter
            filter_vals[filter_num] += np.abs(
                0.5
                * (
                    0.5
                    * (spectrum[i] + spectrum[i + 1])
                    * (lam_obs[i + 1] - lam_obs[i])
                )
            )

            if (filter_num) < (len(filter_edges) - 1):
                # If the current filter is not the last filter (or current
                # upper filter edge is less than the highest filter edge) add
                # the other half of the flux contribution to the next filter.
                # (This prevents trying to add flux to a filter that does not
                # exist if the code is currently on the last filter)

                filter_vals[filter_num + 1] += np.abs(
                    0.5
                    * (
                        0.5
                        * (spectrum[i] + spectrum[i + 1])
                        * (lam_obs[i + 1] - lam_obs[i])
                    )
                )

            # Since the current wavelength range is now moving into the next
            # filter up, switch to the next filter so that the process is
            # repeated
            filter_num += 1
        elif lam_obs[i] == filter_edges[filter_num + 1]:
            # If by some miracle the transition between to wavelengths for the
            # integration happens to exactly match the transition to the next
            # filter range and prevents it from falling into any of the other
            # cases, this will take care of it by adding the current flux
            # contribution to the next filter up, and then switch to the next
            # filter

            filter_vals[filter_num + 1] += np.abs(
                0.5 * (spectrum[i] + spectrum[i + 1]) * (lam_obs[i + 1] - lam_obs[i])
            )

            filter_num += 1

    # Initialize array for filter averages
    filter_aves = np.zeros(len(filter_vals))

    # Determine filter averages
    for i in range(len(filter_aves)):
        filter_aves[i] = filter_vals[i] / (filter_edges[i][-1] - filter_edges[i][0])

    # Initialize array for perturbed values
    perturbed_vals = np.zeros(len(filter_aves))

    # Determine total uncertainty by adding statisitical and instrumental
    # uncertainties in quadrature
    quadrature = np.sqrt((sig_noise) ** 2 + (uncertainty * filter_vals) ** 2)

    # If we want to simulate actual data, include noise by perturbing filter
    # values according to the total uncertainty determined by the quadrature
    if noisy == True:
        for i in range(len(filter_vals)):
            perturbed_vals[i] = np.random.normal(filter_aves[i], quadrature[i])

    return filter_aves, perturbed_vals, quadrature

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 09:41:16 2022

@author: hmfausey

The power_law module contains functions for a power law spectrum. All 
functions return F_nu, the spectral flux density in ergs/cm^2/s/Hz as a 
function of wavelength, for either emitted or observed ranges
"""

##############################################################################
###################################FUNCTIONS##################################
##############################################################################


def fv_emit(lam_emit, lam_0, f_0, beta):
    ##Returns F_nu as a function wavelength.
    # Inputs:
    # lam_emit -- array, emitted wavelengths
    # lam_0 -- scalar, central wavelength of reddest wavelength band to be
    # used as the normalization point for the flux
    # f_0 -- scalar, flux at lam_0
    # beta -- scalar, spectral index
    # Returns:
    # flux -- array, emitted flux at each corresponding wavelength in
    # lam_emit

    flux = f_0 * (lam_0 / lam_emit) ** (-beta)

    return flux


def fv_obs(lam_obs, lam_0, f_0, beta, z):
    ##Function takes an array of observed wavelengths, determines their
    # wavelengths when emitted, and determines F_nu as a function wavelength by
    # using the fv_emit function. This function does not account for diminished
    # flux due to the inverse square law or any redshift affects, but assumes
    # the function shape when observed will be the same as the function shape
    # when emitted, and relies on f_0 for correct normalization/flux when
    # observed.
    # Inputs:
    # lam_obs -- array, emitted wavelengths
    # lam_0 -- scalar, central wavelength of reddest wavelength band to be
    # used as the normalization point for the flux
    # f_0 -- scalar, flux at lam_0
    # beta -- scalar, spectral index
    # z -- scalar, redshift
    # Returns:
    # lam_emit -- array, set of wavelengths when they are emitted
    # flux -- array, emitted flux at each corresponding wavelength in
    # lam_obs

    lam_emit = lam_obs / (1 + z)  # determine corresponding emitted wavelengths
    lam_0_emit = lam_0 / (1 + z)

    flux = fv_emit(lam_emit, lam_0_emit, f_0, beta)

    return lam_emit, flux

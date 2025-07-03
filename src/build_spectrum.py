# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 10:04:30 2022

@author: hmfausey

build_spectrum code combines all models and assuptions pertaining to the GRB
spectrum and creates a spectrum across the relevant wavelengths based on a set 
of input parameters
"""

##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

import power_law as pl

import attenuation

import extinction

##############################################################################
###################################FUNCTIONS##################################
##############################################################################


def build(filter_edges, f_0, beta, z, E_bv, extinction_law="smc"):
    ##Function takes a set of wavelengths and input parameters, and builds the
    # spectrum according to attenuation and extinction models.
    # Inputs:
    # filter_edges -- 2D numpy array, contains the upper and lower edges of
    # each filter in order
    # f_0 -- scalar, flux at lam_0
    # beta -- scalar, spectral index
    # z -- scalar, redshift
    # E_bv -- scalar, extinction, E_{b-v} value
    # extinction_law -- string, indicates desired extinction law for
    # spectrum. Choices are 'smc' for small magellenic cloud, 'lmc' for
    # large magellenic cloud, and 'mw' for milky way. (default 'smc')
    # Returns:
    # lam_obs -- array, observed wavelengths
    # spectrum -- array, observed spectrum when accounting for attenuation
    # and extinction

    # Calculate wavelength range for spectrum
    lam_obs = np.linspace(filter_edges[0][0], filter_edges[-1][-1], 1000)

    # Calculate normalization wavelength with filter edges (set to center of
    # reddest wavelength band)
    lam_0 = (filter_edges[-1][0] + filter_edges[-1][-1]) / 2

    # Get base power law
    lam_emit, base_spectrum = pl.fv_obs(lam_obs, lam_0, f_0, beta, z)

    # Get extincted spectrum at host galaxy
    lam_emit, dinosaur = extinction.get_extincted_curve(
        lam_emit, base_spectrum, E_bv, extinction_law
    )

    # Get attenuated and redshifted spectrum
    lam_obs, spectrum = attenuation.attenuate(lam_emit, dinosaur, z)

    return lam_obs, spectrum

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 15:10:00 2022

@author: hmfausey

This is an example main for how to run the code
"""

import os

import phozzy

import numpy as np

##If you want to use parallelization, the main function needs to be written as
# it is below, or run from the command line.

if __name__ == "__main__":
    ##Define parameters for the simulation

    # Edges of photometric bands
    '''
    filter_edges = np.array(
        [[0.5, 0.64], [0.64, 0.87], [0.87, 1.2], [1.2, 1.7], [1.7, 2.4]]
    )
    '''
    
    filter_edges = np.array([[0.4 - 0.56], [0.56 - 0.81], [0.81 - 1.16], [1.16 - 1.68], [1.68 - 2.4]])
    
    # string that will proceed all output files for this run
    save_string = "lower5bands_EEEN"

    # Number of GRBs to be simulated and fit
    nGRBs = 500

    # If you want to run the code in parallel, set parallel=True in the phozzy
    # function call
    phozzy.phozzy(
        nGRBs,
        filter_edges,
        save_string,
        parallel=True,
        nwalkers=50,
        burnin=250,
        produc=500,
        z_input='expected',
        Ebv_input='evolving',
        Ebv_prior= 'evolving'
        )

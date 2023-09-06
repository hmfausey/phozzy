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
    filter_edges = np.array([[0.5, 0.64], [0.64, 0.87], [0.87, 1.2], [1.2, 1.7], [1.7, 2.4]])


    save_string = 'test_run'
    
    phozzy.phozzy(5, filter_edges, save_string, parallel = True)







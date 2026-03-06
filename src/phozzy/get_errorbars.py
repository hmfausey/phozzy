#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 15:48:18 2025

@author: hfausey
"""

import numpy as np
from matplotlib import pyplot as plt
import random

###############################################################################
def import_one_result(num, walkers, save_string, Ebv_fitting=True):
    ##Imports the results of all walkers from all runs and combines them into
    # arrays of input and output values
    # Inputs:
    # num -- int, number of runs
    # walkers -- int, number of walkers used in the MCMC fitting method
    # save_string -- string, desired string for all input and output data
    # Returns:
    # F_in -- numpy array, array of all input flux values
    # F_out --numpy array, array of all output flux values
    # beta_in -- numpy array, array of all input spectral index values
    # beta_out -- numpy array, array of all output spectral index values
    # z_in -- numpy array, array of all input redshift values
    # z_out -- numpy array, array of all output redshift values
    # Ebv_in -- numpy array, array of all E_{B-V} input values
    # Ebv_out -- numpy array, array of all E_{B-V} output values

    # Initialize parameter arrays with enough space for every walkers from every
    # run
    F_in = np.zeros(walkers)
    F_out = np.zeros(walkers)
    beta_in = np.zeros(walkers)
    beta_out = np.zeros(walkers)
    z_in = np.zeros(walkers)
    z_out = np.zeros(walkers)
    Ebv_in = np.zeros(walkers)
    Ebv_out = np.zeros(walkers)

    filestring = save_string + "_" + str(num) + "_datastore.txt"

    datastore = np.loadtxt(filestring, delimiter=" ")
    
    F_in, F_out, beta_in, beta_out, z_in, z_out, Ebv_in, Ebv_out = datastore[:,0], datastore[:,1], datastore[:,2], datastore[:,3], datastore[:,4], datastore[:,5], datastore[:,6], datastore[:,7]
    
    '''
    F_in[walkers * i : walkers * (i + 1)] = datastore[:, 0]
    F_out[walkers * i : walkers * (i + 1)] = datastore[:, 1]

    beta_in[walkers * i : walkers * (i + 1)] = datastore[:, 2]
    beta_out[walkers * i : walkers * (i + 1)] = datastore[:, 3]

    z_in[walkers * i : walkers * (i + 1)] = datastore[:, 4]
    z_out[walkers * i : walkers * (i + 1)] = datastore[:, 5]

    if Ebv_fitting:
        Ebv_in[walkers * i : walkers * (i + 1)] = datastore[:, 6]
        Ebv_out[walkers * i : walkers * (i + 1)] = datastore[:, 7]

    if Ebv_fitting:
        return F_in, F_out, beta_in, beta_out, z_in, z_out, 
    '''
    return F_in, F_out, beta_in, beta_out, z_in, z_out, Ebv_in, Ebv_out
###############################################################################


stringy_boi = 'ME_6bands_UEEN'
save_string = stringy_boi+'/'+stringy_boi+ "_results"

len_find_error = 250
times = 1000
z_bound = 5

print(stringy_boi, z_bound)

completeness, purity, accuracy10, accuracy20 = np.zeros(times), np.zeros(times), np.zeros(times), np.zeros(times)

for i in range(times):
    comp = 0
    pur = 0
    acc10 = 0
    acc20 = 0
    high_tot = 0
    low_tot = 0
    
    pulled_indices = np.ones(len_find_error)*(-1)
    
    ind = 0
    while ind < len_find_error:
        
        pulled_index = -1
        
        while pulled_index in pulled_indices:
            pulled_index = random.randint(0,499)
        
        pulled_indices[ind] = pulled_index
        
        _, _, _, _, z_in, z_out, _, _ = import_one_result(pulled_index, 50, save_string)
        
        for j in range(len(z_in)):
            if z_in[j] >= z_bound:
                if z_out[j] >= z_bound:
                    comp += 1
                if z_in[j]*0.8 <= z_out[j] <= z_in[j]*1.2:
                    acc20 += 1
                    if z_in[j]*0.9 <= z_out[j] <= z_in[j]*1.1:
                        acc10 += 1
                high_tot += 1
                
            elif z_in[j] < z_bound:
                if z_out[j] < z_bound:
                    pur += 1
                low_tot += 1
                
            

        ind += 1
    
    completeness[i], purity[i], accuracy10[i], accuracy20[i] = (comp/high_tot), (pur/low_tot), (acc10/high_tot), (acc20/high_tot)


print("Completeness: ", np.quantile(completeness, [0.16, 0.5, 0.84]), np.quantile(completeness, [0.5]) - np.quantile(completeness, [0.16]), np.quantile(completeness, [0.84]) - np.quantile(completeness, [0.5]))
print("Purity: ", np.quantile(purity, [0.16, 0.5, 0.84]), np.quantile(purity, [0.5]) - np.quantile(purity, [0.16]), np.quantile(purity, [0.84]) - np.quantile(purity, [0.5]))
print("10% Accuracy: ", np.quantile(accuracy10, [0.16, 0.5, 0.84]), np.quantile(accuracy10, [0.5]) - np.quantile(accuracy10, [0.16]), np.quantile(accuracy10, [0.84]) - np.quantile(accuracy10, [0.5]))
print("20% Accuracy: ", np.quantile(accuracy20, [0.16, 0.5, 0.84]), np.quantile(accuracy20, [0.5]) - np.quantile(accuracy20, [0.16]), np.quantile(accuracy20, [0.84]) - np.quantile(accuracy20, [0.5]))

plt.hist(completeness)
plt.show()
plt.close()

plt.hist(purity)
plt.show()
plt.close()

plt.hist(accuracy10)
plt.show()
plt.close()

plt.hist(accuracy20)
plt.show()
plt.close()


################################################################################

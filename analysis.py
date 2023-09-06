# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 23:56:55 2022

@author: hmfausey

Analyze creates a density plot of input vs output redshifts
"""
##############################################################################
####################################IMPORTS###################################
##############################################################################

import numpy as np

from matplotlib import pyplot as plt

##############################################################################
###################################CONSTANTS##################################
##############################################################################

LAM_ALPHA = 0.121567 #in MICROmeters

##############################################################################
###################################FUNCTIONS##################################
##############################################################################

def import_results(num, walkers, save_string): 
    ##Imports the results of all walkers from all runs and combines them into
     #arrays of input and output values
     #Inputs:
         #num -- int, number of runs
         #walkers -- int, number of walkers used in the MCMC fitting method
         #save_string -- string, desired string for all input and output data
     #Returns:
         #F_in -- numpy array, array of all input flux values
         #F_out --numpy array, array of all output flux values 
         #beta_in -- numpy array, array of all input spectral index values
         #beta_out -- numpy array, array of all output spectral index values
         #z_in -- numpy array, array of all input redshift values
         #z_out -- numpy array, array of all output redshift values
         #Ebv_in -- numpy array, array of all E_{B-V} input values
         #Ebv_out -- numpy array, array of all E_{B-V} output values
    
    #Initialize parameter arrays with enough space for every walkers from every
     #run
    F_in = np.zeros(walkers*num)
    F_out = np.zeros(walkers*num)
    beta_in = np.zeros(walkers*num)
    beta_out = np.zeros(walkers*num)
    z_in = np.zeros(walkers*num)
    z_out = np.zeros(walkers*num)
    Ebv_in = np.zeros(walkers*num)
    Ebv_out = np.zeros(walkers*num)
    
    #Loop through all output results files
    for i in range(num):
        #Using savestring to identify the correct set of runs, and i to pull
         #data from each individual run and add it to each corresponding 
         #parameter array
        filestring = save_string+"_"+str(i)+"_datastore.txt"
        
        datastore = np.loadtxt(filestring, delimiter = ' ')
        
        F_in[walkers*i: walkers*(i+1)] = datastore[:,0]
        F_out[walkers*i: walkers*(i+1)] = datastore[:,1]
        
        beta_in[walkers*i: walkers*(i+1)] = datastore[:,2]
        beta_out[walkers*i: walkers*(i+1)] = datastore[:,3]
        
        z_in[walkers*i: walkers*(i+1)] = datastore[:,4]
        z_out[walkers*i: walkers*(i+1)] = datastore[:,5]
        
        Ebv_in[walkers*i: walkers*(i+1)] = datastore[:,6]
        Ebv_out[walkers*i: walkers*(i+1)] = datastore[:,7]
    
    return F_in, F_out, beta_in, beta_out, z_in, z_out, Ebv_in, Ebv_out

def input_output_density_plot(num, walkers, filter_edges, save_string, highz_cutoff, acc):
    ##Method creates an input vs output redshift density plot to visualize 
     #results
     #Inputs:
         #num -- int, number of runs
         #walkers -- int, number of walkers used by the mcmc method
         #filter_edges -- 2D numpy array, contains the upper and lower edges of 
          #each filter in order 
          #(eg. filter_edges = [[lower_1, upper_1], [lower_2, upper_2], [lower_3, upper_3]])
         #save_string -- string, desired string for all input and output data
         #highz_cutoff -- int, cutoff between low and high redshift
         #acc -- float, desired accuracy (eg. 0.10 for GRBs within 10% of the 
          #input reshift, 0.20 for GRBs within 20% of the input redshift)
     #Returns:
         #None
    
    #Import relevant results (for this plot we are only interested in the input
     #and output redshift)
    _, _, _, _, z_in, z_out, _, _ = import_results(num, walkers, save_string)
    
    #Determine highest redshift input or output to determine redshift range
    
    zin_max = np.amax(z_in)
    zout_max = np.amax(z_out)
    
    z_range = int(max(zin_max, zout_max))+1
    
    #flatten array of filter edges
    filters_flattened = np.reshape(filter_edges,(filter_edges.size))
    #Remove repeat
    unique_filter_edges = np.unique(filters_flattened)
    #Determine redshift associated with band edges
    associated_redshift = (unique_filter_edges/LAM_ALPHA) - 1
    
    #Create 2D histogram with z_in vs z_out
    plt.hist2d(z_in, z_out,  bins=(50,50) , range=[[0,z_range],[0,z_range]], cmap = plt.cm.jet)
    plt.colorbar()
    plt.xlabel("Input redshift")
    plt.ylabel("Output redshift")
    
    #Create extra lines to make plot easier to read
    
    #Show where input redshift matches output redshift
    plt.plot([0,z_range+1],[0,z_range+1], 'w-.')
    #Show accuracy ranges
    plt.plot([0,z_range+1],[0, z_range+1*(1-acc)], 'w:', [0,z_range+1],[0,z_range+1*(1 + acc)], 'w:')
    #Show band edges
    for i in range(len(associated_redshift)):
        plt.plot([associated_redshift[i], associated_redshift[i]], [0,z_range+1], 'w--')
    #Show high vs low redshift cutoff
    plt.plot([0,z_range+1], [highz_cutoff ,highz_cutoff], 'w-', [highz_cutoff, highz_cutoff],[0,z_range+1],'w-')

    plt.savefig(save_string+'_density_plot.pdf')
    plt.show()
    plt.close()

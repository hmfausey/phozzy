# phozzy
'phozzy' is a package designed for simulation and fitting of photometric band measurements to estimate an instruments expected performance for accurately retrieving the redshift of a GRB. It allows the user to input any number of bands, a desired uncertainty and instrument noise, extinction law, parameter input and prior distributions, high-redshift cutoffs, and desired redshift accuracies, so it can be applied to a range of instruments and scenarios. Individual modules can be used for implementing intergalactic attenuation and extinction, creating spectra, and determining photometric band measurements. The key functions of each module and their uses are explained below.

If you make use of this package please cite Fausey et al., 2023 (https://academic.oup.com/mnras/article/526/3/4599/7289245#421711247)

## power_law
The power_law module contains two functions. Both return a power law based on a normalization and spectral index, but one takes emitted wavelengths to create the spectrum, and the other takes observed wavelengths.

## extinction
The extinction module contains multiple functions to determine the transmission function due to host galaxy extinction, and apply it to a spectrum.
The get_extincted_curve function applies the transmission function directly to a spectrum and returns it to the user. To just get the transmission function, one can use the transmission function instead. The other functions in the module are designed to determine parameter values for the selected extinction law. Extinction models are based on work done by Pei et al., 1992 (https://ui.adsabs.harvard.edu/link_gateway/1992ApJ...395..130P/doi:10.1086/171637)

## attenuation
The attenuation module employs the model created by Meiksin, 2006 (https://ui.adsabs.harvard.edu/link_gateway/2006MNRAS.365..807M/doi:10.1111/j.1365-2966.2005.09756.x), and strongly based on Madau, 1996 (https://ui.adsabs.harvard.edu/link_gateway/1995ApJ...441...18M/doi:10.1086/175332). It determines the transmission function due to absorption and scattering by hydrogen in the intergalactic medium based on the input redshift.
The attenuate function can be used to apply the transmission function  directly to a spectrum and return it to the user. To just get the transmission function, one can use teff_total instead.

## build_spectrum
The build_spectrum module pulls together the power_law, extinction, and attenuation modules to create a spectrum that accounts for both effects.
The build function allows the user to input a desired wavelegth range and normalization wavelength, as well as parameter inputs for the normalization, spectral index, wavelength, and extinction, to create a spectrum and return it to the user.

## filters
The filters module determines the flux measurement for a set of photometric bands given by the user.
In particular, the filter_observations function take a spectrum, the set of filter edges (in wavelength), whether noise should be included, and the statistical uncertainty and instrumental noise values desired by the user.

## create_data_set
The create_data_set module creates a sample of GRBs by randomly generating GRB parameters and initial guesses, and determining the corresponding photometric band fluxes and uncertianties. It then saves these values to text files that can be pulled from by the fitting and analysis modules to compare with the fitting results to esimate performance.

## d_L_calc
the d_L_calc module estimates the luminosity distance of a GRB at redshift z assuming a flat Universe. It is used by the create_data_set module to estimate the flux of a GRB seen from different redshifts. See https://ui.adsabs.harvard.edu/abs/1999astro.ph..5116H/abstract for overview.

## mcmc
The mcmc module uses the emcee package (Foreman-Mackey et al., 2013; https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F) to fit a set of photometric band measurements. It uses a Bayesian log-likelihood function and has multiple sets of priors, including uniform, basic, and evolving priors for the extinction, E_{B-V}, as well as the option for upper limits for the evolving extinction prior.

## analysis
The analysis module reads in the results of the mcmc fitting method and creates a 2D histogram of the input vs output redshift for all walkers of all runs of a simulation.

## phozzy
The phozzy module runs the entire simulation from start to finish using the create_data, mcmc, and analysis modules to create a set of simulated photometric band measurements, fit them, and interpret the results.

## example_main
phozzy come with an example main to show how to set up and run the code

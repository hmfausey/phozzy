# phozzy
phozzy is a package designed for simulation and fitting of photometric band measurements to estimate an instruments expected performance for extimating the redshift of a detected GRB. It allows the user to input some number of bands, a desired uncertainty and instrument noise, disired extinction law, parameter input and prior distributions, high-redshift cutoffs, and desired redshift accuracies, so it can be applied to a range of instruments and scenarios. Individual modules can be used for implementing intergalactic attenuation and extinction, creating spectra, and determining photometric band measurements. The key functions of each module and their uses are explained below.

If you make use of this package please cite (INSERT PAPER HERE)

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



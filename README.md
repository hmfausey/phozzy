# phozzy
phozzy is a package designed for simulation and fitting of photometric band measurements to estimate an instruments expected performance for extimating the redshift of a detected GRB. It allows the user to input some number of bands, a desired uncertainty and instrument noise, disired extinction law, parameter input and prior distributions, high-redshift cutoffs, and desired redshift accuracies, so it can be applied to a range of instruments and scenarios. Individual modules can be used for implementing intergalactic attenuation and extinction, creating spectra, and determining photometric band measurements. The key functions of each module and their uses are explained below.

## power_law
The power_law module contains two functions. Both return a power law based on a normalization and spectral index, but one takes emitted wavelengths to create the spectrum, and the other takes observed wavelengths.

## extinction
The extinction module contains multiple functions to determine the transmission function due to host galaxy extinction, and apply it to a spectrum.
The get_extincted_curve function applies the transmission function directly to a spectrum and returns it to the user. To just get the transmission function, one can use the transmission function instead. The other functions in the module are designed to determine parameter values for the selected extinction law.

## attenuation

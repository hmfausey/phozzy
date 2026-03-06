#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 14:29:29 2025

@author: hmfausey
"""

import numpy as np
from matplotlib import pyplot as plt

from build_spectrum import build
from build_spectrum_special import build as build_special

filter_edges = np.array([[0.4, 0.56], [0.56, 0.81], [0.81, 1.16], [1.16, 1.68], [1.68, 2.4]])


lam_obs_meiksin, meiksin = build(filter_edges, 100, 0.7, 8, 0)
lam_obs_totani20, totani20 = build_special(filter_edges, 100, 0.7, 8, 0, 20)

lam_obs_totani, totani = build_special(filter_edges, 100, 0.7, 8, 0, 22)


plt.plot(lam_obs_meiksin, meiksin, 'r-', label= "Madau/Meiksin model")
plt.plot(lam_obs_totani20, totani20, 'g--', label = r"Miralda-Escude/Totani model, $\log(\frac{N_{\rm H}}{\text{cm}^{-2}}) = 20$")
plt.plot(lam_obs_totani, totani, 'b:', label = r"Miralda-Escude/Totani model, $\log(\frac{N_{\rm H}}{\text{cm}^{-2}}) = 22$")
plt.legend()
plt.xlabel(r"Wavelength ($\mu$m)")
plt.ylabel(r"Flux ($\mu$Jy)")
plt.xlim([0.75, 1.5])
plt.savefig("model_compare.pdf")
plt.show()
plt.close()



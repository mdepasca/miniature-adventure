from matplotlib import pyplot as plt
import numpy as np
import supernova
import supernova_fit
import utilities as util

def plot_k(snObj, save=False):
    lambda_obs = [479.66, 638.26, 776.90, 910.82]

    z = (snObj.zSpec if snObj.zSpec else snObj.snPhotHost)
    lambda_rest = [el/(1+z) for el in lambda_obs]




"""
either get the observation file or get a supernova object...

- get a supernova object (independt on storing format!)
- needs:
    .filter's central wavelength
    .to check for observations in all bands at the same *integer* mjd
- blueshifts flux measures to restfram wavelength (either max flux or max SNR)
- plot lines with dots
- plot observer rest frame central wavelength of each filter
- plot shift direction arrow
- title containing source redshift
- axes labels
"""

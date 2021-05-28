#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy import constants as const
from astropy.io import ascii
from astropy.table import Table
from matplotlib import pyplot as plt
#from pathlib import Path

# Path for Throughput files
#path = '/Users/nflagey/WORK/ETC_SVN/etc/cgi-bin/mse/THROUGHPUT/'
#throughput_folder = Path("THROUGHPUT/")
path = 'THROUGHPUT/'

def mse_etc_wavegrid():
    """Makes the wavelength grid for each spectrograph.
    """

    # Read CSV files (digitize curves from CoDR report)
    # wav in micron
    lrblue = Table.read(path + 'LMR_resol_LRblue.csv', format='csv')
    lrgreen = Table.read(path + 'LMR_resol_LRgreen.csv', format='csv')
    lrred = Table.read(path + 'LMR_resol_LRred.csv', format='csv')
    lrj = Table.read(path + 'LMR_resol_LRj.csv', format='csv')
    mrblue = Table.read(path + 'LMR_resol_MRblue.csv', format='csv')
    mrgreen = Table.read(path + 'LMR_resol_MRgreen.csv', format='csv')
    mrred = Table.read(path + 'LMR_resol_MRred.csv', format='csv')
    mrh = Table.read(path + 'LMR_resol_MRh.csv', format='csv')

    # Resolution element as function of wavelength assuming spectral coverage from CoDR report
    # resel = fwhm = lambda / R
    resel_lrblue = lrblue['wav'] / lrblue['resol']
    resel_lrgreen = lrgreen['wav'] / lrgreen['resol']
    resel_lrred = lrred['wav'] / lrred['resol']
    resel_lrj = lrj['wav'] / lrj['resol']
    resel_mrblue = mrblue['wav'] / mrblue['resol']
    resel_mrgreen = mrgreen['wav'] / mrgreen['resol']
    resel_mrred = mrred['wav'] / mrred['resol']
    resel_mrh = mrh['wav'] / mrh['resol']

    # Now we can built the wavelength grid
    # LR - blue
    wgrid_lrblue = np.array([np.min(lrblue['wav'])])
    while True:
        wgrid_lrblue = np.append(wgrid_lrblue, wgrid_lrblue[-1] * (
                    1 + 1 / np.interp(wgrid_lrblue[-1], lrblue['wav'], lrblue['resol'])))
        if wgrid_lrblue[-1] > np.max(lrblue['wav']):
            break
    # LR - green
    wgrid_lrgreen = np.array([np.min(lrgreen['wav'])])
    while True:
        wgrid_lrgreen = np.append(wgrid_lrgreen, wgrid_lrgreen[-1] * (
                    1 + 1 / np.interp(wgrid_lrgreen[-1], lrgreen['wav'], lrgreen['resol'])))
        if wgrid_lrgreen[-1] > np.max(lrgreen['wav']):
            break
    # LR - red
    wgrid_lrred = np.array([np.min(lrred['wav'])])
    while True:
        wgrid_lrred = np.append(wgrid_lrred,
                                wgrid_lrred[-1] * (1 + 1 / np.interp(wgrid_lrred[-1], lrred['wav'], lrred['resol'])))
        if wgrid_lrred[-1] > np.max(lrred['wav']):
            break
    # LR - J
    wgrid_lrj = np.array([np.min(lrj['wav'])])
    while True:
        wgrid_lrj = np.append(wgrid_lrj,
                                wgrid_lrj[-1] * (1 + 1 / np.interp(wgrid_lrj[-1], lrj['wav'], lrj['resol'])))
        if wgrid_lrj[-1] > np.max(lrj['wav']):
            break
    # MR - blue
    wgrid_mrblue = np.array([np.min(mrblue['wav'])])
    while True:
        wgrid_mrblue = np.append(wgrid_mrblue, wgrid_mrblue[-1] * (
                1 + 1 / np.interp(wgrid_mrblue[-1], mrblue['wav'], mrblue['resol'])))
        if wgrid_mrblue[-1] > np.max(mrblue['wav']):
            break
    # MR - green
    wgrid_mrgreen = np.array([np.min(mrgreen['wav'])])
    while True:
        wgrid_mrgreen = np.append(wgrid_mrgreen, wgrid_mrgreen[-1] * (
                1 + 1 / np.interp(wgrid_mrgreen[-1], mrgreen['wav'], mrgreen['resol'])))
        if wgrid_mrgreen[-1] > np.max(mrgreen['wav']):
            break
    # MR - red
    wgrid_mrred = np.array([np.min(mrred['wav'])])
    while True:
        wgrid_mrred = np.append(wgrid_mrred,
                                wgrid_mrred[-1] * (1 + 1 / np.interp(wgrid_mrred[-1], mrred['wav'], mrred['resol'])))
        if wgrid_mrred[-1] > np.max(mrred['wav']):
            break
    # MR - H
    wgrid_mrh = np.array([np.min(mrh['wav'])])
    while True:
        wgrid_mrh = np.append(wgrid_mrh,
                              wgrid_mrh[-1] * (1 + 1 / np.interp(wgrid_mrh[-1], mrh['wav'], mrh['resol'])))
        if wgrid_mrh[-1] > np.max(mrh['wav']):
            break

    # Convert wavelength into Angstroms
    wgrid_lrblue *= 1e4
    wgrid_lrgreen *= 1e4
    wgrid_lrred *= 1e4
    wgrid_lrj *= 1e4
    wgrid_mrblue *= 1e4
    wgrid_mrgreen *= 1e4
    wgrid_mrred *= 1e4
    wgrid_mrh *= 1e4

    print('hello')

mse_etc_wavegrid()

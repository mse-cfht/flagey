#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

# Imports
import numpy as np
from astropy.io import fits
import pickle
import os
#from pathlib import Path

# Path for Throughput files
#path = '/Users/nflagey/WORK/ETC_SVN/etc/cgi-bin/mse/THROUGHPUT/'
#path = os.path.abspath(os.getcwd())
#throughput_folder = Path("THROUGHPUT/")
path = 'THROUGHPUT/'

# Wavelengths in Angstroms (3500 - 18500)
lam = np.arange(1500, dtype=float) * 10 + 3500.
lam0 = np.array([360, 370, 400, 445, 551, 658, 806, 1000, 1214, 1477, 1779]) * 10.
lam1 = np.array([360, 370, 400, 482, 626, 767, 900, 910, 950, 962, 1235, 1300, 1500, 1662, 1800]) * 10.

def mse_etc_injeff():
    """Writes the injection efficiency files for the ITC
    """

    # Open FITS file from IDL
    # JS: This file is not currently available
    hdu = fits.open('/Users/nflagey/WORK/MSE/MSE_Injection/PSFs/segments/6_inj_curves/inj_curves_mapmean_aao.fits')

    # IQ values (FWHM, in arcsecs, at 550 nm)
    iq_vals = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
              0.85, 0.90, 0.95, 1.00]

    # ZD values
    zd_vals = [0, 30, 50, 60]

    # Fiber size values
    lr_vals = [0.9, 1.0]
    mr_vals = [0.9, 1.0]
    hr_vals = [0.7, 0.8]

    # Create Python variables
    ie_lr_mean = hdu[0].data
    ie_mr_mean = hdu[1].data
    ie_hr_mean = hdu[2].data
    # ie_lr_sdev = hdu[3].data
    # ie_mr_sdev = hdu[4].data
    # ie_hr_sdev = hdu[5].data
    # ie_mr_mm = hdu[7].data
    # ie_hr_mm = hdu[8].data
    # ie_lr_sdevmean = hdu[9].data
    # ie_mr_sdevmean = hdu[10].data
    # ie_hr_sdevmean = hdu[11].data
    # ie_lr_sdevsdev = hdu[12].data
    # ie_mr_sdevsdev = hdu[13].data
    # ie_hr_sdevsdev = hdu[14].data
    # ie_lr_sdevmm = hdu[15].data
    # ie_mr_sdevmm = hdu[16].data
    # ie_hr_sdevmm = hdu[17].data

    # Structure the curves into dictionaries of dictionaries ...
    # LR
    ie_lr = {}
    for q in range(len(iq_vals)):
        iq_str = "iq" + "{:.2f}".format(iq_vals[q])
        ie_lr[iq_str] = {}
        for z in range(len(zd_vals)):
            zd_str = "zd" + "{:0^2}".format(zd_vals[z])
            ie_lr[iq_str][zd_str] = {}
            for f in range(len(lr_vals)):
                lr_str = "fib" + "{:.2}".format(lr_vals[f])
                ie_lr[iq_str][zd_str][lr_str] = ie_lr_mean[f, q, z, :]
    # MR
    ie_mr = {}
    for q in range(len(iq_vals)):
        iq_str = "iq" + "{:.2f}".format(iq_vals[q])
        ie_mr[iq_str] = {}
        for z in range(len(zd_vals)):
            zd_str = "zd" + "{:0^2}".format(zd_vals[z])
            ie_mr[iq_str][zd_str] = {}
            for f in range(len(mr_vals)):
                mr_str = "fib" + "{:.2}".format(mr_vals[f])
                ie_mr[iq_str][zd_str][mr_str] = ie_mr_mean[f, q, z, :]
    # HR
    ie_hr = {}
    for q in range(len(iq_vals)):
        iq_str = "iq" + "{:.2f}".format(iq_vals[q])
        ie_hr[iq_str] = {}
        for z in range(len(zd_vals)):
            zd_str = "zd" + "{:0^2}".format(zd_vals[z])
            ie_hr[iq_str][zd_str] = {}
            for f in range(len(hr_vals)):
                hr_str = "fib" + "{:.2}".format(hr_vals[f])
                ie_hr[iq_str][zd_str][hr_str] = ie_hr_mean[f, q, z, :]

    np.save(path + 'injeff_curve_lr.npy', ie_lr)
    np.save(path + 'injeff_curve_mr.npy', ie_mr)
    np.save(path + 'injeff_curve_hr.npy', ie_hr)


mse_etc_injeff()

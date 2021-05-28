#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Flagey"
__email__ = "flagey@cfht.hawaii.edu"
__status__ = "Development"

#2021 Modifications
# J. Sobeck; sobeck@cfht.hawaii.edu

# Imports
from mse_etc import MseSpectrum

# Create spectrum
s0 = MseSpectrum(tgtmag=24, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.5,
                 coating='ZeCoat', fibdiam=1.0, spectro='LR', spatbin=1, specbin=1, meth='getSNR', snr=1, etime=3600)

# s0 = MseSpectrum(tgtmag=20, band='g', template='flat', redshift=0, airmass=1.2, skymag=20.7, seeing=0.5,
#                  coating='ZeCoat', fibdiam=0.8, spectro='HR', spatbin=1, specbin=1, meth='getTime', snr=25, etime=3600)

# Check intermediate steps
# Apply extinction
s1 = s0.apply_atmos_ext()
# Apply frontend throughput
s2 = s1.apply_throughput_front()
# Apply injection efficiency
s3 = s2.apply_injeff()
# Apply backend throughput
s4 = s3.apply_throughput_back()

# Compute SNR at once
s4b, script, div = s0.compute_snr(doplot='offline')

print('Done')

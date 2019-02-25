# Import
import vegapy
import os
import astropy.units as u
import numpy as np

# Settings
verbose = 0
visual = 0

# Initialize
PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, '../source/')
print('>>> TEST: {}'.format(FILE))

# Initilization tests
det = vegapy.Detector((64, 64), pixel_size=0.01*u.arcsec, randomkeywortwithoutmeaning=3)

det_RON = vegapy.Detector((64, 64), pixel_size=0.01*u.arcsec, readout_noise=35*u.electron/u.pix)

exposure_array = np.ones((64, 64)) * u.ph / u.s
RO = det_RON(exposure_array, integration_time=200*u.ms, target_FoV=det_RON.FoV)
if verbose:
	print(RO)

try:
	det.__dict__
	print(">>> TEST DETECTOR: Successful")
except:
	print(">>> TEST DETECTOR: NOT successful")

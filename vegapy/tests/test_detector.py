# Import
import vegapy
import os
import astropy.units as u

# Initialize
PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, '../source/')
print('>>> TEST: {}'.format(FILE))

# Initilization tests
det = vegapy.Detector((64, 64), pixel_size=0.01*u.arcsec, randomkeywortwithoutmeaning=3)

try:
	det.__dict__
	print(">>> TEST DETECTOR: Successful")
except:
	print(">>> TEST DETECTOR: NOT successful")

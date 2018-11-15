# Import
import vegapy
import os
import matplotlib.pyplot as plt
import astropy.units as u

# Settings
visual = 0
verbose = 0

# Initialize
PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, '../source/')
print('>>> TEST: {}'.format(FILE))


# Initilization tests
print(">>> TEST TARGET> tar1")
tar = vegapy.Target(band='H', pixel_scale=0.1*u.arcsec, shape=(64, 64), sky_background=13)

print(">>> TEST TARGET> tar2")
tar2 = vegapy.Target(band='H', pixel_scale=0.1*u.arcsec, FoV=(2*u.arcmin, 2*u.arcmin), star_table=SOURCE_PATH+'example/smalltarget4_stars.dat', sky_background=14)
if verbose > 0:
	print(tar2.shape)
if visual > 0:
	vegapy.imshow(tar2.data, scale='log', colorbar={'pad': 0.0})
# if False:
# 	from astropy.io import fits
# 	fits.writeto('target_data.fits', tar2.data.value, overwrite=True)

print(">>> TEST TARGET> tar3")
tar3 = vegapy.Target(band='H', FoV=(2*u.arcmin, 2*u.arcmin), shape=(64, 64), sky_background=13)



try:
	tar.data.unit
	tar.__dict__
	print(">>> TEST TARGET: Successful")
except:
	print(">>> TEST TARGET: NOT successful")

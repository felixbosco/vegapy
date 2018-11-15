# Import
import vegapy
import os
import astropy.units as u

# Settings
verbose = 0
visual = 0

# Initialize
PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, '../source/')
print('>>> TEST: {}'.format(FILE))



tel = vegapy.Telescope(8.0*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=SOURCE_PATH+'example/scao_2ms_imav.fits')
#tel = Telescope(8.2*u.m, psf_source='seeing', seeing_fwhm=0.8*u.arcsec, psf_resolution=0.022*u.arcsec, size=128)

if visual > 0:
	vegapy.imshow(tel.psf, colorbar={'pad': 0.0})
if verbose > 0:
	print(tel)

try:
	tel.__dict__
	print(">>> TEST TELESCOPE: Successful")
except:
	print(">>> TEST TELESCOPE: NOT successful")

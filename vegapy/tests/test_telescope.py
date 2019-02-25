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


tel_static = vegapy.Telescope(8.0*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=SOURCE_PATH+'example/scao_correction_10s_longexposure.fits')
assert np.abs(np.sum(tel_static.psf) - 1.0) < 1e-6

tel_compute = vegapy.Telescope(8.2*u.m, psf_source='seeing', seeing_fwhm=0.8*u.arcsec, psf_resolution=0.022*u.arcsec, size=128)
assert np.abs(np.sum(tel_compute.psf) - 1.0) < 1e-6

if visual > 0:
	vegapy.imshow(tel_compute.psf, colorbar={'pad': 0.0})
if verbose > 0:
	print(tel_compute)

tel_nonstatic = vegapy.Telescope(8.2*u.m, central_obscuration=0.14, psf_source=SOURCE_PATH+'example/seeing_2ms_shortexposures.fits')
tel_nonstatic(np.ones((64, 64)), 20*u.mas, integration_time=0.2*u.s)
assert np.abs(np.sum(tel_nonstatic.psf) - 1.0) < 1e-6

tel_airy = vegapy.Telescope(8.2*u.m, psf_source='airy_model', wavelength=1.63*u.micron)
# print('Airy model:', tel_airy)

try:
	tel_static.__dict__
	tel_compute.__dict__
	tel_nonstatic.__dict__
	print(">>> TEST TELESCOPE: Successful")
except:
	print(">>> TEST TELESCOPE: NOT successful")

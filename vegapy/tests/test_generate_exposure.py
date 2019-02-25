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

target = vegapy.Target(band='H', FoV=11*u.arcsec*2, shape=(1024, 1024), star_table=SOURCE_PATH+'example/star_table_example.dat', sky_background=14.4)
telescope = vegapy.Telescope(8.2*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=SOURCE_PATH+'example/seeing_2ms_shortexposures.fits')
detector = vegapy.Detector((1024, 1024),
						pixel_size=0.0106*u.arcsec,
						readout_noise=35*u.electron/u.pix,
						system_gain=17*u.electron/u.adu,
						optics_transmission=0.9,
                 		quantum_efficiency=0.9*u.electron/u.ph,
						saturation_level=7200*u.adu)
DIT = 0.2*u.s

if verbose > 0:
	print(target)
	print(telescope)
	print(detector)
if visual > 0:
	imshow(telescope.psf, scale='log')
if verbose > 0:
	print('Resolutions')
	print('Target', target.resolution)
	print('PSF', telescope.psf_resolution)
	print('Detector', detector.resolution)


vegapy.generate_exposure(target, telescope, detector, DIT, number_frames=10, filename=None, verbose=0, maximum_number_frames_per_file=5, randomkeywordwithoutmeaning=3)

print('>>> TEST: Successful')

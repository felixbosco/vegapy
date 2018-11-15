import vegapy
import os
import astropy.units as u


PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, '../source/')


target = vegapy.Target(band='H', FoV=22*u.arcsec, shape=(1024, 1024), star_table=SOURCE_PATH+'example/smalltarget4_stars.dat', sky_background=13)
telescope = vegapy.Telescope(8.2*u.m, central_obscuration=0.14, name="VLT Unit Telescope", psf_source=SOURCE_PATH+'example/noao_2ms_psf.fits')
detector = vegapy.Detector((1024, 1024),
						pixel_size=0.0106*u.arcsec,
						readout_noise=35*u.electron/u.pix,
						system_gain=17*u.electron/u.adu,
						optics_transmission=0.9,
                 		quantum_efficiency=0.9*u.electron/u.ph,
						saturation_level=7200.*u.adu)


vegapy.generate_exposure(target, telescope, detector, DIT=0.2*u.s, number_frames=5, filename='my_exposure.fits')

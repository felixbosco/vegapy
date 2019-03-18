from vegapy import Detector
import astropy.units as u

Hawaii2RG = Detector((2048, 2048),
						pixel_size=0.0106*u.arcsec,
						readout_noise=35*u.electron/u.pix,
						system_gain=17*u.electron/u.adu,
						optics_transmission=0.9,#.5,
                 		quantum_efficiency=0.9*u.electron/u.ph,
						saturation_level=7200.*u.adu)

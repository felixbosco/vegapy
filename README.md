# VEGAPy: A Virtual Exposure Generator for Astronomy in Python

This package is designed for creating realistic (short exposure) images, including e.g. read noise, photon noise and short exposure point spread functions (PSFs). For an example session, see the file example_session.py in the examples sub-package or the example in the _Basic Usage_ section.

## Migration to Specklepy
The package was migrated into [Specklepy](https://github.com/felixbosco/specklepy), see [generate module](https://github.com/felixbosco/specklepy#synthetic-image-generator). You can now use it with a simple `specklepy generate parfile.par` sequence from the terminal, where all your parameters are stored in that parameter file.

## Basic usage
In principle, you just need to define instances of the three package intern classes for the science target (Target), the telescope optics (Telescope) and the detector readout electronics (Detector). The classes are explained in detail below, but also see the documentary strings in the class definitions. Once defined, the class instances are arguments to the function generate_exposure(), which computes the synthetic exposures based on the class parameters.

```python

import vegapy
import astropy.units as u

target = vegapy.Target(band='H',
                        FoV=22*u.arcsec,
                        shape=(1024, 1024),
                        star_table='my_star_table_file.dat',
                        sky_background=13)
telescope = vegapy.Telescope(8.2*u.m,
                             central_obscuration=0.14,
                             name="VLT Unit Telescope",
                             psf_source='my_psf_file.fits')
detector = vegapy.Detector((1024, 1024),
                           pixel_size=0.01*u.arcsec,
                           readout_noise=35*u.electron/u.pix,
                           system_gain=17*u.electron/u.adu,
                           optics_transmission=0.9,
                           quantum_efficiency=0.9*u.electron/u.ph,
                           saturation_level=7200*u.adu)

vegapy.generate_exposure(target, telescope, detector, DIT=0.2*u.s, number_frames=5, filename='my_exposure.fits')
```

## About the classes and principle function
Please find an overview of the three object classes and the principle function below.

### Target class

### Telescope class

### Detector class

### generate_exposure() function

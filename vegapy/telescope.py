import os
import sys
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.modeling import models
import warnings
from scipy.signal import fftconvolve
from scipy.ndimage import zoom


class Telescope(object):

	"""
	Attributes are:
		'diameter':
		'psf_source':
		'psf_plane':

	Optional attributes are:
		'central_obscuration':

	Future features are:
		Mode 'airy_model': Shall compute the psf with the Fourier transform of the aperture instead of a file.
		Mode 'seeing': Shall compute the psf as a Gaussian seeing disk instead of a file.
	"""

	TIME_STEP_KEYS = ['TIMESTEP', 'INTTIME', 'CDELT3']
	RESOLUTION_KEYS = ['PIXSIZE', 'CDELT1']


	def __init__(self, diameter, psf_source, psf_plane=0, **kwargs):
		# Read input parameters
		self.diameter = diameter
		self.psf_source = psf_source
		self.psf_plane = psf_plane
		for key in kwargs:
			self.__setattr__(key, kwargs[key])

		# Compute secondary parameters
		if hasattr(self, 'central_obscuration'):
			self.area = (1. - self.central_obscuration**2) * np.pi * (self.diameter / 2)**2
		else:
			self.area = np.pi * (self.diameter / 2)**2

		if isinstance(psf_source, str):
			if psf_source == 'airy_model':
				self.compute_airy_model()
			elif psf_source == 'seeing':
				self.compute_seeing_disk(kwargs)
			else:
				self.read_psf_file(psf_source)
		else:
			raise TypeError('psf_source must be str-type, but is given as {}'.format(type(psf_source)))


	def __call__(self, flux_array, flux_array_resolution, integration_time=None, verbose=0, **kwargs):
		"""Requires the integration_time only if the PSF is non-static."""
		tmp = flux_array * self.area
		total_flux = np.sum(tmp)
		tmp_unit = tmp.unit

		# Resample flux_array to psf resolution
		try:
			ratio = float(flux_array_resolution / self.psf_resolution)
			#if ratio < 1.0:
			#	print('')
			#	print('Be cautious, the resolution of the PSF is worse than of the target data.')
			#	#raise ValueError('Be cautious, the resolution of the PSF  ({}) is worse than of the target data ({}).'.format(self.psf_resolutio, flux_array_resolution))
		except UnitConversionError as e:
			raise UnitConversionError("The resolution values of the image ({}) and PSF ({}) have different units!".format(flux_array_resolution, self.psf_resolution))

		# Prepare PSF if non-static
		if hasattr(self, 'timestep'):
			if integration_time is None:
				raise ValueError("If the PSF source of Telescope is non-static, the call function requires the integration_time.")
			self.integrate_psf(integration_time=integration_time)

		#convolved = np.zeros(tmp.shape)
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			if ratio < 1.0:
				self.psf = zoom(self.psf, 1/ratio, order=1) / ratio**2
				self.psf = self._normalize(self.psf)
			else:
				tmp = zoom(tmp, ratio, order=1) / ratio**2
		if tmp.shape[0] > 2048+512 or self.psf.shape[0] > 512+256:
			print('With these sizes (image: {} and PSF: {}), the computation will be very expensive. It is suggested to adapt the resolution of the objects.'.format(tmp.shape, self.psf.shape))
			user_input = raw_input('Do you still want to continue? [Y/N]')
			if user_input in ['Y', 'y', 'Yes', 'yes']:
				pass
			else:
				raise Exception('Program aborted, re-define the resolution of the objects.')
		convolved = fftconvolve(tmp, self.psf, mode='same') * tmp_unit
		if verbose > 0:
			print('Check of flux conservation during convolution:')
			print('Before: ', total_flux)
			print('After:  ', np.sum(convolved))
		return convolved


	def __str__(self):
		tmp = "Telescope:\n"
		for key in self.__dict__:
			if key == 'psf':
				continue
			tmp += "{}: {}\n".format(key, self.__dict__[key])
		return tmp


	def read_psf_file(self, filename, hdu_entry=0):
		with fits.open(filename) as hdulist:
			header = hdulist[hdu_entry].header
		if header['NAXIS'] == 2:
			with fits.open(self.psf_source) as hdulist:
				self.psf = self._normalize(hdulist[hdu_entry].data)
		else:
			for key in self.TIME_STEP_KEYS:
				try:
					self.timestep = self._get_value(header, key)
					break
				except KeyError as e:
					continue
		for key in self.RESOLUTION_KEYS:
			try:
				self.psf_resolution = self._get_value(header, key)
				break
			except KeyError as e:
				continue
			raise IOError("No key from {} was found in file for the psf resolution.".format(self.RESOLUTION_KEYS))


	def _get_value(self, header, key, alias_dict={'sec': 's', 'milliarcsec': 'mas', 'microns': 'micron'}, verbose=False):
		"""
		The alias_dict dictionary is used as a mapping from
		unvalid unit strings for u.Unit(str). Feel free to
		add new aliases.
		"""

		value = header[key]
		unit_str = header.comments[key]

		if unit_str == '':
		    if verbose:
		        print("Function 'get_value()' did not find a unit in the comment string.")
		    return value
		else:
		    try:
		        unit = u.Unit(unit_str)
		    except ValueError as e:
		        if verbose:
		            print("ValueError:", e)
		            print("Trying aliases from alias_dict...")
		        try:
		            unit = u.Unit(alias_dict[unit_str])
		        except:
		            raise IOError("Found no matching key in the alias_dict. You may add the corresponding entry.")
		    return value * unit


	def _normalize(self, array, mode='unity_circular'):
		"""Normalizes the array to either have a sum of 1 ('unity' mode) or that the peak value is 1 ('max' mode)."""
		if mode == 'unity':
		    return array / np.sum(array)
		elif mode == 'max':
		    return array / np.max(array)
		elif mode == 'unity_circular':
		    x, y = array.shape
		    low_cut = array[0, int(y/2)]
		    array = np.maximum(array - low_cut, 0)
		    return array / np.sum(array)


	def integrate_psf(self, integration_time, hdu_entry=0):
		number_planes = int(integration_time / self.timestep)
		with fits.open(self.psf_source) as hdulist:
			data = hdulist[hdu_entry].data

			self.psf_plane += 1
			if self.psf_plane + number_planes < data.shape[0]:
				self.psf = np.sum(data[self.psf_plane : self.psf_plane+number_planes], axis=0)
			else:
				self.psf = np.sum(data[self.psf_plane : ], axis=0)
				self.psf += np.sum(data[ : (self.psf_plane+number_planes) % data.shape[0]], axis=0)
			self.psf_plane += number_planes - 1
			self.psf_plane = self.psf_plane % data.shape[0]
			#Normalization
			self.psf = self._normalize(self.psf)


	def _distance(self, index, reference):
	    return np.sqrt(np.square(index[0]-reference[0]) + np.square(index[1]-reference[1]))


	def compute_airy_model(self, size=256, wavelength=None):
		"""Not implemented properly yet."""
		raise NameError("Function 'compute_airy_model()' is not implemented yet.")
		physical_size = self.diameter#8.2*u.m
		#central_obscuration = 1.0*u.m

		shape = (size, size)
		radius = size / 2
		center = (radius, radius)
		pixel_scale = physical_size / size
		array = np.ones(shape=shape)
		mask = np.zeros(shape=shape, dtype=bool)

		for index, value in np.ndenumerate(mask):
		    if self._distance(index, center) > radius or self._distance(index, center) < radius * self.central_obscuration:
		        mask[index] = 1

		aperture = np.ma.masked_array(array, mask=mask)
		psf = np.square(np.abs(np.fft.fftshift(np.fft.fft2(aperture))))
		#self.psf = psf
		#self.psf_resolution = pixel_scale


	def _gaussian_2d(self, index, center, fwhm, amplitude=1):
		if not isinstance(fwhm, tuple):
			std = 2 * np.sqrt(2 * np.log(2)) * fwhm
			std = (std, std)
		else:
			std = (2 * np.sqrt(2 * np.log(2)) * fwhm[0], 2 * np.sqrt(2 * np.log(2)) * fwhm[1])
		return amplitude * np.exp(- np.square((index[0] - center[0]) / std[0]) / 2 - np.square((index[1] - center[1]) / std[1]) / 2)


	def compute_seeing_disk(self, kwargs={}):
		"""function compute_seeing_disk()

		Approximates the seeing disk with a 2D Gaussian of with given FWHM.
		Yet the function is very slow.
		Caution: Not tested extensively yet."""

		print('Caution: The function compute_seeing_disk() has not been tested extensively yet.')

		# Handling input
		for key in ['seeing', 'fwhm', 'seeing_fwhm']:
			if key in kwargs:
				fwhm = kwargs[key]
				break
				# if not isinstance(fwhm, u.Quantity):
				# 	print('Treating the FWHM of the seeing PSF as in units of arcsec.')
				# 	fwhm *= u.arcsec
		try:
			fwhm
		except NameError as e:
			print('Using default FWHM of 1.0arcsec.')
			fwhm = 1.0*u.arcsec

		for key in ['psf_resolution', 'resolution']:
			if key in kwargs:
				fwhm = kwargs[key]
				break
		if not hasattr(self, 'psf_resolution'):
			print('Using default psf_resolution of 20mas per pix.')
			self.psf_resolution = 0.02*u.arcsec

		for key in ['size', 'psf_size']:
			if key in kwargs:
				size = kwargs[key]
				break
		try:
			size
		except NameError as e:
			print('Using default size for PSF of 256.')
			size = 256

		shape = (size, size)
		center = ((size + 1) / 2, (size + 1) / 2)
		xdata, ydata = np.mgrid[ : shape[0], : shape[1]]
		stddev = fwhm / self.psf_resolution / (2 * np.sqrt(2 * np.log(2)))
		Gaussian = models.Gaussian2D(amplitude=1.0, x_mean=center[0], y_mean=center[1], x_stddev=stddev, y_stddev=stddev)
		psf = Gaussian(xdata, ydata)
		#psf = np.zeros(shape)
		#for index, value in np.ndenumerate(psf):
		#	psf[index] = self._gaussian_2d(index, center, fwhm/self.psf_resolution)
		self.psf = self._normalize(psf)


if __name__ == '__main__' and '-notest' not in sys.argv:
	# By giving the '-notest' option, these tests will not be executed.
	print("TEST> Entering tests...")
	import tests.test_telescope

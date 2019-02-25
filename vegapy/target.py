import os
import sys
import numpy as np
import astropy.units as u
from astropy import constants as const
from astropy.table import Table


PATH, FILE = os.path.split(__file__)
SOURCE_PATH = os.path.join(PATH, './source/data/')


class Target(object):

	def __init__(self, band, **kwargs):
		"""
		This function accepts two out of the three keywords:
			Tuple(int): 'shape'
			Tuple(Astropy.Unit): 'FoV'
			Astropy.Unit: 'pixel_scale'

		Optional attributes are:
			'sky_background':
			'config_file':
			'star_table':
			'number_stars':

		Future features:
			Replace the FoV attribute by a single value instead of a tuple, in case that the FoV is a square.
			The object shall obtain a phase center in ICRS, which is set to RA=00:00:00.00 and DEC=00:00:00.00
		"""

		# Read input parameters
		self.band = band
		for key in kwargs:
			self.__setattr__(key, kwargs[key])

		# Compute secondary parameters
		if 'shape' in kwargs and 'FoV' in kwargs and 'pixel_scale' in kwargs:
			# assert that the three values match
			raise NotImplementedError('Handling of all three parameters shape, FoV, and pixel_scale is not implemented yet. Give only two of them.')
		elif 'shape' in kwargs and 'FoV' in kwargs:
			# assert that the two values match
			shape = kwargs['shape']
			FoV = kwargs['FoV']
			if not isinstance(FoV, tuple):
				FoV = (FoV, FoV)
			try:
				assert ( FoV[0] / shape[0] - FoV[1] / shape[1]).value < 1e-6
			except AssertionError as e:
				raise ValueError("The field of view (FoV) and shape of Target must have the same relative size along both axes.")
			self.shape = shape
			self.FoV = FoV
			self.data = np.zeros(self.shape)
			self.pixel_scale = self.FoV[0] / self.shape[0]
		elif 'pixel_scale' in kwargs and 'FoV' in kwargs:
			self.FoV = kwargs['FoV']
			self.pixel_scale = kwargs['pixel_scale']
			self.shape = (int(self.FoV[0] / self.pixel_scale), int(self.FoV[1] / self.pixel_scale))
			self.data = np.zeros(shape=self.shape)
		elif 'shape' in kwargs and 'pixel_scale' in kwargs:
			self.shape = kwargs['shape']
			self.data = np.zeros(self.shape)
			self.pixel_scale = kwargs['pixel_scale']
			self.FoV = (self.data.shape[0] * self.pixel_scale, self.data.shape[1] * self.pixel_scale)
		else:
			print(kwargs)
			if not 'config_file' in kwargs:
				raise ValueError("Target() requires exactly two out of the three keywords 'shape', 'FoV', and 'pixel_scale'.")
		self.data = self.data * u.ph / u.m**2 / u.s
		self.resolution = self.pixel_scale
		self.band_reference_flux = self._get_reference_flux()

		# Handle optional parameters
		if hasattr(self, 'config_file'):
			self._read_config_file()
		if hasattr(self, 'sky_background'):
			self._initialize_sky_background()
		if hasattr(self, 'number_stars'):
			self._generate_stars()
		if hasattr(self, 'star_table'):
			if hasattr(self, 'flux_unit'):
				self._read_star_table(self.flux_unit)
			else:
				self._read_star_table()


	def __call__(self):
		return self.data.decompose()


	def __str__(self):
		tmp = "Target:\n"
		for key in self.__dict__:
			if key == 'stars' or key == 'data':
				continue
			tmp += "{}: {}\n".format(key, self.__dict__[key])
		return tmp


	@property
	def resolution(self):
		return self.pixel_scale


	@resolution.setter
	def resolution(self, value):
		self.pixel_scale = value


	def _get_reference_flux(self, filename=SOURCE_PATH+'photometric_bands.dat', format='ascii'):
		#print("Caution: The function 'Target._get_reference_flux()' uses a problematic relative path.")
		table = Table.read(filename, format=format)
		row_index = np.where(table["Band"] == self.band)
		fwhm = table['FWHM'][row_index][0]
		flux = table['Flux'][row_index][0] * u.Jy
		return (flux / const.h * fwhm * u.ph).decompose()


	def _get_flux(self, magnitude):
		return 10**(magnitude/-2.5) * self.band_reference_flux


	def _initialize_sky_background(self):
		"""
		This function computes the sky background flux per pixel. There is no handler for other
		data types than int, float, or Astropy.unit.Quantity and the latter does only draw the
		value attribute from the quantity.
		"""
		if isinstance(self.sky_background, int) or isinstance(self.sky_background, float):
			# Interpreting as mag / arcsec**2
			self.sky_background_flux = self._get_flux(self.sky_background) / u.arcsec**2
		elif isinstance(self.sky_background, u.Quantity):
			# Interpreting as mag / arcsec**2
			print("Caution: This function interpretes sky_background as in units of mag per arcsec**2.")
			self.sky_background_flux = self._get_flux(self.sky_background.value) / u.arcsec**2
		else:
			raise TypeError("Function 'Target._initialize_sky_background()' does not accept a magnitude of type {}.".format(type(self.sky_background)))
		flux_per_pixel = (self.sky_background_flux * self.pixel_scale**2).decompose()
		# print(type(flux_per_pixel), type(self.data))
		# print(flux_per_pixel.unit, self.data.unit)
		# if isinstance(self.data, np.ndarray):
		# 	self.data = np.maximum(self.data, flux_per_pixel.value) * flux_per_pixel.unit
		# elif isinstance(self.data, u.quantity.Quantity):
		# 	self.data = np.maximum(self.data, flux_per_pixel)
		# else:
		# 	self.data = np.maximum(self.data, flux_per_pixel)
		self.data = np.maximum(self.data, flux_per_pixel)


	def _generate_stars(self, saveto='star_table_latest.dat'):
		pass
		pos_x = np.zeros((self.number_stars))
		pos_y = np.zeros((self.number_stars))
		mag = np.zeros((self.number_stars))
		flux = np.zeros((self.number_stars))
		star_table = Table([pos_x, pos_y, mag, flux], names=('pos_x', 'pos_y', 'mag', 'flux'))
		for n in range(self.number_stars):
			pass
		try:
			star_table.write(saveto, format='ascii')
		except:
			pass


	def _to_map_unit(self, unit, filename=SOURCE_PATH+'photometric_bands.dat', format='ascii'):
		print('Interpreting flux values in units of {}. To change this, please provide a star_table_dict={"flux_unit": "desired unit"}.'.format(unit))
		if isinstance(unit, str):
			unit = u.Unit(unit)
		try:
			tmp =  unit.to(self.data.unit)
		except UnitConversionError as e:
			table = Table.read(filename, format=format)
			row_index = np.where(table["Band"] == self.band)
			fwhm = table['FWHM'][row_index][0]
			wavelength = table['Wavelength'][row_index][0] * u.micron
			tmp = (unit / const.h * u.ph * fwhm * wavelength).decompose()
			tmp = unit.to(self.data.unit)
		return tmp


	def _read_star_table(self, keyword_x='pos_x', keyword_y='pos_y', keyword_flux='flux'):
		self.stars = Table.read(self.star_table, format='ascii')
		for row in self.stars:
			position = (int(row[keyword_x]), int(row[keyword_y]))
			try:
				self.data[position] = np.maximum(self.data[position], row[keyword_flux])
			except KeyError as e:
				flux = self._get_flux(row['mag'])
				self.data[position] = np.maximum(self.data[position], flux)
			except:
				self.data[position] = np.maximum(self.data[position].value, row['flux']) * self.data.unit
		#print(self.data)


	def _read_config_file(self):
		pass


if __name__ == '__main__' and '-notest' not in sys.argv:
	# By giving the '-notest' option, these tests will not be executed.
	print("TEST> Entering tests...")
	import tests.test_target

from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(name='VEGAPy',
      version='0.1.dev1',
      description='A Python generator for synthetic astronomical exposures.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Data Processing :: Astronomy',
        'Operating System :: OS Independent',
      ],
      keywords='astronomy synthetic images short exposure',
      url='https://github.com/felixbosco/vegapy',
      author='Felix Bosco',
      author_email='bosco@mpia.de',
      license='MIT',
      packages=['vegapy'],
      install_requires=[
          'astropy',
      ],
      include_package_data=True,
      zip_safe=False)

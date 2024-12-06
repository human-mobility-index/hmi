# coding: utf-8
import os
import sys
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='hmi',
      version='0.1',
      description="The ``hmi`` package is a python module that provides an interface to compute distances between locations based on [Özak's (2018, ](http://rdcu.be/I4YI)[2010)](http://omerozak.com/pdf/Ozak_voyage.pdf) [Human Mobility Index](https://human-mobility-index.github.io/).",
      url='http://github.com/ozak/hmi',
      keywords="distances mobility migratory gis geospatial geographic raster vector zonal statistics spatial analysis",
      author='Ömer Özak',
      author_email='omer@omerozak.com',
      license='GPLv3',
      #package_dir={'': 'src'},
      packages=['hmi'],
      long_description=read('README.md'),
      install_requires=[
                'numpy',
                'pandas',
                'docopt',
                'GDAL',
                'pyproj',
                'scikit-image',
                'matplotlib',
                'coverage',
                'fiona',
                'geopandas',
                'pysal',
                'affine',
                'rasterstats',
                'georasters',
                'shapely',
                'pysal',
            ],
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Economics',
          'Intended Audience :: History',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: GIS',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
      ],
      zip_safe=False,
      )

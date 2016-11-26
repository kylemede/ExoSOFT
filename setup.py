#!/usr/bin/env python
import numpy as np
import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup
    
from distutils.extension import Extension

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# Handle encoding
major, minor1, minor2, release, serial = sys.version_info
if major >= 3:
    def rd(filename):
        f = open(filename, encoding="utf-8")
        r = f.read()
        f.close()
        return r
else:
    def rd(filename):
        f = open(filename)
        r = f.read()
        f.close()
        return r

setup(    
    name='ExoSOFT', 
    version="1.0.0", 
    author='Kyle Mede',
    author_email = 'kylemede@gmail.com',
    url = 'https://github.com/kylemede/ExoSOFT',
    packages =['ExoSOFT'],
    license = ['GNU GPLv3'],
    description ='Exoplanet Simple Orbit Fitting Toolbox',
    long_description="For further details, please read the documentation at\nhttp://ExoSOFT.readthedocs.io/en/latest/",
    package_data={"": ["LICENSE", "AUTHORS.rst"]},
    include_package_data=True,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    #include_dirs = [np.get_include()],
    install_requires = ['numpy','pyfits','scipy','pylab','matplotlib',\
                        'psutil','astropy','KMlogger','emcee'],
    #ext_modules=[]
)

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
    
#from distutils.extension import Extension
from Cython.Build import cythonize

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
    version="1.0.66", 
    author='Kyle Mede',
    author_email = 'kylemede@gmail.com',
    url = 'https://github.com/kylemede/ExoSOFT',
    packages =['ExoSOFT','ExoSOFT.tools'],
    license = ['GNU GPLv3'],
    description ='Exoplanet Simple Orbit Fitting Toolbox',
    long_description="For further details, please read the documentation at\nhttp://ExoSOFT.readthedocs.io/en/latest/",
    package_data={"": ["LICENSE", "AUTHORS.rst","ExoSOFT/tools/cytools.c","ExoSOFT/tools/cytools.pyx"]},
    package_dir = {"ExoSOFT":'ExoSOFT'},
    scripts = ['scripts/ExoSOFT'],
    include_package_data=True,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    include_dirs = ['ExoSOFT','ExoSOFT/tools','examples'],
    install_requires = ['numpy','pyfits','scipy','matplotlib',\
                        'psutil','astropy','KMlogger','emcee'],
    ext_modules = cythonize(["ExoSOFT/tools/cytools.pyx"]),
    #ext_modules=[]
    #ext_modules=[ Extension("cytools",["ExoSOFT/tools/cytools.c"]),Extension("cytools",["ExoSOFT/tools/cytools.pyx"]) ]
)

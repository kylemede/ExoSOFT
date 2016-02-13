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
    author_email = 'kylemede@astron.s.u-tokyo.ac.jp',
    url = 'https://github.com/kylemede/??????????',
    packages =['exosoft'],
    license = ['GNU GPLv3'],
    description ='Exoplanet Simple Orbit Fitting Toolbox',
    long_description=rd("README.rst") + "\n\n"
                    + "Changelog\n"
                    + "---------\n\n"
                    + rd("HISTORY.rst"),
    package_data={"": ["LICENSE", "AUTHORS.rst"]},
    include_package_data=True,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    #include_dirs = [np.get_include()],
    install_requires = ['numpy''pyfits','scipy','pylab','matplotlib','psutil',\
                        'swig','astropy','PyAstronomy'],
    #ext_modules=[]
)

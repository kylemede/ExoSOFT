#!/usr/bin/env python
#import numpy as np
import os
import sys
#import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup
    
from distutils.extension import Extension
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
    version="0.1.2", 
    author='Kyle Mede',
    author_email = 'kylemede@gmail.com',
    url = 'https://github.com/kylemede/ExoSOFT',
    packages =['ExoSOFT','ExoSOFT.tools'],
    license = ['GNU GPLv3'],
    description ='The Exoplanet Simple Orbit Fitting Toolbox (ExoSOFT)',
    long_description="For further details, please read the documentation at http://ExoSOFT.readthedocs.io/en/latest/,"+\
    " and visit the GitHub repository https://github.com/kylemede/ExoSOFT",
    # Forcefully adding the non-python files to the package to ensure they get in
    data_files = [('ExoSOFT/tools',["./ExoSOFT/tools/cytools.c"]) , ('ExoSOFT/tools',["./ExoSOFT/tools/cytools.pyx"])],
    package_data = {'ExoSOFT/tools':["./ExoSOFT/tools/cytools.c"] , 'ExoSOFT/tools':["./ExoSOFT/tools/cytools.pyx"]},
    package_dir = {"ExoSOFT":'ExoSOFT', "ExoSOFT.tools":'ExoSOFT/tools'},
    scripts = ['scripts/ExoSOFT'],
    include_package_data=True,
    classifiers = [
        'Development Status :: 3 - Alpha',#5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    include_dirs = ['ExoSOFT','ExoSOFT/tools'],
    install_requires = ['cython','numpy','scipy','matplotlib','jdcal',\
                        'astropy','KMlogger','emcee','pathos','corner','pyyaml'],
    # use cythonize to compile cython module cytools.pyx to both a .c and .so
    ext_modules = cythonize([Extension("ExoSOFT.tools.cytools",["./ExoSOFT/tools/cytools.pyx"])]),
)

ExoSOFT
=======

[![Build Status](https://travis-ci.org/kylemede/ExoSOFT.svg?branch=master)](https://travis-ci.org/kylemede/ExoSOFT)
[![DOI](https://zenodo.org/badge/74822648.svg)](https://zenodo.org/badge/latestdoi/74822648)
[![PyPI version](https://badge.fury.io/py/ExoSOFT.svg)](https://badge.fury.io/py/ExoSOFT)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://github.com/kylemede/ExoSOFT/blob/master/LICENSE)
<!--[![Coverage Status](https://coveralls.io/repos/github/kylemede/ExoSOFT/badge.svg?branch=master)](https://coveralls.io/github/kylemede/ExoSOFT?branch=master)-->

*(Build tests on Travis CI currently being done for Python versions 2.7, 3.3 and 3.4.)*

Documentation of the code can be found [HERE](http://exosoft.readthedocs.io/en/latest/index.html) *(Note: Documentation is still a work in progress)*

**The Exoplanet Simple Orbit Fitting Toolbox**

ExoSOFT is a toolbox for the orbital analysis of exoplanets and binary star
systems.  It is free to compile, open source, fits any combination of
astrometric and radial velocity data, and offers 4 parameter space exploration
techniques, including MCMC.  Additionally, it is packaged together with an
automated set of post-processing and plotting routines to summarize the results
.  Verifications were performed by fitting noisy synthetic data produced with
an independent Keplerian model for a variety of systems ranging the entire
parameter space.  

The release paper reviewing ExoSOFT and its verification are available [here](http://iopscience.iop.org/article/10.3847/1538-3881/aa5e4a/meta;jsessionid=9AABE88A2AA6F5D49DFDEDA7E23B9824.c4.iopscience.cld.iop.org).

Examples of its use on real data include:

 - An investigation of the V450And system [Helminiak & Kuzuhara & Mede et. al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...832...33H),

 - and the Tau AB binary [Mede & Brandt (2014)](http://adsabs.harvard.edu/abs/2014IAUS..299...52M>).

ExoSOFT would be a suitable choice for performing orbital analysis during surveys with new RV and direct imaging instruments.

Installation
============

**Option 1:**


 $ pip install ExoSOFT

 An example of a full installation from scratch on a Mac, using conda:

 $ conda create --name ExoSOFTcondaEnv python Cython gcc scipy numpy matplotlib

 $ source activate ExoSOFTcondaEnv

 $ pip install ExoSOFT

**Option 2:**
 $ conda create --name ExoSOFTcondaEnv python Cython gcc scipy numpy matplotlib

 $ source activate ExoSOFTcondaEnv

 - clone this repository

 - from its directory on your local machine, run:

 $ python setup.py install


How to run ExoSOFT
==================

Make a valid settings.yaml file, matching the format in the [repository's examples directory](https://github.com/kylemede/ExoSOFT/tree/master/examples).
  Depending on the system you are investigating, you will also need to have data files of the correct format ([data format examples here](https://github.com/kylemede/ExoSOFT/tree/master/examples)).

** *Make sure to update the directory where ExoSOFT should write all output files to, as indicated by the 'outDir' setting in your local settings.py file.* **

Then, from the directory containing the settings.yaml, simply start ExoSOFT with:

 $ ExoSOFT

Should you wish to provide custom priors, a proper priors.py file must also be in the same directory.  To make one, copy the [default file](https://github.com/kylemede/ExoSOFT/blob/master/examples/priors.py) and edit accordingly.

There are two working examples, one for that runs the stages Simulated Annealing then emcee (settings_SAemcee.py).  The other runs Simulated Annealing, Sigma Tuning, then MCMC (settings_SASTMCMC.py).  Either copy their settings file in the examples directory to a new file or rename them to settings.py and run it the same as above.


Additionally, there are two basic examples of how to interact with the model more directly.  These examples are given in:

 - model\_example\_basic.py
 - model\_example\_emcee.py

 These will show the user how to write their own scripts that use the ExoSOFT model instead of using ExoSOFT itself.

Solutions to some install or runtime errors
===========================================

**Solution to 'binary incompatible' error:**

Depending on the install of the scientific python stack on your machine, you may get a 'binary incompatibility' error.
[This arrises from the copy of scipy currently installed having being built against a version of numpy<1.8, which does not catch the error thrown by Cython.](http://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility)  
The solution is to uninstall and re-build numpy, scipy.  It takes a while to build scipy especially, but this worked for me.  

The commands in order are:

 $ pip uninstall -y numpy   

 $ pip uninstall -y scipy

 <!--$pip uninstall -y scikit-learn-->

 $ pip install numpy --ignore-installed --no-cache-dir --no-binary :all:

 $ pip install scipy --ignore-installed --no-cache-dir --no-binary :all:

 <!--$pip install scikit-learn --ignore-installed --no-cache-dir --no-binary :all:-->

**Solution to problems with matplotlib**

Again, depending on your installation of matplotlib, some errors can arise.  
These include missing __init__.py in one of the matplotlib directories (mpl_toolkits in particular).
Another is the inability to import matplotlib.pyplot.

To fix this:

 - Make sure there is now locally installed version of matplotlib in your home directory that could confuse the python path.

 - re-install matplotlib from scratch

 $ pip uninstall -y matplotlib

 $ pip install matplotlib

Additional Modules
==================

There are some additional modules available in the git repository that are not part of the base ExoSOFT package.
These can be found in the [extra](https://github.com/kylemede/ExoSOFT/tree/master/extra) directory and include:

    - customPost.py
    - artificialDataMaker.py

customPost.py is for re-performing the post analysis of ExoSOFT with the ability to customize the plotting settings.  
This is in order to produce the exact versions of the plots and statistical calculations necessary for publication.

artificialDataMaker.py is to produce simulated/artificial data with pre-determined orbital elements.

For more details, please see the documentation in each module.

Attribution
===========

Please cite our [paper](http://iopscience.iop.org/article/10.3847/1538-3881/aa5e4a/meta;jsessionid=9AABE88A2AA6F5D49DFDEDA7E23B9824.c4.iopscience.cld.iop.org) if you find this code useful in your
research.  

The Bibtex entry for this paper is::

@article{1538-3881-153-3-135,
  author={Kyle Mede and Timothy D. Brandt},
  title={The Exoplanet Simple Orbit Fitting Toolbox (ExoSOFT): An Open-source Tool for Efficient Fitting of Astrometric and Radial Velocity Data},
  journal={The Astronomical Journal},
  volume={153},
  number={3},
  pages={135},
  url={http://stacks.iop.org/1538-3881/153/i=3/a=135},
  year={2017},
}

License
=======

Copyright 2016 Kyle Mede and contributors.

ExoSOFT is free software made available under the GNU GPLv3 license.
For details see the license.txt file.

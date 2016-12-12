ExoSOFT
=======

**The Exoplanet Simple Orbit Fitting Toolbox**

ExoSOFT is a toolbox for the orbital analysis of exoplanets and binary star 
systems.  It is free to compile, open source, fits any combination of 
astrometric and radial velocity data, and offers 4 parameter space exploration 
techniques, including MCMC.  Additionally, it is packaged together with an 
automated set of post-processing and plotting routines to summarize the results
.  Verifications were performed by fitting noisy synthetic data produced with 
an independent Keplerian model for a variety of systems ranging the entire 
parameter space.  

The release paper reviewing ExoSOFT and its verification is currently under review.

Examples of its use on real data include:

 -An investigation of the V450And system `Helminiak & Kuzuhara & Mede et. al. (2016)<http://adsabs.harvard.edu/abs/2016ApJ...832...33H`,
 
 -and the Tau AB binary `Mede & Brandt (2014)<http://adsabs.harvard.edu/abs/2014IAUS..299...52M>`.

ExoSOFT would be a suitable choice for performing orbital analysis during surveys with new RV and direct imaging instruments.


For further details, please read the `documentation <http://ExoSOFT.readthedocs.io/en/latest/>',
and visit the `GitHub repository <https://github.com/kylemede/ExoSOFT>`

Dependencies:
-------------
Note: Installing python packages with pip is best as it handles the version and 
dependency issues for you.  Here are the dependencies for ExoSOFT.

 ## Those NOT available on pip
 
Python, gcc, gcc-c++, texlive-epstopdf

 #Available for install with pip
 
cython, numpy, scipy, pylab, matplotlib, psutil, astropy, corner, KMlogger, emcee, pathos


Attribution
-----------

Please cite our soon to be publish paper if you find this code useful in your
research.  The Bibtex entry for this paper is::

 ?? STILL UNDER SUBMISSION AND NOT ON ARXIV ??


License
-------

Copyright 2016 Kyle Mede and contributors.

ExoSOFT is free software made available under the GNU GPLv3 license. 
For details see the license.txt file.

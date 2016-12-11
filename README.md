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
parameter space.  Examples of its use on real data include 
[Helminiak & Kuzuhara & Mede et. al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...832...33H).  
, for the V450And system, and the tau AB binary in 
[Mede & Brandt (2014)](http://adsabs.harvard.edu/abs/2014IAUS..299...52M>).  
ExoSOFT would be a suitable choice for performing orbital analysis during 
surveys with new RV and direct imaging instruments.

Installation
============

Option 1:

 $ pip install ExoSOFT
 
Option 2: 

 -clone this repository
 -from its directory on your local machine, run:
 
 $ python setup.py install
 
How to run ExoSOFT
==================

Make a valid settings.py file, matching the format in the [repository's examples directory](https://github.com/kylemede/ExoSOFT/tree/master/examples). 
  Depending on the system you are investigating, you will also need to have data files of the corect format ([data format examples here](https://github.com/kylemede/ExoSOFT/tree/master/examples)].

Then, from the directory containing the settings.py, simply start ExoSOFT with:

 $ExoSOFT
 
Should you wish to provide custom priors, a proper priors.py file must also be in the same directory.  To make one, copy the [default file](https://github.com/kylemede/ExoSOFT/blob/master/ExoSOFT/tools/priors.py) and edit accordingly.


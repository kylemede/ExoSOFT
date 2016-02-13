Installation Notes
==================

The package uses SWIG to interface the C++ code with Python.  To compile the 
user will also need the gcc-c++ compiler.  Aditionally, after making the plots
to sumarize the results, the eps files are converted to pdf for easier viewing.
To make this work the user will also need the texlive-epstopdf or the verison 
for your OS to allow the epstopdf command from the cmd.

Following install either by downloading from this repository or pip, make 
sure the directory to ExoSOFT/exosoft/ is added to your PATH in your .bashrc or
 .bash_profile.  
 
 Then, to run:

$python ExoSOFT.py /full/path/to/*settings.py

If ExoSOFT.py is not executable after install, you can use $chmod +x ExoSOFT.py
to make it so and avoid typing 'python' when starting it.

A set of example settings and data files are in ExoSOFT/examples.  This can be 
ran with:

$ExoSOFT.py 

as it will just assume that without a provided path, you want to play with the 
example.

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
`Helminiak & Kuzuhara & Mede et. al. (2016) <????????>`_.  
, for the V450And system, and the tau AB binary in 
`Mede & Brandt (2014) <http://adsabs.harvard.edu/abs/2014IAUS..299...52M>`_.  
ExoSOFT would be a suitable choice for performing orbital analysis during 
surveys with new RV and direct imaging instruments.

Attribution
-----------

Please cite our soon to be publish paper if you find this code useful in your
research.  The Bibtex entry for this paper is::

????


License
-------

Copyright 2016 Kyle Mede and contributors.

ExoSOFT is free software made available under the GNU GPLv3 license. 
For details see the license.txt file.
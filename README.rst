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



Dependencies:
-------------
Note: Installing python packages with pip is best as it handles the version and 
dependency issues for you.  Here are the dependencies for ExoSOFT.

 ## Those NOT available on pip
 
Python, gcc, gcc-c++, texlive-epstopdf

 #Available for install with pip
 
numpy, scipy, pylab, matplotlib, psutil, astropy, PyAstronomy, corner, KMlogger


--------------
Install & Run:
--------------
To get around the occiasional difficulties installing SWIG on a Mac or 
Windows machine, you can use a platform independent installation with conda.

FIRST SET UP AN ANACONDA ENVIRONMENT
------------------------------------
Download and install miniconda from: http://conda.pydata.org/docs/install/quick.html

1. $ conda create -n ExoSOFTcondaEnv pip python ipython numpy scipy matplotlib pylab swig astropy jdcal
2. $ source activate ExoSOFTcondaEnv
 ## add 'source activate ExoSOFTcondaEnv' to end of .bashrc to have it load at beginning of terminal load every time.
4. $ conda install -c asmeurer pango
5. $ pip install PyAstronomy corner psutil KMlogger
 ## matplotlib.pyplot doesn't exists till frontend is built, to do this:
7. $ ipython
8. $ import pylab 



Installation Notes
==================

The package uses SWIG to interface the C++ code with Python.  To compile the 
user will also need the gcc-c++ compiler.  Aditionally, after making the plots
to sumarize the results, the eps files are converted to pdf for easier viewing.
To make this work the user will also need the texlive-epstopdf or the verison 
for your OS to allow the epstopdf command from the cmd.

Following install either by downloading from this repository or pip, make 
sure the directory to ExoSOFT/exosoft/ is added to your PATH in your .bashrc or .bash_profile.  
 
Then, to run:

 $python ExoSOFT.py /full/path/to/*settings.py

If ExoSOFT.py is not executable after install, you can use $chmod +x ExoSOFT.py
to make it so and avoid typing 'python' when starting it.

A set of example settings and data files are in ExoSOFT/examples.  This can be 
ran with:

 $ExoSOFT.py 

as it will just assume that without a provided path, you want to play with the 
example.


NEXT UNPACK AND SET UP ExoSOFT
------------------------------
1. unzip in directory of choice 
   (lets say '/home/ExoSOFT/' for demonstration purposes.)
2. open '/home/ExoSOFT/exosoft/exosoftpath.py' and update the directory string
   for rootDir.
3. open the settings.py file in '/home/ExoSOFT/examples/' and update the 
   directories in the keys 'outDir', 'DIdataFile' and 'RVdataFile' 
   on lines 39-47.
4. Possibly not be necessary, but the cpp files might require being compiled 
   again on your machine to run.
   If so, from a bash terminal:
    $cd '/home/ExoSOFT/exosoft/tools/cppTools/'
    
    $make clean
    
    $make
   NOTE: If you have difficulties compiling, make sure SWIG is installed correctly.  The documentation for this is provided here:
    http://www.swig.org/Doc3.0/Preface.html#Preface_osx_installation
5. If the directories are updated to match the location on your machine and the 
   cpp code compiled, let's try and run ExoSOFT by:
    $cd '/home/ExoSOFT/exosoft/'
    
    $python ExoSOFT.py
6. If it runs properly, then check the outputs when finished in the directory 
   you set 'outDir' to.  Else, the errors are most likely dependancy based, so 
   please check the traceback to solve.  Setting the 'logLevel' on line 25 of 
   settings.py to 10 will give you all the debug messages to help track down 
   the problem.
7. The current settings are the minimum to converge to a single posteriors peak
   and perform all three stages of ExoSOFT in a couple minutes.  Running it for 
   more samples by increasing the 'nSamples' parameter at the top or increasing 
   the 'nChains' and 'nMCMCcns' to higher matching values would be produce more 
   well sampled posteriors.  For example, 7 chains each of 5e7 were used to 
   produce the results in the ExoSOFT release paper, which took our computer 
   ~5hrs to complete.



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

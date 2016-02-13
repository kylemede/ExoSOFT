Packages installed on top of standard Fedora 22 libraries.
-gcc compiler for both c and c++ (although we only use c++...)
($sudo yum install gcc) 
($sudo yum install gcc-c++) 
epstopdf
($sudo yum install texlive-epstopdf)
-psutil
available at (https://pypi.python.org/pypi/psutil)
($sudo pip install psutil)
-scipy 
($sudo pip install scipy)
-pylab (usuall installed as a part of scipy or python, but this will ensure it is complete)
($sudo pip install pylab)
-matplotlib (installed with pylab, but again just to make sure)
($sudo pip install matplotlib)
-swig 
($sudo pip install swig)
-pyfits
($sudo pip install pyfits)
-astropy
($sudo pip install astropy)

NOTE: Need sudo for swig compiling --> '$ sudo make' or '$ sudo make clean' for cpp swig code!!!  I guess this will come into the setup.py area.
      NOT TRUE if the user has read/write privilages for the SMODT dirs/files!!!

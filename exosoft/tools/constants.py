import numpy as np

Grav = 6.67384e-11 #from physics.nist.gov
Gcgs = 6.67259e-8 #Grav in cgs converted from NIST value
pi = np.pi #=3.141592653589793
KGperMsun = 1.9884e30 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
daysPerYear = 365.2422 #For Gregorian calendar.  For Julian Calendar it is 365.25. In Astronomy we always use JD values for epochs, T and Tc. which to use...?? http://pumas.jpl.nasa.gov/files/04_21_97_1.pdf
secPerYear = 60*60*24*daysPerYear  #=31556926.080000006 (for Gregorian days/year)
MperAU = 149597870700.0 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
KGperMjupiter = 1.8983e27 #from http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
metersPerParsec = 3.08567758e16 #google
AUperParsec = 206264.806 #google
from __future__ import absolute_import
import numpy as np
##Note: for the moment, these are doubly defined until camel case/underscore inconsistancies are cleaned up
grav = Grav = 6.67384e-11 #from physics.nist.gov
grav_cgs = Gcgs = 6.67259e-8 #Grav in cgs converted from NIST value
pi = np.pi #=3.141592653589793
kg_per_msun = KGperMsun = 1.9884e30 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
days_per_year = daysPerYear = 365.2422 #For Gregorian calendar.  For Julian Calendar it is 365.25. In Astronomy we always use JD values for epochs, T and Tc. which to use...?? http://pumas.jpl.nasa.gov/files/04_21_97_1.pdf
sec_per_year = secPerYear = 60*60*24*daysPerYear  #=31556926.080000006 (for Gregorian days/year)
m_per_au = MperAU = 149597870700.0 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
kg_per_mjup = KGperMjupiter = 1.8983e27 #from http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
m_per_parsec = metersPerParsec = 3.08567758e16 #google
au_per_parsec = AUperParsec = 206264.806 #google
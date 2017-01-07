import numpy as np
import os
from PyAstronomy import pyasl
from astropy import constants as const

days_per_year = 365.2422
sec_per_year = 60*60*24*days_per_year

def calcOrbit(outDir='',outBaseName='mockdata_'):
    """
    This is a tool to produce a Keplerian orbit in RA,Dec,RV to verify exosoft.
    It just computes Kepler orbit positions and velocities then rotates
    them into the observed values in the plane of the sky.

    NOTE: the accuracy of the output values is to 'one part in 10^5',
          due to the simple center differencing used to calculate the velocities
          from the positions and times.

    Outputs are 2 files matching formats used in exosoft DI and RV.
    outBaseName+'DIdata.dat' has columns:
    #1. JD
    #2. RA (x) ["]
    #3. RA ERROR ["]
    #4. Dec (y) ["]
    #5. Dec ERROR ["]
    outBaseName+'RVdata.dat' has colunns:
    #1. JD
    #6. RV of primary (or secondary) rel to CofM [m/s]
    #7. RV ERROR [m/s]
    """
    quiet = False
    #Computer Directory
    if outDir=='':
        outDir='/mnt/HOME/MEGA/ExoSOFT-outputCopies/after-Dec25-2016/'#'/mnt/HOME/MEGA/MEGA-ExoSOFT-outputCopies-ForV2'#$$$$$$$$$$$$$$$$$$$$ MAKE SURE THIS IS SET TO MACH YOUR COMPUTER!!!
    #baseSaveDir = '/run/media/kmede/SharedData/Data/data_SMODT/'
    NumDataPointsOutRV = 25 #must be much less than 10000.  values between 10-500 are suitable.
    NumDataPointsOutDI = 10 #must be much less than 10000.  values between 10-500 are suitable.
    storePrimaryRVs = True
    percentError = 5 #error is set to a percentage of the median
    realizeErrors = True
    percentCoverage = 100.00 #percent of total orbit for data to span.  Over 100% is ok if you want overlapping data.

    ## System settings
    ##################
    M_secondary =  1.0*(const.M_jup.value/const.M_sun.value)    # [Solar masses]
    M_primary = 1.0 # [Solar masses]
    distance = 20.0 # [parsecs]
    #Orbital Elements
    TimeLastPeri =2450639.5  #JD  #2450817.5
    e =0.3 #
    period = 11.9 # [years]
    Omega =100.6*np.pi/180# Longitude of ascending node [deg]
    omega = 14.8*np.pi/180 # Argument of periastron [deg]
    i = 45.0*np.pi/180 # Inclination [deg]

    km_to_arcsec = 1.0/(const.au.value/1000.0)/distance # convert km to arcsecond
    massratio=M_primary/M_secondary
    mu = const.G.cgs.value*M_primary*(const.M_sun.value*1000.0)*(1.0 + 1.0/massratio) #gravitational parameter
    a = (mu*(period*sec_per_year)**2/4.0/np.pi**2)**(1.0/3.0) #in cm
    a_km = a/1.0e5 #to km
    a_AU = a_km/(const.au.value/1000.0) #to AU
    a2 = a_km/(massratio + 1.0)
    a1 = a_km - a2 # Semimajor axis of the low-mass component (in km)

    # print input orbital elements
    if quiet==False:
        print "\n\nOrbital Elements Used:\ne = "+str(e)
        print "period = "+str(period)+" Years"
        print "LongAN = "+str(Omega*180.0/np.pi)+" deg"
        print "ArgPeri = "+str(omega*180.0/np.pi)+" deg"
        print "a_total = "+str(a_AU)+" AU"
        print "inclination = "+str(i*180.0/np.pi)+" deg"
        print "Time of Last Periapsis = "+str(TimeLastPeri)+" JD"
        print "Mass 1 = "+str(M_primary)+" Msun"
        print "Mass 2 = "+str(M_secondary)+" Msun"
        print "Mass 2 = "+str(M_secondary*(const.M_sun.value/const.M_jup.value))+" Mjupiter"
        print "System distance = "+str(distance)+" PC, or "+str(1.0/(distance/1000.0))+' [mas]'
        #settings prints
        if storePrimaryRVs:
            print "Saving RVs of primary star relative to Center of Mass\n"
        else:
            print "Saving RVs of companion relative to Center of Mass\n"
        print 'Errors were calculated as '+str(percentError)+"% of the median value in the observables"
        if realizeErrors:
            print 'Data values were realized from the errors'
        else:
            print 'Data values are perfect with NO realization of the errors'
        print str(NumDataPointsOutRV)+" RV, and "+str(NumDataPointsOutDI)+" DI epochs will be calculated and stored"
        print 'The data will cover '+str(percentCoverage)+'% of the total orbit.\n'

    # Positions of both components in km relative to center of mass
    ke = pyasl.KeplerEllipse(a1, period, e=e, Omega=0.)
    NptsBIG = 10000
    t = (np.arange(NptsBIG) - 1)/(NptsBIG - 2.)*period

    ## update t to cover the percent of total orbit requested.
    #print 'len(t) before = '+str(len(t))
    t2=np.empty(( int((percentCoverage/100)*len(t)) ))
    if percentCoverage<=100.0:
        t2[0:t2.size]=t[0:t2.size]
    elif percentCoverage<200.0:
        t2[0:t.size]=t
        t2[t.size:]=t[0:( int(((percentCoverage-100.0)/100)*len(t)) )]+period
    t=t2
    #print 'len(t) after = '+str(len(t))
    pos_A = ke.xyzPos(t)
    pos_B = -pos_A/massratio

    # Velocities in km/s using centered differencing
    vel_A = (pos_A[2:] - pos_A[:-2])/(t[2] - t[0])/(86400*365.2422)
    pos_A = pos_A[1:-1]
    vel_B = (pos_B[2:] - pos_B[:-2])/(t[2] - t[0])/(86400*365.2422)
    pos_B = pos_B[1:-1]
    t = t[1:-1]

    # Construct rotation matrix (from wikipedia [http://en.wikipedia.org/wiki/Orbital_elements#Euler_angle_transformations])
    x1 = np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.cos(i)*np.sin(omega)
    x2 = np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.cos(i)*np.sin(omega)
    x3 = np.sin(i)*np.sin(omega)

    y1 = -np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(i)*np.cos(omega)
    y2 = -np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(i)*np.cos(omega)
    y3 = np.sin(i)*np.cos(omega)

    z1 = np.sin(i)*np.sin(Omega)
    z2 = -np.sin(i)*np.cos(Omega)
    z3 = np.cos(i)

    rotmat = np.asarray([[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]])

    # Rotate positions, velocities
    pos_A = np.dot(pos_A, rotmat)
    vel_A = np.dot(vel_A, rotmat)
    pos_B = np.dot(pos_B, rotmat)
    vel_B = np.dot(vel_B, rotmat)

    ## Randomly re-sample position, vel and time arrays to requested number of samples
    pos_Anew = []
    pos_Bnew = []
    vel_Anew = []
    vel_Bnew = []
    tnew = []
    i=0
    js=[]
    while i<NumDataPointsOutRV:
        j = np.random.randint(0,len(t))
        if j not in js:
            js.append(j)
            pos_Anew.append(pos_A[j])
            pos_Bnew.append(pos_B[j])
            vel_Anew.append(vel_A[j])
            vel_Bnew.append(vel_B[j])
            tnew.append(t[j])
            i+=1

    pos_A = np.array(pos_Anew)
    pos_B = np.array(pos_Bnew)
    vel_A = np.array(vel_Anew)
    vel_B = np.array(vel_Bnew)
    t=np.array(tnew)

    #make an ary for raw forms of the data.
    data = np.zeros((pos_A.shape[0], 8))
    data[:, 0] = t*2*np.pi/period #1. phase
    data[:, 1] = t # 2. time (years)
    data[:, 2] = pos_A[:, 0]*km_to_arcsec #3. x position of secondary (arcsec)
    data[:, 3] = pos_A[:, 1]*km_to_arcsec #4. y position of secondary (arcsec)
    data[:, 4] = vel_A[:, 2] #5. radial velocity of secondary (km/s)
    data[:, 5] = pos_B[:, 0]*km_to_arcsec #6. x position of primary (arcsec)
    data[:, 6] = pos_B[:, 1]*km_to_arcsec #7. y position of primary (arcsec)
    data[:, 7] = vel_B[:, 2] #8. radial velocity of primary (km/s)

    #update raw forms to initial NewBEAT versions
    data2 = np.zeros((pos_A.shape[0],5))
    data2[:,0] = data[:, 1]*days_per_year+TimeLastPeri #JD
    data2[:,1] = pos_A[:, 1]*km_to_arcsec - pos_B[:, 1]*km_to_arcsec #Ythi=Xplot=RA  separation between two bodies based on primary being at 0,0 ["]
    data2[:,2] = pos_A[:, 0]*km_to_arcsec - pos_B[:, 0]*km_to_arcsec #Xthi=Yplot=Dec  separation between two bodies based on primary being at 0,0 ["]
    data2[:,3] = vel_B[:, 2]*1000.0 # RV of primary compared to center of mass origin[ m/s]
    data2[:,4] = vel_A[:, 2]*1000.0 # RV of secondary compared to center of mass origin[ m/s]

    #calculate error and use it to realize the errors in the DI data if requested
    errorRA = np.median(np.abs(data2[:,1]))*(percentError/100.0)
    errorDec = np.median(np.abs(data2[:,2]))*(percentError/100.0)
    errorRA = np.max([errorRA,errorDec])
    errorDec = np.max([errorRA,errorDec])
    #print 'Using larger of two DI errors for both = '+str(errorDec)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,1]+=np.random.normal(0,errorRA)
            data2[i,2] += np.random.normal(0,errorDec)
    #calculate error and use it to realize the errors in the RV data if requested
    errorRVprimary = np.median(np.abs(data2[:,3]))*(percentError/100.0)
    errorRVsecondary = np.median(np.abs(data2[:,4]))*(percentError/100.0)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,3] += np.random.normal(0,errorRVprimary)
            data2[i,4] += np.random.normal(0,errorRVsecondary)

    #########################################################
    #load up data for into arys for exosoft DI, exosoft RV.
    #########################################################
    #dataDI2 has columns:
    #1. JD
    #2. RA (x) ["]
    #3. RA ERROR ["]
    #4. Dec (y) ["]
    #5. Dec ERROR ["]
    #dataRV has colunns:
    #1. JD
    #6. RV of primary (or secondary) rel to CofM [m/s]
    #7. RV ERROR [m/s]
    dataDI2 = np.empty((pos_A.shape[0],5))
    dataDI2[:,0] = data2[:, 0]#1. JD
    dataDI2[:,1] = data2[:,1]#2. RA (x) ["]
    dataDI2[:,2] = errorRA #3. RA ERROR ["]
    dataDI2[:,3] = data2[:,2]#4. Dec (y) ["]
    dataDI2[:,4] = errorDec#5. Dec ERROR ["]
    dataRV = np.empty((pos_A.shape[0],3))
    dataRV[:,0] = data2[:, 0]#1. JD
    if storePrimaryRVs:
        dataRV[:,1] = data2[:,3] #RV primary [m/s]
        dataRV[:,2] = errorRVprimary#RV primary error [m/s]
    else:
        dataRV[:,1] = data2[:,4]#RV secondary [m/s]
        dataRV[:,2] = errorRVsecondary#RV secondary error [m/s]

    if True:
        ## Randomly re-sample the DI data to half that of the RV to mimick the fact that there is usually more much RV data than DI data
        dataDI3=[]
        i=0
        js=[]
        while i<NumDataPointsOutDI:
            j = np.random.randint(0,len(dataDI2[:,0]))
            if j not in js:
                js.append(j)
                dataDI3.append(dataDI2[j,:])
                i+=1
        dataDI3 = np.array(dataDI3)
    else:
        dataDI3 = dataDI2
    #print "resulting RV data files have "+str(len(dataRV[:,0]))+" epochs"
    #print "resulting DI data files have "+str(len(dataDI3[:,0]))+" epochs"

    ##write files to disk
    if False:
        # raw all-in-one format, NOT for ExoSOFT use.
        np.savetxt(os.path.join(outDir,outBaseName+'.dat'), data, fmt="%.10g")
    if True:
        # 2 files for ExoSOFT use
        np.savetxt(os.path.join(outDir,outBaseName+'RVdata.dat'), dataRV, fmt="%.10g")
        np.savetxt(os.path.join(outDir,outBaseName+'DIdata.dat'), dataDI3, fmt="%.10g")
    if quiet==False:
        print '\nOutput data files written to:\n'+outDir

if __name__ == "__main__":
    calcOrbit()

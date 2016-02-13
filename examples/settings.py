#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
from priors import ePriorRatio,pPriorRatio,incPriorRatio,mass1PriorRatio,mass2PriorRatio,paraPriorRatio,chabrierPrior

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (240000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
# For MCMC mode this is the number of SA and ST chains.
'nChains' : (25,"Number MC/SA/ST of chains"),
# Number of MCMC chains to run in parallel. ONLY available in 'MCMC' mode. [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nMCMCcns' : (7,"Number MCMC of chains"),
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 30,
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# Run in Automatic mode? This will perform checks and select the stages to run automatically. [bool]
'autoMode' : (True, 'Run in Automatic mode?'),
# mode to run simulation in, choices {'MC','SA','ST','SAST','SASTMCMC,'MCMC'} [string]
# NOTE: 'ST' and 'MCMC' modes need a full list of parameters for startParams, else they fail!
#       'MCMC' also needs a full list of sigmas in startSigmas.
'stages' : 'MCMC',
# If in autoMode, how strict should the initialization (SA & ST) be? [string]
# choices ('loose','enough','tight')
'initCrit' : 'enough',
############################################
# Starting parameters and sigmas for MCMC  #
# Can be found with prior run in SAST mode #
############################################
# if unknown, set to False!! else [comma separated list of doubles]
'startParams' : False,
# if unknown, set to False!! else [comma separated list of doubles]
'startSigmas' : False,
# If better startParams are found during exosoft, push them and the newest sigmas into this file? [bool]
# NOTE: this can be helpful, but use caution if you do not wish to overwrite the values in here.
"pushToSettFiles":False,
}

directoriesDict = {
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/mnt/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "TEST-ArtificialJupiter-5percentError",
#*************************************************************************************************************************               
# full path to input astrometry data file. [string]
'DIdataFile': '/mnt/HOME/Dropbox/EclipseWorkspaceDB/ExoSOFT/examples/DIdata.dat',                
# full path to input radial velocity data file. [string]
'RVdataFile': '/mnt/HOME/Dropbox/EclipseWorkspaceDB/ExoSOFT/examples/RVdata.dat',
}

advancedSettingsDict = {
#NOTE: key max = 8characters, value+comment max = 68 characters, comment Max=47 it seems in testing.
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file during MC mode. [double]
'chiMAX' : (300.0,"Max reduced chiSquared during MC"),
# maximum allowed reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(5,'Max reduced chiSquared to enter ST.'),
# maximum allowed reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(1.0,'Max reduced chiSquared to enter MCMC.'),
#number of times to produce a summary log msg during a stage's progress [int]
'nSumry'  :10,
# make plot of posterior distributions? [bool]
'pltDists' :True,
# make plots of RV and DI/AM orbit fits [bool]
'pltOrbit' :True,
# Delete chain files after simulation is complete? [bool]
'delChains' :True,
# Delete combined data files after simulation is complete? [bool]
'delCombined' :False,
# run 'make' on C++/SWIG code to make sure it is up-to-date [bool]
'remake' :False,
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  
'CopyToDB' :False,
'dbFolder' : '/mnt/HOME/Dropbox/SMODT-outputCopies/',
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
'CalcBurn' :True,
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) [bool]
# (should already be handled by SimAnneal stage2 though...)?
'rmBurn' : (True,"Remove Burn-in?"),
# Calculate the Correlation lengths and number of effective points of each chain (must be more than 1 chain)? [bool]
# NOTE: CAUTION, can take a long time for long runs.  Still needs to be sped up somehow.
'calcCL' :False,
# number of samples to draw for simulated annealing stage [int] 
'nSAsamp' :(3000000,"Num SA samples"),
# Simulated Annealing starting temperature [double]
'strtTemp' : (50.0,"SA start temp."),
# Maximum unitless bias-corrected standard deviation allowed between best reduced chi squareds of SA results. [double]
'maxUstd': 0.02,
# Starting sigma size, % of parameter range, recommend [0.05,0.25].  [double]
# After first trial of SA and ST, take ST output and use here.
'strtSig' : (0.01,"start percent param range for SA"),
# Number of samples till temperature drop. [int]
# Allowed vals [1,nSAsamp), Ideal is ~50.
'tempInt'  : (50,"Num steps till temp drops in SA."),
# number of samples to draw for sigma tuning stage [int].
'nSTsamp' :(1000000,"Num ST samples"),
# number of steps per varying parameter until calculating the acceptance rate and tuning sigmas. [int]
# Allowed vals [1,nSTsamp), testing shows a value of ~200 works well.
'sigInt': (200,"Num steps/par till calc acc rate/tune sigs."),
# Maximum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMax' :(1.0,'Max ratio of params range,for step size.'),
# Minimum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMin' :(0.02,'Min ratio of params range,for step size.'),
# interval of accepted values between storing in output array (for SA,ST,MCMC, not MC) [int]
# Make sure to save enough that R~1.0 at max, posteriors look smooth, BUT not too much data is saved that you are just wasting disk space.
'saveInt' : (10,"Int between saving params, for all but MC."),
# Interval of saved values before write/dump the data to disk to avoid consuming too much RAM during long runs. They take 11MB/100000.
'dmpInt'   : 100000,
## NOTE: progress plots have no code yet, so MUST be False!!!
# Make plots of MCMC progress plots? [bool]  
'pltMCMCprog' :False,
# Make plots of Simulated Annealing progress plots? [bool]
'pltSAprog' :False,
# Calculate the Gelman-Rubin statistic? [bool]
'CalcGR' :True,
# How many times do you want the Gelman-Rubin statistic calculated [int]
'nGRcalc' :10,
#####################################
# Special Settings for the models ###
#####################################
# Operate in low eccenctricity mode? [bool]
# Then step through in sqrt(e)sin(omega) and sqrt(e)cos(omega) instead of e & omega directly
'lowEcc'   : (True,"low eccentricty stepping?"),
# fit to the primary's RV orbit [bool]
'fitPrime' : (False,"Fit primary's orbit?"),
# Are the RVs in the RVdata.dat for the Primary star? [bool]
'primeRVs' : (True,"RVs measured from Primary?"),
# Draw values for K directly, do NOT calculate it [bool]. Kills varying of Inclination.  Only possible in RV only mode.
'Kdirect'  : (False,'Vary K direct, do not calc it'),
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'TcStep' : (False,"Step in Tc not T?"),
# take the time of center transit (inferior conjunction) into account? [bool]
'TcEqualT' : (True,"Fix Tc=T?"),
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omegaPrv' : (0.0,"Custom fixed val added to RV omega in model"),
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omegaPdi' : (0.0,"Custom fixed val added to DI omega in model"),
# Is the data in the DIdata.dat in PA,SA format? else, it is in E,N (x,y) format [bool]
'pasa'     : (False,"Is astrometry data in PA,SA format?"),
}


rangesDict={
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
# NOTE: For DI only cases, use mass1 values as total mass and set mass2 values to zero.
'mass1MIN' : 0.2,
'mass1MAX' : 2.35,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.0001,
'mass2MAX' : 0.005,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 30.0,
'paraMAX' : 70.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 80.0,
'OmegaMAX' : 130.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.00000001,
'eMAX' : 0.15,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2449550,
'TMAX' : 2452500,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 8.0,
'PMAX' : 16.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 20,
'incMAX' : 65,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : -70,
'omegaMAX' : 170,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[-3],
'vMAXs' :[3],
}

priorsDict={
############################
#    System Information    #
# ONLY FOR GAUSSIAN PRIORS #
############################
#best estimate of primary's mass, and error [double][Msun]
'mass1Est' : (0.0,"Primary's estimated mass"),
'mass1Err' : (0.0,"Primary's estimated mass error"),
#best estimate of secondary's mass, and error [double][Msun]
'mass2Est' : (0.0,"Secondary's estimated mass"),
'mass2Err' : (0.0,"Secondary's estimated mass error"),
#best estimate of parallax, and error [double][mas]
'paraEst'  : (50,"Estimated parallax"),
'paraErr'  : (2.5,"Estimated parallax error"),
##################################
# Push prior functions into dict #
##################################
'ePrior'    :(True,'Use prior for eccentricity?',ePriorRatio),
'pPrior'    :(True,'Use prior for period?',pPriorRatio),
'incPrior'  :(True,'Use prior for inclination?',incPriorRatio),
'M1Prior':(True,'Use prior for m1?',mass1PriorRatio),
'M2Prior':(True,'Use prior for m2?',mass2PriorRatio),
'parPrior' :(True,'Use prior for parallax?',paraPriorRatio),
}

######################
# Merge All dicts#
######################
settingsDict = {}
for key in simpleSettingsDict:
    settingsDict[key]=simpleSettingsDict[key]
for key in directoriesDict:
    settingsDict[key]=directoriesDict[key]
for key in advancedSettingsDict:
    settingsDict[key]=advancedSettingsDict[key]
for key in rangesDict:
    settingsDict[key]=rangesDict[key]
for key in priorsDict:
    settingsDict[key]=priorsDict[key]

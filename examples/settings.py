#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (10000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
# For MCMC mode this is the number of SA and ST chains.
'nChains' : (3,"Number MC/SA/ST of chains"),
# Number of MCMC chains to run in parallel. ONLY available in 'MCMC' mode. [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nMCMCcns' : (3,"Number MCMC of chains"),
# Number of emcee walkers to run.  ONLY available in 'emcee' mode. [1,500] [int]
'n_wlkrs' :  (200,"Number of emcee walkers"),
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'IMPORTANTINFO'=25,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 25,
# data mode, choices {'RV','DI','3D'} [string]
'data_mode' : ('3D',"Data Mode (RV,DI,3D)"),
# Run in Automatic mode? This will perform checks and select the stages to run automatically. [bool]
'autoMode' : (False, 'Run in Automatic mode?'),
# mode to run simulation in, choices {'MC','SA','ST','SAST','SASTMCMC,'MCMC'} [string]
# NOTE: 'ST' and 'MCMC' modes need a full list of parameters for startParams, else they fail!
#       'MCMC' also needs a full list of sigmas in startSigmas.
'stages' : 'emcee',
# If in autoMode, how strict should the initialization (SA & ST) be? [string]
# choices ('loose','enough','tight')
'initCrit' : 'tight',
############################################
# Starting parameters and sigmas for MCMC  #
# Can be found with prior run in SAST mode #
############################################
# if unknown, set to False!! else [comma separated list of doubles]
'startParams' : [0.989258173441,0.000969607646338,49.4824152464,101.824621571,0.0570635830548,2450656.6178,2450656.6178,12.0122193117,44.5819933928,0.194526776635,5.22763017736,46.9083271473,8.91689222077,-0.0413666288362],
# if unknown, set to False!! else [comma separated list of doubles]
'startSigmas' : [0.0025495,4.99e-06,0.099,0.36,0.000101373206355,4.5,0.049,0.089,0.00098488578018,0.006],
}

directoriesDict = {
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/mnt/Data1/Todai_Work/Data/data_ExoSOFT',
# full path to input astrometry data file. [string]
'di_dataFile': '/mnt/HOME/MEGA/Dropbox/EclipseWorkspaceDB/ExoSOFT/examples/DIdata.dat',      
# full path to input radial velocity data file. [string]
'rv_dataFile': '/mnt/HOME/MEGA/Dropbox/EclipseWorkspaceDB/ExoSOFT/examples/RVdata.dat',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "TEST-ArtificialJupiter-5percentError-ExoSOFT2-tst-tst",
#*************************************************************************************************************************               
}

advancedSettingsDict = {
#NOTE: key max = 8characters, value+comment max = 68 characters, comment Max=47 it seems in testing.
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file during MC mode. [double]
'chiMAX' : (50.0,"Max reduced chiSquared during MC"),
# maximum allowed reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(10,'Max reduced chiSquared to enter ST.'),
# maximum allowed reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(2.0,'Max reduced chiSquared to enter MCMC.'),
# Start all MCMC chains from the same best found during the SAST stages? [bool]
'strtMCMCatBest':True,
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
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
'CalcBurn' :True,
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) [bool]
# (should already be handled by SimAnneal stage2 though...)?
'rmBurn' : (True,"Remove Burn-in?"),
# number of burn-in samples to remove from beginning of emcee walker chains. [int]
'n_emcee_burn' : 1000,
# Calculate the Correlation lengths and number of effective points of each chain (must be more than 1 chain)? [bool]
# NOTE: CAUTION, can take a long time for long runs.  Still needs to be sped up somehow.
'calcCL' :True,
# Calculate the integrated autocorrelation time for each parameter?
# Wraps the 'emcee.autocorr.integrated_time function. [bool]
'calcIAC':True,
# Calculate the Gelman-Rubin statistic? [bool]
'CalcGR' :True,
# How many times do you want the Gelman-Rubin statistic calculated [int]
'nGRcalc' :10,
# number of samples to draw for simulated annealing stage [int] 
'nSAsamp' :(50000,"Num SA samples"),
# Simulated Annealing starting temperature [double]
'strtTemp' : (100.0,"SA start temp."),
# Number of samples till temperature drop. [int]
# Allowed vals [1,nSAsamp), Ideal is ~50.
'tempInt'  : (50,"Num steps till temp drops in SA."),
# Maximum unitless bias-corrected standard deviation allowed between best reduced chi squareds of SA results. [double]
'maxUstd': 0.5,
# number of samples to draw for sigma tuning stage [int].
'nSTsamp' :(100000,"Num ST samples"),
# Starting sigma size, ratio of parameter range, recommend [0.05,0.25].  [double]
# After first trial of SA and ST, take ST output and use here.
'strtSig' : (0.01,"start percent param range for SA"),
# number of steps per varying parameter until calculating the acceptance rate and tuning sigmas. [int]
# Allowed vals [1,nSTsamp), testing shows a value of ~200 works well.
'sigInt': (200,"Num steps/par till calc acc rate/tune sigs."),
# Minimum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMin' :(0.001,'Min ratio of params range,for step size.'),
# Minimum and maximum allowed acceptance rates to tune parameter step sizes to meet.  [double,double]
'accRates' :([0.25,0.35],'[min,max] acceptance rates for ST/MCMC.'),
# interval of accepted values between storing in output array (for SA,ST,MCMC, not MC) [int]
# Make sure to save enough that R~1.0 at max, posteriors look smooth, BUT not too much data is saved that you are just wasting disk space.  ### rename this "thinning"??
'saveInt' : (10,"Int between saving params, for all but MC."),
# thin emcee samples? [bool]
'thin_emcee': True,
# emcee thinning rate. [int]
#EX if =10, then save every 10th sample.  #$$$$$$$$$$$ check this isn't greater or even close to nsamp/nwalkers
'thin_rate': 50,
# Interval of saved values before write/dump the data to disk to avoid consuming too much RAM during long runs. They take 11MB/100000.
'dmpInt'   : 100000,
#####################################
# Special Settings for the models ###
#####################################
# Operate in low eccenctricity mode? [bool]
# Then step through in sqrt(e)sin(omega) and sqrt(e)cos(omega) instead of e & omega directly
'low_ecc'   : (True,"low eccentricty stepping?"),
# Draw values for K directly, do NOT calculate it [bool]. Kills varying of Inclination.  Only possible in RV only mode.
'Kdirect'  : (False,'Vary K direct, do not calc it'),
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'vary_tc' : (False,"Step in Tc not T?"),
# take the time of center transit (inferior conjunction) into account? [bool]
'tc_equal_to' : (True,"Fix Tc=T?"),
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omega_offset_rv' : (0.0,"Custom fixed val added to RV omega in model"),
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omega_offset_di' : (0.0,"Custom fixed val added to DI omega in model"),
# Is the data in the DIdata.dat in PA,SA format? else, it is in E,N (x,y) format [bool]
'pasa'     : (False,"Is astrometry data in PA,SA format?"),
}


rangesDict={
###################################################
# Ranges for acceptable random number inputs ######
###################################################
## For Omega, e, T, inc and omega, None indicates to use default ranges.
## For Omega and omega, values can vary outside ranges, but are shifted by 
## +/-360 befire being stored.
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
# NOTE: For DI only cases, use mass1 values as total mass and set mass2 values to zero.
#       This will be done during start up otherwise.
'm1_min' : 0.0005,
'm1_max' : 2.55,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
# Note: 1Mj ~ 0.00095 Msun
'm2_min' : 0.00001,
'm2_max' : 0.005,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'para_min' : 1.0,
'para_max' : 100.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg] OR None
# Default: [0,360] if 'dataMode' is '3D', else [0,180].
'long_an_min' : None,
'long_an_max' : None,
# Minimum/Maximum allowed value for the Eccentricity [double] OR None
# Default: [0,0.98]
'ecc_min' : None,
'ecc_max' : None,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD] OR None
# Default: [earliestsEpoch-period,earliestEpoch]
't_min' : 2449000,
't_max' : 2453500,
# Minimum/Maximum allowed value for the Period [double][yrs]
'p_min' : 1.0,
'p_max' : 50.0,
# Minimum/Maximum allowed value for the Inclination [double][deg] OR None
# Default: [0,180].  [0,90] is another popular choice.
'inc_min' : 1,
'inc_max' : 90,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg] OR None
# Default: [0,360]
'arg_peri_min' : -50,
'arg_peri_max' : 90,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# NONE   THEY ARE GONE  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'offset_mins' :[-3],
'offset_maxs' :[3],
}

priorsDict={
############################
#    System Information    #
# ONLY FOR GAUSSIAN PRIORS #
############################
#best estimate of primary's mass, and error [double][Msun]
'm1_est' : (0.0,"Primary's estimated mass"),
'm1_err' : (0.0,"Primary's estimated mass error"),
#best estimate of secondary's mass, and error [double][Msun]
'm2_est' : (0.0,"Secondary's estimated mass"),
'm2_err' : (0.0,"Secondary's estimated mass error"),
#best estimate of parallax, and error [double][mas]
'para_est'  : (50,"Estimated parallax"),
'para_err'  : (2.5,"Estimated parallax error"),
##################################
# Push prior functions into dict #
##################################
# For ALL: False indicates a flat prior. True indicates to use defaults.
'ecc_prior'    :(True,'Use prior for eccentricity?'),
'p_prior'    :(True,'Use prior for period?'),
# For the inclination prior, use strings or booleans to inducate the specific 
# function to use.  Either sin(i), cos(i) or flat.  True indicates sin(i).
'inc_prior'  :('cos',"inclination prior ['cos','sin',True or False]"),
# For m1 and m2: use strings to indicate specific prior function.
# m1 default is 'PDMF', m2 default is 'CMF'.
'm1_prior':(True,"m1 prior ['PDMF', 'IMF', True or False]"),
'm2_prior':(True,"m2 prior ['PDMF', 'IMF', 'CMF', True or False]"),
'para_prior' :(True,'Use prior for parallax?'),
}

######################
# Merge All dicts#
######################
settings = {}
for key in simpleSettingsDict:
    settings[key]=simpleSettingsDict[key]
for key in directoriesDict:
    settings[key]=directoriesDict[key]
for key in advancedSettingsDict:
    settings[key]=advancedSettingsDict[key]
for key in rangesDict:
    settings[key]=rangesDict[key]
for key in priorsDict:
    settings[key]=priorsDict[key]

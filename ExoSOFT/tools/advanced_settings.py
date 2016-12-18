#######################
## Advanced settings ##
#######################
advanced_settings_dict = {
# full path to input astrometry data file. [string]
'di_dataFile': './DIdata.dat',
# full path to input radial velocity data file. [string]
'rv_dataFile': './RVdata.dat',

# Run in Automatic mode? This will perform checks and select the stages to run automatically. [bool]
'autoMode' : False,
# If in autoMode, how strict should the initialization (SA & ST) be? [string]
# choices ('loose','enough','tight')
'initCrit' : 'tight',

# Start all MCMC chains from the same best found during the SAST stages? [bool]
'strtMCMCatBest': True,
#number of times to produce a summary log msg during a stage's progress [int]
'nSumry'  : 10,

# How many times do you want the Gelman-Rubin statistic calculated [int]
'nGRcalc' : 10,

# number of samples to draw for simulated annealing stage [int] 
'nSAsamp' : 50000,
# Simulated Annealing starting temperature [double]
'strtTemp' :  100.0,
# Number of samples till temperature drop. [int]
# Allowed vals [1,nSAsamp), Ideal is ~50.
'tempInt'  : 50,
# Maximum unitless bias-corrected standard deviation allowed between best reduced chi squareds of SA results. [double]
'maxUstd': 0.5,
# number of samples to draw for sigma tuning stage [int].
'nSTsamp' : 100000,
# Starting sigma size, ratio of parameter range, recommend [0.05,0.25].  [double]
# After first trial of SA and ST, take ST output and use here.
'strtSig' : 0.01,
# number of steps per varying parameter until calculating the acceptance rate and tuning sigmas. [int]
# Allowed vals [1,nSTsamp), testing shows a value of ~200 works well.
'sigInt': 200,
# Minimum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMin' : 0.001,
# Minimum and maximum allowed acceptance rates to tune parameter step sizes to meet.  [double,double]
'accRates' : [0.25,0.35],
# Interval of saved values before write/dump the data to disk to avoid consuming too much RAM during long runs. They take 11MB/100000.
'dmpInt'   : 100000,

#####################################
# Special Settings for the models ###
#####################################
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'vary_tc' : False,
# take the time of center transit (inferior conjunction) into account? [bool]
'tc_equal_to' : True,
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omega_offset_rv' : 0.0,
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omega_offset_di' : 0.0,

# Advanced settings
KMAX: 0,
KMIN: 0,
# Draw values for K directly, do NOT calculate it [bool]. Kills varying of Inclination.  Only possible in RV only mode.
'Kdirect'  : False,
}


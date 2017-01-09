from __future__ import absolute_import
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from astropy import constants as const
import emcee
import yaml
import os
import corner
from six.moves import range
import KMlogger
from ExoSOFT import tools

log = KMlogger.getLogger('main',lvl=20,addFH=False)

def main(settings_in,advanced_settings_in, priors_in):
    """
    An example of using ExoSOFTmodel with emcee (https://github.com/dfm/emcee).
    This will show:
    1. How to intantiate the necessary objects with their required input 
    parameters
    2. How prepare the inputs to emcee
    3. Then run emcee and finally make some simple plots of the results.
    """
    ## Simple boolean flags to control the basic steps of this script.
    # Just do a quick test run? 
    # Else, run a longer one to give smooth posteriors
    quick = True
    # What plots should be made?
    show_burnin = False
    show_chains = True
    show_posteriors = True
    show_corner = False
    
    sd = tools.startup(settings_in,advanced_settings_in,priors_in,rePlot=True)
    log = KMlogger.getLogger('main',lvl=sd['logLevel'],addFH=False)
    
    ## Instantiate main objects/classes: 
    #  ExoSOFTpriors, ExoSOFTdata and ExoSOFTparams.  
    #  These are all instantiated as member variables of ExoSOFTmodel class.
    Model = tools.ExoSOFTmodel(sd)
    
    ## define a set of starting parameters
    # The user can use any reasonable guess here.  For the 5% Jupiter analogue 
    # used in this example, we will use expected values for simplicity. 
    m2 = const.M_jup.value/const.M_sun.value
    sqrte_sinomega = np.sqrt(0.048)*np.sin((np.pi/180.0)*14.8)
    sqrte_cosomega = np.sqrt(0.048)*np.cos((np.pi/180.0)*14.8)
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    start_params = [1.0,m2,50,100.6,sqrte_sinomega,2450639.0,11.9,45.0,sqrte_cosomega,0.1]
    
    ncpu = multiprocessing.cpu_count()
    ndim = len(sd['range_maxs']) # number of parameters in the model 
    if quick:
        nwalkers = 50 # number of MCMC walkers 
        nburn = 50 # "burn-in" to stabilize chains 
        nsteps = 1000 # number of MCMC steps to take 
    else:
        nwalkers = 500 
        nburn = 500 
        nsteps = 50000 
    ## NOTE: starting_guesses array must be a numpy array with dtype=np.dtype('d').
    starting_guesses = tools.make_starting_params(start_params,nwalkers,scale=0.01)
    
    ## Call emcee to explore the parameter space
    sampler = emcee.EnsembleSampler(nwalkers, ndim, tools.ln_posterior, 
                                    args=[Model], threads=ncpu)
    sampler.run_mcmc(starting_guesses, nsteps)
    
    # chain is of shape (nwalkers, nsteps, ndim)
    # discard burn-in points and reshape
    trace = sampler.chain[:, nburn:, :] 
    trace = trace.reshape(-1, ndim)

    labels = ['m2', 'period', 'inclination']
    
    ## Show walkers during burn-in 
    if show_burnin:
        fig = plt.figure(figsize=(10,5))
        j=0
        for i, chain in enumerate(sampler.chain[:, :nburn, :].T): 
            if i in [1,6,7]:
                plt.subplot(3, 1, j+1)
                if i==1:
                    # convert m2 units to Mjup instead of Msun
                    chain*= const.M_sun.value/const.M_jup.value
                plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
                plt.ylabel(labels[j])
                j+=1
        plt.show()
    
    ## Show walkers after burn-in
    if show_chains:
        fig = plt.figure(figsize=(10,5))
        j=0
        for i, chain in enumerate(sampler.chain[:, nburn:, :].T):
            if i in [1,6,7]:
                plt.subplot(3, 1, j+1)
                if i==1:
                    # convert m2 units to Mjup instead of Msun
                    chain*= const.M_sun.value/const.M_jup.value
                plt.plot(chain, drawstyle='steps', color='k', alpha=0.2)
                log.info('\nparameter # '+str(i))
                log.info('mean = '+str(np.mean(chain)))
                log.info('median = '+str(np.median(chain)))
                log.info('variance = '+str(np.var(chain)))
                #print('chain ',repr(chain))
                plt.ylabel(labels[j])
                j+=1
        plt.show()
        
    ## inspect posteriors
    if show_posteriors:
        fig = plt.figure(figsize=(12,3))
        j=0
        for i in range(ndim):
            if i in [1,6,7]:
                plt.subplot(1,3,j+1)
                true_val = start_params[i]
                trace_use = trace[:,i]
                if i==1:
                    # convert m2 units to Mjup instead of Msun
                    trace_use*= const.M_sun.value/const.M_jup.value
                    true_val*= const.M_sun.value/const.M_jup.value
                plt.hist(trace_use, 100, color="k", histtype="step")
                yl = plt.ylim()
                plt.vlines(true_val, yl[0], yl[1], color='blue', lw=3, alpha=0.25, label='true')
                plt.title("{}".format(labels[j]))
                plt.legend()
                j+=1
        plt.show()
    
    ## make a corner plot
    if show_corner:
        fig = corner.corner(trace, labels=labels, quantiles=[0.16, 0.5, 0.84], truths=start_params)
        plt.show()
    
    log.critical("\nFinished emcee example :D \n")

if __name__ == '__main__':
    ## Load in any available settings, advanced settings and priors available in the CWD.
    ## Then push the path to those available to dictionary keys for later use.
    settings_in = None
    priors_in = None
    advanced_settings_in = None
    pf = None
    if os.path.exists('./settings.yaml'):
        f = open('./settings.yaml','r')
        settings_in = yaml.load(f)
        f.close()
        sf = os.path.join(os.getcwd(),'settings.yaml')
        settings_in['settings_in_path'] = sf
    if os.path.exists('./advanced_settings.yaml'):
        f = open('./advanced_settings.yaml','r')
        advanced_settings_in = yaml.load(f)
        f.close()
        asf = os.path.join(os.getcwd(),'advanced_settings.yaml')
        advanced_settings_in['advanced_settings_in_path'] = asf
    if os.path.exists('./priors.py'):
        from priors import ExoSOFTpriors as priors_in
        pf = os.path.join(os.getcwd(),'priors.py')
        if settings_in != None:
            settings_in['priors_in_path'] = pf
    main(settings_in,advanced_settings_in, priors_in)

#EOF
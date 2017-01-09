from __future__ import absolute_import
import numpy as np
#import multiprocessing
#import matplotlib.pyplot as plt
from astropy import constants as const
from six.moves import range
import KMlogger
from ExoSOFT import tools

log = KMlogger.getLogger('main',lvl=20,addFH=False)

def main():
    """
    A simple example of how to instantiate the objects used by ExoSOFTmodel 
    with their required input parameters.  This example uses the 5% Jupiter 
    analogue discussed in the release paper and calculates the chi squared
    of the expected parameters for the simulated data.
    """
    ## load up default settings dictionary
    sd = tools.load_settings('./settings.py')
    
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
    start_params = [1.0,m2,50,100.6,sqrte_sinomega,2450639.0,11.9,45.0,sqrte_cosomega,0.01]
    ## NOTE: this array must be a numpy array with dtype=np.dtype('d').
    start_params = np.array(start_params,dtype=np.dtype('d'))
    
    ## calculate the log posterior for the provided data and input parameters
    ln_post = tools.ln_posterior(start_params, Model)
    
    ## print a few basic results of the fit
    log.info('chi_squared_3d '+str(Model.chi_squared_3d))
    log.info('reduced chi_squared_3d '+str(Model.chi_squared_3d/35.0))
    log.info('prior '+str(Model.prior))
    log.info('ln_post '+str(ln_post))

if __name__ == '__main__':
    main()
    
#EOF
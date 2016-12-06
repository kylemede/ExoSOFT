#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import sys

## import from modules in ExoSOFT ##
from . import constants as const
#import constants as const

class ExoSOFTpriors(object):
    """
    Choices:
    inclination: (None,True,'sin','cos') True indicates default of 'sin'
    eccentricity: (None,True) True indicates default of '2e'
    period: (None,True) True indicates default of 'log'
    mass of primary: (None,True,'IMF','PDMF') True indicates default of 'PDMF'
    mass of companion: (None,True, 'IMF','PDMF','CMF') True indicates default of 'CMF'
    parallax: (None,True) True indicates default of 'gauss'
    """
    def __init__(self, ecc_prior=True, p_prior=True, inc_prior=True, 
                 m1_prior=True, m2_prior=True, para_prior=True, inc_min=0.0,
                 inc_max=180.0, p_min=0.0, p_max=300.0, para_est=0, 
                 para_err=0, m1_est=0, m1_err=0, m2_est=0, m2_err=0):   
        ## choices   
        # choices:2e  
        self.e_prior = ecc_prior
        # choices:log
        self.p_prior = p_prior
        # choices:sin, cos
        self.inc_prior = inc_prior
        # choices:IMF, PDMF
        self.m1_prior = m1_prior
        # choices:CMF, IMF, PDMF
        self.m2_prior = m2_prior
        # choices:gauss
        self.para_prior = para_prior
        ## values necessary to calculate priors
        self.inc_min = inc_min
        self.inc_max = inc_max
        self.p_min = p_min
        self.p_max = p_max
        self.para_est = para_est
        self.para_err = para_err
        self.m1_est = m1_est
        self.m1_err = m1_err
        self.m2_est = m2_est
        self.m2_err = m2_err
        
    def test_priors(self, pars_prop, pars_last):
        """Just for testing no errors occur while trying to caculate the
        priors ratio.  Nothing returned. """
        val = self.combinedPriors(pars_prop, pars_last)
        if False:
            print('test priors ratio value = '+str(val))
                
    def priors_ratio(self, pars_prop, pars_last):
        """
        Calculates the priors ratio.
        
        Input arrays need the first elements to be:
        [m1, m2, parallax, long_an, e, to, tc, p, inc, arg_peri]
        This can be the full 'model_in_pars' or a truncated version with 
        just these parameters.
        """
        try:
            priorsRatio = self.priors(pars_prop) / self.priors(pars_last)
            return priorsRatio
        except:
            print("An error occured while trying to calculate the priors ratio.")
            sys.exit(0)
        
    def priors(self, pars):
        """
        Combined priors for a single step in the chain.  This can be called 
        to form the priors ratio OR when accounting for the priors when 
        producing the posteriors during post-processing.
        
        Input array needs the first elements to be:
        [m1, m2, parallax, long_an, e, to, tc, p, inc, arg_peri]
        This can be the full 'model_in_pars', 'stored_pars' or a truncated 
        version with just these parameters.
        """
        # model_in_params: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,arg_peri_di,arg_peri_rv,a_tot_au,K]
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        comboPriors = 1.0
        try:
            if self.e_prior:
                #print 'ePrior before'#$$$$$$$$$$
                comboPriors*=self.ecc_prior_fn(pars[4])
                #print('e prior = ',self.ecc_prior_fn(pars[4]))
            if self.p_prior:
                #print 'pPrior before'#$$$$$$$$$$
                comboPriors*=self.p_prior_fn(pars[7])
                #print('p prior = ',self.p_prior_fn(pars[7]))
            if self.inc_prior:
                #print 'incPrior before'#$$$$$$$$$$
                comboPriors*=self.inc_prior_fn(pars[8])
                #print('inc prior = ',self.inc_prior_fn(pars[8]))
            if self.m1_prior:
                #print 'M1Prior before'#$$$$$$$$$$
                comboPriors*=self.m1_prior_fn(pars[0])
                #print('m1 prior = ',self.m1_prior_fn(pars[0]))
                #print 'M1Prior after'#$$$$$$$$$$
            if self.m2_prior:
                #print 'M2Prior before'#$$$$$$$$$$
                comboPriors*=self.m2_prior_fn(pars[1],pars[0])
                #print('m2 prior = ',self.m2_prior_fn(pars[1],pars[0]))
                #print 'M2Prior after'#$$$$$$$$$$
            if self.para_prior:
                #print 'parPrior before'#$$$$$$$$$$
                comboPriors*=self.para_prior_fn(pars[2])
                #print('para prior = ',self.para_prior_fn(pars[2]))
                #print 'parPrior after'#$$$$$$$$$$
            return comboPriors
        except:
            print("An error occured while trying to calculate the combined sigle priors.")
            sys.exit(0)
        
    ##########################################################################
    ####### Not to those who wish to write their own prior functions #########
    # Only change the code and not the name of the functions or their inputs.#  
    ##########################################################################          
    def ecc_prior_fn(self, ecc):
        ## UPGRADE TO INCLUDE FURTHER STRINGS TO INDICATE ECC PRIOR OPTIONS FROM KEPLER AND OTHER SURVEY RESULTS #$$$$$$$$$$$$$$$$$$$$
        ret = 1.0
        if ecc!=0:
            if (self.e_prior == True) or (self.e_prior == '2e'):
                if (self.p_min*const.daysPerYear)>1000.0:
                    ret = 2.0*ecc
        return ret
                  
    def p_prior_fn(self, p):
        ret = 1.0
        if p!=0.0:
            if (self.p_prior == True) or (self.p_prior == 'log'):
                ret = 1.0 / ( p * np.log( self.p_max / self.p_min ) )
        return ret
        
    def inc_prior_fn(self, inc):
        ret = 1.0
        if inc not in [0.0,90.0,180.0]:
            mn = np.radians(self.inc_min)
            mx = np.radians(self.inc_max)
            inc_rad = np.radians(inc)
            if (self.inc_prior == True) or (self.inc_prior == 'sin'):
                ret = np.sin(inc_rad) / np.abs(np.cos(mn)-np.cos(mx))
            elif self.inc_prior == 'cos':
                ret =  np.cos(inc_rad) / np.abs(np.cos(mn)-np.cos(mx))
                
        return ret
        
    def m1_prior_fn(self, mass):
        ret = 1.0
        if mass!=0:
            # First check if a gauss prob should be also calculated.
            if 0 not in [self.m1_est,self.m1_err]:
                ret*=self.gaussian(mass, self.m1_est, self.m1_err)
            # Then caculate requested mass function prior.
            if (self.m1_prior == True) or (self.m1_prior == "PDMF"):
                ret*=self.pdmf_prior(mass)
            elif self.m1_prior == "IMF":
                ret*=self.imf_prior(mass)
                
        return ret
        
    def m2_prior_fn(self, m2, m1):
        ret = 1.0
        if 0.0 not in [m2,m1]:
            # First check if a gauss prob should be also calculated.
            if 0 not in [self.m2_est,self.m2_err]:
                ret*=self.gaussian(m2, self.m2_est, self.m2_err)
            # Then caculate requested mass function prior.
            if (self.m2_prior == True) or (self.m2_prior == "CMF"):
                ret*=self.cmf_prior(m2, m1)
            elif self.m2_prior == "PDMF":
                ret*=self.pdmf_prior(m2)
            elif self.m2_prior == "IMF":
                ret*=self.imf_prior(m2)
                
        return ret
            
    def para_prior_fn(self, para):
        ret = 1.0
        if para!=0.0:
            ## this needs parallax in arc sec, but standard in ExoSOFT is mas, so /1000
            ret = 1.0/((para/1000.0)**4.0)
            #print('before gauss '+repr(ret))
            #print('para 1/**4 ',1.0/(para**4.0))
            if 0 not in [self.para_est,self.para_err]:
                ## a Gaussian prior centered on hipparcos and width of 
                #  hipparcos estimated error
                
                ret*=self.gaussian(para, self.para_est, self.para_err)
                #print('para gauss ',self.gaussian(para, self.para_est, self.para_err))
        return ret
            
    def imf_prior(self, m):
        """From table 1 of Chabrier et al. 2003"""
        a = 0.068618528140713786
        b = 1.1023729087095586
        c = 0.9521999999999998
        # calculate IMF specific to mass of object
        if m<1.0:
            d = (a/m) * np.exp( (-(np.log10(m)+b)**2) / c )
        else:
            d = 0.019239245548314052*(m**(-2.3))
        return d
    
    def pdmf_prior(self, m):
        """From table 1 of Chabrier et al. 2003"""
        a = 0.068618528140713786
        b = 1.1023729087095586
        c = 0.9521999999999998
        # calculate PDMF specific to mass of object
        if m<1.0:
            d = (a/m)*np.exp((-(np.log10(m)+b)**2)/c)
        elif m<3.47:
            d = 0.019108957203743077*(m**(-5.37))
        elif m<18.20:
            d = 0.0065144172285487769*(m**(-4.53))
        else:
            d = 0.00010857362047581295*(m**(-3.11))
        return d
    
    def cmf_prior(self, m2, m1):
        """from equation 8 of Metchev & Hillenbrand 2009"""
        beta = -0.39
        d = (m2**(beta)) * (m1**(-1.0*beta-1.0))
        return d
    
    def gaussian(self,x,mu,sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

#END OF FILE
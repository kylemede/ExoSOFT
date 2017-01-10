#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from scipy import stats
#from astropy import constants as const
import sys
#############################################################################
##  Do not import any non-standard modules here, only use those universaly ##
##  available to minimize any errors that could occur.                     ##
#############################################################################

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
                 m1_prior=True, m2_prior=True, para_prior=True,
                 para_est=0, para_err=0, m1_est=0, m1_err=0, m2_est=0, m2_err=0,
                 ecc_beta_a=0.867, ecc_beta_b=3.03,
                 ecc_J08_sig=0.3, ecc_Rexp_lamda=5.12, ecc_Rexp_a=0.781,
                 ecc_Rexp_sig=0.272, ecc_ST08_a=4.33, ecc_ST08_k=0.2431,p_gamma=-0.7,
                 mins_ary=[],maxs_ary=[]):
        ## check min and max range arrays
        if (len(mins_ary)>1) and (len(maxs_ary)>1):
            self.mins_ary = mins_ary
            self.maxs_ary = maxs_ary
        else:
            raise IOError('\n\n No min/max ranges were provided to the priors object!!')
        ## push in two manual constants
        self.days_per_year = 365.2422
        self.sec_per_year = 60*60*24*self.days_per_year
        ## choices
        ## choices:'2e', 'ST08','J08', 'RayExp', 'beta', 'uniform'. Default is 'beta'.
        self.e_prior = ecc_prior
        ## `best-fit' alpha and beta from Kipping+2013
        self.ecc_beta = stats.beta(ecc_beta_a,ecc_beta_b) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html
        ## 'best-fit' for the basic Rayleigh from Juric+2008
        self.ecc_J08_sig = ecc_J08_sig      ### Put this into the advanced settings ???
        ## `best-fit' alpha, lambda and sig from Kipping+2013
        self.ecc_Rexp_lamda = ecc_Rexp_lamda  ### Put this into the advanced settings ???
        self.ecc_Rexp_a = ecc_Rexp_a     ### Put this into the advanced settings ???
        self.ecc_Rexp_sig = ecc_Rexp_sig ### Put this into the advanced settings ???
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rayleigh.html#scipy.stats.rayleigh
        self.ecc_R = stats.rayleigh() 
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.expon.html
        self.ecc_exp = stats.expon()  
        ## `best-fit' for the Shen&Turner 2008 pdf
        self.ecc_ST08_a = ecc_ST08_a     ### Put this into the advanced settings ???
        self.ecc_ST08_k = ecc_ST08_k     ### Put this into the advanced settings ???

        #self.ecc_norm = stats.norm
        #self.ecc_norm_mean = ##
        #self.ecc_norm_sig = ##
        #self.ecc_norm.pdf(ecc,loc=self.ecc_norm_mean, scale=self.ecc_norm_sig)

        ## choices: power-law, Jeffrey's
        self.p_prior = p_prior
        self.p_gamma = p_gamma
        ## choices:sin, cos
        self.inc_prior = inc_prior
        ## choices:IMF, PDMF
        self.m1_prior = m1_prior
        ## choices:CMF, IMF, PDMF
        self.m2_prior = m2_prior
        ## choices:gauss
        self.para_prior = para_prior
        ## values necessary to calculate priors
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

    def priors_ratio(self, pars_prop, pars_last,ranges=False):
        """
        Calculates the priors ratio.

        Input arrays need the first elements to be:
        [m1, m2, parallax, long_an, e, to, tc, p, inc, arg_peri]
        This can be the full 'model_in_pars' or a truncated version with
        just these parameters.
        """
        try:
            priorsRatio = self.priors(pars_prop,ranges=ranges) / self.priors(pars_last,ranges=ranges)
            return priorsRatio
        except:
            print("An error occured while trying to calculate the priors ratio.")
            sys.exit(0)

    def priors(self, pars, ranges=False):
        """
        Combined priors for a single step in the chain.  This can be called
        to form the priors ratio OR when accounting for the priors when
        producing the posteriors during post-processing.

        Input array needs the first elements to be:
        [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        This can be the full 'stored_pars' parameters.
        """
        # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
        comboPriors = 1.0
        try:
            #print("\nEach prior was:\n")

            #print('M1Prior before')#$$$$$$$$$$
            comboPriors*=self.m1_prior_fn(pars[0])
            #print('m1 prior = ',self.m1_prior_fn(pars[0]))

            #print('M2Prior before')#$$$$$$$$$$
            comboPriors*=self.m2_prior_fn(pars[1],pars[0])
            #print('m2 prior = ',self.m2_prior_fn(pars[1],pars[0]))

            #print('parPrior before')#$$$$$$$$$$
            comboPriors*=self.para_prior_fn(pars[2])
            #print('para prior = ',self.para_prior_fn(pars[2]))

            ## for long_an
            comboPriors*=self.uniform_fn(pars[3],3)
            #print('long_an prior = ',self.uniform_fn(pars[3],3))

            #print('ePrior before')#$$$$$$$$$$
            comboPriors*=self.ecc_prior_fn(pars[4])
            #print('e prior = ',self.ecc_prior_fn(pars[4]))

            ## for To or Tc
            comboPriors*=self.uniform_fn(pars[5],5)
            #print('To/Tc prior = ',self.uniform_fn(pars[5],5))

            ## up to here the indices are the same for ranges and stored version of parameter arrays
            if ranges:
                #print('pPrior before')#$$$$$$$$$$
                comboPriors*=self.p_prior_fn(pars[6])
                #print('p prior = ',self.p_prior_fn(pars[6]))

                #print('incPrior before')#$$$$$$$$$$
                comboPriors*=self.inc_prior_fn(pars[7])
                #print('inc prior = ',self.inc_prior_fn(pars[7]))

                ## for arg_peri
                comboPriors*=self.uniform_fn(pars[8],8)
                #print('arg_peri prior = ',self.uniform_fn(pars[8],8))
                ## for the RV velocity offsets
                #print('\nlen(pars)', len(pars)-9)
                for i in range(len(pars)-9):
                    #print('V '+str(i))
                    comboPriors*=self.uniform_fn(pars[9+i],9+i)
                    #print('V '+str(i)+' prior = ',self.uniform_fn(pars[9+i],9+i))
            else:
                #print('pPrior before')#$$$$$$$$$$
                comboPriors*=self.p_prior_fn(pars[7])
                #print('p prior = ',self.p_prior_fn(pars[7]))

                #print('incPrior before')#$$$$$$$$$$
                comboPriors*=self.inc_prior_fn(pars[8])
                #print('inc prior = ',self.inc_prior_fn(pars[8]))

                ## for arg_peri
                comboPriors*=self.uniform_fn(pars[9],8)
                #print('arg_peri prior = ',self.uniform_fn(pars[9],8))
                ## for the RV velocity offsets
                #print('\nlen(pars)', len(pars)-13)
                for i in range(len(pars)-13):
                    comboPriors*=self.uniform_fn(pars[13+i],9+i)
                    #print('V '+str(i)+' prior = ',self.uniform_fn(pars[13+i],9+i))
            #print('returning comboPriors ',comboPriors)
            return comboPriors
        except:
            #print("An error occured while trying to calculate the combined sigle priors.")
            sys.exit(0)

    ##########################################################################
    ####### Not to those who wish to write their own prior functions #########
    # Only change the code and not the name of the functions or their inputs.#
    ##########################################################################
    def ecc_prior_fn(self, ecc):
        ## UPGRADE TO INCLUDE FURTHER STRINGS TO INDICATE ECC PRIOR OPTIONS FROM KEPLER AND OTHER SURVEY RESULTS #$$$$$$$$$$$$$$$$$$$$
        ret = 1.0
        if (self.e_prior==False) or (self.e_prior=="uniform"):
            ret = self.uniform_fn(ecc,4)
        else:
            if ecc!=0:
                if (self.e_prior == True) or (self.e_prior=='beta'):
                    ret = self.ecc_beta.pdf(ecc)
                elif self.e_prior == '2e':
                    if (self.mins_ary[6]*self.days_per_year)>1000.0:
                        ret = 2.0*ecc
                elif self.e_prior=='STO8':
                    ret =(1.0/self.ecc_ST08_k)*(1.0/(1.0+ecc)**self.ecc_ST08_a)-(ecc/2.0**self.ecc_ST08_a)
                elif self.e_prior=='J08':
                    ret = self.ecc_R.pdf(ecc,scale=self.ecc_R_sig)
                elif self.e_prior=='RayExp':
                    A = self.ecc_Rexp_a*self.ecc_exp.pdf(ecc, scale=1.0/self.ecc_Rexp_lamda)
                    B = (1.0-self.ecc_Rexp_a)*self.ecc_R.pdf(ecc,scale=self.ecc_Rexp_sig)/self.ecc_Rexp_sig
                    ret = A+B
        if ret==0: ret=-np.inf
        return ret

    def p_prior_fn(self, p):
        ret = 1.0
        if (self.p_prior==False) or (self.p_prior=="uniform"):
            ret = self.uniform_fn(p,6)
        else:
            if p!=0.0:
                if (self.p_prior == True) or (self.p_prior == 'jeffrey'):
                    # A Jeffrey's prior based on Gregory 2005.
                    ret = 1.0 / ( p * np.log( self.maxs_ary[6] / self.mins_ary[6] ) )
                elif self.p_prior == 'power-law':
                    ret = p **(self.p_gamma)
        if ret==0: ret=-np.inf
        return ret

    def inc_prior_fn(self, inc):
        ret = 1.0
        if (self.inc_prior==False) or (self.inc_prior=="uniform"):
            ret = self.uniform_fn(inc,7)
        else:
            if inc not in [0.0,90.0,180.0]:
                mn = np.radians(self.mins_ary[7])
                mx = np.radians(self.maxs_ary[7])
                inc_rad = np.radians(inc)
                if (self.inc_prior == True) or (self.inc_prior == 'sin'):
                    ret = np.sin(inc_rad) / np.abs(np.cos(mn)-np.cos(mx))
                elif self.inc_prior == 'cos':
                    ret =  np.cos(inc_rad) / np.abs(np.cos(mn)-np.cos(mx))
        if ret==0: ret=-np.inf
        return ret

    def m1_prior_fn(self, mass):
        ret = 1.0
        if (self.m1_prior==False) or (self.m1_prior=="uniform"):
            ret = self.uniform_fn(mass,0)
        else:
            if mass!=0:
                # First check if a gauss prob should be also calculated.
                if 0 not in [self.m1_est,self.m1_err]:
                    ret*=self.gaussian(mass, self.m1_est, self.m1_err)
                # Then caculate requested mass function prior.
                if (self.m1_prior == True) or (self.m1_prior == "PDMF"):
                    ret*=self.pdmf_prior(mass)
                elif self.m1_prior == "IMF":
                    ret*=self.imf_prior(mass)
        if ret==0: ret=-np.inf
        return ret

    def m2_prior_fn(self, m2, m1):
        ret = 1.0
        if (self.m2_prior==False) or (self.m2_prior=="uniform"):
            ret = self.uniform_fn(m2,1)
        else:
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
        if ret==0: ret=-np.inf
        return ret

    def para_prior_fn(self, para):
        ret = 1.0
        if (self.para_prior==False) or (self.para_prior=="uniform"):
            ret = self.uniform_fn(para,2)
        else:
            if para!=0.0:
                ## this needs parallax in arc sec, but standard in ExoSOFT is mas, so /1000
                ret = 1.0/((para/1000.0)**4.0)
                ##print('before gauss '+repr(ret))
                if 0 not in [self.para_est,self.para_err]:
                    ## a Gaussian prior centered on hipparcos and width of
                    #  hipparcos estimated error
                    ret*=self.gaussian(para, self.para_est, self.para_err)
                    ##print('para gauss ',self.gaussian(para, self.para_est, self.para_err))
        if ret==0: ret=-np.inf
        return ret

    def imf_prior(self, m):
        """From table 1 of Chabrier et al. 2003"""
        a = 0.068618528140713786
        b = 1.1023729087095586
        c = 0.9521999999999998
        # calculate IMF specific to mass of object
        if m<1.0:
            ret = (a/m) * np.exp( (-(np.log10(m)+b)**2) / c )
        else:
            ret = 0.019239245548314052*(m**(-2.3))
        if ret==0: ret=-np.inf
        return ret

    def pdmf_prior(self, m):
        """From table 1 of Chabrier et al. 2003"""
        a = 0.068618528140713786
        b = 1.1023729087095586
        c = 0.9521999999999998
        # calculate PDMF specific to mass of object
        if m<1.0:
            ret = (a/m)*np.exp((-(np.log10(m)+b)**2)/c)
        elif m<3.47:
            ret = 0.019108957203743077*(m**(-5.37))
        elif m<18.20:
            ret = 0.0065144172285487769*(m**(-4.53))
        else:
            ret = 0.00010857362047581295*(m**(-3.11))
        if ret==0: ret=-np.inf
        return ret

    def cmf_prior(self, m2, m1):
        """from equation 8 of Metchev & Hillenbrand 2009"""
        beta = -0.39
        ret = (m2**(beta)) * (m1**(-1.0*beta-1.0))
        if ret==0: ret=-np.inf
        return ret

    def gaussian(self,x,mu,sig):
        return np.exp( (-(x - mu)*(x - mu)) / (2.0*(sig*sig)) )

    def uniform_fn(self,val, i):
        """ a uniform prior to be used by many parameters """
        ret = 1.0
        if (len(self.mins_ary)>=i) and (len(self.maxs_ary)>=i):
            if self.maxs_ary[i]!=0:
                if (self.maxs_ary[i] - self.mins_ary[i]) > 0:
                    ## For all with uniform priors!!
                    mn = self.mins_ary[i]
                    mx = self.maxs_ary[i]
                    if mn<0:
                        mx = mx-mn
                    ret = stats.uniform.pdf(val,loc=mn, scale=mx)
                    #return 1.0 / (self.maxs_ary[i] - self.mins_ary[i])
        if ret==0: ret=-np.inf
        return ret
#END OF FILE

#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import constants
import sys

class Priors(object):
    def __init__(self,settings,log):
        self.settings = settings
        self.log = log
        
    def testPriors(self,last,proposed):
        val = self.combinedPriors(proposed,last)
        if False:
            print 'test priors ratio value = '+str(val)
                
    def combinedPriors(self,parsCurr,parsLast):
        """
        A function to combine priors in the settings dict.
        This can be used at the Simulator's instantiation to make sure the
        priors have no errors before starting a run.  This will be done 
        using the minimum range values.
        Else, it is just called during accept to calc the priors ratio.
        
        NOTE: -priors in the Advanced settings dict must be a tuple 
              of (bool, comment string, function).
              -Also, remember that non of the range values are allowed to be zero
              as it breaks this and a few other functions.
        """
        try:
            priorsRatio = self.combinedPriorsSingle(parsCurr)/self.combinedPriorsSingle(parsLast)
            return priorsRatio
        except:
            self.log.critical("An error occured while trying to calculate the priors ratio.")
            sys.exit(0)
        
    def combinedPriorsSingle(self,pars):
        """
        Combined priors for a single step in the chain.  This can be called 
        to form the priors ratio OR when accounting for the priors when 
        producing the posteriors during post-processing.
        """
        comboPriors = 1.0
        try:
            if self.settings['ePrior']:
                #print 'ePrior'
                comboPriors*=self.ePrior(pars[4])
            if self.settings['pPrior']:
                #print 'pPrior'
                comboPriors*=self.pPrior(pars[7])
            if self.settings['incPrior']:
                #print 'incPrior'
                comboPriors*=self.incPrior(pars[8])
            if self.settings['M1Prior']:
                #print 'M1Prior'
                comboPriors*=self.mass1Prior(pars[0])
                #print 'M1Prior'
            if self.settings['M2Prior']:
                #print 'M2Prior'
                comboPriors*=self.mass2(pars[1],pars[0])
                #print 'M2Prior'
            if self.settings['parPrior']:
                #print 'parPrior'
                comboPriors*=self.paraPrior(pars[2])
                #print 'parPrior out'
            return comboPriors
        except:
            self.log.critical("An error occured while trying to calculate the combined sigle priors.")
            sys.exit(0)
        
    #NOTE: only change the code and not the name of the functions or their inputs.            
    def ePrior(self,ecc):
        ## UPGRADE TO INCLUDE STRING INDICATED ECC PRIOR OPTIONS FROM KEPLER AND OTHER SURVEY RESULTS
        ret = 1.0
        if ecc!=0:
            if (self.settings['lowEcc']==False)and(self.settings['eMAX']!=0):
                if (self.settings['PMIN']*constants.daysPerYear)>1000.0:
                    ret = 2.0*ecc
        return ret
                  
    def pPrior(self,p):
        ret = 1.0
        if self.settings['PMAX']!=0:
            if p!=0.0:
                ret = 1.0/(p*(np.log(self.settings['PMAX']/self.settings['PMIN'])))
        return ret
        
    def incPrior(self,inc):
        ret = 1.0
        if self.settings['incMAX']!=0:
            if inc not in [0.0,90.0,180.0]:
                mn = self.settings['incMIN']*(constants.pi/180.0)
                mx = self.settings['incMAX']*(constants.pi/180.0)
                if self.settings['incPrior'] is 'sin':
                    ret = np.sin(inc*(constants.pi/180.0))/np.abs(np.cos(mn)-np.cos(mx))
                elif self.settings['incPrior'] is 'cos':
                    ret =  np.cos(inc*(constants.pi/180.0))/np.abs(np.cos(mn)-np.cos(mx))
        return ret
        
    def mass1Prior(self,mass):
        ret = 1.0
        if (self.settings['mass1MAX']!=0):
            if mass!=0:
                if self.settings['mass1Est']!=self.settings['mass1Err']!=0:
                    ret*=self.gaussian(mass, self.settings['mass1Est'], self.settings['mass1Err'])
                if (self.settings['M1Prior']=="PDMF")or(self.settings['M1Prior']==True):
                    ret*=self.pdmfPrior(mass)
                elif self.settings['M1Prior']=="IMF":
                    ret*=self.imfPrior(mass)
        return ret
        
    def mass2Prior(self,m2,m1):
        ret = 1.0
        if (self.settings['mass2MAX']!=0):
            if 0.0 not in [m2,m1]:
                if self.settings['mass2Est']!=self.settings['mass2Err']!=0:
                    ret*=self.gaussian(m2, self.settings['mass2Est'], self.settings['mass2Err'])
                if (self.settings['M2Prior']=="CMF")or(self.settings['M2Prior']==True):
                    ret*=self.cmfPrior(m2,m1)
                elif self.settings['M2Prior']=="PDMF":
                    ret*=self.pdmfPrior(m2)
                elif self.settings['M2Prior']=="IMF":
                    ret*=self.imfPrior(m2)
        return ret
            
    def paraPrior(self,para):
        ret = 1.0
        if 0.0 not in [para,self.settings['paraMAX']]:
            ret = 1.0/(para**4.0)
            if self.settings['paraEst']!=self.settings['paraErr']!=0:
                ## a Gaussian prior centered on hipparcos and width of hipparcos estimated error
                ret*=self.gaussian(para, self.settings['paraEst'], self.settings['paraErr'])
        return ret
            
    def imfPrior(self,m):
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        else:
            d = 0.019239245548314052*(m**(-2.3))
        return d
    
    def pdmfPrior(self,m):
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        elif m<3.47:
            d = 0.019108957203743077*(m**(-5.37))
        elif m<18.20:
            d = 0.0065144172285487769*(m**(-4.53))
        else:
            d = 0.00010857362047581295*(m**(-3.11))
        return d
    
    def cmfPrior(self,m2,m1):
        beta = -0.39
        d = (m2**(beta))*(m1**(-1.0*beta-1.0))
        return d
    
    def gaussian(self,x,mu,sig):
        return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

#END OF FILE
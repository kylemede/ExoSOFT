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
        priorsRatio = 1.0
        try:
            if self.settings['ePrior']:
                #print 'ePrior'
                priorsRatio*=self.ePriorRatio(parsCurr[4],parsLast[4])
            if self.settings['pPrior']:
                #print 'pPrior'
                priorsRatio*=self.pPriorRatio(parsCurr[7],parsLast[7])
            if self.settings['incPrior']:
                #print 'incPrior'
                priorsRatio*=self.incPriorRatio(parsCurr[8],parsLast[8])
            if self.settings['M1Prior']:
                #print 'M1Prior'
                priorsRatio*=self.mass1PriorRatio(parsCurr[0],parsLast[0])
                #print 'M1Prior'
            if self.settings['M2Prior']:
                #print 'M2Prior'
                priorsRatio*=self.mass2PriorRatio(parsCurr[1],parsLast[1],parsCurr[0],parsLast[0])
                #print 'M2Prior'
            if self.settings['parPrior']:
                #print 'parPrior'
                priorsRatio*=self.paraPriorRatio(parsCurr[2],parsLast[2])
                #print 'parPrior out'
            return priorsRatio
        except:
            self.log.critical("An error occured while trying to calculate the priors.")
            sys.exit(0)
        
    #NOTE: only change the code and not the name of the functions or their inputs.
    def ePriorRatio(self,eProposed,eLast):
        if (self.settings['lowEcc']==False)and(self.settings['eMAX']!=0):
            if eProposed!=eLast!=0:
                if (self.settings['PMIN']*constants.daysPerYear)>1000.0:
                    return eProposed/eLast
                else:
                    return 1.0
            else:
                return 1.0
        else:
            return 1.0
        
    def pPriorRatio(self,Pproposed,Plast):
        if self.settings['PMAX']!=0:
            if Pproposed!=0:
                return Plast/Pproposed
            else:
                return 1.0
        else:
            return 1.0
        
    def incPriorRatio(self,incProposed,incLast):
        if self.settings['incMAX']!=0:
            if (incLast%90.0)!=0:
                return np.sin(incProposed*(constants.pi/180.0))/np.sin(incLast*(constants.pi/180.0))
            else:
                return 1.0
        else:
            return 1.0
        
    def mass1PriorRatio(self,MProposed,MLast):
        if (self.settings['mass1MAX']!=0):
            if MProposed!=MLast!=0:
                gaussRatio = 1.0
                if self.settings['mass1Est']!=self.settings['mass1Err']!=0:
                    ## a Gaussian prior
                    top = self.gaussian(MProposed, self.settings['mass1Est'], self.settings['mass1Err'])
                    btm = self.gaussian(MLast, self.settings['mass1Est'], self.settings['mass1Err'])
                    gaussRatio = top/btm
                prop=1.0
                lst=1.0
                if (self.settings['M1Prior']=="PDMF")or(self.settings['M1Prior']==True):
                    prop = self.pdmfPrior(MProposed)
                    lst = self.pdmfPrior(MLast)
                elif self.settings['M1Prior']=="IMF":
                    prop = self.imfPrior(MProposed)
                    lst = self.imfPrior(MLast)
                return (prop/lst)*gaussRatio
            else:
                return 1.0
        else:
            return 1.0
        
    def mass2PriorRatio(self,m2Prop,m2Last,m1Prop,m1Last):
        if (self.settings['mass2MAX']!=0):
            if 0.0 not in [m2Prop,m2Last,m1Prop,m1Last]:
                gaussRatio = 1.0
                if self.settings['mass2Est']!=self.settings['mass2Err']!=0:
                    ## a Gaussian prior
                    top = self.gaussian(m2Prop, self.settings['mass2Est'], self.settings['mass2Err'])
                    btm = self.gaussian(m2Last, self.settings['mass2Est'], self.settings['mass2Err'])
                    gaussRatio = top/btm
                prop=1.0
                lst=1.0
                if (self.settings['M2Prior']=="CMF")or(self.settings['M2Prior']==True):
                    prop = self.cmfPrior(m2Prop,m1Prop)
                    lst = self.cmfPrior(m2Last,m1Last)
                elif self.settings['M2Prior']=="PDMF":
                    prop = self.pdmfPrior(m2Prop)
                    lst = self.pdmfPrior(m2Last)
                elif self.settings['M2Prior']=="IMF":
                    prop = self.imfPrior(m2Prop)
                    lst = self.imfPrior(m2Last)
                return (prop/lst)*gaussRatio
            else:
                return 1.0
        else:
            return 1.0
        
    def paraPriorRatio(self,paraProposed,paraLast):
        if paraProposed!=paraLast!=self.settings['paraMAX']!=0:
            ratioA = (paraLast**4.0)/(paraProposed**4.0)
            ratioB = 1.0
            if self.settings['paraEst']!=self.settings['paraErr']!=0:
                ## a Gaussian prior centered on hipparcos and width of hipparcos estimated error
                top = self.gaussian(paraProposed, self.settings['paraEst'], self.settings['paraErr'])
                btm = self.gaussian(paraLast, self.settings['paraEst'], self.settings['paraErr'])
                ratioB = top/btm
            return ratioA*ratioB
        else:
            return 1.0
        
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
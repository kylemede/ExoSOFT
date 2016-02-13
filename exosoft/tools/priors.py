#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import constants
import settings as sett

########################################
#Define the priors as python functions #
########################################
#NOTE: only change the code and not the name of the functions or their inputs.
def ePriorRatio(eProposed,eLast):
    if (sett.settingsDict['lowEcc'][0]==False)and(sett.settingsDict['eMAX']!=0):
        if eProposed!=eLast!=0:
            if (sett.settingsDict['PMIN']*constants.daysPerYear)>1000.0:
                return eProposed/eLast
            else:
                return 1.0
        else:
            return 1.0
    else:
        return 1.0
    
def pPriorRatio(Pproposed,Plast):
    if sett.settingsDict['PMAX']!=0:
        if Pproposed!=0:
            return Plast/Pproposed
        else:
            return 1.0
    else:
        return 1.0
    
def incPriorRatio(incProposed,incLast):
    if sett.settingsDict['incMAX']!=0:
        if (incLast%90.0)!=0:
            return np.sin(incProposed*(constants.pi/180.0))/np.sin(incLast*(constants.pi/180.0))
        else:
            return 1.0
    else:
        return 1.0
    
def mass1PriorRatio(MProposed,MLast):
    if (sett.settingsDict['mass1MAX']!=0)and True:
        if MProposed!=MLast!=0:
            prop = chabrierPrior(MProposed,IMF=False)
            lst = chabrierPrior(MLast,IMF=False)
            return prop/lst
        else:
            return 1.0
    else:
        return 1.0
    
def mass2PriorRatio(MProposed,MLast):
    if (sett.settingsDict['mass2MAX']!=0)and True:
        if MProposed!=MLast!=0:
            prop = chabrierPrior(MProposed,IMF=False)
            lst = chabrierPrior(MLast,IMF=False)
            return prop/lst
        else:
            return 1.0
    else:
        return 1.0
    
def paraPriorRatio(paraProposed,paraLast):
    if paraProposed!=paraLast!=sett.settingsDict['paraMAX']!=0:
        ratioA = (paraLast**4.0)/(paraProposed**4.0)
        ratioB = 1.0
        if sett.settingsDict['paraEst'][0]!=0:
            ## a Gaussian prior centered on hipparcos and width of hipparcos estimated error
            top = gaussian(paraProposed, sett.settingsDict['paraEst'][0], sett.settingsDict['paraErr'][0])
            btm = gaussian(paraLast, sett.settingsDict['paraEst'][0], sett.settingsDict['paraErr'][0])
            ratioB = top/btm
        return ratioA*ratioB
    else:
        return 1.0
    
def chabrierPrior(m,IMF=False):
    if IMF==False:
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        elif m<3.47:
            d = 0.019108957203743077*(m**(-5.37))
        elif m<18.20:
            d = 0.0065144172285487769*(m**(-4.53))
        else:
            d = 0.00010857362047581295*(m**(-3.11))
    else:
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        else:
            d = 0.019239245548314052*(m**(-2.3))
    return d

def gaussian(x,mu,sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

#END OF FILE
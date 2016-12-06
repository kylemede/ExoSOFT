#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
import gc
import copy
#import sys
#from scipy.constants.codata import precision
import timeit
import datetime
import emcee
from . import tools
import KMlogger
from six.moves import range
#from tools.model import constants as const

class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settings):
        self.paramsLast = 0
        self.paramsBestRaw = 0
        self.paramsBestStored = 0
        self.nSaved = 0
        self.nSavedPeriodic = 0
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.bestRedChiSqr = np.inf
        self.stgNsampDict = {'SA':'nSAsamp','ST':'nSTsamp','MC':'nSamples','MCMC':'nSamples','emcee':'nSamples'}
        self.acceptCount = 0
        self.acceptStr = ''
        self.shiftStr = ''
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.sigmasInRangeCounter = 0
        self.chainNum =0
        self.settings = settings
        self.log = KMlogger.getLogger('main.simulator',lvl=100,addFH=False)
        self.log.logSystemInfo()
        #$$$$$$$$$$$$$$$$$$$$$
        self.Model = tools.ExoSOFTmodel(self.settings)
        self.priors_last = 1
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        (self.realData,self.rangeMaxsRaw,self.rangeMinsRaw,self.rangeMaxs,self.rangeMins,self.starterSigmas,self.paramInts,self.nu,self.nuDI,self.nuRV) = self.starter() 
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #self.Orbit = tools.cppTools.Orbit()
        #self.Orbit.loadStaticVars(self.settings['omegaFdi'],self.settings['omegaFrv'],self.settings['lowEcc'],self.settings['pasa'])
        #self.Orbit.loadRealData(self.realData)
        #self.Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #Just initial seed val, reset in resetTracked() to be unique for each chain.
        self.seed = int(timeit.default_timer())
        np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.settings['finalFolder'],"tmpOutdata-"+str(self.chainNum)+".npy")
        
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        ## Recover parameter related items in the settings put there during startup
        realData = self.settings['realData']
        sigmas = self.settings['startSigmas']
        rangeMinsRaw = self.settings['rangeMinsRaw']
        rangeMaxsRaw = self.settings['rangeMaxsRaw'] 
        rangeMins = self.settings['range_mins']
        rangeMaxs = self.settings['range_maxs']
        paramInts = self.settings['paramInts']

        ##find total number of RV and DI epochs in real data
        nDIepochs = np.sum(np.where(realData[:,2]<1e6,1,0))
        nRVepochs = np.sum(np.where(realData[:,6]<1e6,1,0))
        nEpochs = len(realData[:,0])
        ##Take mass1, dist, inc and period from those include in nu calcs
        ##as they have clear priors.
        paramIntsClean = copy.deepcopy(paramInts)
        notInNuInts = [2,6,7] #ie. parallax, period and inclination as they always have well defined priors.
        for val in notInNuInts:
            paramIntsClean=paramIntsClean[np.where(paramIntsClean!=val)]
        diVars = paramIntsClean[np.where(paramIntsClean<9)] #all, but rv offsets
        rvVars = paramIntsClean[np.where(paramIntsClean!=3)] #all but long_an
        self.log.debug('DIvars cleaned = '+repr(diVars))
        self.log.debug('RVvars cleaned = '+repr(rvVars))
        nDIvars = len(diVars)
        nRVvars = len(rvVars)
        allVars = paramInts
        if nDIepochs==0:
            nVars = nRVvars
            allVars = diVars
        elif nRVepochs==0:
            nVars = nDIvars
            allVars = rvVars
        else:
            nVars = len(paramInts)
        self.log.debug("vars = "+repr(paramInts))
        self.log.debug('[nEpochs, nDIepochs, nRVepochs] = ['+str(nEpochs)+', '+str(nDIepochs)+', '+str(nRVepochs)+']')
        self.log.debug('[nVars, nDIvars, nRVvars] = ['+str(nVars)+', '+str(nDIvars)+', '+str(nRVvars)+']')
        nuDI = 1
        nuRV = 1
        nu = 1
        if nDIepochs*2>nDIvars:
            nuDI = nDIepochs*2-nDIvars
        else:
            self.log.debug("nDIepochs*2>nDIvars is False so setting nuDI=1")
        if nRVepochs>nRVvars:
            nuRV = nRVepochs-nRVvars
        else:
            self.log.debug("nRVepochs>nRVvars is False so setting nuRV=1")
        if (nDIepochs*2+nRVepochs)>nVars:
            nu = nDIepochs*2+nRVepochs-nVars
        else:
            self.log.debug("(nDIepochs*2+nRVepochs)>nVars is False so setting nu=1")
        self.log.debug('[nu, nuDI, nuRV] = ['+str(nu)+', '+str(nuDI)+', '+str(nuRV)+']')
        #load these into settings dict
        self.settings["nRVdsets"] = len(self.settings['offset_mins'])
        self.settings['commentsDict']['nRVdsets'] = "Number of RV data sets"
        self.settings['nDIepoch'] = nDIepochs
        self.settings['commentsDict']['nDIepoch'] = "Number of DI epochs"
        self.settings['nRVepoch'] = nRVepochs
        self.settings['commentsDict']['nRVepoch'] = "Number of RV epochs"
        self.settings['n3Depoch'] = nEpochs
        self.settings['commentsDict']['n3Depoch'] = "Number of 3D epochs"
        self.settings['nu'] = nu
        self.settings['commentsDict']['nu'] = "Total nu"
        self.settings['nuDI'] = nuDI
        self.settings['commentsDict']['nuDI'] = "nu for DI"
        self.settings['nuRV'] = nuRV
        self.settings['commentsDict']['nuRV'] = "nu for RV"
        paramIntsStr = repr(paramInts).replace(' ','')
        self.settings['parInts'] = paramIntsStr
        self.settings['commentsDict']['parInts'] = "Varried params"
        self.settings['chainNum'] = self.chainNum
        self.settings['commentsDict']['chainNum'] = "chain number"
        self.settings['3DVars'] = allVars
        self.settings['DIvars'] = diVars
        self.settings['RVvars'] = rvVars
        
        sigMaxs = np.zeros(rangeMinsRaw.shape)
        sigMins = np.zeros(rangeMinsRaw.shape)
        for i in paramInts:
            sigMaxs[i]=(rangeMaxsRaw[i]-rangeMinsRaw[i])*self.settings['sigMax']
            sigMins[i]=(rangeMaxsRaw[i]-rangeMinsRaw[i])*self.settings['sigMin']
        self.settings['sigMaxs'] = sigMaxs
        self.settings['sigMins'] = sigMins
        #self.log.debug('sigMaxs = '+repr(sigMaxs))
        #self.log.debug('sigMins = '+repr(sigMins))
        ## check priors are ok with range mins
        ## will return the value if all good, or call for sys.exit() if not.
        _ = self.Model.Priors.priors_ratio(rangeMins,rangeMins)
        self.log.debug("priors ratio function completed with no errors.")
        #Priors.testPriors(rangeMins,rangeMins)
        
        return (realData,rangeMaxsRaw,rangeMinsRaw,rangeMaxs,rangeMins,sigmas,paramInts,nu,nuDI,nuRV)
    
    def increment(self,pars=[],sigs=[],stage=''):
        """
        Increment all varying parameters if MC.
        Else, just increment one of them at random.
        """
        parsOut = copy.deepcopy(pars)
        varyInt=0
        sig = 0
        ## vary all the params if MC, 
        ##(or special cases at beginning of SA where this func is called pretending to be MC.)
        if  (('MCMC' not in stage) and (stage=='MC')) or (stage=='SA'and(self.acceptCount==0)):
            for i in range(0,len(pars)):
                if i in self.paramInts:
                    parsOut[i]=np.random.uniform(self.rangeMinsRaw[i],self.rangeMaxsRaw[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            self.parIntVaryAry.append(varyInt)
            sig = sigs[varyInt]
            parsOut[varyInt]=np.random.uniform(pars[varyInt]-sig,pars[varyInt]+sig)        
        
        ## if Kdirect not set, then inclination varys.
        ## then K=0 going into Orbit so Orbit will calc it
        #if 8 in self.paramInts:
        #    parsOut[12] = 0
        
        return parsOut
    
    #def rangeCheck(self,pars,sample,stage=''):
    #    """
    #    Check if values inside allowed ranges.  For those that are periodic, 
    #    wrap them into the allowed ranges ([0,360] for Omega and omega).
    #    """
    #    inRange=True
    #    #convert from raw (directly varied) params, to storable versions
    #    
    #    if (sample>=10)and((self.acceptCount==0)and(stage=='SA')):
    #        ##Jump as starting position after first 10 tries was in poor part of param space. for SA only.
    #        paramsOut = self.increment(self.rangeMinsRaw,np.zeros(pars.shape),stage='MC')
    #        ## convert from Raw form if in lowEcc mode
    #        self.Orbit.convertParsFromRaw(paramsOut)
    #        inRange=True
    #    return (paramsOut,inRange)
    
    def accept(self,sample,pars,ln_post=1,temp=1.0,stage=''):
        """
        First this will calculate chi squared for model vs real data.
        
        For mcOnly it performs simple chisquared cut-off acceptance 
        based on 'chiMAX' value in settings.
        Else, it will calculate the priors and accept based on 
        the Metropolis-Hastings algorithm. The temp factor will 
        be set to 1.0 for MCMC and Sigma Tuning, and should be provided 
        for Simulated Annealing.
        """
        #print "ln210"
        paramsOut = copy.deepcopy(pars) #stored versions
        ## Calculate chi squareds for 3D,DI,RV and update bestPars and bestSumStr if this is better than the best
        raw3D = self.Model.chi_squared_3d
        reduced3D = raw3D/self.nu
        #print("\n\nreduced3D = "+str(reduced3D))#$$$$$$$$
        reducedDI = self.Model.chi_squared_di/self.nuDI
        reducedRV = self.Model.chi_squared_rv/self.nuRV
        #print "ln216"
        
        #(raw3D, reducedDI, reducedRV, reduced3D) = tools.chiSquaredCalc3D(self.realData,modelData,self.nuDI,self.nuRV,self.nu)
        #paramsOut[11] = raw3D
        # Store posterior probability (instead of chi squared)
        #print("ln_post = "+repr(ln_post))
        #print("paramsOut[11] = "+repr(paramsOut[11]))
        if self.bestSumStr=='':
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+' Nothing accepted yet below chi squared max = '+str(self.settings['chiMAX'])
            self.latestSumStr="Latest reduced chiSquared : [total,DI,RV] = ["+str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
        #print "ln223"
        if (reduced3D)<self.bestRedChiSqr:
            self.bestRedChiSqr=(reduced3D)
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+\
            ' BEST reduced chiSquareds so far: [total,DI,RV] = ['\
            +str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
            #bestPars = copy.deepcopy(paramsOut)
            ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
            #self.Orbit.convertParsToRaw(bestPars)
            #self.paramsBestRaw = bestPars
            self.paramsBestRaw = copy.deepcopy(self.Model.Params.direct_pars)
            self.paramsBestStored = paramsOut
            #print "ln234"
        ## check if this step is accepted
        accept = False
        if (stage=='MC')or(stage=='SA'and(self.acceptCount==0)):
            ## for MC or first step of SA
            if (raw3D/self.nu)<self.settings['chiMAX']:
                accept=True
            #print "ln241"
        else:
            ## For SA after first sample, MCMC, and ST
            try:               
                likelihoodRatio = np.exp((self.paramsLast[11] - raw3D)/(2.0*temp))
                #print("\nlikelihoodRatio  = "+repr(likelihoodRatio))
                priorsRatio = self.Model.prior/self.priors_last
                #print("priorsRatio  = "+repr(priorsRatio))
                #print("priorsRatio*likelihoodRatio  = "+repr(priorsRatio*likelihoodRatio))
                #prob = priorsRatio*likelihoodRatio
                #print "ln249"
                ############################################################
                ## Decide using Metropolis-Hastings basic rejection function
                if np.random.uniform(0.0, 1.0)<=(priorsRatio*likelihoodRatio):
                    accept = True
                #print "ln253"
            except:
                accept = False
        ## check for all modes to make sure m2 is never >m1
        if paramsOut[1]>paramsOut[0]:
            accept = False
        ## Update counters and strings
        if accept:
            self.acceptCount+=1
            self.acceptBoolAry.append(1)
            self.paramsLast=paramsOut
            self.priors_last = self.Model.prior
            self.latestSumStr = stage+" chain #"+str(self.chainNum)+\
            ' Latest accepted reduced chiSquareds: [total,DI,RV] = ['+\
            str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
            #if (stage=='SA')and(self.acceptCount>10):
            #    self.Orbit.NewtonWarningsOn(True)
            
        else:
            self.acceptBoolAry.append(0)
        ##log a status summary?
        if self.settings['nSumry']>0:
            if sample%(self.settings[self.stgNsampDict[stage]]//self.settings['nSumry'])==0:
                perc = sample*100//self.settings[self.stgNsampDict[stage]]
                ####str(self.nSaved)+' (curr '+str(self.nSavedPeriodic)+"), Finished: "+str(sample)+"/"+\
                sumStr = "below\n"+stage+" chain #"+str(self.chainNum)+", # Accepted: "+str(self.acceptCount)+", # Saved: "+\
                str(self.nSaved)+", Finished: "+str(sample)+"/"+\
                str(self.settings[self.stgNsampDict[stage]])+" = "+str(perc)+"%, Current T: "+str(temp)+"\n"
                sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
                #print "Accepted: "+str(self.acceptCount)+'\n'+self.bestSumStr #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                self.log.debug(sumStr)
        
        
        paramsOut2 = copy.deepcopy(paramsOut)
        ## enabling this would make the 11th column of both all algorithm outputs be the same.
        if False:
            ## Replace 3D chi squared with posterior probability
            paramsOut2[11] = ln_post#np.exp(ln_post)
        
        return (accept,paramsOut2)
    
    def tempDrop(self,sample,strtTemp,temp,stage=''):
        """
        Determine if it is time to drop the temp, and drop if it is.
        Total temperature range is [strtTemp,0.01), so the minimum 
        temperature is actually <1.0 meaning the last few temperature drops 
        will really push the currently found minimum towards its peak.
        There will be a fixed number of temperature steps = 'nTmpStps'.
        """
        if stage=='SA':
            if sample%self.settings['tempInt']==0:
                temp-=(strtTemp-0.01)*(float(self.settings['tempInt'])/float(self.settings['nSAsamp']))
        return temp
    
    def sigTune(self,sample,sigs=[],stage=''):
        """
        Check if it is time to calculate the acceptance rate.
        If stage is ST, then it will also tune the sigmas.
        If MCMC then it will just calculate the acceptance 
        rate.  In both cases, 'nSigStps' throughout the simulation 
        a summary message will also be written to the log.
        """
        sigmasOut = copy.deepcopy(sigs)
        if (stage=='ST')or(stage=='MCMC'):
            if (sample%(len(self.paramInts)*self.settings['sigInt'])==0)and(self.acceptCount>1):
                self.acceptStr = '\n'+stage+" chain #"+str(self.chainNum)+'\n'
                self.shiftStr = ''
                self.parIntVaryAry = np.array(self.parIntVaryAry)
                self.acceptBoolAry = np.array(self.acceptBoolAry)
                self.acceptStr+="Number of steps used to calculate acceptance rate = "+repr(len(self.acceptBoolAry))+'\n'
                acceptArray = []
                for i in self.paramInts:
                    ##calculate acceptance rate for each param
                    nAcc = np.sum(np.where(self.parIntVaryAry==i,self.acceptBoolAry,0))
                    nTot = len(np.where(self.parIntVaryAry==i)[0])
                    self.acceptStr+= 'parameter # '+str(i)+' acceptance = '+str(float(nAcc)/float(nTot))+'\n'
                    acceptArray.append(float(nAcc)/float(nTot))
                    if stage=='ST':
                        ##check each rate to choose up/down shift and do so and update shiftStr
                        self.shiftStr+= '\n'+stage+" chain #"+str(self.chainNum)+'\nparameter # '+str(i)+" shifting sigma "+str(sigs[i])+" -> "
                        if ((float(nAcc)/float(nTot))>self.settings["accRates"][1])and(sigs[i]<self.settings['sigMaxs'][i]):
                            sigmasOut[i]+=self.settings['sigMins'][i]
                        elif ((float(nAcc)/float(nTot))<self.settings["accRates"][0])and(sigs[i]>self.settings['sigMins'][i]):
                            sigmasOut[i]-=self.settings['sigMins'][i]
                        #fix chances that sigs can get out of their allowed range
                        if sigmasOut[i]<self.settings['sigMins'][i]:
                            sigmasOut[i]=self.settings['sigMins'][i]
                            self.log.debug('chain# '+str(self.chainNum)+', sigma for par # '+str(i)+' went below allowed range, so setting to max '+str(sigmasOut[i]))
                        if sigmasOut[i]>self.settings['sigMaxs'][i]:
                            sigmasOut[i]=self.settings['sigMaxs'][i]
                            self.log.debug('chain# '+str(self.chainNum)+', sigma for par # '+str(i)+' went above allowed range, so setting to max '+str(sigmasOut[i]))
                        self.shiftStr+=str(sigmasOut[i])+"\n"
                self.acceptBoolAry = []
                self.parIntVaryAry = []
                ##log a status summary?
                if self.settings['nSumry']>0:
                    if sample%(self.settings[self.stgNsampDict[stage]]//self.settings['nSumry'])==0:
                        self.log.debug(self.acceptStr+self.shiftStr)
                ## check if sigmas are same as last time and update tracker
                accRatesInRange = True
                accMin = 0.9*self.settings["accRates"][0]
                accMax = 1.1*self.settings["accRates"][1]
                for acRt in acceptArray:
                    if (acRt<accMin) or (acRt>accMax):
                        accRatesInRange = False
                if accRatesInRange:
                    self.sigmasInRangeCounter+=1
                    #print 'sigmas were the same '+str(self.sigmasInRangeCounter)+" times"
                else:
                    self.sigmasInRangeCounter=0                
        else:
            #No need to track these, so reset them to empty to save any RAM they are using
            self.acceptBoolAry = []
            self.parIntVaryAry = []
        return sigmasOut
    
    def endSummary(self,temp,sigmas,stage=''):
        """
        Make a final summary of important statistics for the chain.
        """
        sumStr = '\n'+"="*70+"\nEND OF "+stage+" CHAIN #"+str(self.chainNum)+" SUMMARY:\nFinalTemp = "
        sumStr+= str(temp)+"\nTotal number of steps accepted = "+str(self.acceptCount)+"\n"
        sumStr+= "Average acceptance rate = "
        sumStr+=str(float(self.acceptCount)/float(self.settings[self.stgNsampDict[stage]]))+"\n"
        sumStr+= "Total number of steps stored = "+str(self.nSaved)+"\n"
        sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
        sumStr+="Last params  = "+repr(self.paramsLast)+'\n'
        sumStr+="Best params (Raw) = "+repr(self.paramsBestRaw)+'\n'
        if (stage=="ST")or(stage=="MCMC"):
            if stage=='ST':
                sumStr+="Final Sigmas = "+repr(sigmas)+'\n'
            sumStr+=self.acceptStr+self.shiftStr
        sumStr+='\n'+'='*70+'\n'
        self.log.info(sumStr)
        return float(self.acceptCount)/float(self.settings[self.stgNsampDict[stage]])
    
    def resetTracked(self,stage):
        """
        Reset the internal strings, arys and counters.
        """
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.acceptStr = ''
        self.shiftStr = ''
        self.bestRedChiSqr = 1e6
        self.nSaved = 0
        self.acceptCount = 0
        self.nSavedPeriodic = 0
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.settings['chainNum'] = self.chainNum
        self.settings['commentsDict']['nDIepoch'] = "chain number"
        self.settings['curStg']= stage
        self.settings['commentsDict']['curStg'] = 'Current stage either [SA,ST,MCMC or MC]'
        #if stage=='SA':
        #    self.Orbit.NewtonWarningsOn(False)
        #else:
        #    self.Orbit.NewtonWarningsOn(True)
        # make a very random seed value to ensure each chain is different.  
        # Should we make this value an optional input and pass on as a return value to keep a process number using the same seed?? $$$
        # if so, it needs to be pushed into the results file as well.
        t = np.random.uniform(1,1e6)
        self.seed = int((timeit.default_timer()/(self.chainNum+1))/t)
        self.log.debug("Chain# "+str(self.chainNum)+" has random number seed = "+str(self.seed))
        np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.settings['finalFolder'],"tmpOutdata-"+str(self.chainNum)+".npy")
        if os.path.exists(self.tmpDataFile):
            os.remove(self.tmpDataFile)
            self.log.debug("just removed data file from disk:\n"+self.tmpDataFile)
    
    def startSummary(self,pars,sigs,stage):
        startStr=stage+" chain #"+str(self.chainNum)+' VALS AT START OF '+stage+' SIM:\n'
        startStr+= 'params = '+repr(pars)+'\n'   
        startStr+= 'rangeMins = '+repr(self.rangeMins)+'\n'
        startStr+= 'rangeMaxs = '+repr(self.rangeMaxs)+'\n'
        if self.settings['low_ecc']:
            startStr+= 'rangeMinsRaw = '+repr(self.rangeMinsRaw)+'\n'
            startStr+= 'rangeMaxsRaw = '+repr(self.rangeMaxsRaw)+'\n' 
        startStr+= 'sigmas = '+repr(sigs)+'\n'
        startStr+= 'paramInts = '+repr(self.paramInts)+'\n'  
        startStr+= '[nu,nuDI,nuRV] = ['+str(self.nu)+', '+str(self.nuDI)+', '+str(self.nuRV)+']\n'
        self.log.debug(startStr)
        
    def simulatorFunc(self,stage='',chainNum=1,startParams=[],startSigmas=[],temp=1.0):
        """
        The core function to perform the requested stage of the simulation ('MC','SA','ST','MCMC').
        If stage is SA or ST: final (params,sigmas) are returned, else nothing.
        """
        tic=timeit.default_timer()
        lastTic=tic
        timesAry = []
        self.log.debug("Trying "+str(self.settings[self.stgNsampDict[stage]])+" samples for chain #"+str(chainNum)+" in "+stage+" mode.")
        self.chainNum = chainNum
        self.resetTracked(stage)
        bar = tools.ProgBar(total=100,barLength=0)
        modelData = np.zeros((len(self.realData),3))
        acceptedParams = []
        strtTemp = temp      
        endDatetime = ''
        sigmas = copy.deepcopy(self.starterSigmas)
        ## if valid startSigmas provided, start with them, else use defaults.
        if (type(startSigmas)==list)or(type(startSigmas)==np.ndarray):
            if len(startSigmas)>0:
                sigmas = copy.deepcopy(startSigmas)
        ## if valid startParams provided, start there, else start at random point.
        if (type(startParams)==list)or(type(startParams)==np.ndarray):
            if len(startParams)>0:
                #convert 'stored' to 'direct/raw' versions
                paramsLast = copy.deepcopy(self.Model.Params.stored_to_direct(startParams))
                #paramsLast = copy.deepcopy(startParams)
                self.log.debug('initial/latest pars have reduced chi sqr of '+str(startParams[11]/self.nu))
            else: 
                paramsLast = self.increment(self.rangeMinsRaw,sigmas,stage='MC')
        else: 
            paramsLast = self.increment(self.rangeMinsRaw,sigmas,stage='MC')    
        ## load up starting params as 'latest' and perform first increment from these to start loop with.
        latestParsRaw = copy.deepcopy(paramsLast)
        proposedParsRaw = self.increment(latestParsRaw,sigmas,stage)
        ## convert from Raw form if in lowEcc mode     
        ln_post = tools.ln_posterior(latestParsRaw,self.Model)       
        paramsLast = copy.deepcopy(self.Model.Params.stored_pars)
        #print('proposed pars have reduced chi sqr of '+str(paramsLast[11]/self.nu))#$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #self.Orbit.convertParsFromRaw()
        self.paramsLast = paramsLast
        self.startSummary(paramsLast,sigmas,stage)
        self.acceptBoolAry.append(0)
        #if stage =='MCMC':
        #    #load starting point to make accepted numbers match perfect to requested samples
        #    acceptedParams.append(self.paramsLast)
        #    self.nSaved += 1 
        #    self.nSavedPeriodic+=1
        #    self.acceptBoolAry.append(0)
    
        
        ############################################################################
        ##  All the general ExoSOFT inputs are now ready to go.  Next, either      #
        ##  prepare, call, and wrap-up the use of emcee, OR do a loop over the     #
        ##  samples to perform ExoSOFT implimentations of MC, SA, ST or M-H MCMC.  #
        ############################################################################
        avgAcceptRate = 0
        if stage=='emcee':
            #print("preparing to call emcee")
            ###$$$$$$ Make lots of these part of the settings dict $$$$$$
            ncpu = self.settings['nMCMCcns']
            ndim = len(startParams) # number of parameters in 'stored' format
            ndim_raw = len(self.settings['rangeMinsRaw']) # number of parameters in the model, in 'raw'/'direct' format
            nwalkers = self.settings['n_wlkrs']
            nsteps = int(float(self.settings['nSamples'])/float(nwalkers))  ##$$$ leave this way, or add to settings dict?
            #print('\n'*3+'Calling make_start_params')
            #print('\n'+repr(latestParsRaw)+'\n')
            ## NOTE: starting_guesses array must be a numpy array with dtype=np.dtype('d').
            starting_guesses = tools.make_starting_params(latestParsRaw,nwalkers,scale=0.01)###$$$$$$$$$$$$$$$$$$$$$$$$$ maybe move this over to the built in starting samples func
            #print('\n'*3+'back from make_start_params')
            ## Call emcee to explore the parameter space
            sampler = emcee.EnsembleSampler(nwalkers, ndim_raw, tools.ln_posterior, 
                                            args=[self.Model], threads=ncpu)
            self.log.importantinfo("\nObjects and inputs prepared, now calling emcee.")
            sampler.run_mcmc(starting_guesses, nsteps)
            self.log.importantinfo("ensamble sampling with emcee complete.\n")
            
            ################################################################
            ### Refactor emcee trace/chain into ExoSOFT stored params format
            ################################################################
            #t1 = timeit.default_timer()
            # chain is of shape (nwalkers, nsteps, ndim)
            # discard burn-in points and reshape
            
            trace = sampler.chain
            probs = sampler.lnprobability
            print("\nBefore burn-in strip")
            print("\n"+repr(trace.shape))
            print('\n'+repr(probs.shape))               
            if self.settings['rmBurn']:
                trace = trace[:, self.settings['n_emcee_burn']:, :]
                probs = probs[:, self.settings['n_emcee_burn']:]
            print("\nAfter burn-in strip")
            print("\n"+repr(trace.shape))
            print('\n'+repr(probs.shape))
            if self.settings['thin_emcee']:
                trace = trace[:, ::self.settings['thin_rate'], :]
                probs = probs[:, ::self.settings['thin_rate']]
            print("\nAfter burn-in strip AND thinning")
            print("\n"+repr(trace.shape))
            print('\n'+repr(probs.shape))
            probs = probs.reshape(-1)
            #print('\n'+repr(probs.shape))
            trace = trace.reshape(-1, ndim_raw)
            # Now the trace shape is (nwalkers*nsteps, ndim)
            ## Go through every sample output by emcee and convert to 
            ## 'stored' format for later use by ExoSOFT
            #print('ndim = '+str(ndim))
            #print('ndim_raw = '+str(ndim_raw))
            #print('\ntrace.shape = '+repr(trace.shape))
            acceptedParams = np.ones((trace.shape[0],ndim))
            for i in range(trace.shape[0]):
                acceptedParams[i] = self.Model.Params.direct_to_stored(trace[i])
                #print(repr(a)+'\n')
            ## push probability = np.exp(lnprob) into chi squared column
            acceptedParams[:,11] = probs#np.exp(probs)
            #print('\nrepr accepted: '+repr(acceptedParams))
            #print("Converting emcee sample outputs to ExoSOFT format took: "+tools.timeStrMaker(timeit.default_timer()-t1))
            
            self.acceptStr = "Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction))
            #print(repr(sampler.acceptance_fraction))
            avgAcceptRate = np.mean(sampler.acceptance_fraction)
            tools.periodicDataDump(self.tmpDataFile,acceptedParams)
        else:
            ###################################################################
            ##loop through each sample 
            ##Follows these steps: in Range?,calc model,accept?,Store?,increment,lower temp?,tune sigmas? dump data to disk? DONE ->write output data
            sample=0
            #last_ln_post=0#$$$$$$$$
            while sample<(self.settings[self.stgNsampDict[stage]]+1):
                sample+=1
                # call ln_posterior as it will check the ranges first.
                # if in range, then it will run the model, else return -np.inf
                #if np.array_equal(proposedParsRaw, latestParsRaw):
                #    print 'proposed and latest are same'
                #else:
                #    print 'proposed and latest NOT same'
                ln_post = tools.ln_posterior(proposedParsRaw,self.Model)
                #print 'ln_post = '+repr(ln_post)
                # check if parameters in range (out of range ln_post = -np.inf)
                if ln_post>(-np.inf):
                    proposedPars = copy.deepcopy(self.Model.Params.stored_pars)
                    #if np.array_equal(self.paramsLast,proposedPars):
                    #    print 'same params'
                    (accept,paramsOut) = self.accept(sample,proposedPars,ln_post,temp,stage)
                    #print 'self.Model.chi_squared_3d = '+repr(self.Model.chi_squared_3d)
                    if accept:
                        #print "self.Model.Params.stored_pars = "+repr(self.Model.Params.stored_pars)
                        
                        #print 'self.Model.chi_squared_3d = '+repr(self.Model.chi_squared_3d)
                        #last_ln_post = ln_post
                        latestParsRaw = copy.deepcopy(self.Model.Params.direct_pars)
                        ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
                        #self.Orbit.convertParsToRaw(latestParsRaw)
                        if ('MCMC' not in stage)and(stage in ['MC','SA','ST']):
                            acceptedParams.append(paramsOut)#proposedPars)
                            #print 'proposedPars :\n'+repr(proposedPars)
                            #print "acceptedParams= "+repr(acceptedParams)
                            
                            self.nSaved+=1
                            self.nSavedPeriodic+=1                 
                else:
                    self.acceptBoolAry.append(0)
                if (stage in ['MCMC']) and (sample%self.settings['saveInt']) == 0:
                    acceptedParams.append(self.paramsLast)
                    self.nSaved += 1 
                    self.nSavedPeriodic+=1
                proposedParsRaw = self.increment(latestParsRaw,sigmas,stage)
                temp = self.tempDrop(sample,strtTemp,temp,stage)
                sigmas = self.sigTune(sample,sigmas,stage)
                if (self.nSavedPeriodic>0)and(self.nSavedPeriodic==self.settings['dmpInt']):
                    ## dump acceptedParams array to disk and collect garbage
                    self.log.debug('Dumping data to filename:\n'+self.tmpDataFile+\
                                   '\nThe acceptedParams Ary had '+str(np.array(acceptedParams).shape[0])+\
                                   ' param sets, and size = '+str(np.array(acceptedParams).nbytes/1.0e6)+' MB')
                    tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
                    acceptedParams = []
                    self.nSavedPeriodic = 0
                    self.log.debug('about to collect the garbage')
                    gc.collect()
                if sample%(self.settings[self.stgNsampDict[stage]]//100)==0:
                    #update predicted completion time every 1%
                    timesAry.append(timeit.default_timer()-lastTic)
                    lastTic = timeit.default_timer()
                    timeRemSec = np.mean(timesAry)*(100.0-(float(sample)*100.0)/float(self.settings[self.stgNsampDict[stage]]))
                    #timeStr = ' about '+tools.timeStrMaker(timeRemSec)+' remaining.'
                    endDatetime = " Will be done"+tools.dateStrMaker(datetime.datetime.now(),timeRemSec)
                    if self.settings['logLevel']<30:
                        perc = int(sample*100//self.settings[self.stgNsampDict[stage]])
                        bar.render(perc, suffix=stage+str(chainNum)+' '+endDatetime)#timeStr)
                if sample%(self.settings[self.stgNsampDict[stage]]//10)==0:
                    #push predicted completion time log file
                    self.log.fileonly(stage+str(chainNum)+' '+endDatetime)            
                if (stage=='ST') and ((self.sigmasInRangeCounter>10)and(self.nSaved>1)):
                    self.log.debug("sigmasInRangeCounter>10 so breaking sample loop.")
                    break
            if self.settings['logLevel']<30:
                bar.render(100,suffix=stage+str(chainNum)+' Complete!\n')
            self.log.debug(stage+" took: "+tools.timeStrMaker(timeit.default_timer()-tic))
            avgAcceptRate = self.endSummary(temp,sigmas,stage)
            #print str(avgAcceptRate)+'\n'
            tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
            
            if stage=='SA':
                ##start ST at the best location with tight sigmas, and it will tune to ideal sigmas
                sigmas = self.settings['sigMins']
                
        clobber = False
        if self.settings['curStg']=="SA":
            clobber = True
        outFname = tools.writeFits('outputData'+stage+str(chainNum)+'.fits',self.tmpDataFile,self.settings,clob=clobber)
        ## A couple extra wrap-up tasks needed to adapt emcee outputs to match ExoSOFT later stages
        if stage=='emcee':
            self.paramsBestStored = tools.findBestOrbit(outFname,bestToFile=True,findAgain=False,by_ln_prob=True)
            #self.bestRedChiSqr
            #self.acceptStr
                
        return (outFname,self.paramsBestStored,sigmas,self.bestRedChiSqr,avgAcceptRate,self.acceptStr)
#END OF FILE      
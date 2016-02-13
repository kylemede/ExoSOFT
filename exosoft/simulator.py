#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import gc
import copy
import sys
#from scipy.constants.codata import precision
import tools
import timeit
from tools import constants as const

class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settingsDict):
        self.paramsLast = 0
        self.paramsBest = 0
        self.nSaved = 0
        self.nSavedPeriodic = 0
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.bestRedChiSqr = 1e6
        self.stgNsampDict = {'SA':'nSAsamp','ST':'nSTsamp','MC':'nSamples','MCMC':'nSamples'}
        self.acceptCount = 0
        self.acceptStr = ''
        self.shiftStr = ''
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.chainNum =0
        self.settingsDict = settingsDict
        self.log = tools.getLogger('main.simulator',lvl=100,addFH=False)
        tools.logSystemInfo(self.log)
        (self.realData,self.rangeMaxsRaw,self.rangeMinsRaw,self.rangeMaxs,self.rangeMins,self.starterSigmas,self.paramInts,self.nu,self.nuDI,self.nuRV) = self.starter() 
        self.Orbit = tools.cppTools.Orbit()
        self.Orbit.loadStaticVars(self.dictVal('omegaFdi'),self.dictVal('omegaFrv'),self.dictVal('lowEcc'),self.dictVal('pasa'))
        self.Orbit.loadRealData(self.realData)
        self.Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
        #Just initial seed val, reset in resetTracked() to be unique for each chain.
        self.seed = int(timeit.default_timer())
        np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.dictVal('finalFolder'),"tmpOutdata-"+str(self.chainNum)+".npy")
        
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        ## Recover parameter related items in the settingsDict put there during startup
        realData = self.dictVal('realData')
        sigmas = self.dictVal('startSigmas')
        rangeMinsRaw = self.dictVal('rangeMinsRaw')
        rangeMaxsRaw = self.dictVal('rangeMaxsRaw') 
        rangeMins = self.dictVal('rangeMins')
        rangeMaxs = self.dictVal('rangeMaxs')
        paramInts = self.dictVal('paramInts')

        ##find total number of RV and DI epochs in real data
        nDIepochs = np.sum(np.where(realData[:,2]<1e6,1,0))
        nRVepochs = np.sum(np.where(realData[:,6]<1e6,1,0))
        nEpochs = len(realData[:,0])
        ##Take mass1, dist, inc and period from those include in nu calcs
        ##as they have clear priors.
        paramIntsClean = copy.deepcopy(paramInts)
        notInNuInts = [2,7,8]      
        for val in notInNuInts:
            paramIntsClean=paramIntsClean[np.where(paramIntsClean!=val)]
        nDIvars = np.sum(np.where(paramIntsClean<10,1,0))
        self.log.debug('DIvars cleaned = '+repr(paramIntsClean[np.where(paramIntsClean<10)]))
        self.log.debug('RVvars cleaned = '+repr(paramIntsClean[np.where(paramIntsClean!=3)]))
        nRVvars = np.sum(np.where(paramIntsClean!=3,1,0))
        if nDIepochs==0:
            nVars = nRVvars
        elif nRVepochs==0:
            nVars = nDIvars
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
        self.settingsDict["nRVdsets"] = (len(self.dictVal('vMINs')),"Number of RV data sets")
        self.settingsDict['nDIepoch'] = (nDIepochs,"Number of DI epochs")
        self.settingsDict['nRVepoch'] = (nRVepochs,"Number of RV epochs")
        self.settingsDict['n3Depoch'] = (nEpochs,"Number of 3D epochs")
        self.settingsDict['nu'] = (nu,"Total nu")
        self.settingsDict['nuDI'] = (nuDI,"nu for DI")
        self.settingsDict['nuRV'] = (nuRV,"nu for RV")
        paramIntsStr = repr(paramInts).replace(' ','')
        self.settingsDict['parInts'] = (paramIntsStr,"Varried params")
        self.settingsDict['chainNum'] = (self.chainNum,"chain number")
        ## check priors are ok with range mins
        self.combinedPriors(rangeMins,rangeMins,True)
        
        return (realData,rangeMaxsRaw,rangeMinsRaw,rangeMaxs,rangeMins,sigmas,paramInts,nu,nuDI,nuRV)
    
    def combinedPriors(self,parsCurr,parsLast,test=False):
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
            if self.dictVal('ePrior'):
                #print 'ePrior'
                priorsRatio*=self.settingsDict['ePrior'][2](parsCurr[4],parsLast[4])
            if self.dictVal('pPrior'):
                #print 'pPrior'
                priorsRatio*=self.settingsDict['pPrior'][2](parsCurr[7],parsLast[7])
            if self.dictVal('incPrior'):
                #print 'incPrior'
                priorsRatio*=self.settingsDict['incPrior'][2](parsCurr[8],parsLast[8])
            if self.dictVal('M1Prior'):
                #print 'M1Prior'
                priorsRatio*=self.settingsDict['M1Prior'][2](parsCurr[0],parsLast[0])
                #print 'M1Prior'
            if self.dictVal('M2Prior'):
                #print 'M2Prior'
                priorsRatio*=self.settingsDict['M2Prior'][2](parsCurr[1],parsLast[1])
                #print 'M2Prior'
            if self.dictVal('parPrior'):
                #print 'parPrior'
                priorsRatio*=self.settingsDict['parPrior'][2](parsCurr[2],parsLast[2])
                #print 'parPrior out'
            if test==False:
                return priorsRatio
        except:
            self.log.critical("An error occured while trying to calculate the priors.")
            sys.exit(0)
            
    def dictVal(self,key):
        """
        Get the value for a key in the settings dictionary.
        This will handle the values that are tuples and not
        returning the value.
        """
        try:
            if type(self.settingsDict[key])==tuple:
                return self.settingsDict[key][0]
            else:
                return self.settingsDict[key]
        except:
            return False
    
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
        if  ('MCMC' not in stage) and (stage=='MC'):
            for i in range(0,len(pars)):
                if i in self.paramInts:
                    parsOut[i]=np.random.uniform(self.rangeMinsRaw[i],self.rangeMaxsRaw[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            self.parIntVaryAry.append(varyInt)
            sig = sigs[varyInt]*(self.rangeMaxsRaw[varyInt]-self.rangeMinsRaw[varyInt])
            parsOut[varyInt]=np.random.uniform(pars[varyInt]-sig,pars[varyInt]+sig)        
        ## if TcEqualT, push the varied one into the other
        if self.dictVal('TcEqualT'):
            if self.dictVal('TcStep'):
                parsOut[5]=parsOut[6]
            else:
                parsOut[6]=parsOut[5]
        ## if Kdirect not set, then inclination varys.
        ## then K=0 going into Orbit so Orbit will calc it
        if 8 in self.paramInts:
            parsOut[12] = 0
        return parsOut
    
    def rangeCheck(self,pars,sample,stage=''):
        """
        Check if values inside allowed range
        """
        inRange=True
        paramsOut = copy.deepcopy(pars)
        ## convert from Raw form if in lowEcc mode
        self.Orbit.convertParsFromRaw(paramsOut)
        for i in range(0,len(pars)):
            if i in self.paramInts:
                if (i==6)and(self.dictVal('TcStep')):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i==5)and(self.dictVal('TcStep')==False):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i!=5) and (i!=6):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
        if (sample>=10)and((self.acceptCount==0)and(stage=='SA')):
            ##Jump as starting position after first 10 tries was in poor part of param space. for SA only.
            paramsOut = self.increment(self.rangeMinsRaw,np.zeros(pars.shape),stage='MC')
            ## convert from Raw form if in lowEcc mode
            self.Orbit.convertParsFromRaw(paramsOut)
            inRange=True
        return (paramsOut,inRange)
    
    def accept(self,sample,pars,modelData,temp=1.0,stage=''):
        """
        First this will calculate chi squared for model vs real data.
        
        For mcOnly it performs simple chisquared cut-off acceptance 
        based on 'chiMAX' value in settingsDict.
        Else, it will calculate the priors and accept based on 
        the Metropolis-Hastings algorithm. The temp factor will 
        be set to 1.0 for MCMC and Sigma Tuning, and should be provided 
        for Simulated Annealing.
        """
        paramsOut = copy.deepcopy(pars)
        ## Calculate chi squareds for 3D,DI,RV and update bestPars and bestSumStr if this is better than the best
        (raw3D, reducedDI, reducedRV, reduced3D) = tools.chiSquaredCalc3D(self.realData,modelData,self.nuDI,self.nuRV,self.nu)
        paramsOut[11] = raw3D
        if self.bestSumStr=='':
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+' Nothing accepted yet below chi squared max = '+str(self.dictVal('chiMAX'))
            self.latestSumStr="Latest reduced chiSquared : [total,DI,RV] = ["+str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
        if (reduced3D)<self.bestRedChiSqr:
            self.bestRedChiSqr=(reduced3D)
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+\
            ' BEST reduced chiSquareds so far: [total,DI,RV] = ['\
            +str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
            bestPars = copy.deepcopy(paramsOut)
            ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
            self.Orbit.convertParsToRaw(bestPars)
            self.paramsBest = bestPars
        ## check if this step is accepted
        accept = False
        if (stage=='MC')or(stage=='SA'and(self.acceptCount==0)):
            ## for MC or first step of SA
            if (paramsOut[11]/self.nu)<self.dictVal('chiMAX'):
                accept=True
        else:
            ## For SA after first sample, MCMC, and ST
            try:
                likelihoodRatio = np.e**((self.paramsLast[11] - raw3D)/(2.0*temp))
                priorsRatio = self.combinedPriors(paramsOut,self.paramsLast)
                if np.random.uniform(0.0, 1.0)<=(priorsRatio*likelihoodRatio):
                    accept = True
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
            self.latestSumStr = stage+" chain #"+str(self.chainNum)+\
            ' Latest accepted reduced chiSquareds: [total,DI,RV] = ['+\
            str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]"
        else:
            self.acceptBoolAry.append(0)
        ##log a status summary?
        if self.dictVal('nSumry')>0:
            if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
                perc = sample*100//self.dictVal(self.stgNsampDict[stage])
                ####str(self.nSaved)+' (curr '+str(self.nSavedPeriodic)+"), Finished: "+str(sample)+"/"+\
                sumStr = "below\n"+stage+" chain #"+str(self.chainNum)+", # Accepted: "+str(self.acceptCount)+", # Saved: "+\
                str(self.nSaved)+", Finished: "+str(sample)+"/"+\
                str(self.dictVal(self.stgNsampDict[stage]))+" = "+str(perc)+"%, Current T: "+str(temp)+"\n"
                sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
                self.log.debug(sumStr)
        return (paramsOut,accept)
    
    def tempDrop(self,sample,strtTemp,temp,stage=''):
        """
        Determine if it is time to drop the temp, and drop if it is.
        Total temperature range is [strtTemp,0.01), so the minimum 
        temperature is actually <1.0 meaning the last few temperature drops 
        will really push the currently found minimum towards its peak.
        There will be a fixed number of temperature steps = 'nTmpStps'.
        """
        if stage=='SA':
            if sample%self.dictVal('tempInt')==0:
                temp-=(strtTemp-0.01)*(float(self.dictVal('tempInt'))/float(self.dictVal('nSAsamp')))
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
            if (sample%(len(self.paramInts)*self.dictVal('sigInt'))==0)and(self.acceptCount>1):
                self.acceptStr = '\n'+stage+" chain #"+str(self.chainNum)+'\n'
                self.shiftStr = ''
                self.parIntVaryAry = np.array(self.parIntVaryAry)
                self.acceptBoolAry = np.array(self.acceptBoolAry)
                self.acceptStr+="Number of steps used to calculate acceptance rate = "+repr(len(self.acceptBoolAry))+'\n'
                for i in self.paramInts:
                    ##calculate acceptance rate for each param
                    nAcc = np.sum(np.where(self.parIntVaryAry==i,self.acceptBoolAry,0))
                    nTot = len(np.where(self.parIntVaryAry==i)[0])
                    self.acceptStr+= 'parameter # '+str(i)+' acceptance = '+str(float(nAcc)/float(nTot))+'\n'
                    if stage=='ST':
                        ##check each rate to choose up/down shift and do so and update shiftStr
                        self.shiftStr+= '\n'+stage+" chain #"+str(self.chainNum)+'\nparameter # '+str(i)+" shifting sigma "+str(sigs[i])+" -> "
                        if ((float(nAcc)/float(nTot))>0.35)and(sigs[i]<self.dictVal('sigMax')):
                            sigmasOut[i]+=self.dictVal('sigMin')
                        elif ((float(nAcc)/float(nTot))<0.25)and(sigs[i]>self.dictVal('sigMin')):
                            sigmasOut[i]-=self.dictVal('sigMin')
                        self.shiftStr+=str(sigmasOut[i])+"\n"
                self.acceptBoolAry = []
                self.parIntVaryAry = []
                ##log a status summary?
                if self.dictVal('nSumry')>0:
                    if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
                        self.log.debug(self.acceptStr+self.shiftStr)
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
        sumStr+=str(float(self.acceptCount)/float(self.dictVal(self.stgNsampDict[stage])))+"\n"
        sumStr+= "Total number of steps stored = "+str(self.nSaved)+"\n"
        sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
        sumStr+="Last params  = "+repr(self.paramsLast)+'\n'
        sumStr+="Best params (Raw) = "+repr(self.paramsBest)+'\n'
        if (stage=="ST")or(stage=="MCMC"):
            if stage=='ST':
                sumStr+="Final Sigmas = "+repr(sigmas)+'\n'
            sumStr+=self.acceptStr+self.shiftStr
        sumStr+='\n'+'='*70+'\n'
        self.log.info(sumStr)
    
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
        self.settingsDict['chainNum'] = (self.chainNum,"chain number")
        # make a very random seed value to ensure each chain is different.  
        # Should we make this value an optional input and pass on as a return value to keep a process number using the same seed?? $$$
        # if so, it needs to be pushed into the results file as well.
        t = np.random.uniform(1,1e6)
        self.seed = int((timeit.default_timer()/(self.chainNum+1))/t)
        self.log.debug("Chain# "+str(self.chainNum)+" has random number seed = "+str(self.seed))
        np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.dictVal('finalFolder'),"tmpOutdata-"+str(self.chainNum)+".npy")
        if os.path.exists(self.tmpDataFile):
            os.remove(self.tmpDataFile)
            self.log.debug("just removed data file from disk:\n"+self.tmpDataFile)
    
    def startSummary(self,pars,sigs,stage):
        startStr=stage+" chain #"+str(self.chainNum)+' VALS AT START OF '+stage+' SIM:\n'
        startStr+= 'params = '+repr(pars)+'\n'
        startStr+= 'rangeMins = '+repr(self.rangeMins)+'\n'
        startStr+= 'rangeMaxs = '+repr(self.rangeMaxs)+'\n'
        if self.dictVal('lowEcc'):
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
        self.log.debug("Trying "+str(self.dictVal(self.stgNsampDict[stage]))+" samples for chain #"+str(chainNum)+" in "+stage+" mode.")
        self.chainNum = chainNum
        self.resetTracked(stage)
        bar = tools.ProgressBar('green',width=10,block='=',empty='-',lastblock='>')
        modelData = np.zeros((len(self.realData),3))
        acceptedParams = []
        self.settingsDict['curStg']=(stage,'Current stage either [SA,ST,MCMC or MC]')
        strtTemp = temp      
        sigmas = copy.deepcopy(self.starterSigmas)
        ## if valid startSigmas provided, start with them, else use defaults.
        if (type(startSigmas)==list)or(type(startSigmas)==np.ndarray):
            if len(startSigmas)>0:
                sigmas = copy.deepcopy(startSigmas)
        ## if valid startParams provided, start there, else start at random point.
        if (type(startParams)==list)or(type(startParams)==np.ndarray):
            if len(startParams)>0:
                paramsLast = copy.deepcopy(startParams)
                self.log.info('initial/latest pars have reduced chi sqr of '+str(paramsLast[11]/self.nu))
            else: 
                paramsLast = self.increment(self.rangeMinsRaw,sigmas,stage='MC')
        else: 
            paramsLast = self.increment(self.rangeMinsRaw,sigmas,stage='MC')           
        ## load up starting params as 'latest' and perform first increment from these to start loop with.
        latestParsRaw = copy.deepcopy(paramsLast)
        proposedParsRaw = self.increment(latestParsRaw,sigmas,stage)
        self.acceptBoolAry.append(0)
        #print 'proposed pars have reduced chi sqr of '+str(proposedParsRaw[11]/self.nu)#$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ## convert from Raw form if in lowEcc mode
        self.Orbit.convertParsFromRaw(paramsLast)
        self.paramsLast = paramsLast
        #self.startSummary(proposedPars,sigmas,stage)
        ##loop through each sample 
        ##Follows these steps: in Range?,calc model,accept?,Store?,increment,lower temp?,tune sigmas? dump data to disk? DONE ->write output data
        sample=0
        while sample<(self.dictVal(self.stgNsampDict[stage])+1):
            sample+=1
            (proposedPars,inRange)=self.rangeCheck(proposedParsRaw,sample,stage)
            if inRange:
                self.Orbit.calculate(modelData,proposedPars)
                (params,accept) = self.accept(sample,proposedPars,modelData,temp,stage)
                if accept:
                    latestParsRaw = copy.deepcopy(params)
                    ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
                    self.Orbit.convertParsToRaw(latestParsRaw)
                    if ('MCMC' not in stage)and(stage=='MC'):
                        acceptedParams.append(params)
                        self.nSaved+=1
                        self.nSavedPeriodic+=1  
                    elif (self.acceptCount%self.dictVal('saveInt'))==0:
                        acceptedParams.append(params)  
                        self.nSaved+=1   
                        self.nSavedPeriodic+=1                
            else:
                self.acceptBoolAry.append(0)
            proposedParsRaw = self.increment(latestParsRaw,sigmas,stage)
            temp = self.tempDrop(sample,strtTemp,temp,stage)
            sigmas = self.sigTune(sample,sigmas,stage)
            if (self.nSavedPeriodic>0)and((self.nSaved%self.dictVal('dmpInt'))==0):
                ## dump acceptedParams array to disk and collect garbage
                self.log.debug('Dumping data to filename:\n'+self.tmpDataFile+\
                               '\nThe acceptedParams Ary had '+str(np.array(acceptedParams).shape[0])+\
                               ' param sets, and size = '+str(np.array(acceptedParams).nbytes/1.0e6)+' MB')
                tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
                acceptedParams = []
                self.nSavedPeriodic = 0
                self.log.debug('about to collect the garbage')
                gc.collect()
            if sample%(self.dictVal(self.stgNsampDict[stage])//100)==0:
                timesAry.append(timeit.default_timer()-lastTic)
                lastTic = timeit.default_timer()
            if (self.dictVal('logLevel')<30)and(sample%(self.dictVal(self.stgNsampDict[stage])//100)==0):
                timeRemSec = np.mean(timesAry)*(100.0-(float(sample)*100.0)/float(self.dictVal(self.stgNsampDict[stage])))
                timeStr = ' about '+tools.timeStrMaker(timeRemSec)+' remaining.'
                bar.render(sample*100//self.dictVal(self.stgNsampDict[stage]), stage+str(chainNum)+' Completed,'+timeStr)
        if self.dictVal('logLevel')<30:
            bar.render(100,stage+str(chainNum)+' Complete!\n')
        self.log.debug(stage+" took: "+tools.timeStrMaker(timeit.default_timer()-tic))
        self.endSummary(temp,sigmas,stage)
        tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
        outFname = tools.writeFits('outputData'+stage+str(chainNum)+'.fits',self.tmpDataFile,self.settingsDict)
        if stage=='SA':
            ##start ST at the best location with tight sigmas, and it will tune to ideal sigmas
            sigmas = np.ones(np.array(sigmas).shape)*self.dictVal('sigMin')   
        return (outFname,self.paramsBest,sigmas,self.bestRedChiSqr)
#END OF FILE      
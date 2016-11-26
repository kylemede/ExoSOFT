from __future__ import absolute_import
from __future__ import print_function
import KMlogger
from . import generalTools as genTools
from . import readWriteTools as rwTools
from multiprocessing import Process
import pickle
import timeit
import os
import numpy as np
from six.moves import range

log = KMlogger.getLogger('main.chainTools',lvl=100,addFH=False) 

class singleProc(Process):
    """
    This is the Manager object that controls the a single processes for a 
    ExoSOFT simulation run.  It is called by the multiProcessStarter once for 
    each chain/process requested by the user through the simulation settings 
    file.
    
    :param str settings: settings Dictionary
    :param str fNameBase: File name, including the full path, for the output 
        data files.
    :param list stageList: List of stages to run ex.['MC','SA','ST','MCMC'] 
        lives.
    :param int chainNum: number of this chain
    """
    def __init__(self, settings, SimObj, stage, chainNum, pklFilename = '', params=[],sigmas=[],strtTemp=1.0):
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = log
        self.settings = settings 
        self.stage = stage
        self.params = params
        self.sigmas = sigmas
        self.strtTemp = strtTemp
        self.Sim = SimObj
        self.pklFilename = pklFilename
        
    def run(self):        
        ## run the requested stage and [ickle its return values
        self.log.debug('Starting to run process #'+str(self.chainNum))
        (outFname,params,sigmas,bestRedChiSqr,avgAcceptRate,acceptStr) = self.Sim.simulatorFunc(self.stage,self.chainNum,self.params,self.sigmas,self.strtTemp)
        self.log.debug('chain #'+str(self.chainNum)+" of "+self.stage+' stage  OUTFILE :\n'+outFname)
        pickle.dump([outFname,params,sigmas,bestRedChiSqr,avgAcceptRate,acceptStr], open(self.pklFilename,'wb'))
        
    def loadResult(self):
        return pickle.load(open(self.pklFilename,'rb'))

class multiProcObjResults(object):
    def __init__(self,bestRedChiSqrs,avgAcceptRates,acceptStrs,stage,retStr,latestRetStr):
        self.bestRedChiSqrs = bestRedChiSqrs
        self.avgAcceptRates = avgAcceptRates
        self.acceptStrs = acceptStrs
        self.stage = stage
        self.retStr = retStr
        self.latestRetStr = latestRetStr
        
class multiProcObj(object):
    def __init__(self,settings,Sim,stage):
        self.outFnames = []
        self.params = []
        self.sigmas = []
        self.bestRedChiSqrs = []
        self.avgAcceptRates = []
        self.acceptStrs = []
        self.settings = settings
        self.Sim = Sim
        self.stage = stage
        if stage=="emcee":
            self.numProcs = 1
        elif stage=='MCMC':
            self.numProcs = settings['nMCMCcns']
        else:
            self.numProcs = settings['nChains']
        self.retStr = ''
        self.latestRetStr = ''
        
    def _loadUpArys(self,master):
        for procNumber in range(len(master)):
            ret = master[procNumber].loadResult()
            if os.path.exists(ret[0]):
                self.outFnames.append(ret[0])
                self.params.append(ret[1])
                self.sigmas.append(ret[2])
                self.bestRedChiSqrs.append(ret[3])
                self.avgAcceptRates.append(ret[4])
                self.acceptStrs.append(ret[5])
            elif ret[0]!='':
                log.error("Resulting file from MPO.run does not exist:\n"+\
                          ret[0]+'\nit had a reduced chi of '+str(ret[3])+'\n')
        self.sortResults()
        
    def resultsOnly(self):
        accRts = self.avgAcceptRates
        stg = self.stage
        chis = self.bestRedChiSqrs
        accStrs = self.acceptStrs
        retStr = self.retStr
        latRetStr = self.latestRetStr 
        resObj = multiProcObjResults(chis,accRts,accStrs,stg,retStr,latRetStr)
        return resObj
        
    def sortResults(self):
        """
        Sort all result arrays by ascending reduced chi squared values.
        """
        if len(self.bestRedChiSqrs)>0:
            chis = np.array(self.bestRedChiSqrs)
            chisSorted = np.sort(self.bestRedChiSqrs)
            outFnames = []
            params = []
            sigmas = []
            bestRedChiSqrs = []
            avgAcceptRates = []
            acceptStrs = []
            for chi in chisSorted:
                if chi not in bestRedChiSqrs:
                    try:
                        goodPts = np.where(chis==chi)[0].tolist()
                        outFnames.append(self.outFnames[goodPts[0]])
                        params.append(self.params[goodPts[0]])
                        sigmas.append(self.sigmas[goodPts[0]])
                        bestRedChiSqrs.append(chi)
                        avgAcceptRates.append(self.avgAcceptRates[goodPts[0]])
                        acceptStrs.append(self.acceptStrs[goodPts[0]])
                    except:
                        log.error("A problem occured while trying to sort "+\
                                  "results.  For chi = "+str(chi)+" the "+\
                                  "resulting gootPts were "+repr(goodPts))
                else:
                    log.debug("already a chain with a reduced chi of "+str()+\
                              " in the sorted MPO array.")
            self.outFnames = outFnames
            self.params = params
            self.sigmas = sigmas
            self.bestRedChiSqrs = bestRedChiSqrs
            self.avgAcceptRates = avgAcceptRates
            self.acceptStrs = acceptStrs
    
    def getTopProcs(self,maxRedChiSqr,fillToNumProc=False,nProcs=None,allBest=False):
        if nProcs==None:
            nProcs = self.numProcs
        (params,sigmas,chis,outFnames,avgAcceptRates,acceptStrs)=self.findGoodBadOnes(maxRedChiSqr) 
        if chis!=None:        
            if allBest:
                params = [params[0]]
                sigmas = [sigmas[0]]
                chis = [chis[0]]
                outFnames = [outFnames[0]]
                while len(params)<nProcs:
                    params.append(params[0])
                    sigmas.append(sigmas[0])
                    chis.append(chis[0])
                    outFnames.append(outFnames[0])
            elif (len(params)<nProcs) and fillToNumProc:
                while len(params)<nProcs:
                    rndVal = np.random.randint(0,len(params))
                    params.append(params[rndVal])
                    sigmas.append(sigmas[rndVal])
                    chis.append(chis[rndVal])
                    outFnames.append(outFnames[rndVal])
            elif len(params)>nProcs:
                params = params[:nProcs]
                sigmas = sigmas[:nProcs]
                chis = chis[:nProcs]
                outFnames = outFnames[:nProcs]
        return (params,sigmas,chis,outFnames)      
    
    def findGoodBadOnes(self,maxRedChiSqr,findGoodOnes=True):
        self.sortResults()
        outFnames = []
        params = []
        sigmas = []
        chis = []
        avgAcceptRates = []
        acceptStrs = []
        chis2 = np.array(self.bestRedChiSqrs)
        if findGoodOnes:
            goodBads = np.where(chis2<maxRedChiSqr)[0].tolist()
        else:
            goodBads = np.where(chis2>maxRedChiSqr)[0].tolist()
        if len(goodBads)>0:            
            for gb in goodBads:
                outFnames.append(self.outFnames[gb])
                params.append(self.params[gb])
                sigmas.append(self.sigmas[gb])
                chis.append(self.bestRedChiSqrs[gb])
                avgAcceptRates.append(self.avgAcceptRates[gb])
                acceptStrs.append(self.acceptStrs[gb])
        else:
            params=sigmas=chis=outFnames=avgAcceptRates=acceptStrs=None
        return (params,sigmas,chis,outFnames,avgAcceptRates,acceptStrs)
    
    def killBadOnes(self,maxRedChiSqr,limitToNumProcs=True,nProcs=None):
        if nProcs==None:
            nProcs = self.numProcs
        (paramsA,sigmasA,chisA,outFnamesA,avgAcceptRatesA,acceptStrsA)=self.findGoodBadOnes(maxRedChiSqr) 
        (paramsB,sigmasB,chisB,outFnamesB,avgAcceptRatesB,acceptStrsB)=self.findGoodBadOnes(maxRedChiSqr,findGoodOnes=False)
        rwTools.rmFiles(outFnamesB) 
        if paramsA!=None:
            if limitToNumProcs and (len(paramsA)>nProcs):
                self.outFnames = outFnamesA[:nProcs]
                self.bestRedChiSqrs = chisA[:nProcs]
                self.params = paramsA[:nProcs]
                self.sigmas = sigmasA[:nProcs]
                self.avgAcceptRates = avgAcceptRatesA[:nProcs]
                self.acceptStrs = acceptStrsA[:nProcs]
            else:
                self.outFnames = outFnamesA
                self.bestRedChiSqrs = chisA
                self.params = paramsA
                self.sigmas = sigmasA
                self.avgAcceptRates = avgAcceptRatesA
                self.acceptStrs = acceptStrsA
        else:
            self.outFnames = []
            self.params = []
            self.sigmas = []
            self.bestRedChiSqrs = []
            self.avgAcceptRates = []
            self.acceptStrs = []
            
    def writeBest(self):
        try:
            (bstChi,bstInt) = self._best()
            rwTools.writeBestsFile(self.settings,self.params[bstInt],self.sigmas[bstInt],bstChi,self.stage)
            return True
        except:
            log.critical("No parameters were accepted!!!")
            return False
        
    def _best(self):
        if len(self.acceptStrs)>0:
            bstChi = 1e9
            bstInt = 0
            for i in range(len(self.acceptStrs)):
                if self.bestRedChiSqrs[i]<bstChi:
                    bstChi = self.bestRedChiSqrs[i]
                    bstInt = i
            return (bstChi,bstInt)
        else:
            log.debug("No chains in this MPO yet.")
    def getBest(self):
        (bstChi,bstInt) = self._best()
        return (self.outFnames[bstInt],self.params[bstInt],self.sigmas[bstInt],bstChi,self.avgAcceptRates[bstInt],self.acceptStrs[bstInt])
            
    def run(self,params=[],sigmas=[],strtTemp=1.0):
        """
        The function to run multiple processes of the simulator and absorb the 
        results into the object.
        """
        if self.stage != 'MC':
            if len(params)==0:
                log.error("params and sigmas parameter inputs to multiproc must have a set of values for each requested proc to start for SA, ST and MCMC modes!")
        else:
            #load up basically empty arrays for params and sigmas for MC case
            sigmas = params = list(range(self.numProcs))
        master = []
        tic=timeit.default_timer()
        extra = ''
        if self.stage=='SA':
            extra+=" with a starting temperature of "+str(strtTemp)
        log.info("Going to start "+str(self.numProcs)+" chains for the "+self.stage+" stage"+extra)
        for procNumber in range(self.numProcs):
            pklFilename = os.path.join(self.settings['finalFolder'],'pklTemp'+"-"+self.stage+'-'+str(procNumber)+".p")
            master.append(singleProc(self.settings,self.Sim,self.stage,procNumber,pklFilename=pklFilename,params=params[procNumber],sigmas=sigmas[procNumber],strtTemp=strtTemp))
            master[procNumber].start()
        for procNumber in range(self.numProcs):
            master[procNumber].join()  
        toc=timeit.default_timer()
        s = "ALL "+str(self.numProcs)+" chains of the "+self.stage+" stage took a total of "+genTools.timeStrMaker(int(toc-tic))
        retStr =s+"\n"
        log.info(s)
        self._loadUpArys(master)
        self.latestRetStr = retStr
        ############################################################
        ## Note to me:
        ## a Cool way to start up an interactive terminal to see the value of all current local variables
        ## great for devel and debugging.
        ############################################################
        #import code; code.interact(local=locals())
        ############################################################   

def iterativeSA(settings,Sim,internalTemp=None):
    """
    Perform SA with multiProc nSAiters times, droping the starting temperature each time by strtTemp/nSAiters.
    """
    testing=False
    tic=timeit.default_timer()
    numProcs = settings['nChains']
    nSAiters = 7
    max_nSAiters = 10
    nSAstrtIters = 3
    chisForCalc = None
    strtPars = list(range(numProcs))
    strtsigmas = list(range(numProcs))
    SAmultiProc = multiProcObj(settings,Sim,'SA')
    uSTD = 1e6
    iter = -1
    if internalTemp==None:
        internalTemp = settings['strtTemp']
    temp = internalTemp
    while uSTD>settings['maxUstd']:
        iter+=1
        if iter>0:
            if (iter>nSAstrtIters)and(chisForCalc==None):
                ## Try SA again from the start
                rwTools.rmFiles(SAmultiProc.outFnames)
                s = "3 itterations of SA failed to find a suitable starting point."
                s+= " So starting it again at double the initial temperature."
                log.critical(s)
                iterativeSA(settings,Sim,internalTemp=internalTemp*2.0)
            elif iter<=nSAiters:
                temp -= internalTemp/nSAiters
        if iter>max_nSAiters:
            s = "iterativeSA took more than "+str(max_nSAiters)+", so process was terminated!!!"
            log.raisemsg(s)
            raise  RuntimeError('\n\n'+s)
            
        log.info("\nIteration #"+str(iter+1))
        SAmultiProc.retStr +="Iteration #"+str(iter+1)+"\n"
        ## run multiProc for this temperature, kill bad chains 
        ## and get best as start positions of next iteration.
        SAmultiProc.run(params=strtPars,sigmas=strtsigmas,strtTemp=temp)
        SAmultiProc.retStr +=SAmultiProc.latestRetStr
        SAmultiProc.killBadOnes(settings['chiMaxST'])
        (pars,sigs,chisForCalc,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST']) 
        #print 'so far top chis are: '+repr(chisForCalc) 
        if testing:
            print('so far top chis are: '+repr(chisForCalc))   
            if outFnms!=None:     
                print('matching file names are: ')
                for n in outFnms:
                    print(n)
        if pars==None:
            strtPars = list(range(numProcs))
        else:
            ## shift through the names of the top procs and rename 
            ## their files to new temp names.
            # NOTE: This code is a bit of overkill and can trim down after it 
            #       proves working for a month or so (dated:April 20 2016).
            for i in range(len(outFnms)):
                worked=False
                try:
                    curNm = outFnms[i]
                    outNm = os.path.join(os.path.dirname(curNm),"SAtempData-"+str(iter)+"-"+str(i)+".fits")
                    rwTools.renameFits(curNm,outNm,killInput=True)
                    worked=True
                except:
                    log.critical("\nfailed to rename file during iteration "+str(iter)+" of SA")
                    s="\ncurNm = "+curNm+"\noutNm = "+outNm+"\nSAmultiProc.outFnames = "
                    for n in SAmultiProc.outFnames:
                        s+="\n"+n
                    log.critical(s)
                if worked:
                    for j in range(len(SAmultiProc.outFnames)):
                        if SAmultiProc.outFnames[j]==curNm:
                            SAmultiProc.outFnames[j] = outNm
                            if testing:
                                print('outfile renamed to '+SAmultiProc.outFnames[j])
            (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'],fillToNumProc=True)
            strtPars = pars
            ## Calc 'uniform' STD and make final msgs for this iteration          
            log.debug(str(len(strtPars))+' sets of starting parameters being passed to next iteration.')
            if len(chisForCalc)==numProcs:
                uSTD = genTools.unitlessSTD(chisForCalc)
            log.importantinfo("After iteration #"+str(iter+1)+" the top "+str(len(chisForCalc))+" solutions with reduced chiSquared < "+str(settings['chiMaxST'])+" have a unitless STD of "+str(uSTD))
            SAmultiProc.retStr +="The latest top "+str(len(chisForCalc))+" reduced chiSquareds had a unitless STD of "+str(uSTD)+'\n'
    #############
    ## wrap up ##
    #############
    (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'])   
    #$$$$$$$$$$$ Debug $$$$$$$$$$$$$$$$
    if testing:
        (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'])
        print('\nFINAL before renaming\nchis = '+repr(chis)+'\noutFnms = ')
        for n in SAmultiProc.outFnames:
            print(n)
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #rename final data files to standard SA convention
    if len(outFnms)>1:
        for i in range(len(outFnms)):
            curNm = outFnms[i]
            outNm = os.path.join(os.path.dirname(curNm),'outputDataSA'+str(i)+'.fits')
            rwTools.renameFits(curNm,outNm)
            for j in range(len(SAmultiProc.outFnames)):
                if SAmultiProc.outFnames[j]==curNm:
                    SAmultiProc.outFnames[j] = outNm    
    #$$$$$$$$$$$ Debug $$$$$$$$$$$$$$$$
    if testing:
        (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'])
        print('\nFINAL after renaming\noutFnms = ')
        for n in SAmultiProc.outFnames:
            print(n)
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
    SAmultiProc.writeBest()
    toc=timeit.default_timer()
    s = "ALL "+str(iter+1)+" iterations of SA took a total of "+genTools.timeStrMaker(int(toc-tic))
    SAmultiProc.retStr +=s+"\n"
    log.warning(s)
    return SAmultiProc            

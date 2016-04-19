import exoSOFTlogger
import generalTools as genTools
import readWriteTools as rwTools
from multiprocessing import Process
import pickle
import timeit
import os
import numpy as np

log = exoSOFTlogger.getLogger('main.chainTools',lvl=100,addFH=False) 

class singleProc(Process):
    """
    This is the Manager object that controls the a single processes for a 
    exosoft simulation run.  It is called by the multiProcessStarter once for 
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
        if stage=='MCMC':
            self.numProcs = settings['nMCMCcns']
        else:
            self.numProcs = settings['nChains']
        self.retStr = ''
        self.latestRetStr = ''
        
    def _loadUpArys(self,master):
        for procNumber in range(len(master)):
            ret = master[procNumber].loadResult()
            self.outFnames.append(ret[0])
            self.params.append(ret[1])
            self.sigmas.append(ret[2])
            self.bestRedChiSqrs.append(ret[3])
            self.avgAcceptRates.append(ret[4])
            self.acceptStrs.append(ret[5])
        self.sortResults()
        
    def sortResults(self):
        """
        Sort all result arrays by ascending reduced chi squared values.
        """
        chis = np.array(self.bestRedChiSqrs)
        chisSorted = np.sort(self.bestRedChiSqrs)
        outFnames = []
        params = []
        sigmas = []
        bestRedChiSqrs = []
        avgAcceptRates = []
        acceptStrs = []
        for chi in chisSorted:
            try:
                goodPts = np.where(chis==chi)[0].tolist()
                outFnames.append(self.outFnames[goodPts[0]])
                params.append(self.params[goodPts[0]])
                sigmas.append(self.sigmas[goodPts[0]])
                bestRedChiSqrs.append(chi)
                avgAcceptRates.append(self.avgAcceptRates[goodPts[0]])
                acceptStrs.append(self.acceptStrs[goodPts[0]])
            except:
                log.error("A problem occured while trying to sort results.  For chi="+\
                          str(chi)+" the resulting gootPts were "+repr(goodPts))
        self.outFnames = outFnames
        self.params = params
        self.sigmas = sigmas
        self.bestRedChiSqrs = bestRedChiSqrs
        self.avgAcceptRates = avgAcceptRates
        self.acceptStrs = acceptStrs
    
    def getTopProcs(self,maxRedChiSqr,killBadOnes=True,fillToNumProc=False,nProcs=None):
        self.sortResults()
        if nProcs==None:
            nProcs = self.numProcs
        outFnames = []
        params = []
        sigmas = []
        chis = []
        avgAcceptRates = []
        acceptStrs = []
        goods = np.where(np.array(self.bestRedChiSqrs)<maxRedChiSqr)[0].tolist()
        if len(goods)>0:
            for good in goods:
                outFnames.append(self.outFnames[good])
                params.append(self.params[good])
                sigmas.append(self.sigmas[good])
                chis.append(self.bestRedChiSqrs[good])
                avgAcceptRates.append(self.avgAcceptRates[good])
                acceptStrs.append(self.acceptStrs[good])
            if (len(params)<nProcs) and fillToNumProc:
                while len(params)<nProcs:
                    rndVal = np.random.randint(0,len(params))
                    params.append(params[rndVal])
                    sigmas.append(sigmas[rndVal])
                    chis.append(chis[rndVal])
            elif len(params)>nProcs:
                params = params[:nProcs]
                sigmas = sigmas[:nProcs]
                chis = chis[:nProcs]
                if killBadOnes:
                    rwTools.rmFiles(outFnames[nProcs:])
                    self.outFnames = outFnames[:nProcs]
                    self.bestRedChiSqrs = chis
                    self.params = params
                    self.sigmas = sigmas
                    self.avgAcceptRates = avgAcceptRates[:nProcs]
                    self.acceptStrs = acceptStrs[:nProcs]
                outFnames = outFnames[:nProcs]
        else:
            params=sigmas=chis=outFnames=None
        return (params,sigmas,chis,outFnames)
    
    def writeBest(self):
        (bstChi,bstInt) = self._best()
        rwTools.writeBestsFile(self.settings,self.params[bstInt],self.sigmas[bstInt],bstChi,self.stage)
        
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
            sigmas = params = range(self.numProcs)
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

def iterativeSA(settings,Sim):
    """
    Perform SA with multiProc nSAiters times, droping the starting temperature each time by strtTemp/nSAiters.
    """
    tic=timeit.default_timer()
    numProcs = settings['nChains']
    nSAiters = 7.0
    strtPars = range(numProcs)
    strtsigmas = range(numProcs)
    SAmultiProc = multiProcObj(settings,Sim,'SA')
    uSTD = 1e6
    iter = -1
    temp = settings['strtTemp']
    while uSTD>settings['maxUstd']:
        iter+=1
        if iter>0:
            if iter<=nSAiters:
                temp -= settings['strtTemp']/nSAiters
            elif iter>nSAiters:
                ## Try SA again from the start
                rwTools.rmFiles(SAmultiProc.outFnames)
                log.critical("Nothing found on this round of iterativeSA, so trying again from the top.")
                iterativeSA(settings,Sim)
        log.info("\nIteration #"+str(iter+1))
        SAmultiProc.retStr +="Iteration #"+str(iter+1)+"\n"
        ## run multiProc
        SAmultiProc.run(params=strtPars,sigmas=strtsigmas,strtTemp=temp)
        SAmultiProc.retStr +=SAmultiProc.latestRetStr
        (pars,sigs,chisForCalc,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'],killBadOnes=False)            
        if pars==None:
            strtPars = range(numProcs)
        else:
            ## sift through the names of the top procs and rename their files to new temp names.
            if False:
                print '\n0\nchis = '+repr(chisForCalc)+'\noutFnms = '+repr(outFnms)
            for i in range(len(outFnms)):
                curNm = outFnms[i]
                outNm = os.path.join(os.path.dirname(curNm),"SAtempData-"+str(iter)+"-"+str(i)+".fits")
                rwTools.renameFits(curNm,outNm,killInput=True)
                for j in range(len(SAmultiProc.outFnames)):
                    if SAmultiProc.outFnames[j]==curNm:
                        SAmultiProc.outFnames[j] = outNm
            #print '\nSAmultiProc.outFnames = '+repr(SAmultiProc.outFnames)
            (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'],killBadOnes=True,fillToNumProc=True)
            strtPars = pars
            ## Calc 'uniform' STD and make final msgs for this iteration          
            log.debug(str(len(strtPars))+' sets of starting parameters being passed to next iteration.')
            if len(chisForCalc)==numProcs:
                uSTD = genTools.unitlessSTD(chisForCalc)
            log.info("After iteration #"+str(iter+1)+" the top "+str(len(chisForCalc))+" solutions with reduced chiSquared < "+str(settings['chiMaxST'])+" have a unitless STD of "+str(uSTD))
            SAmultiProc.retStr +="The latest top "+str(len(chisForCalc))+" reduced chiSquareds had a unitless STD of "+str(uSTD)+'\n'
    ############
    ## wrap up #
    ############
    (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'],killBadOnes=False)   
    if False: 
        print '\n1\nchis = '+repr(chis)+'\noutFnms = '+repr(outFnms)
    #rename final data files to standard SA convention
    if len(outFnms)>1:
        for i in range(len(outFnms)):
            curNm = outFnms[i]
            outNm = os.path.join(os.path.dirname(curNm),'outputDataSA'+str(i)+'.fits')
            rwTools.renameFits(curNm,outNm)
            for j in range(len(SAmultiProc.outFnames)):
                if SAmultiProc.outFnames[j]==curNm:
                    SAmultiProc.outFnames[j] = outNm    
    #Debug
    if False:
        (pars,sigs,chis,outFnms) = SAmultiProc.getTopProcs(settings['chiMaxST'],killBadOnes=False)
        print '\n2\nchis = '+repr(chis)+'\noutFnms = '+repr(outFnms)
        
    SAmultiProc.writeBest()
    toc=timeit.default_timer()
    s = "ALL "+str(iter+1)+" iterations of SA took a total of "+genTools.timeStrMaker(int(toc-tic))
    SAmultiProc.retStr +=s+"\n"
    log.warning(s)
    return SAmultiProc            
                        
                        
                                
        
#         if topProcs!=None:
#             #print 'filtering'
#             chisSorted = [] 
#             goodParams = []           
#             #Filter inputs if more than max num MCMC proc available to use the best ones
#             chisSorted = np.sort(retAry[3])
#             chisSorted = chisSorted[np.where(chisSorted<settings['chiMaxST'])]
#             if len(chisSorted)==0:
#                 strtPars = range(numProcs)
#             elif (len(chisSorted)==1)and(numProcs==1):
#                 bestRetAry=retAry
#                 uSTD=1e-6
#             else:
#                 #first updated bestRetAry
#                 for i in range(len(retAry[0])):
#                     if retAry[3][i] in chisSorted:                       
#                         if len(bestRetAry[0])<numProcs:
#                             for j in range(6):
#                                 bestRetAry[j].append(retAry[j][i])
#                         else:
#                             bestChis = np.sort(bestRetAry[3])
#                             if retAry[3][i]<bestChis[-1]:
#                                 for j in range(6):
#                                     bestRetAry[j].append(retAry[j][i])
#                 #now make list of best ones to use in next round
#                 if len(bestRetAry[0])>numProcs:
#                     bestChis = np.sort(bestRetAry[3])
#                     log.debug("len best before filtering = "+str(len(bestRetAry[0]))+", worst was "+str(bestChis[-1]))
#                     bestRetAry2 = [[],[],[],[],[],[]]
#                     #trim best lists down to size
#                     for i in range(0,len(bestChis)):
#                         if bestRetAry[3][i] in bestChis[:numProcs]:
#                             for j in range(6):
#                                 bestRetAry2[j].append(bestRetAry[j][i])
#                     bestRetAry = bestRetAry2
#                 #copy resulting data files to new temp names.
#                 for i in range(0,len(bestRetAry[0])):
#                     curNm = bestRetAry[0][i]
#                     outNm = os.path.join(os.path.dirname(curNm),"SAtempData-"+str(iter)+"-"+str(i)+".fits")
#                     rwTools.renameFits(curNm,outNm,killInput=True)
#                     bestRetAry[0][i] = outNm
#                 #kill those that were not good enough
#                 for i in range(len(retAry[0])):
#                     if os.path.exists(retAry[0][i]):
#                         rwTools.rmFiles([retAry[0][i]])
#                 ## Now fill out an array of starting parameter sets from the best above.
#                 ## first load up with one set of goodParams, then randomly from it till full.
#                 log.debug('best chis:\n' +repr(np.sort(bestRetAry[3])))
#                 goodParams = bestRetAry[1]
#                 strtPars=[]
#                 if len(goodParams)>0:
#                     for i in range(0,len(goodParams)):
#                         strtPars.append(goodParams[i])
#                     while len(strtPars)<numProcs:
#                         rndVal = np.random.randint(0,len(goodParams))
#                         strtPars.append(goodParams[rndVal])
#                 log.info(str(len(strtPars))+' sets of starting parameters being passed to next iteration.')
#                 #print 'STD = '+str(np.std(bestRetAry[3]))
#                 if len(bestRetAry[3])==numProcs:
#                     uSTD = genTools.unitlessSTD(bestRetAry[3])
#                 log.info("After iteration #"+str(iter+1)+" the top "+str(len(bestRetAry[3]))+" solutions with reduced chiSquared < "+str(settings['chiMaxST'])+" have a unitless STD of "+str(uSTD))
#                 retStr2 +="The latest top "+str(len(bestRetAry[3]))+" reduced chiSquareds had a unitless STD of "+str(uSTD)+'\n'
#     ## wrap up
#     if len(bestRetAry[0])>1:
#         #rename final data files to standard SA convention
#         for i in range(0,len(bestRetAry[0])):
#             curNm = bestRetAry[0][i]
#             outNm = os.path.join(os.path.dirname(curNm),'outputDataSA'+str(i)+'.fits')
#             rwTools.renameFits(curNm,outNm)
#             bestRetAry[0][i] = outNm
#     #Find best fit, write to file, maybe push into settings files if better than one in there already.
#     if len(bestRetAry[0])>0:
#         bstChiSqr = np.sort(bestRetAry[3])[0]
#         for i in range(len(bestRetAry[0])):
#             if bestRetAry[3][i] == bstChiSqr:
#                 bestpars = bestRetAry[1][i]
#                 bestsigs = bestRetAry[2][i]
#         rwTools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'SA')
#     toc=timeit.default_timer()
#     s = "ALL "+str(iter+1)+" iterations of SA took a total of "+genTools.timeStrMaker(int(toc-tic))
#     retStr2 +=s+"\n"
#     log.warning(s)
#     return (bestRetAry,retStr2)
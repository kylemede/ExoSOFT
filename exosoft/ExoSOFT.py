#!/usr/bin/python python
#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import matplotlib
# Force matplotlib to not use any Xwindows backend, to further avoid Display issues or when ExoSOFT is ran through ssh without -X.
matplotlib.use('Agg') 
import tools
import simulator
from exosoftpath import rootDir as ExoSOFTdir
import sys
import os
import timeit
import numpy as np
from multiprocessing import Process
import pickle

"""
    This is the 'main' of exosoft. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
""" 
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
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
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
        
def multiProc(settings,Sim,stage,numProcs,params=[],sigmas=[],strtTemp=1.0):
    log = tools.getLogger('main.multiProc',lvl=100,addFH=False)
    if stage != 'MC':
        if len(params)==0:
            log.error("params and sigmas parameter inputs to multiproc must have a set of values for each requested proc to start for SA, ST and MCMC modes!")
    else:
        #load up basically empty arrays for params and sigmas for MC case
        sigmas = params = range(numProcs)
    master = []
    tic=timeit.default_timer()
    extra = ''
    if stage=='SA':
        extra+=" with a starting temperature of "+str(strtTemp)
    log.info("Going to start "+str(numProcs)+" chains for the "+stage+" stage"+extra)
    for procNumber in range(numProcs):
        pklFilename = os.path.join(settings['finalFolder'],'pklTemp'+"-"+stage+'-'+str(procNumber)+".p")
        master.append(singleProc(settings,Sim,stage,procNumber,pklFilename=pklFilename,params=params[procNumber],sigmas=sigmas[procNumber],strtTemp=strtTemp))
        master[procNumber].start()
    for procNumber in range(numProcs):
        master[procNumber].join()    
    toc=timeit.default_timer()
    s = "ALL "+str(numProcs)+" chains of the "+stage+" stage took a total of "+tools.timeStrMaker(int(toc-tic))
    retStr =s+"\n"
    log.info(s)
    retAry = [[],[],[],[],[],[]]
    for procNumber in range(numProcs):
        ret = pickle.load(open(master[procNumber].pklFilename,'rb'))
        for i in range(6):
            retAry[i].append(ret[i])
    return (retAry,retStr)

def iterativeSA(settings,Sim):
    """
    Perform SA with multiProc nSAiters times, droping the starting temperature each time by strtTemp/nSAiters.
    """
    tic=timeit.default_timer()
    log = tools.getLogger('main.iterativeSA',lvl=100,addFH=False)
    maxNumMCMCprocs = settings['nMCMCcns']
    numProcs = settings['nChains']
    nSAiters = 7.0
    strtPars = range(numProcs)
    strtsigmas = range(numProcs)
    bestRetAry = [[],[],[],[],[],[]]
    uSTD = 1e6
    iter = -1
    retStr2 = ''
    temp = settings['strtTemp']
    while uSTD>settings['maxUstd']:
        iter+=1
        if iter>0:
            #clean up previous SA data files on disk to avoid clash
            #tools.rmFiles(bestRetAry[0][:])
            if iter<nSAiters:
                temp -= settings['strtTemp']/nSAiters
        log.info("\nIteration #"+str(iter+1))
        retStr2 +="Iteration #"+str(iter+1)+"\n"
        (retAry,retStr) = multiProc(settings,Sim,'SA',numProcs,params=strtPars,sigmas=strtsigmas,strtTemp=temp)
        #print '\n'*5
        retStr2 +=retStr
        if len(retAry)>0:
            #print 'filtering'
            chisSorted = [] 
            goodParams = []           
            #Filter inputs if more than max num MCMC proc available to use the best ones
            chisSorted = np.sort(retAry[3])
            chisSorted = chisSorted[np.where(chisSorted<settings['chiMaxST'])]
            if len(chisSorted)==0:
                strtPars = range(0,numProcs)
            elif (len(chisSorted)==1)and(numProcs==1):
                bestRetAry=retAry
                uSTD=1e-6
            else:
                #first updated bestRetAry
                for i in range(len(retAry[0])):
                    if retAry[3][i] in chisSorted:                       
                        if len(bestRetAry[0])<numProcs:
                            for j in range(6):
                                bestRetAry[j].append(retAry[j][i])
                        else:
                            bestChis = np.sort(bestRetAry[3])
                            if retAry[3][i]<bestChis[-1]:
                                for j in range(6):
                                    bestRetAry[j].append(retAry[j][i])
                #now make list of best ones to use in next round
                if len(bestRetAry[0])>numProcs:
                    bestChis = np.sort(bestRetAry[3])
                    log.debug("len best before filtering = "+str(len(bestRetAry[0]))+", worst was "+str(bestChis[-1]))
                    bestRetAry2 = [[],[],[],[],[],[]]
                    #trim best lists down to size
                    for i in range(0,len(bestChis)):
                        if bestRetAry[3][i] in bestChis[:numProcs]:
                            for j in range(6):
                                bestRetAry2[j].append(bestRetAry[j][i])
                    bestRetAry = bestRetAry2
                #copy resulting data files to new temp names.
                for i in range(0,len(bestRetAry[0])):
                    curNm = bestRetAry[0][i]
                    outNm = os.path.join(os.path.dirname(curNm),"SAtempData-"+str(iter)+"-"+str(i)+".fits")
                    tools.renameFits(curNm,outNm,killInput=True)
                    bestRetAry[0][i] = outNm
                #kill those that were not good enough
                for i in range(len(retAry[0])):
                    if os.path.exists(retAry[0][i]):
                        tools.rmFiles([retAry[0][i]])
                ## Now fill out an array of starting parameter sets from the best above.
                ## first load up with one set of goodParams, then randomly from it till full.
                log.debug('best chis:\n' +repr(np.sort(bestRetAry[3])))
                goodParams = bestRetAry[1]
                strtPars=[]
                if len(goodParams)>0:
                    for i in range(0,len(goodParams)):
                        strtPars.append(goodParams[i])
                    while len(strtPars)<numProcs:
                        rndVal = np.random.randint(0,len(goodParams))
                        strtPars.append(goodParams[rndVal])
                log.info(str(len(strtPars))+' sets of starting parameters being passed to next iteration.')
                #print 'STD = '+str(np.std(bestRetAry[3]))
                if len(bestRetAry[3])==numProcs:
                    uSTD = tools.unitlessSTD(bestRetAry[3])
                log.info("After iteration #"+str(iter+1)+" the top "+str(len(bestRetAry[3]))+" solutions with reduced chiSquared < "+str(settings['chiMaxST'])+" have a unitless STD of "+str(uSTD))
                retStr2 +="The latest top "+str(len(bestRetAry[3]))+" reduced chiSquareds had a unitless STD of "+str(uSTD)+'\n'
    ## wrap up
    if len(bestRetAry[0])>1:
        #rename final data files to standard SA convention
        for i in range(0,len(bestRetAry[0])):
            curNm = bestRetAry[0][i]
            outNm = os.path.join(os.path.dirname(curNm),'outputDataSA'+str(i)+'.fits')
            tools.renameFits(curNm,outNm)
            bestRetAry[0][i] = outNm
    #Find best fit, write to file, maybe push into settings files if better than one in there already.
    if len(bestRetAry[0])>0:
        bstChiSqr = np.sort(bestRetAry[3])[0]
        for i in range(len(bestRetAry[0])):
            if bestRetAry[3][i] == bstChiSqr:
                bestpars = bestRetAry[1][i]
                bestsigs = bestRetAry[2][i]
        tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'SA')
    toc=timeit.default_timer()
    s = "ALL "+str(iter+1)+" iterations of SA took a total of "+tools.timeStrMaker(int(toc-tic))
    retStr2 +=s+"\n"
    log.warning(s)
    return (bestRetAry,retStr2)

def exoSOFT():
    """
    'main'
    """
    ## Call startup to get dict and load up final directories into it.
    settings = tools.startup(sys.argv,ExoSOFTdir)
    log = tools.getLogger('main',dir=settings['finalFolder'],lvl=settings['logLevel'])
    tools.logDict(log,settings)
    #log.debug("Prepend string passed in was '"+settings['prepend']+"'")
    Sim = simulator.Simulator(settings)
       
    ###########################################
    # Run nChains for MC/SA/ST mode requested #
    #  Then up to nMCMCcns if MCMC requested  #
    ###########################################     
    tic=timeit.default_timer()
    stageList = settings['stageList']
    durationStrings = ''
    log.warning("Starting requested stages")
    if 'MC' in stageList:
        log.warning("Starting MC stage")
        (returns,b) = (returnsMC,durStr) = multiProc(settings,Sim,'MC',settings['nChains'])
        if len(returnsMC[0])>0:
            bstChiSqr = np.sort(returnsMC[3])[0]
            for i in range(len(returnsMC[0])):
                if returnsMC[3][i] == bstChiSqr:
                    bestpars = returnsMC[1][i]
                    bestsigs = []
            tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'MC')
        durationStrings+='** MC stage **\n'+durStr
        toc=timeit.default_timer()
        log.warning(durStr)
    if 'SA' in stageList:
        log.warning("Starting SA stage")
        (returns,b) = (returnsSA,durStr) = iterativeSA(settings,Sim)
        durationStrings+='** Iterative SA stage **\n'+durStr
    if 'ST' in stageList:
        log.warning("Starting ST stage")
        tic=timeit.default_timer()
        startParams = []
        startSigmas = []
        if settings['stages'] in ['ST','STMCMC']:
            for i in range(0,settings['nChains']):
                startParams.append(settings['startParams'])
                startSigmas.append(settings['startSigmas'])
        elif len(returnsSA)>0:
            for i in range(len(returnsSA[0])):
                if returnsSA[3][i]<settings['chiMaxST']:
                    startParams.append(returnsSA[1][i])
                    startSigmas.append(returnsSA[2][i])
        else:
            log.critical("No SA results available to start the ST chains with.")
        if len(startSigmas)>0:
            (returns,b) = (returnsST,durStr) = multiProc(settings,Sim,'ST',len(startSigmas),startParams,startSigmas)
            durationStrings+='** ST stage **\n'+durStr
        # check best results of ST and store to a file.
        # Maybe replace pars and sigs in original settings files?
        if len(returnsST[0])>0:
            bstChiSqr = np.sort(returnsST[3])[0]
            for i in range(len(returnsST[0])):
                if returnsST[3][i] == bstChiSqr:
                    bestpars = returnsST[1][i]
                    bestsigs = returnsST[2][i]
            tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'ST')
        toc=timeit.default_timer()
        s = "ST took a total of "+tools.timeStrMaker(int(toc-tic))
        log.warning(s)
    if 'MCMC' in stageList:
        log.warning("Starting MCMC stage")
        startParams = []
        startSigmas = []
        if settings['stages']=='MCMC':
            chisSorted = range(0,settings['nMCMCcns'])
            for i in range(0,settings['nMCMCcns']):
                startParams.append(settings['startParams'])
                startSigmas.append(settings['startSigmas'])
        elif len(returnsST)>0:
            chisSorted = []            
            #Filter inputs if more than max num MCMC proc available to use the best ones
            chisSorted = np.sort(returnsST[3])
            chisSorted = chisSorted[np.where(chisSorted<settings['cMaxMCMC'])]
            if len(chisSorted)>settings['nMCMCcns']:
                chisSorted = chisSorted[:settings['nMCMCcns']]
            for i in range(len(returnsST[0])):
                if returnsST[3][i] in chisSorted:
                    startParams.append(returnsST[1][i])
                    startSigmas.append(returnsST[2][i])
            ##$$$$$$$$$$$$$$$$$ MAKE THIS WORK IF nMCMCchains<nSTchains $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        else:
            log.critical("No ST results available to start the MCMC chains with.")
        if len(chisSorted)>0:
            (returns,b) = (returnsMCMC,durStr) = multiProc(settings,Sim,'MCMC',len(chisSorted),startParams,startSigmas)
            durationStrings+='** MCMC stage **\n'+durStr
            # Maybe replace pars in original settings files?
            bstChiSqr = np.sort(returnsMCMC[3])[0]
            for i in range(len(returnsMCMC[0])):
                if returnsMCMC[3][i] == bstChiSqr:
                    bestpars = returnsMCMC[1][i]
                    bestsigs = returnsMCMC[2][i]
            tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'MCMC')
        log.warning(durStr)
    outFiles = returns[0]
    toc=tic2=timeit.default_timer()
    s = "ALL stages took a total of "+tools.timeStrMaker(int(toc-tic))
    durationStrings+=s+'\n'
    log.info(s)
    
    ###################
    # Post-processing # 
    ###################
    log.warning("Starting Post-Processing")  
    
    #figure out
    ## combine the data files
    allFname = ''
    if len(outFiles)>0:
        allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+stageList[-1]+"data.fits")
        tools.combineFits(outFiles,allFname)
    
    ## calc and strip burn-in?
    burnInStr = ''
    if (len(outFiles)>1)and(settings['CalcBurn'] and ('MCMC' in stageList)):
        if 'MCMC' in outFiles[0]:
            (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
            if settings['rmBurn']:
                strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
                outFiles = strippedFnames
                ## combine stripped files to make final file?
                if len(strippedFnames)>0:
                    strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                    tools.combineFits(strippedFnames,strippedAllFname)
                    ## replace final combined filename with new stripped version
                    allFname = strippedAllFname
                
    ## find best fit
    if os.path.exists(allFname):
        bestFit = tools.findBestOrbit(allFname)
            
    ## orbit plots?
    if settings['pltOrbit'] and os.path.exists(allFname):
        plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+stageList[-1])
        tools.orbitPlotter(bestFit,settings,plotFnameBase,format='eps')
    
    ## plot posteriors?
    clStr = ''
    if settings['pltDists'] and os.path.exists(allFname):
        plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+stageList[-1])
        clStr = tools.summaryPlotter(allFname,plotFilename,bestVals=bestFit,stage=settings['stageList'][-1],shadeConfLevels=True,plotALLpars=True)
    
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settings['CalcGR'] and ('MCMC' in stageList)):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settings['nSamples'])
    
    ## progress plots?  INCLUDE?? maybe kill this one. Function exists, but not decided how to use it here.
    
    ## calc correlation length & number effective points? 
    effPtsStr = ''
    if ((len(outFiles)>1)and('MCMC' in stageList))and (settings['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)

    ## Make a summary file of results 
    toc=timeit.default_timer()
    postTime = toc-tic2
    allTime = toc-tic
    if os.path.exists(allFname):
        tools.summaryFile(settings,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings)
    
    ##clean up files (move to folders or delete them)
    tools.cleanUp(settings,stageList,allFname)
    try:
        if settings['CopyToDB']:
            tools.copyToDB(settings)
    except:
        s = " \nAn error occured while trying to copy files to dropbox."
        s+= " \nThis could have happend if the key doesn't exist in your settings dictionary."
        s+= " \nThis could have happend if the key doesn't exist in your settings dictionary.  No biggie as this function will most likely be removed prior to going public."
        log.debug(s)
        
    ## Final log messages and end
    log.debug("Post-processing took a total of "+tools.timeStrMaker(postTime))
    log.warning(" ExoSOFT is Done :-)\n EVERYTHING took a total of "+tools.timeStrMaker(allTime)+'\n')
    log.debug("End of exosoft main")
    ##END MAIN 

if __name__ == '__main__':
    exoSOFT()
#!/usr/bin/env python
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

"""
    This is the 'main' of exosoft. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
""" 
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
    MCmpo= None
    SAmpo = None
    STmpo = None
    MCMCmpo = None
    log.warning("Starting requested stages")
    if 'MC' in stageList:
        log.warning("Starting MC stage")
        MCmpo = tools.multiProcObj(settings,Sim,'MC')
        MCmpo.run()
        MCmpo.writeBest()
        (pars,sigs,chis,outFiles) = MCmpo.getTopProcs(settings['chiMAX'])
        durationStrings+='** MC stage **\n'+MCmpo.latestRetStr
        toc=timeit.default_timer()
        log.warning(MCmpo.latestRetStr)
    if 'SA' in stageList:
        log.warning("Starting SA stage")
        SAmpo = tools.iterativeSA(settings,Sim)
        durationStrings+='** Iterative SA stage **\n'+SAmpo.latestRetStr
    if 'ST' in stageList:
        log.warning("Starting ST stage")
        tic=timeit.default_timer()
        startParams = []
        startSigmas = []
        if settings['stages'] in ['ST','STMCMC']:
            for i in range(0,settings['nChains']):
                startParams.append(settings['startParams'])
                startSigmas.append(settings['startSigmas'])
        elif SAmpo!=None:
            (startParams,startSigmas,chis,outFiles) = SAmpo.getTopProcs(settings['chiMaxST'],fillToNumProc=True)
        else:
            log.critical("No SA results available to start the ST chains with.")
        if len(startParams)>0:
            STmpo = tools.multiProcObj(settings,Sim,'ST')
            STmpo.run(params=startParams,sigmas=startSigmas)
            STmpo.writeBest()
            durationStrings+='** ST stage **\n'+STmpo.latestRetStr
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
        elif STmpo!=None:
            try:
                (startParams,startSigmas,chisSorted,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=True,nProcs=settings['nMCMCcns'],allBest=settings['strtMCMCatBest'])
                log.debug('Starting all MCMC chain at the same top fit '+\
                          'found during the ST stage.')
            except:
                (startParams,startSigmas,chisSorted,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=True,nProcs=settings['nMCMCcns'])
        else:
            log.critical("No ST results available to start the MCMC chains with.")
        if len(startParams)>0:
            MCMCmpo = tools.multiProcObj(settings,Sim,'MCMC')
            MCMCmpo.run(params=startParams,sigmas=startSigmas)
            MCMCmpo.writeBest()
            durationStrings+='** MCMC stage **\n'+MCMCmpo.latestRetStr
        log.warning(MCMCmpo.latestRetStr)
    
    ## Done all stages 
    toc=tic2=timeit.default_timer()
    s = "ALL stages took a total of "+tools.timeStrMaker(int(toc-tic))
    durationStrings+=s+'\n'
    log.info(s)
    
    ###################
    # Post-processing # 
    ###################
    log.warning("Starting Post-Processing")  
    
    FINALmpo = None
    if MCMCmpo!=None:
        FINALmpo = MCMCmpo
        tools.pklIt(settings,MCMCmpo.resultsOnly(),'MCMCmpoRO')
        tools.pklIt(settings, FINALmpo.resultsOnly(),'FINALmpoRO')
    elif STmpo!=None:
        FINALmpo = STmpo
        tools.pklIt(settings,STmpo.resultsOnly(),'STmpoRO')
        tools.pklIt(settings, FINALmpo.resultsOnly(),'FINALmpoRO')
    elif SAmpo!=None:
        FINALmpo = SAmpo
        tools.pklIt(settings,SAmpo.resultsOnly(),'SAmpoRO')
        tools.pklIt(settings, FINALmpo.resultsOnly(),'FINALmpoRO')
    elif MCmpo!=None:
        FINALmpo = MCmpo
        tools.pklIt(settings, MCmpo.resultsOnly(),'MCmpoRO')
        tools.pklIt(settings, FINALmpo.resultsOnly(),'FINALmpoRO')
    else:
        log.critical("\nNo FINALmpo exists!!! \nExoSOFT failed to complete any of the requested stages!!")
    
    if FINALmpo!=None:
        outFiles = FINALmpo.outFnames
        allFname = ''
        bestFit = []
        burnInStr = ''
        clStr = ''
        grStr = ''
        effPtsStr = ''
        postTime = ''
        allTime = ''
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        ## combine the data files
        if len(outFiles)>0:
            allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+FINALmpo.stage+"data.fits")
            tools.combineFits(outFiles,allFname)
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        ## calc and strip burn-in?
        if (len(outFiles)>1)and(settings['CalcBurn'] and (FINALmpo.stage=='MCMC')):
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
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        ## find best fit
        if os.path.exists(allFname):
            bestFit = tools.findBestOrbit(allFname)
        ## orbit plots?
        if settings['pltOrbit'] and os.path.exists(allFname):
            plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+FINALmpo.stage)
            tools.orbitPlotter(bestFit,settings,plotFnameBase,format='eps')
        ## plot posteriors?
        if settings['pltDists'] and os.path.exists(allFname):
            plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+FINALmpo.stage)
            clStr = tools.summaryPlotter(allFname,plotFilename,bestVals=bestFit,stage=FINALmpo.stage,shadeConfLevels=True,plotALLpars=True)
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        ##calc Gelman-Rubin convergence statistics?
        if (len(outFiles)>1) and (settings['CalcGR'] and (FINALmpo.stage=='MCMC')):
            grStr = tools.gelmanRubinCalc(outFiles,settings['nSamples'])
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        
        ## progress plots?  INCLUDE?? maybe kill this one. Function exists, but not decided how to use it here.
        
        ## calc correlation length & number effective points? 
        if ((len(outFiles)>1)and(FINALmpo.stage=='MCMC'))and (settings['calcCL'] and os.path.exists(allFname)):
            effPtsStr = tools.mcmcEffPtsCalc(allFname)
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        
        ## Make a summary file of results 
        toc=timeit.default_timer()
        postTime = toc-tic2
        allTime = toc-tic
        if os.path.exists(allFname):
            tools.summaryFile(settings,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings,MCmpo,SAmpo,STmpo,MCMCmpo)
        ## pickle final versions of all the results
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        tools.pklIt(settings,settings,'settings')
        
        ##clean up files (move to folders or delete them)
        tools.cleanUp(settings,stageList,allFname)
        try:
            if settings['CopyToDB']:
                tools.copyToDB(settings)
        except:
            s = " \nAn error occured while trying to copy files to dropbox."
            s+= " \nThis could have happend if the key doesn't exist in your settings dictionary."
            s+= " \nNo biggie as this function will most likely be removed prior to going public."
            log.debug(s)
            
        ## Final log messages and end
        log.debug("Post-processing took a total of "+tools.timeStrMaker(postTime))
        log.warning(" ExoSOFT is Done :-)\n EVERYTHING took a total of "+tools.timeStrMaker(allTime)+'\n')
        log.debug("End of exosoft main")
        ##END MAIN 

if __name__ == '__main__':
    exoSOFT()
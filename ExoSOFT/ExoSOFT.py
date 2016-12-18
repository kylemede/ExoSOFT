#@Author: Kyle Mede, kylemede@gmail.com
from __future__ import absolute_import
#import matplotlib
# Force matplotlib to not use any Xwindows backend, to further avoid Display issues or when ExoSOFT is ran through ssh without -X.
#matplotlib.use('Agg') 
#import sys
import os
import timeit
import KMlogger
#import numpy as np
from six.moves import range

from . import tools
from . import simulator

"""
    This is the 'main' of ExoSOFT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
""" 
def exoSOFT(settings_in, advanced_settings_in, priors_in):
    ## Call startup to get dict and load up final directories into it.
    #settings = tools.startup(sett_file_path)
    settings = tools.startup(settings_in, advanced_settings_in, priors_in)
    log = KMlogger.getLogger('main',dr=settings['finalFolder'],lvl=settings['logLevel'])
    log.logDict(settings)
    #log.debug("Prepend string passed in was '"+settings['prepend']+"'")
    Sim = simulator.Simulator(settings)
       
    ###########################################
    # Run nChains for MC/SA/ST mode requested #
    #  Then up to nMCMCcns if MCMC requested  #
    ###########################################     
    tic=timeit.default_timer()
    stageList = settings['stageList']
    durationStrings = ''
    stgsPassed = True
    [MCMCmpo,SAmpo,STmpo,MCmpo,emcee_mpo] = [None,None,None,None,None]
    log.warning("Starting requested stages")
    if 'MC' in stageList:
        log.warning("Starting MC stage")
        MCmpo = tools.multiProcObj(settings,Sim,'MC')
        MCmpo.run()
        stgsPassed = MCmpo.writeBest()
        if stgsPassed:
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
            (startParams,startSigmas,_,outFiles) = SAmpo.getTopProcs(settings['chiMaxST'],fillToNumProc=True)
        else:
            log.critical("No SA results available to start the ST chains with.")
        if len(startParams)>0:
            STmpo = tools.multiProcObj(settings,Sim,'ST')
            STmpo.run(params=startParams,sigmas=startSigmas)
            stgsPassed = STmpo.writeBest()
            if stgsPassed:
                durationStrings+='** ST stage **\n'+STmpo.latestRetStr
        toc=timeit.default_timer()
        s = "ST took a total of "+tools.timeStrMaker(int(toc-tic))
        log.warning(s)
    if 'MCMC' in stageList:
        log.warning("Starting MCMC stage")
        startParams = []
        startSigmas = []
        if settings['stages']=='MCMC':
            _ = list(range(0,settings['nMCMCcns']))
            for _ in range(0,settings['nMCMCcns']):
                startParams.append(settings['startParams'])
                startSigmas.append(settings['startSigmas'])
        elif STmpo!=None:
            try:
                (startParams,startSigmas,_,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=True,nProcs=settings['nMCMCcns'],allBest=settings['strtMCMCatBest'])
                log.debug('Starting all MCMC chain at the same top fit '+\
                          'found during the ST stage.')
            except:
                (startParams,startSigmas,_,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=True,nProcs=settings['nMCMCcns'])
        else:
            log.critical("No ST results available to start the MCMC chains with.")
        if len(startParams)>0:
            MCMCmpo = tools.multiProcObj(settings,Sim,'MCMC')
            MCMCmpo.run(params=startParams,sigmas=startSigmas)
            stgsPassed = MCMCmpo.writeBest()
            if stgsPassed:
                durationStrings+='** MCMC stage **\n'+MCMCmpo.latestRetStr
            log.warning(MCMCmpo.latestRetStr)
    if 'emcee' in stageList:
        log.warning("Starting emcee stage")
        startParams = []
        startSigmas = []
        if STmpo!=None:
            (startParams,startSigmas,_,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=False,nProcs=1)
            log.debug('Use best fit from ST stage as the center of the ball for the emcee walkers')
        elif SAmpo!=None:
            (startParams,startSigmas,_,outFiles) = SAmpo.getTopProcs(settings['cMaxMCMC'],fillToNumProc=False,nProcs=1)
            log.debug('Use best fit from SA stage as the center of the ball for the emcee walkers')
        else:
            startParams = [settings['startParams']]
            startSigmas = [settings['startSigmas']]
            log.debug('Use startParams in settings dictionary as the center of the ball for emcee walkers')
        emcee_mpo = tools.multiProcObj(settings,Sim,'emcee')
        emcee_mpo.run(params=startParams,sigmas=startSigmas)
        stgsPassed = emcee_mpo.writeBest()
        if stgsPassed:
            durationStrings+='** emcee stage **\n'+emcee_mpo.latestRetStr
        log.warning(emcee_mpo.latestRetStr)
    
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
    if emcee_mpo!=None:
        FINALmpo = emcee_mpo
        tools.pklIt(settings,emcee_mpo.resultsOnly(),'emcee_mpo')
        tools.pklIt(settings, FINALmpo.resultsOnly(),'FINALmpoRO')
    elif MCMCmpo!=None:
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
        [allFname,burnInStr,clStr,grStr,effPtsStr,iacStr,postTime,allTime] = ['','','','','','','','']
        bestFit = []
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
            bestFit = tools.findBestOrbit(allFname,by_ln_prob= FINALmpo.stage=="emcee")
            #bestFit = tools.findBestOrbit(allFname)
            
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
        
        ## calc mean correlation length & number effective points? 
        if settings['calcCL'] and os.path.exists(allFname):
            if ((len(outFiles)>1)and(FINALmpo.stage=='MCMC'))or(FINALmpo.stage=='emcee'):
                effPtsStr = tools.mcmcEffPtsCalc(allFname)
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings],'finalSummaryStrs')
        
        ## calc auto-correlation time using tools from the emcee package
        if ((FINALmpo.stage=='emcee')or(FINALmpo.stage=='MCMC')) and (settings['calcIAC'] and os.path.exists(allFname)):
            #print("\n about to call autocorr")
            iacStr = tools.autocorr(allFname)
        #print("iacStr = "+iacStr)
        
        ## Make a summary file of results 
        toc=timeit.default_timer()
        postTime = toc-tic2
        allTime = toc-tic
        if os.path.exists(allFname):
            tools.summaryFile(settings,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,allTime,postTime,durationStrings,MCmpo,SAmpo,STmpo,MCMCmpo,emcee_mpo)
        ## pickle final versions of all the results
        tools.pklIt(settings,[allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,allTime,postTime,durationStrings],'finalSummaryStrs')
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
        log.debug("End of ExoSOFT main")
# EOF
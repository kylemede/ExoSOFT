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
#         (returns,b) = (returnsMC,durStr) = tools.multiProc(settings,Sim,'MC',settings['nChains'])
#         if len(returnsMC[0])>0:
#             bstChiSqr = np.sort(returnsMC[3])[0]
#             for i in range(len(returnsMC[0])):
#                 if returnsMC[3][i] == bstChiSqr:
#                     bestpars = returnsMC[1][i]
#                     bestsigs = []
#             tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'MC')
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
            (startParams,startSigmas,chis,outFiles) = SAmpo.getTopProcs(settings['chiMaxST'],killBadOnes=False,fillToNumProc=True)
        else:
            log.critical("No SA results available to start the ST chains with.")
        if len(chis)>0:
            STmpo = tools.multiProcObj(settings,Sim,'ST')
            STmpo.run(params=startParams,sigmas=startSigmas)
            STmpo.writeBest()
            #(returns,b) = (returnsST,durStr) = tools.multiProc(settings,Sim,'ST',len(startSigmas),startParams,startSigmas)
            durationStrings+='** ST stage **\n'+STmpo.latestRetStr
#         # check best results of ST and store to a file.
#         # Maybe replace pars and sigs in original settings files?
#         if len(returnsST[0])>0:
#             bstChiSqr = np.sort(returnsST[3])[0]
#             for i in range(len(returnsST[0])):
#                 if returnsST[3][i] == bstChiSqr:
#                     bestpars = returnsST[1][i]
#                     bestsigs = returnsST[2][i]
#             tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'ST')
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
            (startParams,startSigmas,chisSorted,outFiles) = STmpo.getTopProcs(settings['cMaxMCMC'],killBadOnes=False,fillToNumProc=True,nProcs=settings['nMCMCcns'])
#             chisSorted = []            
#             #Filter inputs if more than max num MCMC proc available to use the best ones
#             chisSorted = np.sort(returnsST[3])
#             chisSorted = chisSorted[np.where(chisSorted<settings['cMaxMCMC'])]
#             if len(chisSorted)>settings['nMCMCcns']:
#                 chisSorted = chisSorted[:settings['nMCMCcns']]
#             for i in range(len(returnsST[0])):
#                 if returnsST[3][i] in chisSorted:
#                     startParams.append(returnsST[1][i])
#                     startSigmas.append(returnsST[2][i])
#             ##$$$$$$$$$$$$$$$$$ MAKE THIS WORK IF nMCMCchains<nSTchains $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        else:
            log.critical("No ST results available to start the MCMC chains with.")
        if len(chis)>0:
            MCMCmpo = tools.multiProcObj(settings,Sim,'MCMC')
            MCMCmpo.run(params=startParams,sigmas=startSigmas)
            MCMCmpo.writeBest()
            #(returns,b) = (returnsMCMC,durStr) = tools.multiProc(settings,Sim,'MCMC',len(chisSorted),startParams,startSigmas)
            durationStrings+='** MCMC stage **\n'+MCMCmpo.latestRetStr
            # Maybe replace pars in original settings files?
#             bstChiSqr = np.sort(returnsMCMC[3])[0]
#             for i in range(len(returnsMCMC[0])):
#                 if returnsMCMC[3][i] == bstChiSqr:
#                     bestpars = returnsMCMC[1][i]
#                     bestsigs = returnsMCMC[2][i]
#             tools.writeBestsFile(settings,bestpars,bestsigs,bstChiSqr,'MCMC')
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
    elif STmpo!=None:
        FINALmpo = STmpo
    elif SAmpo!=None:
        FINALmpo = SAmpo
    elif MCmpo!=None:
        FINALmpo = MCmpo
    else:
        log.critical("\nNo FINALmpo exists!!! \nExoSOFT failed to complete any of the requested stages!!")
    
    if FINALmpo!=None:
        outFiles = FINALmpo.outFnames
        ## combine the data files
        allFname = ''
        if len(outFiles)>0:
            allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+FINALmpo.stage+"data.fits")
            tools.combineFits(outFiles,allFname)
        
        ## calc and strip burn-in?
        burnInStr = ''
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
                    
        ## find best fit
        if os.path.exists(allFname):
            bestFit = tools.findBestOrbit(allFname)
                
        ## orbit plots?
        if settings['pltOrbit'] and os.path.exists(allFname):
            plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+FINALmpo.stage)
            tools.orbitPlotter(bestFit,settings,plotFnameBase,format='eps')
        
        ## plot posteriors?
        clStr = ''
        if settings['pltDists'] and os.path.exists(allFname):
            plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+FINALmpo.stage)
            clStr = tools.summaryPlotter(allFname,plotFilename,bestVals=bestFit,stage=FINALmpo.stage,shadeConfLevels=True,plotALLpars=True)
        
        ##calc Gelman-Rubin convergence statistics?
        grStr = ''
        if (len(outFiles)>1) and (settings['CalcGR'] and (FINALmpo.stage=='MCMC')):
            grStr = tools.gelmanRubinCalc(outFiles,settings['nSamples'])
        
        ## progress plots?  INCLUDE?? maybe kill this one. Function exists, but not decided how to use it here.
        
        ## calc correlation length & number effective points? 
        effPtsStr = ''
        if ((len(outFiles)>1)and(FINALmpo.stage=='MCMC'))and (settings['calcCL'] and os.path.exists(allFname)):
            effPtsStr = tools.mcmcEffPtsCalc(allFname)
    
        ## Make a summary file of results 
        toc=timeit.default_timer()
        postTime = toc-tic2
        allTime = toc-tic
        if os.path.exists(allFname):
            #MCmpo= None
            #SAmpo = None
            #STmpo = None
            #MCMCmpo = None
            tools.summaryFile(settings,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings,MCmpo,SAmpo,STmpo,MCMCmpo)
        
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
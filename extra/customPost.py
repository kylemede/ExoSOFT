#@Author: Kyle Mede, kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
# Force matplotlib to not use any Xwindows backend, to further avoid Display issues or when ExoSOFT is ran through ssh without -X.
matplotlib.use('Agg')
from ExoSOFT import tools
import sys
import os
import numpy as np
import timeit
import yaml
import glob
import KMlogger
from six.moves import range

"""
 The function to re-run the post-analysis routines of ExoSOFT, either with
 their default values, or your own custom ones.

 PLACE A COPY OF THIS FILE INTO THE SAME DIRECTORY AS YOUR SETTINGS AND DATA FILES!!

 Make any changes to the input plotting and other post-processing
 functions necessary to get them perfectly ready for your publishing.

 Then, start with :
 $python customPost.py
"""

def custom_post(settings_in,advanced_settings_in, priors_in):
    ## load up settings that were passed in
    tic=timeit.default_timer()
    settings = tools.startup(settings_in,advanced_settings_in,priors_in,rePlot=True)
    log = KMlogger.getLogger('main',dr=settings['finalFolder'],lvl=settings['logLevel'])
    skipBurnInStrip=True

    # try to bring back pickeled fleshed out settings file and final results strings
    # if they exists.
    [allFname,burnInStr,clStr,grStr,effPtsStr,durationStrings,postTime,allTime] = ['','','','','','','','']
    [MCMCmpoRO,SAmpoRO,STmpoRO,MCmpoRO,emcee_mpoRO,FINALmpoRO] = [None,None,None,None,None,None]
    outFiles = []
    if False:
        try:
            settings = tools.unPklIt(settings,'settings')
            finalSummaryStrs = tools.unPklIt(settings,'finalSummaryStrs')
            [allFname,outFiles,stageList,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings] = finalSummaryStrs
            FINALmpoRO = tools.unPklIt(settings,'FINALmpoRO')
            # mpoRO contains members: bestRedChiSqrs,avgAcceptRates,acceptStrs,stage,retStr,latestRetStr
            [MCmpoRO,SAmpoRO,STmpoRO,MCMCmpoRO] = tools.reloadMpoROs(settings)
            'For MCMC stage: '+str(np.mean(MCMCmpoRO.avgAcceptRates))
        except:
             print('AN ERROR OCCURRED WHILE TRYING TO LOAD PKLS BACK IN')

    ## make list of output files
    outFiles = []
    if settings['stageList'][-1]=='MCMC':
        outFiles = np.sort(glob.glob(os.path.join(settings['finalFolder'],\
                                              "outputDataMCMC*_BIstripped.fits")))
        if len(outFiles)==0:
            outFiles = np.sort(glob.glob(os.path.join(settings['finalFolder'],\
                                                  "outputDataMCMC*.fits")))
    else:
        outFiles = np.sort(glob.glob(os.path.join(settings['finalFolder'],\
                                              "outputData"+\
                                              settings['stageList'][-1]+"*.fits")))
    if len(outFiles)==0:
        log.error("No appropriate list of individual chain outputs to "+\
                  "perform burn-in stripping on, or combine.")
    ## get name for combined output data file
    allFname = ''
    if os.path.exists(os.path.join(settings['finalFolder'],\
                                   "combined-BIstripped-MCMCdata.fits")):
        allFname = os.path.join(settings['finalFolder'],\
                                "combined-BIstripped-MCMCdata.fits")
        skipBurnInStrip=True
    elif os.path.exists(os.path.join(settings['finalFolder'],\
                                     "combined"+settings['stageList'][-1]+"data.fits")):
        allFname = os.path.join(settings['finalFolder'],\
                                "combined"+settings['stageList'][-1]+"data.fits")
        skipBurnInStrip=False
    else:
        skipBurnInStrip=True
        log.critical("No combined MCMC file available.  Please check and/or "+\
                     "modify the customPost.py script as needed.")
    ## combine the individual chain data files needed
    if (len(outFiles)>0) and (allFname==''):
        log.importantinfo('about to combine data files together')
        if 'BIstripped' in outFiles[0]:
            allFname = os.path.join(os.path.dirname(outFiles[0]),\
                                    "combined-BIstripped-MCMCdata.fits")
        else:
            allFname = os.path.join(os.path.dirname(outFiles[0]),\
                                    "combined"+settings['stageList'][-1]+"data.fits")
        tools.combineFits(outFiles,allFname)

    ## calc and strip burn-in?
    burnInStr = ''
    if True:
        if skipBurnInStrip==False:
            log.importantinfo('about to strip burn-in')
            if (len(outFiles)>1)and(settings['CalcBurn'] and\
                                    (settings['stageList'][-1]=='MCMC')):
                (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)
                if settings['rmBurn'][0]:
                    strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
                    outFiles = strippedFnames
                    ## combine stripped files to make final file?
                    if len(strippedFnames)>0:
                        strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),\
                                                        "combined-BIstripped-MCMCdata.fits")
                        tools.combineFits(strippedFnames,strippedAllFname)
                        ## replace final combined filename with new stripped version
                        allFname = strippedAllFname

    ## find best fit
    if True:
        log.importantinfo('about to find best orbit')
        bestFit = tools.findBestOrbit(allFname,bestToFile=False,findAgain=False)
    else:
        # for manual input of best values
        bestFit = np.array([0.9892581,0.0009696,49.4824152,101.82462,0.05706358,2450656.61,2450656.61,12.012219,44.5819933,0.1945267,5.227630,46.908327,8.916892,-0.04136662])
    
    ## If using the 'expected' values for synthetic data, they are in 'stored_pars' format.
    # stored_pars: [m1,m2,parallax,long_an,e,to,tc,p,inc,arg_peri,a_tot_au,chi_sqr,K,v1,v2...]
    expectedValues = np.array([1.0,0.0009543,50.0,101.6,0.048,2450639.0,2450639.0,11.9,45.0,14.8,5.21,  1.0,   8.9,0.0])

    if False:
        log.importantinfo('about to make orbit plots')
        plotFnameBase = os.path.join(settings['finalFolder'],'orbPlot-Manual')
        ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
        ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
        tools.orbitPlotter(bestFit,settings,plotFnameBase,fl_format='eps',DIlims=[],RVlims=[])
        ## Simulated jupiter
        #tools.orbitPlotter(bestFit,settings,plotFnameBase,fl_format='eps',DIlims=[],RVlims=[[-9.95,9.9],[-1.5,1.5],[-0.55,0.55]],diErrMult=1,diLnThk=2.5)

    clStr=''
    if True:
        log.importantinfo('about to plot posteriors')
        tic2=timeit.default_timer()
        plotFilename = os.path.join(settings['finalFolder'],'posteriors')
        clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],bestVals=bestFit,stage=settings['stageList'][-1], shadeConfLevels=True,forceRecalc=False,plotALLpars=True)
        #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,8,10],xLims=[[0.5,1.6],[0.5,1.5],[35,52],[4.0,6.3]],bestVals=expectedValues,stage=settings['stageList'][-1], shadeConfLevels=True,forceRecalc=False)
        log.importantinfo("It took a total of "+tools.timeStrMaker(timeit.default_timer()-tic2))

    ## make corner plot?
    if False:
        log.importantinfo('about to make corner plot')
        tic2=timeit.default_timer()
        plotFilename = os.path.join(settings['finalFolder'],'cornerPlot_expectedVals')
        label_kwargs={'fontsize':28,'weight':'heavy','labelpad':1}
        other_kwargs={'labelsize':18,'xy_factor':1.5}
        tools.cornerPlotter(allFname,plotFilename,paramsToPlot=[0,1,8,10],bestVals=expectedValues,smooth=True,label_kwargs=label_kwargs,latex=True,fl_format='png',other_kwargs=other_kwargs)
        log.importantinfo("It took a total of "+tools.timeStrMaker(timeit.default_timer()-tic2))

    if False:
        log.importantinfo('about to make 2D density plot')
        plotFilename = os.path.join(settings['finalFolder'],'m1m2-densPlot')
        #ranges=[[xMin,xMax],[yMin,yMax]]
        tools.densityPlotter2D(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]])

    if False:
        log.importantinfo('about to make progress plots')
        ##make progress plot
        for stgStr in ['SA','ST','MCMC']:
            emcee_stage = False
            if stgStr=='emcee':
                emcee_stage = True
            for procNum in range(0,2):
                dataFname = os.path.join(settings['finalFolder'],'outputData'+stgStr+str(procNum)+'.fits')
                if os.path.exists(dataFname):
                    log.info("About to make progress plots for file: "+dataFname)
                    for parNum in [1,4,8]:
                        plotFilename = os.path.join(settings['finalFolder'],'progressPlot-'+stgStr+str(procNum)+'-'+str(parNum))
                        try:
                            #outputDataFilename,plotFilename,paramToPlot,yLims=[],xLims = [],expectedVal=None,emcee_stage=False,downSample=True
                            tools.progressPlotter(dataFname,plotFilename,parNum,yLims=[],xLims = [],expectedVal=None,emcee_stage=emcee_stage,downSample=False)
                        except:
                            log.importantinfo('could not make plot for proc# '+str(procNum)+', and par# '+str(parNum))
    ##calc R?
    grStr = ''
    if False:
        log.importantinfo('about to calc GR')
        if (len(outFiles)>1) and (settings['CalcGR'] and (settings['stageList'][-1]=='MCMC')):
            tic2=timeit.default_timer()
            grStr = tools.gelmanRubinCalc(outFiles,settings['nSamples'])
            log.info(grStr)
            log.importantinfo("It took a total of "+tools.timeStrMaker(timeit.default_timer()-tic2))

    ## custom re check of the orbit fit
    if False:
        log.importantinfo('about to calculate predicted location for custom date')
        orbParams = bestFit
        epochs=[2457327.500]
        tools.predictLocation(orbParams,settings,epochs)

    ## calc correlation length & number effective points?
    effPtsStr = ''
    if False:
        log.importantinfo('about to calc # effective points')
        if ((len(outFiles)>1)and(settings['stageList'][-1]=='MCMC'))and (settings['calcCL'] and os.path.exists(allFname)):
            tic2=timeit.default_timer()
            effPtsStr = tools.mcmcEffPtsCalc(allFname)
            log.info(effPtsStr)
            log.importantinfo("It took a total of "+tools.timeStrMaker(timeit.default_timer()-tic2))

    ## calc auto-correlation time using tools from the emcee package
    if False:
        log.importantinfo('about to calc integrated autocorr time')
        if ((settings['stageList'][-1]=='emcee')or(settings['stageList'][-1]=='MCMC')) and (settings['calcIAC'] and os.path.exists(allFname)):
            tic2=timeit.default_timer()
            iacStr = tools.autocorr(allFname)
            log.info(iacStr)
            log.importantinfo("It took a total of "+tools.timeStrMaker(timeit.default_timer()-tic2))

    postTime = timeit.default_timer()-tic
    ## following post-processing stages can take a long time, so write the current
    ## summary information to the summary file and add the rest later
    if False:
        if os.path.exists(allFname):
            log.importantinfo('about to make summary file')
            tools.summaryFile(settings,settings['stageList'],allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,allTime,postTime,durationStrings,MCmpoRO,SAmpoRO,STmpoRO,MCMCmpoRO,emcee_mpoRO)

    ## CAUTION!!  This way of cleaning is permenant!!
    ## Double check your settings and that this is what you want first!!
    ## clean up files (move to folders or delete them)
    if False:
        log.importantinfo('about to clean up output directory')
        tools.cleanUp(settings,settings['stageList'],allFname)

    ## copy to dropbox is a legacy function that will be repurposed or killed later
    if False and ('copyToDB' in settings):
        if settings['copyToDB']:
            log.importantinfo('about to copy files to dropbox')
            tools.copyToDB(settings)


    log.importantinfo("It took a total of "+tools.timeStrMaker(postTime))
    ### DONE customPost
    print("\ncustomPost script complete :-D\n")


if __name__ == '__main__':
    settings_in = None
    priors_in = None
    advanced_settings_in = None
    if os.path.exists('./settings.yaml'):
        f = open('./settings.yaml','r')
        settings_in = yaml.load(f)
        f.close()
    if os.path.exists('./advanced_settings.yaml'):
        f = open('./advanced_settings.yaml','r')
        advanced_settings_in = yaml.load(f)
        f.close()
    if os.path.exists('./priors.py'):
        from priors import ExoSOFTpriors as priors_in

    custom_post(settings_in,advanced_settings_in, priors_in)
#END OF FILE

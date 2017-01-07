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
    settings = tools.startup(settings_in, advanced_settings_in,priors_in,rePlot=True)
    log = KMlogger.getLogger('main',dr=settings['finalFolder'],lvl=100)
    skipBurnInStrip=True
    
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
        bestFit = np.array([0.9892581,0.0009696,49.4824152,101.82462,0.05706358,2450656.61,2450656.61,12.012219,44.5819933,0.1945267,5.227630,46.908327,8.916892,-0.04136662])
    if True:
        log.importantinfo('about to make orbit plots')
        ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
        ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
        plotFnameBase = os.path.join(settings['finalFolder'],'orbPlot-Manual')
        tools.orbitPlotter(bestFit,settings,plotFnameBase,fl_format='eps',DIlims=[],RVlims=[])
        
    clStr=''
    if True:
        log.importantinfo('about to plot posteriors')
        plotFilename = os.path.join(settings['finalFolder'],'posteriors')
        clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],bestVals=bestFit,stage=settings['stageList'][-1], shadeConfLevels=True,forceRecalc=True,plotALLpars=True)
        
    if True: 
        log.importantinfo('about to make 2D density plot')
        plotFilename = os.path.join(settings['finalFolder'],'m1m2-densPlot')
        #ranges=[[xMin,xMax],[yMin,yMax]]
        tools.densityPlotter2D(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]])
        
    if True:
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
    if True:
        log.importantinfo('about to calc GR')
        if (len(outFiles)>1) and (settings['CalcGR'] and (settings['stageList'][-1]=='MCMC')):
            grStr = tools.gelmanRubinCalc(outFiles,settings['nSamples'])
            log.info(grStr)
        
    ## custom re check of the orbit fit     
    if False: 
        log.importantinfo('about to calculate predicted location for custom date')
        orbParams = bestFit
        epochs=[2457327.500]
        tools.predictLocation(orbParams,settings,epochs)
    
    ## calc correlation length & number effective points? 
    effPtsStr = ''
    if True:
        log.importantinfo('about to calc # effective points')
        if ((len(outFiles)>1)and(settings['stageList'][-1]=='MCMC'))and (settings['calcCL'] and os.path.exists(allFname)):
            effPtsStr = tools.mcmcEffPtsCalc(allFname)
            log.info(effPtsStr)
            
    ## calc auto-correlation time using tools from the emcee package
    if True:
        log.importantinfo('about to calc integrated autocorr time')
        if ((settings['stageList'][-1]=='emcee')or(settings['stageList'][-1]=='MCMC')) and (settings['calcIAC'] and os.path.exists(allFname)):
            #print("\n about to call autocorr")
            iacStr = tools.autocorr(allFname)
            log.info(iacStr)
    ## following post-processing stages can take a long time, so write the current
    ## summary information to the summary file and add the rest later
    if True:
        if os.path.exists(allFname):
            log.importantinfo('about to make summary file')
            tools.summaryFile(settings,settings['stageList'],allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,1,1,'',None,None,None,None,None)
            ##FULL VERSION IN EXOSOFT, but not all available to customPost
            #tools.summaryFile(settings,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,allTime,postTime,durationStrings,MCmpo,SAmpo,STmpo,MCMCmpo,emcee_mpo)
        
    ## CAUTION!!  This way of cleaning is permenant!!
    ## Double check your settings and that this is what you want first!!
    ##clean up files (move to folders or delete them)
    if False:
        log.importantinfo('about to clean up output directory')
        tools.cleanUp(settings,settings['stageList'],allFname)
        
    ## copy to dropbox is a legacy function that will be repurposed or killed later 
    if False and settings['CopyToDB']:
        log.importantinfo('about to copy files to dropbox')
        tools.copyToDB(settings)
        
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
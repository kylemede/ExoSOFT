import tools
#from tools import constants as const
from exosoftpath import rootDir as ExoSOFTdir
import sys
import os
import numpy as np
import glob

def customPost():
    settingsDict = tools.startup(sys.argv,ExoSOFTdir,rePlot=True)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=100)
    skipBurnInStrip=True
    allFname = ''
    if os.path.exists(os.path.join(settingsDict['finalFolder'],\
                                   "combined-BIstripped-MCMCdata.fits")):
        allFname = os.path.join(settingsDict['finalFolder'],\
                                "combined-BIstripped-MCMCdata.fits")
    elif os.path.exists(os.path.join(settingsDict['finalFolder'],\
                                     "combinedMCMCdata.fits")):
        allFname = os.path.join(settingsDict['finalFolder'],\
                                "combinedMCMCdata.fits")
        skipBurnInStrip=False
    else:
        skipBurnInStrip=True
        log.critical("No combined MCMC file available.  Please check and/or "+\
                     "modify the customPost.py script as needed.")
    ## run make for swig if requested??
    if settingsDict['remake'] and False:
        cwd = os.getcwd()
        log.debug("-"*45+" Starting to remake CPP/SWIG tools "+45*"-")
        os.chdir(os.path.join(settingsDict['ExoSOFTdir'],'tools/cppTools/'))
        os.system('make clean')
        os.system('make')
        os.chdir(cwd)
        log.debug("-"*45+" Done re-making CPP/SWIG tools "+45*"-")

    ##make hack list of output files
    outFiles = []
    if 'BIstripped' in allFname:
        outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],\
                                              "outputDataMCMC*_BIstripped.fits")))
    elif os.path.exists(allFname):
        outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],\
                                              "outputDataMCMC*.fits")))
    else:
        log.error("No appropriate list of output non-combined outputs to "+\
                  "perform burn-in stripping on.")
    
    ## calc and strip burn-in?
    burnInStr = ''
    if True:
        if skipBurnInStrip==False:
            print 'about to strip burn-in'
            if (len(outFiles)>1)and(settingsDict['CalcBurn'] and\
                                    (settingsDict['stageList'][-1]=='MCMC')):
                (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
                if settingsDict['rmBurn'][0]:
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
    if False:
        print 'about to find best orbit'
        bestFit = tools.findBestOrbit(allFname,bestToFile=False,findAgain=False)
    else:
        bestFit = np.array([0.605838520481,0.146024522998,48.0692312367,98.8893962294,0.640666046184,2445709.98532,2445709.98532,8.2093394771,40.2779824043,166.743572004,3.39834214958,43.3616906776,0.0,9730.87607461,11612.1046117,12586.0296274])
    if True:
        print 'about to make orbit plots'
        ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
        ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
        plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-Manual')
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[])
        
    clStr=''
    if False:
        print 'about to plot posteriors'
        plotFilename = os.path.join(settingsDict['finalFolder'],'posteriors')
        clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],bestVals=bestFit,stage=settingsDict['stageList'][-1], shadeConfLevels=True,forceRecalc=True,plotALLpars=True)
        
    if False: 
        print 'about to make 2D density plot'
        plotFilename = os.path.join(settingsDict['finalFolder'],'m1m2-densPlot-Dec23PASAfit-151226-autoRanges')
        #ranges=[[xMin,xMax],[yMin,yMax]]
        tools.densityPlotter2D(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]])
        
    if False:
        print ' about to make progress plots'
        #make progress plot
        for stgStr in ['SA','ST','MCMC']:
            for procNum in range(0,4):
                dataFname = os.path.join(settingsDict['finalFolder'],'outputData'+stgStr+str(procNum)+'.fits')
                for parNum in [1,4,8]:
                    plotFilename = os.path.join(settingsDict['finalFolder'],'progressPlot-'+stgStr+str(procNum)+'-'+str(parNum))
                    try:
                        tools.progressPlotter(dataFname,plotFilename,parNum,yLims=[],bestVals=[])
                    except:
                        print 'could not make plot for proc# '+str(procNum)+', and par# '+str(parNum)
    ##calc R?
    grStr = ''
    if False:
        print 'about to calc GR'
        if (len(outFiles)>1) and (settingsDict['CalcGR'] and (settingsDict['stageList'][-1]=='MCMC')):
            (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
        
    ## custom re check of the orbit fit     
    if False: 
        print 'about to calculate predicted location for custom date'
        orbParams = bestFit
        epochs=[ 2457327.500]
        tools.predictLocation(orbParams,settingsDict,epochs)
    
    ## calc correlation length & number effective points? 
    effPtsStr = ''
    if False:
        print 'about to calc # effective points'
        if ((len(outFiles)>1)and(settingsDict['stageList'][-1]=='MCMC'))and (settingsDict['calcCL'] and os.path.exists(allFname)):
            effPtsStr = tools.mcmcEffPtsCalc(allFname)
            print effPtsStr
    
    ## following post-processing stages can take a long time, so write the current
    ## summary information to the summary file and add the rest later
    if False:
        if os.path.exists(allFname):
            print 'about to make summary file'
            tools.summaryFile(settingsDict,settingsDict['stageList'],allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,1,1,'')
        
    ##clean up files (move to folders or delete them)
    if False:
        print 'about to clean up output directory'
        tools.cleanUp(settingsDict,settingsDict['stageList'],allFname)
    if False and settingsDict['CopyToDB']:
        print 'about to copy files to dropbox'
        tools.copyToDB(settingsDict)

if __name__ == '__main__':
    customPost()    

#END OF FILE
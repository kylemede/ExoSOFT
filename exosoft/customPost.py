import tools
from tools import constants as const
from exosoftpath import rootDir as ExoSOFTdir
import sys
import os
import numpy as np
import glob

def customPost():
    settingsDict = tools.startup(sys.argv,ExoSOFTdir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],"combined-BIstripped-MCMCdata.fits")
    #allFname = os.path.join(settingsDict['finalFolder'],"combinedSubSampledMCMC.fits")
    allFname = os.path.join(settingsDict['finalFolder'],"combinedSTdata.fits")
    skipBurnInStrip=True
    if os.path.exists(allFname)==False:
        allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
        skipBurnInStrip=False
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=100)

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
    outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputDataMCMC*_BIstripped.fits")))
    #outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"subSampledDataMCMC*.fits")))
    
    ## calc and strip burn-in?
    burnInStr = ''
    if True:
        if skipBurnInStrip==False:
            print 'about to strip burn-in'
            if (len(outFiles)>1)and(settingsDict['CalcBurn'] and(settingsDict['stageList'][-1]=='MCMC')):
                (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
                if settingsDict['rmBurn'][0]:
                    strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
                    outFiles = strippedFnames
                    ## combine stripped files to make final file?
                    if len(strippedFnames)>0:
                        strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                        tools.combineFits(strippedFnames,strippedAllFname)
                        ## replace final combined filename with new stripped version
                        allFname = strippedAllFname
        
    ## find best fit
    if False:
        print 'about to find best orbit'
        bestFit = tools.findBestOrbit(allFname,bestToFile=False,findAgain=False)
    else:
        bestFit = np.array([0.605838520481,0.146024522998,48.0692312367,98.8893962294,0.640666046184,2445709.98532,2445709.98532,8.2093394771,40.2779824043,166.743572004,3.39834214958,43.3616906776,0.0,9730.87607461,11612.1046117,12586.0296274])
        #best assuming 3 offsets due possibly to a linear trend of a distant companion
        #bestFit = np.array([0.535838520481,0.0466024522998,48.0692312367,98.8893962294,0.640666046184,2445709.98532,2445709.98532,8.2093394771,40.2779824043,166.743572004,3.39834214958,43.3616906776,0.0,9730.87607461,11612.1046117,12586.0296274])
    if True:
        print 'about to make orbit plots'
        ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
        ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
        plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-Manual2')
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[])
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-9.9,9.5],[-0.8,0.8],[-0.515,0.515]],diErrMult=0,diLnThk=1)
        
        ##super cropped HIP10321 trials
        #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-MANUAL-1mult-3thk-crop1')
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[[[-210,425],[-310,325]],[[346.2,350.6],[236.7,241.9]]],RVlims=[[-1250,600],[-65,65],[-0.515,0.515]],diErrMult=1,diLnThk=3)
        #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-MANUAL-1mult-3thk-crop2')
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[[[-210,425],[-310,325]],[[391.3,395.7],[176.6,181.8]]],RVlims=[[-1250,600],[-65,65],[-0.515,0.515]],diErrMult=1,diLnThk=3)
        #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-MANUAL-1mult-3thk-crop3')
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[[[-210,425],[-310,325]],[[413.3,421.1],[124.7,139.6]]],RVlims=[[-1250,600],[-65,65],[-0.515,0.515]],diErrMult=1,diLnThk=3)
        
        #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-MANUAL-1mult-2thk-')
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-1250,620],[-60,60],[-0.515,0.515]],diErrMult=1,diLnThk=2)
        #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbPlot-MANUAL-5mult-2thk-')
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-1250,620],[-60,60],[-0.515,0.515]],diErrMult=5,diLnThk=2)
        
    clStr=''
    if False:
        print 'about to plot posteriors'
        plotFilename = os.path.join(settingsDict['finalFolder'],'posteriors')
        clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],bestVals=bestFit,stage=settingsDict['stageList'][-1], shadeConfLevels=True,forceRecalc=True,plotALLpars=True)
        #for fake jupiter
        #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,4,8],xLims=[[0.5,2.01],[0.5,1.7],[-0.00,0.1],[39,56]],bestVals=bestFit,stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=True)   [345.9,355]
        #for HIP10321
        #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,4,7,2,8,3,9,5,13,14,15],xLims=[[0.7,1.5],[0.19,0.52],[0.36,0.40],[20.5,21.5],[35,40],[150,171],[242,248],[345.9,355],[2452295,2452393],[6170,6220],[360,405],[6300,6350]],bestVals=bestFit,stage=settingsDict['stages'][-1], shadeConfLevels=True,forceRecalc=False)
        #for HIP10321-SIMPLE
        #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,4,7],xLims=[[0.7,1.5],[0.19,0.5],[0.36,0.40],[20.5,21.5]],bestVals=bestFit,stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
    if False: 
        print 'about to make 2D density plot'
        plotFilename = os.path.join(settingsDict['finalFolder'],'m1m2-densPlot-Dec23PASAfit-151226-autoRanges')
        #ranges=[[xMin,xMax],[yMin,yMax]]
        #tools.densityPlotter2D(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]])
        plotFilename = os.path.join(settingsDict['finalFolder'],'m1m2-densPlot-Dec26fit-151226-manualRanges')
        tools.densityPlotter2D(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]],ranges=[[0.64,1.45],[0.18,0.62]])
        #plotFilename = os.path.join(settingsDict['finalFolder'],'cornerPlot-2')
        #tools.cornerPlotter(allFname, plotFilename,paramsToPlot=[0,1],bestVals=[bestFit[0],bestFit[1]])
        
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
        finalFits=''
        nus = [19, 19, 1]
        epochs=[ 2457327.500]
        tools.predictLocation(orbParams,settingsDict,epochs)
        #tools.recheckFit3D(orbParams,settingsDict,finalFits,nus)
    
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
    
def stackedPosteriorsPlotterHackStarter():
    outputDataFilenames = []
    
    #outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-SyntheticJUPITER-3D-20percent-startAtBest-lowEccTrue/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-MCMC-10percent-lowEcc-PDMFm1m2-newParaPrior-SUPERlong2/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-MCMC-5percent-lowEcc-PDMFm1m2-newParaPrior-SUPERlong2/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-MCMC-1percent-lowEcc-PDMFm1m2-newParaPrior-SUPERlong2/combined-BIstripped-MCMCdata.fits')
    
    plotFilename = os.path.join(os.path.abspath('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT'),'stackedPosterior-lowEccTrue-151118')
    tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[1,4,8],xLims=[[0.65,1.5],[0.00,0.10],[30.0,60]],centersOnly=True)
    #tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[],xLims=[])
    #print 'Final stacked plot file written to:\n'+plotFilename
    if True:
        print 'converted to PDF as well'
        os.system("epstopdf "+plotFilename+'.eps')
        
def latexMatplotlibTest():
    import pylab
    plt = pylab.matplotlib.pyplot
    plotFilename = '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/latexMatplotlibTest'
    
    latex=True
    plotFormat = 'eps'   
    #plt.rcParams['ps.useafm']= True
    #plt.rcParams['pdf.use14corefonts'] = True
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        blah=1
        #plt.rc('text', usetex=True)
        #plt.rcParams['text.latex.unicode']=True
        #plt.rcParams['text.latex.preamble'] = '\usepackage{amssymb}'
    #if latex:
    #    #blah=1
    #    plt.rc('text', usetex=True)
    #    plt.rcParams['text.latex.unicode']=True 
    #    plt.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    #    #plt.rcParams['text.latex.preamble'] = '\usepackage{sfmath}' 
    else:
        plt.rc('font',family='serif')
        plt.rc('text', usetex=False)
    
    f = plt.figure()
    subPlot = plt.subplot(111)
    subPlot.scatter(np.arange(100),np.arange(100))
    subPlot.axes.set_xlabel(r'$\Delta \delta$')
    subPlot.axes.set_ylabel(r'$\varpi$')
    ## Save file if requested.
    if plotFilename!='':
        plt.savefig(plotFilename,format=plotFormat)
        print 'density contour plot saved to: '+plotFilename
    plt.close()
    if (plotFormat=='eps') and True:
        print 'converting to PDF as well'
        try:
            os.system("epstopdf "+plotFilename)
        except:
            print "Seems epstopdf failed.  Check if it is installed properly."
        
def paramConverterTest():
    
    ## Make Orbit cpp obj
    Orbit = tools.cppTools.Orbit()
    Orbit.loadStaticVars(0,0,True,False)
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    bestFit = np.array([ 1.08040127e+00,   9.60493873e-04,   5.02111885e+01,
          1.00808667e+02,   6.29123973e-02,   2.45151739e+06,
          2.45151739e+06,   1.19376399e+01,   4.75040520e+01,
          2.71592539e+02,   5.36101436e+00,   4.15438649e+01,
          8.77774317e+00,  -5.73899482e-02]) 
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in bestFit:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    
    es = [0.1,0.25,0.5,0.75,0.99]
    omegas=[1.,45.,90.,135.,180.,225.,270.,315.,360.]
    
    for i in range(len(omegas)):
        params[4]=0.25
        params[9]=omegas[i]
        print 'before: e= '+str(params[4])+', omega= '+str(params[9])
        Orbit.convertParsToRaw(params)
        print 'raw: par[4]= '+str(params[4])+', par[9]= '+str(params[9])
        Orbit.convertParsFromRaw(params)
        print 'after: e= '+str(params[4])+', omega= '+str(params[9])
    
if __name__ == '__main__':
    customPost()
    #stackedPosteriorsPlotterHackStarter()
    #paramConverterTest()
    #latexMatplotlibTest()
    

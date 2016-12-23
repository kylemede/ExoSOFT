#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp  or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
##matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. to further avoid the Display issue.
plt.ioff() #turns off I/O for matplotlib so it doesn't need to plot to screen, which is impossible during ssh/screen sessions.
Polygon =  matplotlib.patches.Polygon  
gridspec =  matplotlib.gridspec
patches = matplotlib.patches
MultLoc = matplotlib.ticker.MultipleLocator
import copy
#import glob
#import shutil
#from math import fabs
#from PyAstronomy.pyasl.asl import lineWidth
import timeit
import scipy.optimize as so
from scipy import ndimage
import warnings
import KMlogger
from astropy import constants as const
from six.moves import range

## import from modules in ExoSOFT  ##
#from . import constants as const
from .model import ExoSOFTmodel, ln_posterior
from .generalTools import getParStrs, timeStrMaker, findBestOrbit, PASAtoEN, jdToGcal, confLevelFinder
from .readWriteTools import loadFits, loadRealData 

days_per_year = 365.2422

warnings.simplefilter("error")
log = KMlogger.getLogger('main.plotTools',lvl=100,addFH=False)  
colorsList =['Red','Orange','Purple','Fuchsia','Crimson','Green','Aqua','DarkGreen','Gold','DarkCyan','OrangeRed','Plum','Chartreuse','Chocolate','Teal','Salmon','Brown','Blue']


def histMakeAndDump(chiSquareds,data,outFilename='',nbins=100,weight=False, normed=False, parRange=False,retHist=False):
    """
    This will make a matplotlib histogram using the input settings, then writing the resulting  
    centers of the bins and number of data points in said bin values to disk, with '.dat' extension
    if not added already.
    
    This function is designed to work with a follow up like histLoadAndPlot_** to produce publication worthy plots.
    """
    #print('data :\n'+repr(data))
    #print('parRange :\n'+repr(parRange))
    if weight:
        ## use the likelihoods as the weights
        theWeights = np.exp(-chiSquareds/2.0)
    else:
        theWeights = np.ones(len(data))      
    if parRange==False:
        (hst,bin_edges) = np.histogram(data,bins=nbins,normed=False,weights=theWeights,density=None)
    elif len(parRange)==2:
        (hst,bin_edges) = np.histogram(data,bins=nbins,range=(parRange[0],parRange[1]),normed=False,weights=theWeights,density=None)
    else:
        log.critical('the range value provided to histMakeAndDump did not have required length of zero or 2.')
    #find center of bins
    if type(bin_edges)!=np.ndarray:
        bin_edges = np.array(bin_edges)
    binCenters = (bin_edges[1:]+bin_edges[:-1])/2.0
    histData=np.zeros((len(hst),2))
    histData[:,0]=binCenters
    histData[:,1]=hst
    #print('hst :\n'+repr(hst))
    if outFilename!='':
        if outFilename[-4:]!='.dat':
            outFilename=outFilename+'.dat'
        np.savetxt(outFilename,histData)
    if retHist:
        return histData


def histLoadAndPlot_StackedPosteriors(plot,outFilename='',xLabel='X',lineColor='k',xLims=False,latex=False,showYlabel=False,parInt=0,centersOnly=False,trueVal=None,lgndStr=''):
    """
    Loads previously plotted histograms that were written to disk by histPlotAndDump, and plot them up 
    in a way that is ready for publication.  This version is to plot a posterior of the same parameter 
    for multiple simulation runs to see how they differ.
    
    It is foreseen that many versions of this function will exist for different specific publication ready plots.
    NOTE: this is extremely similar to histLoadAndPlot_ShadedPosteriors.  Re-factor to remove doubled code!!
    """
    if outFilename[-4:]!='.dat':
        outFilename=outFilename+'.dat'
    histData = np.loadtxt(outFilename)  
    ys=[]
    xs=[]
    maxN = np.max(histData[:,1])
    minSub = 0
    valRange = np.max(histData[:,0])-np.min(histData[:,0])
    ## check if m2 and if it should be in jupiter masses
    if parInt==1:
        if np.max(histData[:,0])<0.02:
            histData[:,0]=histData[:,0]*(const.M_sun.value/const.M_jup.value)
            trueVal = trueVal*(const.M_sun.value/const.M_jup.value)
            valRange = np.max(histData[:,0])-np.min(histData[:,0])
            xLabel='m2 [Mjupiter]'
            if latex:
                xLabel=r'$m_2$ [$M_{J}$]'
    
    if ('[JD]' in xLabel) and (valRange>(np.min(histData[:,0])/1000.0)):
        #a wide spanning T or Tc hist, so convert to gregorian doubles.
        jdDates = histData[:,0]
        xlabelOrig = xLabel
        try:
            for i in range(len(histData[:,0])):
                histData[i,0] = jdToGcal(histData[i,0])
            if latex:
                xLabel=xLabel[:-5]+'year]}$'
            else:
                xLabel = xLabel[:-3]+'year]'   
            valRange = np.max(histData[:,0])-np.min(histData[:,0])
        except:
            #reset to original values and continue
            xLabel=xlabelOrig
            for i in range(len(histData[:,0])):
                histData[i,0] = jdDates[i]
                
    if ('[JD]' in xLabel) or (valRange<(np.min(histData[:,0])/100.0)):
        #must be the To or Tc, so subtract int(min) and add to x-axis label
        #doing this as it doesn't go well allowing matplotlib to do it itself formatting wise                
        minSub = int(np.min(histData[:,0]))
        histData[:,0]-=minSub
        if latex:
            if xLabel not in [r'$e$',r'$\chi^2$']:
                xLabel = xLabel[:-3]+"+"+str(minSub)+"]}$"
            elif xLabel== r'$m_2$ [$M_{J}$]':
                xLabel = xLabel[:-1]+"+"+str(minSub)+"]"
            else:            
                xLabel = xLabel[:-1]+"+"+str(minSub)+"$"
        else:
            if xLabel not in ['e','chiSquared']:
                xLabel = xLabel[:-1]+"+"+str(minSub)+"]"
            else:
                xLabel = xLabel+"+"+str(minSub)
        
    halfBinWidth = (histData[1][0]-histData[0][0])/2.0
    # load up list of x,y values for tops of bins
    for i in range(0,histData.shape[0]):
        ys.append(histData[i][1]/maxN)
        if centersOnly:
            xs.append(histData[i][0]-halfBinWidth)
        else:
            ys.append(histData[i][1]/maxN)
            xs.append(histData[i][0]-halfBinWidth)
            xs.append(histData[i][0]+halfBinWidth)
        
    plot.plot(xs,ys,color=lineColor,linewidth=2,label=lgndStr)
    plot.axes.set_ylim([0.0,1.02])
    if xLims!=False:
        plot.axes.set_xlim(xLims)
    plot.locator_params(axis='x',nbins=4) # maximum number of x labels
    plot.locator_params(axis='y',nbins=5) # maximum number of y labels
    plot.tick_params(axis='x',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.tick_params(axis='y',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.spines['right'].set_linewidth(0.7)
    plot.spines['bottom'].set_linewidth(0.7)
    plot.spines['top'].set_linewidth(0.7)
    plot.spines['left'].set_linewidth(0.7)
    # add axes label
    if showYlabel:
        if latex:
            plot.axes.set_ylabel(r'$\frac{dp}{dx} \times {\rm constant} $',fontsize=27)
        else:
            plot.axes.set_ylabel('dp/dx(*constant)',fontsize=25)
    else:
        plot.axes.set_yticklabels(['','',''])
    fsize=23
    if xLabel in ['e', r'$e$']:
        fsize=fsize+10
    if 'JD' in xLabel:
        fsize=fsize-3
    plot.axes.set_xlabel(xLabel,fontsize=fsize)
    if trueVal is not None:
        try:
            plot.plot([trueVal-minSub,trueVal-minSub],[0.0,1.02],'-',color='k',linewidth=2.0)
        except:
            log.error("Tried to plot a line on the stacked histogram for the true val, but failed")
    
    return plot

def histLoadAndPlot_ShadedPosteriors(plot,outFilename='',confLevels=False,xLabel='X',xLims=False,bestVal=False,latex=False,showYlabel=False,parInt=0):
    """
    Loads previously plotted histograms that were written to disk by histPlotAndDump, and plot them up 
    in a way that is ready for publication.  This is the standard plotter used for plotting simple posteriors
    with shaded regions matching the 68% and 95% confidence.
    
    It is foreseen that many versions of this function will exist for different specific publication ready plots.
    """
    convertJDover10yrs=True
    if outFilename[-4:]!='.dat':
        outFilename=outFilename+'.dat'
    histData = np.loadtxt(outFilename)
    ys=[]
    xs=[]
    maxN = np.max(histData[:,1])
    if maxN == 0:
        maxN = 1.0
    minSub = 0
    valRange = np.max(histData[:,0])-np.min(histData[:,0])
    toOrtc=False
    if ('[JD]' in xLabel):
        toOrtc=True
    ## check if M2 and if it should be in jupiter masses
    if parInt==1:
        if np.max(histData[:,0])<0.02:
            histData[:,0]=histData[:,0]*(const.M_sun.value/const.M_jup.value)
            bestVal = bestVal*(const.M_sun.value/const.M_jup.value)
            valRange = np.max(histData[:,0])-np.min(histData[:,0])
            xLabel='m2 [Mjupiter]'
            if latex:
                xLabel=r'$m_2$ [$M_{J}$]'
            confLevels=confLevels*(const.M_sun.value/const.M_jup.value)
    
    if convertJDover10yrs:
        #print("\n\n'[JD]' in xLabel = "+repr('[JD]' in xLabel))
        #print('valRange = '+repr(valRange)+'\n\n')
        if toOrtc and (valRange>3650):
            #a wide spanning T or Tc hist, so convert to gregorian doubles.
            #print('about to try and convert to years')
            jdDates = histData[:,0]
            xlabelOrig = xLabel
            bestValOrig = bestVal
            try:
                for i in range(len(histData[:,0])):
                    histData[i,0] = jdToGcal(histData[i,0])
                bestVal = jdToGcal(bestVal)
                for i in range(2):
                    for j in range(2):
                        confLevels[i][j] = jdToGcal(confLevels[i][j])
                if latex:
                    xLabel=xLabel[:-5]+'year]}$'
                else:
                    xLabel = xLabel[:-3]+'year]'   
                valRange = np.max(histData[:,0])-np.min(histData[:,0])
            except:
                #reset to original values and continue
                xLabel=xlabelOrig
                bestVal = bestValOrig
                for i in range(len(histData[:,0])):
                    histData[i,0] = jdDates[i]
            #print('made it to end of convertJDover10yrs block')
            ##print('xLabel = '+xLabel)
            #print('new vals :\n'+repr(histData[:,0]))
            #print('valRange = '+str(valRange))
    if ('[JD]' in xLabel) or (valRange<(np.min(histData[:,0])/200.0)):
        #must be the To or Tc, so subtract int(min) and add to x-axis label
        #doing this as it doesn't go well allowing matplotlib to do it itself formatting wise                
       # print "\n\n'[JD]' in xLabel = "+repr('[JD]' in xLabel)
       # print('np.min(histData[:,0])/100.0 = '+repr(np.min(histData[:,0])/100.0)+'\n\n')
        minSub = int(np.min(histData[:,0]))
        histData[:,0]-=minSub
        if latex:
            if xLabel not in [r'$e$',r'$\chi^2$']:
                xLabel = xLabel[:-3]+"+"+str(minSub)+"]}$"
            elif xLabel== r'$m_2$ [$M_{J}$]':
                xLabel = xLabel[:-1]+"+"+str(minSub)+"]"
            else:            
                xLabel = xLabel[:-1]+"+"+str(minSub)+"$"
        else:
            if xLabel not in ['e','chiSquared']:
                xLabel = xLabel[:-1]+"+"+str(minSub)+"]"
            else:
                xLabel = xLabel+"+"+str(minSub)
        
    halfBinWidth = (histData[1][0]-histData[0][0])/2.0
    # load up list of x,y values for tops of bins
    #print('histData = '+repr(histData))
    #print('maxN = '+repr(maxN))
    #print('histData[0][1] = '+repr(histData[0][1]))
    for i in range(0,histData.shape[0]):
        ys.append(histData[i][1]/maxN)
        ys.append(histData[i][1]/maxN)
        xs.append(histData[i][0]-halfBinWidth)
        xs.append(histData[i][0]+halfBinWidth)
    #print('xs = '+repr(xs))
    # load up list of shaded rectangle objects if confidence levels were provided
    recs = []
    if (type(confLevels)==list)or(type(confLevels)==np.ndarray):
        for i in range(0,histData.shape[0]):
            x=histData[i][0]-halfBinWidth+minSub
            # >95% confidence color
            c = 'w'
            if (x>confLevels[1][0])and(x<confLevels[1][1]):
                # 95% confidence color
                c = '#80b3ff'
            if (x>confLevels[0][0])and(x<confLevels[0][1]):
                # 68% confidence color
                c = '#3366ff'
            recs.append(patches.Rectangle(xy=(histData[i][0]-halfBinWidth,0), width=halfBinWidth*2.0,height=histData[i][1]/maxN,facecolor=c, edgecolor=c))
        # draw updated patches on plot
        for rec in recs:
                plot.add_patch(rec)
            
    # draw the top line of hist
    plot.plot(xs,ys,color='k',linewidth=1)
    plot.axes.set_ylim([0.0,1.02])
    if bestVal is not False:
        try:
            if False:
                plot.plot([bestVal-minSub,bestVal-minSub],[0.0,1.02],'--',color='green',linewidth=2)
            else:
                plot.plot([bestVal-minSub,bestVal-minSub],[0.0,1.02],color='k',linewidth=2)
        except:
            log.error("Tried to plot a line on the shaded histogram for the best val, but failed")
    if xLims!=False:
        plot.axes.set_xlim((xLims[0]-minSub,xLims[1]-minSub))
    plot.locator_params(axis='x',nbins=3) # maximum number of x labels
    plot.locator_params(axis='y',nbins=5) # maximum number of y labels
    plot.tick_params(axis='x',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.tick_params(axis='y',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.spines['right'].set_linewidth(0.7)
    plot.spines['bottom'].set_linewidth(0.7)
    plot.spines['top'].set_linewidth(0.7)
    plot.spines['left'].set_linewidth(0.7)
    # add axes label
    if showYlabel:
        if latex:
            plot.axes.set_ylabel(r'$\frac{dp}{dx} \times {\rm constant} $',fontsize=27)
        else:
            plot.axes.set_ylabel('dp/dx(*constant)',fontsize=25)
    else:
        plot.axes.set_yticklabels(['','',''])
    fsize=23
    if xLabel in ['e', r'$e$']:
        fsize=fsize+5
    #if toOrtc:
    #    fsize=fsize
    log.debug('xlabel = '+repr(xLabel)+", fsize = "+str(fsize))
    plot.axes.set_xlabel(xLabel,fontsize=fsize)
    
    return plot

def addRVdataToPlot(subPlot,epochsORphases,RVs,RVerrs,datasetInts=[],alf=1.0,markersize=9,plotErrorBars=False):
    """
    Add '+' markers for the data locations with respective y axis errors 
    shown as the height of the markers. 
    """
    for i in range(0,RVs.shape[0]):
        if RVerrs[i]<1e3:
            xs = [epochsORphases[i],epochsORphases[i]]
            ys = [RVs[i]-RVerrs[i],RVs[i]+RVerrs[i]]
            #print str(RVerrs[i])+", -> ["+str(epochsORphases[i])+", "+str(RVs[i])+']'
            if plotErrorBars:
                subPlot.plot(xs,ys,c='k',linewidth=2,alpha=alf)
            if len(datasetInts)<len(RVs):
                clr = 'red'
            else:
                clr=colorsList[int(datasetInts[i])]
            subPlot.plot(epochsORphases[i],RVs[i],c=clr,marker='.',markersize=markersize)
    return subPlot

def addDIdataToPlot(subPlot,realData,asConversion,errMult=1.0,thkns=1.0,pasa=False):
    """
    To plot a '+' for each data point with width and height matching the errors converted 
    to x,y coords.
    NOTE:
    errMult is a multiplier of the horizontal and vertical error lengths.  
    A value of '1' would original size the error lengths, '2' would be double.
    """
    ## copy realData and kill off parts where DI errors are 1e6
    diData = copy.deepcopy(realData)
    diData = diData[np.where(diData[:,2]<1e6)[0],:]
    ## plot a cross for either DI data format
    if pasa:
        ## convert PASA data and errors into EN versions, calc max and then plot if errMult>0
        (xcenters, E_error, ycenters, N_error)=PASAtoEN(diData[:,1],0,diData[:,3],0)
        (xas, E_error, yas, N_error)=PASAtoEN(diData[:,1]-diData[:,2]*errMult,0,diData[:,3],0)
        (xbs, E_error, ybs, N_error)=PASAtoEN(diData[:,1]+diData[:,2]*errMult,0,diData[:,3],0)
        (xcs, E_error, ycs, N_error)=PASAtoEN(diData[:,1],0,diData[:,3]-diData[:,4]*errMult,0)
        (xds, E_error, yds, N_error)=PASAtoEN(diData[:,1],0,diData[:,3]+diData[:,4]*errMult,0)
        xALL = np.concatenate((xcenters,xas,xbs,xcs,xds))
        yALL = np.concatenate((ycenters,yas,ybs,ycs,yds))
        xmin = np.min(xALL)*asConversion
        xmax = np.max(xALL)*asConversion
        ymin = np.min(yALL)*asConversion
        ymax = np.max(yALL)*asConversion
        if errMult>0.0:
            for i in range(0,len(xas)):
                #print('plotting DI line1: xs '+repr([xas[i]*asConversion,xbs[i]*asConversion])+', ys '+repr([yas[i]*asConversion,ybs[i]*asConversion]))
                subPlot.plot([xas[i]*asConversion,xbs[i]*asConversion],[yas[i]*asConversion,ybs[i]*asConversion],linewidth=thkns,color='k',alpha=1.0)
                #print('plotting DI line2: xs '+repr([xcs[i]*asConversion,xds[i]*asConversion])+', ys '+repr([ycs[i]*asConversion,yds[i]*asConversion]))
                subPlot.plot([xcs[i]*asConversion,xds[i]*asConversion],[ycs[i]*asConversion,yds[i]*asConversion],linewidth=thkns,color='k',alpha=1.0)
    else:
        xmin = np.min(diData[:,1]-diData[:,2])*asConversion
        xmax = np.max(diData[:,1]+diData[:,2])*asConversion
        ymin = np.min(diData[:,3]-diData[:,4])*asConversion
        ymax = np.max(diData[:,3]+diData[:,4])*asConversion
        if errMult>0.0:
            for i in range(0,diData.shape[0]):
                xCent = diData[i,1]*asConversion
                yCent = diData[i,3]*asConversion
                #print('data [x,y] = ['+str(xCent/asConversion)+', '+str(yCent/asConversion)+']')#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                left = xCent-diData[i,2]*asConversion
                right = xCent+diData[i,2]*asConversion
                top = yCent+diData[i,4]*asConversion
                btm = yCent-diData[i,4]*asConversion
                hfWdth = abs(right-left)*errMult*0.5
                hfHgt = abs(top-btm)*errMult*0.5
                subPlot.plot([left-hfWdth,right+hfWdth],[yCent,yCent],linewidth=thkns,color='k',alpha=1.0)
                subPlot.plot([xCent,xCent],[btm-hfHgt,top+hfHgt],linewidth=thkns,color='k',alpha=1.0)
    return (subPlot,[xmin,xmax,ymin,ymax])

def stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[],xLims=[],stage='MCMC',centersOnly=False,plotALLpars=False,trueVals=None,legendStrs=[]):
    """
    This will plot a simple posterior distribution for each parameter in the data files
    stacked ontop of each other for comparison between different runs.
    It can only be ran on folders which have their histograms already calculated and saved in the 
    plotData subfolder!
    
    It is called by the custom function stackedPosteriorsPlotterHackStarter in rePlot.py.
    
    Very similar to summaryPlotter, but the loop is first in parameter int, then
    file to ensure all are stacked on same subplot properly.  
    NOTE: might be able to extract doubled code to clean things up...
    """
    log.setStreamLevel(lvl=30)
    latex=True
    plotFormat='eps'
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
    
    if type(outputDataFilenames)!=list:
        outputDataFilenames = [outputDataFilenames]
    
    colorsList =['Blue','Red','Black','Chocolate','Purple','Fuchsia','Crimson','Aqua','Gold','OrangeRed','Plum','Chartreuse','Chocolate','Teal','Salmon','Brown']
    colorsList = ['blue','magenta','sandybrown','lime']
    
    #colorsList =['Blue','#ff751a']
    #colorsList =['Green','Blue','Red']
    colorsList2 = []
    while len(outputDataFilenames)>len(colorsList2):
        for color in colorsList:
            colorsList2.append(color)
    colorsList = colorsList2
    
    if os.path.exists(outputDataFilenames[0]):  
        log.debug('\nCreating a simple plot of some key posteriors for files:\n'+repr(outputDataFilenames))
        log.debug("writing resulting figure to:\n"+plotFilename)
        
        ## load first data file to get param lists 
        (head,data) = loadFits(outputDataFilenames[0])
        ## get parameter lists and filter accordingly
        (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=plotALLpars)
        (paramList2,paramStrs2,paramFileStrs2) = getParStrs(head,latex=False,getALLpars=plotALLpars)
        # modify x labels to account for DI only situations where M1=Mtotal
        if (np.var(data[:,1])==0)and (0 in paramList):
            paramStrs2[0] = 'm total [Msun]'
            paramStrs[0] = r'$m_{\rm total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'm-total'
        # check if a subset is to be plotted or the whole set
        # remake lists of params to match subset.
        if len(paramsToPlot)!=0:
            if plotALLpars==False:
                (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=True)
                (paramList2,paramStrs2,paramFileStrs2) = getParStrs(head,latex=False,getALLpars=True)
                paramStrs2Use = []
                paramStrsUse = []
                paramFileStrsUse = []
                paramListUse = []
                for par in paramsToPlot:
                    paramStrs2Use.append(paramStrs2[par])
                    paramStrsUse.append(paramStrs[par])
                    paramFileStrsUse.append(paramFileStrs[par])
                    paramListUse.append(par)
                paramStrs2 = paramStrs2Use
                paramStrs = paramStrsUse
                paramFileStrs = paramFileStrsUse 
                paramList = paramListUse
            else:
                s = "\nSpecific params to plot were provided, yet the plotALLpars flag was set True."
                s+="\nPlease set it to False if you do not want to plot ALL params."
                s+="\nALL params will be plotted."
                log.critical(s)
        
        ## determine appropriate figure size for number of params to plot
        figSizes =  [(7,4),(8,4),(9,4),(10,4), (12,8),(12,11),(12,13),(12,15)]
        gridSizes = [(1,1),(1,2),(1,3),(1,4),  (2,4),(3,4),(4,4),(5,4)]
        sz = 0
        if len(paramStrs2)>20:
            log.critical("summaryPlotter is only capable of handling up to 20 parameters, "+str(len(paramStrs2))+" were passed in!")
        elif len(paramStrs2)>16:
            sz = 7
        elif len(paramStrs2)>12:
            sz = 6
        elif len(paramStrs2)>8:
            sz = 5
        elif len(paramStrs2)>4:
            sz = 4
        elif len(paramStrs2)>3:
            sz = 3
        elif len(paramStrs2)>2:
            sz = 2
        elif len(paramStrs2)>1:
            sz = 1
        ## Create empty figure to be filled up with plots
        stackedFig = plt.figure(figsize=figSizes[sz],tight_layout=True)
        #stackedFig = plt.figure(figsize=figSizes[sz])     
        
        ## Go through params and re-load the hist for each files and plot them 
        for i in range(0,len(paramStrs2)):
            colorInt = 0
            log.debug('Starting to plot stacked hist for '+paramStrs2[i])
            subPlot = stackedFig.add_subplot(gridSizes[sz][0],gridSizes[sz][1],i+1)
            xLim=False
            if len(paramsToPlot)!=0:
                xLim=xLims[i]
            showYlabel=False
            if i in [0,4,8,12]:
                showYlabel = True
            par=0
            try:
                par = paramList[i]
                trueVal = None
                if trueVals is not None:
                    trueVal = trueVals[paramList[i]]
                #print str(i)
                #print str(trueVal)
                #print repr(trueVals)
                #print repr(paramList)
            except:
                log.warning("Parameter "+str(i)+" not in paramList: \n"+repr(paramList))
            ## go through each file and plot this param's hist on same plot
            for outputDataFilename in outputDataFilenames:
                log.debug('Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]+" for file:\n"+outputDataFilename)
                ## check if plot data dir exists
                plotDataDir = os.path.join(os.path.dirname(outputDataFilename),"plotData")
                if os.path.exists(plotDataDir)==False:      
                    log.critical("PlotDataDir doesn't exist!! at:\n"+plotDataDir)
                else:
                    plotDataDir+='/'
                    histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+paramFileStrs[i])
                    if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')):
                        histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i])
                    if os.path.exists(histDataBaseName+'.dat'):
                        log.debug("plotting file:\n"+histDataBaseName)
                        lgndStr=''
                        #print repr(len(legendStrs))
                        if len(legendStrs)>=colorInt:
                            lgndStr=legendStrs[colorInt]
                            #print lgndStr
                        if colorInt>len(colorsList):
                            log.warning("More plots requested than colors available in colorsList!! "+str(len(colorsList))+' < '+str(colorInt))
                        subPlot = histLoadAndPlot_StackedPosteriors(subPlot,outFilename=histDataBaseName,xLabel=paramStrs[i],lineColor=colorsList[colorInt],xLims=xLim,latex=latex,showYlabel=showYlabel,parInt=par,centersOnly=centersOnly,trueVal=trueVal,lgndStr=lgndStr)
                    else:
                        log.debug("Not plotting hist for "+paramStrs2[i]+" as its hist file doesn't exist:\n"+histDataBaseName)
                    colorInt+=1
                    if True:
                        combinedLbl = 'combined'
                        if (colorInt==4) and (paramStrs[i]=='$i{\\rm  [deg]}$'):
                            histDataBaseName = '/mnt/HOME/MEGA/ExoSOFT-outputCopies/after-June14-2016/ClusterRuns/CombinedNormalized-IncHist'
                            #histDataBaseName = '/mnt/HOME/MEGA/ExoSOFT-outputCopies/after-June14-2016/ClusterRuns/chiOne/incHist'
                            #print 'about to try and add combined hist with file '+histDataBaseName
                            subPlot = histLoadAndPlot_StackedPosteriors(subPlot,outFilename=histDataBaseName,xLabel=paramStrs[i],lineColor='black',xLims=xLim,latex=latex,showYlabel=showYlabel,parInt=par,centersOnly=centersOnly,trueVal=trueVal,lgndStr=combinedLbl)
                            #print 'special combined posterior plotted it seems'
                        elif colorInt==4:
                            subPlot.plot([-100,-90],[-100,-90],color='k',lineWidth=2,label=combinedLbl)
        subPlot.legend()
        #plt.tight_layout()        
        ## Save file if requested.
        log.debug('\nStarting to save stacked param hist figure:')
        if plotFilename!='':
            plt.savefig(plotFilename+'.'+plotFormat,fl_format=plotFormat)
            s= 'stacked hist plot saved to: '+plotFilename+'.'+plotFormat
            log.info(s)
        plt.close()
        #print 'figure written to: '+plotFilename
        if plotFormat=='eps':
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename+'.'+plotFormat)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")

        
def summaryPlotter(outputDataFilename,plotFilename,paramsToPlot=[],xLims=[],bestVals=[],stage='MCMC',shadeConfLevels=True,forceRecalc=True,plotALLpars=False,nbins=None):
    """
    This advanced plotting function will plot all the data in a grid on a single figure.  The data will be plotted
    in histograms that will be normalized to a max of 1.0.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    
    NOTE: currently a maximum of 20 parameters can be plotted, ie. for a system with up to 7 RV data sets.
    """
    latex=True
    plotFormat='eps'
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
    ## check if plot data dir exists, else make it
    plotDataDir = os.path.join(os.path.dirname(outputDataFilename),"plotData")
    #print 'plotDataDir = '+plotDataDir
    if os.path.exists(plotDataDir)==False:      
        os.mkdir(plotDataDir)
    plotDataDir+='/'
    
    #print 'about to load fits file'
    (head,data) = loadFits(outputDataFilename)
    if head!=False:  
        log.debug(' Inside summaryPlotter')
        s= '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        log.info(s)
        
        ## check if the passed in value for plotFilename includes fl_format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            log.debug('updating plotFilename to:\n'+plotFilename)
        else:
            plotFilename = plotFilename
        #print "about to try and extract parstrs from header"  #$$$$$$$$$$$$$$$$$$$$$$$$   
        ## get parameter lists and filter accordingly
        (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=plotALLpars)
        (paramList2,paramStrs2,paramFileStrs2) = getParStrs(head,latex=False,getALLpars=plotALLpars)
        #print "got extract parstrs from header"  #$$$$$$$$$$$$$$$$$$
        #print repr((paramList,paramStrs,paramFileStrs))#$$$$$$$$$$$$$$$$$$
        #print repr((paramList2,paramStrs2,paramFileStrs2))#$$$$$$$$$$$$$$$$$$
        # modify x labels to account for DI only situations where M1=Mtotal
        if (np.var(data[:,1])==0)and (0 in paramList):
            paramStrs2[0] = 'm total [Msun]'
            paramStrs[0] = r'$m_{\rm total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'm-total'
        #print 'updated m1 to mtotal'
        # check if a subset is to be plotted or the whole set
        # remake lists of params to match subset.
        if len(paramsToPlot)!=0:
            if plotALLpars==False:
                (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=True)
                (paramList2,paramStrs2,paramFileStrs2) = getParStrs(head,latex=False,getALLpars=True)
                paramStrs2Use = []
                paramStrsUse = []
                paramFileStrsUse = []
                paramListUse = []
                bestValsUse = []
                for par in paramsToPlot:
                    paramStrs2Use.append(paramStrs2[par])
                    paramStrsUse.append(paramStrs[par])
                    paramFileStrsUse.append(paramFileStrs[par])
                    paramListUse.append(par)
                    bestValsUse.append(bestVals[par])
                paramStrs2 = paramStrs2Use
                paramStrs = paramStrsUse
                paramFileStrs = paramFileStrsUse 
                paramList = paramListUse
                bestVals = bestValsUse
            else:
                s = "\nSpecific params to plot were provided, yet the plotALLpars flag was set True."
                s+="\nPlease set it to False if you do not want to plot ALL params."
                s+="\nALL params will be plotted."
                log.critical(s)
            #print 'cleaned up par and str lists'
        # just the varying params are to be plotted, so clean the bestVals list to match
        elif len(bestVals)>len(paramList):
            bestValsUse = []
            for par in paramList:
                bestValsUse.append(bestVals[par])
            bestVals = bestValsUse
        log.debug('\nparamStrs2 = '+repr(paramStrs2)+'\nparamStrs = '+repr(paramStrs)+'\nparamFileStrs = '+repr(paramFileStrs)+'\nparamList = '+repr(paramList)+'\nbestVals = '+repr(bestVals)+'\n')
        ## determine appropriate figure size for number of params to plot
        figSizes =  [(5.5,7),(8,7),(9,7),(10,3.5),(11,8),(11,11),(11,14),(11,16)]
        gridSizes = [(1,1),(1,2),(1,3),(1,4),(2,4),(3,4),(4,4),(5,4)]
        sz = 0
        if len(paramStrs2)>20:
            log.critical("summaryPlotter is only capable of handling up to 20 parameters, "+str(len(paramStrs2))+" were passed in!")
        elif len(paramStrs2)>16:
            sz = 7
        elif len(paramStrs2)>12:
            sz = 6
        elif len(paramStrs2)>8:
            sz = 5
        elif len(paramStrs2)>4:
            sz = 4
        elif len(paramStrs2)>3:
            sz = 3
        elif len(paramStrs2)>2:
            sz = 2
        elif len(paramStrs2)>1:
            sz = 1
        
        ## run through all the data files and parameters requested and make histogram files
        completeCLstr = '-'*22+'\nConfidence Levels are:\n'+'-'*80+'\n'
        for i in range(0,len(paramList)):
            #print paramStrs2[i]
            if (os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat'))==False)or forceRecalc:
                log.debug('Checking parameter has useful data '+str(i)+"/"+str(len(paramStrs2)-1)+": "+paramStrs2[i]+", for file:\n"+outputDataFilename)
                weightHists = False
                if stage in ['SA','MC']:
                    (CLevels,data,chiSquareds,bestDataVal,clStr) = confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=True, returnBestDataVal=True)
                    weightHists = True
                else:
                    (CLevels,data,bestDataVal,clStr) = confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=True)
                    chiSquareds = []
                #print 'ln704'
                if bestDataVal!=0:
                    completeCLstr+=paramStrs2[i]+clStr+'\n'+'-'*80+'\n'
                    log.debug('Making hist file for parameter '+str(i)+"/"+str(len(paramStrs2)-1)+": "+paramStrs2[i]+", for file:\n"+outputDataFilename)
                    histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')
                    #print 'histDataBaseName = '+histDataBaseName
                    if nbins is None:
                        if stage=='MC':
                            nbins = 50
                        else:
                            nbins = 100
                    xLim=False
                    if len(xLims)>0:
                        xLim = xLims[i] 
                    #print 'ln715'
                    histMakeAndDump(chiSquareds,data,outFilename=histDataBaseName,nbins=nbins,weight=weightHists, normed=False,parRange=xLim)
                    #print 'ln721'
                    if (os.path.exists(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat'))==False)or forceRecalc:
                        np.savetxt(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat'),CLevels)
                        log.debug('confidence levels data stored to:\n'+os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat'))
                else:
                    log.debug("Nope! no useful data for "+str(i)+"/"+str(len(paramStrs2)-1)+": "+paramStrs2[i])#+", in file:\n"+outputDataFilename)
        #print 'done making hists files'  #$$$$$$$$$$$$$$$$$$$$$$$$
        ## Create empty figure to be filled up with plots
        sumFig = plt.figure(figsize=figSizes[sz],tight_layout=True)       
        ## make shaded posterior for each param
        for i in range(0,len(paramStrs2)):
            histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')
            #print '\n\n'+os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')+'\n'
            if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')):
                histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')
            if os.path.exists(histDataBaseName):
                log.debug('Starting to plot shaded hist for '+paramStrs2[i])
                subPlot = sumFig.add_subplot(gridSizes[sz][0],gridSizes[sz][1],i+1)#gs[i])
                #print '\nLoading and re-plotting parameter '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]
                log.debug('Loading and re-plotting parameter '+str(i)+"/"+str(len(paramStrs2)-1)+": "+paramStrs2[i])#+" for file:\n"+outputDataFilename)
                CLevels=False
                xLim=False
                if len(xLims)>0:
                    xLim = xLims[i]   
                bestVal = False
                if len(bestVals)>0:
                        bestVal = bestVals[i]                        
                if shadeConfLevels:
                    clFile = os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')
                    if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')):
                        clFile = os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')
                    CLevels=np.loadtxt(clFile)
                    #print 'CLevels = '+repr(CLevels)
                showYlabel=False
                if i in [0,4,8,12]:
                    showYlabel = True
                par=0
                try:
                    par = paramList[i]
                except:
                    log.debug("Parameter "+str(i)+" not in paramList: \n"+repr(paramList))
                #print 'about to make hist plot for file base '+histDataBaseName
                subPlot = histLoadAndPlot_ShadedPosteriors(subPlot,outFilename=histDataBaseName,confLevels=CLevels,xLabel=paramStrs[i],xLims=xLim,bestVal=bestVal,latex=latex,showYlabel=showYlabel,parInt=par)         
                log.debug('Done to plot shaded hist for '+paramStrs2[i])
            else:
                log.debug("Not plotting shaded hist for "+paramStrs2[i]+" as its hist file doesn't exist:\n"+histDataBaseName)
        #print 'done making posterior plots files'  #$$$$$$$$$$$$$$$$$$$$$$$$
        #plt.tight_layout()
        ## Save file if requested.
        log.debug('\nStarting to save param hist figure:')
        if plotFilename!='':
            plt.savefig(plotFilename,fl_format=plotFormat)
            s= 'Summary plot saved to: '+plotFilename
            log.info(s)
        plt.close()
        if (plotFormat=='eps') and True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
        
        return completeCLstr
            
def star(R, x0, y0, color='w', N=5, thin = 0.5):
    """
    Returns an N-pointed star of size R at (x0, y0) (matplotlib patch).
    NOTE: code base taken from star_patch.py in 'Beginning-Python-Visualization'
    """
    polystar = np.zeros((2*N, 2))
    for i in range(2*N):
        angle = i*np.pi/N
        r = R*(1-thin*(i%2))
        polystar[i] = [r*np.cos(angle)+x0, r*np.sin(angle)+y0]
    return Polygon(polystar, fc=color, ec='black',linewidth=1.5)  

def epochsToPhases(epochs,Tc,P_yrs, halfOrbit=False):
    """
    Convert the epochs (from a realData ary) into phase values 
    (ratio of how far from Tc it is), shifted to lie inside [0,1].
    if 'halfOrbit'=True, the vals will lie inside [-0.5,0.5].
    """    
    verbose=False         
    phases = []
    P_days = P_yrs*days_per_year
    for epoch in epochs:
        phaseTimeDiff = epoch - int((epoch-Tc)/P_days)*P_days-Tc #phase shifted into [Tc,Tc+P]
        if verbose:
            print(str(epoch)+" - "+str(int((epoch-Tc)/P_days)*P_days)+" - "+str(Tc)+" = "+str(phaseTimeDiff))
        phase = phaseTimeDiff/P_days#phase shifted into [0,1]
        if halfOrbit:
            if phase>0.5:
                phase = phase-1.0#phase shifted into [-0.5,0.5]
            elif phase<-0.5:
                phase = phase+1.0#phase shifted into [-0.5,0.5]
        phases.append(phase)
        if verbose:
            print('\nepoch = ',epoch)
            print('period [days] = ',P_days)
            print('phase = ',phase)  
    return phases

def orbitPlotter(orbParams,settings,plotFnameBase="",fl_format='png',DIlims=[],RVlims=[],diErrMult=1,diLnThk=1.0,legendStrs=[]):
    """
    Make both the DI and RV plots.
    '-DI.png' and/or '-RV.png' will be added to end of plotFnameBase 
    to make the filenames for each type of plot.
    
    Optional tweaks:
    DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]
    RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
    """
    # There is some custom code to place 'zoom-in' inserts that some could 
    # modify to use with their work.  Review current code and tweak to match 
    # your situaiton accordingly.
    plotCustomInsets = False
    savePlotDataToFile = True
    autoUnits = True
    latex=True
    colorsList = ['blue','magenta','sandybrown','lime']
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
    log.debug("Starting to make orbit plots")
    

    ## check if plot data dir exists, else make it
    plotDataDir = os.path.join(os.path.dirname(plotFnameBase),"plotData")
    if os.path.exists(plotDataDir)==False:      
        os.mkdir(plotDataDir)
    plotDataDir+='/'
    
    # prepare variables and lists for multi-planet plots
    if type(orbParams[0]) not in [list,np.ndarray]:
        orbParams = [orbParams]
    if type(orbParams[0][0]) not in [list,np.ndarray]:
        orbParams = [orbParams]
    if type(settings)==dict:
        settings = [settings]
        
    
    try:
        pasa = settings[0]["pasa"]
    except:
        pasa = False
    
    if len(settings)==len(orbParams):
        log.debug( 'all good, both settings and orbParams have length '+str(len(settings)))
    else:
        print('PROBLEM! number of settings files were '+str(len(settings))+', while number of orbParam sets were '+str(len(orbParams)))
        #print repr(orbParams)
        
    ################
    # Make DI plot #
    ################
    #settings are in a list, orbit params are in a double list
    #[settings] and [[[par1,par2,...]]]
    #[settings1,settings2] and [[[par1,par2,...],[par1,par2,...],...],[[par1,par2,...],[par1,par2,...],...]]
    if settings[0]['data_mode']!='RV':
        diFig = plt.figure(2,figsize=(10,9))
        diMain = diFig.add_subplot(111)
        for settInt in range(0,len(settings)):            
            #print('plotting for settings file # '+str(settInt))
            settingsDI = settings[settInt]
            orbParamsDI = orbParams[settInt]
            
            ## instantiate Model class for this set of settings
            Model = ExoSOFTmodel(settingsDI)
            real_epochs_di = copy.deepcopy(Model.Data.epochs_di)
            real_decsa = copy.deepcopy(Model.Data.decsa)
            real_rapa = copy.deepcopy(Model.Data.rapa)
            
            ##get the real data
            realData = loadRealData(diFilename=settingsDI['di_dataFile'],rvFilename=settingsDI['rv_dataFile'],dataMode=settingsDI['data_mode'])
            for parSetInt in range(0,len(orbParamsDI)):
                #print('plotting orbit parameter set # '+str(parSetInt))
                #convert 'stored' to 'direct/raw' versions
                paramsDIraw = copy.deepcopy(Model.Params.stored_to_direct(orbParamsDI[parSetInt]))
                paramsDI = orbParamsDI[parSetInt]
                realDataDI = copy.deepcopy(realData)
                realDataDI = realDataDI[np.where(realDataDI[:,2]<1e6)[0],:]
                
                # ExoSOFTmodel will calculate predicted/ExoSOFTmodel data for each epoch in Model.Data.epochs_di
                ## calculate the fit locations for the DI epochs to calculate 0-C
                _ = ln_posterior(paramsDIraw, Model)
                predicted_decsa_model = copy.deepcopy(Model.Data.decsa_model)
                predicted_rapa_model = copy.deepcopy(Model.Data.rapa_model)
                
                ##Make ExoSOFTmodel data for 100~1000 points for plotting fit
                nPts = 1000
                fakeEpochs = np.zeros((nPts),dtype=np.dtype('d'))
                Model.Data.rapa = np.ones((nPts),dtype=np.dtype('d'))
                Model.Data.rapa_err = np.ones((nPts),dtype=np.dtype('d'))
                Model.Data.decsa = np.ones((nPts),dtype=np.dtype('d'))
                Model.Data.decsa_err = np.ones((nPts),dtype=np.dtype('d'))
                Model.Data.rapa_model = np.ones((nPts),dtype=np.dtype('d'))
                Model.Data.decsa_model = np.ones((nPts),dtype=np.dtype('d'))
                for i in range(0,nPts-1):
                    fakeEpochs[i] = paramsDI[5]+(days_per_year*paramsDI[7]*(i/float(nPts)))
                fakeEpochs[nPts-1]  = fakeEpochs[0]+days_per_year*paramsDI[7]
                Model.Data.epochs_di = fakeEpochs
                _ = ln_posterior(paramsDIraw, Model)
                fit_epochs = copy.deepcopy(Model.Data.epochs_di)
                fit_decsa_model = copy.deepcopy(Model.Data.decsa_model)
                fit_rapa_model = copy.deepcopy(Model.Data.rapa_model)
                
                ## Get locations of start/end for semi-major axis or COM, and AN/DN for line-of-nodes
                ## Get 1/4 locations (useful for drawing semi-major axis, and finding loc of COM)
                fakeQuarterEpochs = np.zeros((4),dtype=np.dtype('d'))
                Model.Data.rapa = np.ones((4),dtype=np.dtype('d'))
                Model.Data.rapa_err = np.ones((4),dtype=np.dtype('d'))
                Model.Data.decsa = np.ones((4),dtype=np.dtype('d'))
                Model.Data.decsa_err = np.ones((4),dtype=np.dtype('d'))
                Model.Data.rapa_model = np.ones((4),dtype=np.dtype('d'))
                Model.Data.decsa_model = np.ones((4),dtype=np.dtype('d'))
                for i in range(0,4):
                    fakeQuarterEpochs[i] = paramsDI[5]+(days_per_year*paramsDI[7]*(i/4.0))
                # Calculate orbit for quarter epochs
                Model.Data.epochs_di = fakeQuarterEpochs
                _ = ln_posterior(paramsDIraw, Model)
                quarter_epochs = copy.deepcopy(Model.Data.epochs_di)
                quarter_decsa_model = copy.deepcopy(Model.Data.decsa_model)
                quarter_rapa_model = copy.deepcopy(Model.Data.rapa_model)
            
            
                ## make semi-major locs
                semiMajorLocs = np.array([[quarter_rapa_model[0],quarter_decsa_model[0]] , [quarter_rapa_model[2],quarter_decsa_model[2]]])
                ## find loc of COM for possible use
                xCOM = (quarter_rapa_model[3]+quarter_rapa_model[0])/2.0
                yCOM = (quarter_decsa_model[3]+quarter_decsa_model[0])/2.0
                
                ## Find Ascending and Descending Node locations
                nodeEpochs = nodeEpochsCalc(paramsDI,settingsDI["omegaFdi"]) 
                #print('period/2 = '+repr(days_per_year*paramsDI[7]*(1.0/2.0)))
                Model.Data.epochs_di = np.array(nodeEpochs,dtype=np.dtype('d'))
                Model.Data.rapa = np.ones((2),dtype=np.dtype('d'))
                Model.Data.rapa_err = np.ones((2),dtype=np.dtype('d'))
                Model.Data.decsa = np.ones((2),dtype=np.dtype('d'))
                Model.Data.decsa_err = np.ones((2),dtype=np.dtype('d'))
                Model.Data.rapa_model = np.ones((2),dtype=np.dtype('d'))
                Model.Data.decsa_model = np.ones((2),dtype=np.dtype('d'))
                _ = ln_posterior(paramsDIraw, Model)
                lon_decsa_model = copy.deepcopy(Model.Data.decsa_model)
                lon_rapa_model = copy.deepcopy(Model.Data.rapa_model)
                lonXYs = np.array([[lon_rapa_model[0],lon_decsa_model[0]],[lon_rapa_model[1],lon_decsa_model[1]]])
                
                ##load resulting data to file for re-plotting by others, along with calculating and storing O-C values
                #real [x,xerr,y,yerr] OR [PA,PAerr,SA,SAerr] depending on 'pasa' bool in settingsDI dict.
                outDIdataReal = realDataDI[:,1:5]
                #fit [x,y]
                outDIdataFit = []
                for i in range(0,len(fit_epochs)):
                    outDIdataFit.append([fit_rapa_model[i],fit_decsa_model[i]])
                outPredictedDIdata = []
                for i in range(0,len(real_epochs_di)):
                    outPredictedDIdata.append([predicted_rapa_model[i],predicted_decsa_model[i]])
                fnameBase = os.path.join(os.path.dirname(plotDataDir),'DIplotData')
                if pasa:
                    hReal = "[PA,PAerr,SA,SAerr]"
                    hFit = hPredicted = '[PA,SA]'
                else:
                    hReal = "[x,xerr,y,yerr]"
                    hFit = hPredicted = '[x,y]'
                residualDIdata = []
                if pasa:
                    (xcenters, _, ycenters, _)= PASAtoEN(real_rapa[:],0,real_decsa[:],0)
                    for i in range(0,len(real_epochs_di)):
                        residualDIdata.append([xcenters[i]-predicted_rapa_model[i],ycenters[i]-predicted_decsa_model[i]])
                else:
                    for i in range(0,len(real_epochs_di)):
                        residualDIdata.append([real_rapa[i]-predicted_rapa_model[i],real_decsa[i]-predicted_decsa_model[i]])
                if savePlotDataToFile:
                    np.savetxt(fnameBase+'-real'+str(parSetInt)+'.dat',outDIdataReal,header=hReal)
                    np.savetxt(fnameBase+'-predicted'+str(parSetInt)+'.dat',outPredictedDIdata,header=hPredicted)
                    np.savetxt(fnameBase+'-fit'+str(parSetInt)+'.dat',outDIdataFit,header=hFit)
                    #O-C [RAo-c,DECo-c]
                    np.savetxt(fnameBase+'-O-C'+str(parSetInt)+'.dat',residualDIdata,header="[RAo-c,DECo-c]")
                #np.savetxt(fnameBase+'-fit'+str(parSetInt)+'-2.dat',outDIdataFit,header=hFit)
                    
                #determine if to plot [mas] or ["]
                asConversion=1.0
                unitStr = ' ["]'
                if autoUnits:
                    if abs(np.min([np.min(real_rapa[:]),np.min(real_decsa[:])]))<1.0:
                        asConversion = 1000.0
                        unitStr = ' [mas]'
                ## Draw orbit fit
                clr = colorsList[settInt]
                lgndStr =''
                try:
                    lgndStr = legendStrs[settInt]
                except:
                    log.debug('legendStrs param passed in was not cool')
                diMain.plot(fit_rapa_model[:]*asConversion,fit_decsa_model[:]*asConversion,linewidth=diLnThk,color=clr,label=lgndStr)
                if len(orbParamsDI)==1: 
                    ## Draw line-of-nodes
                    diMain.plot(lonXYs[:,0]*asConversion,lonXYs[:,1]*asConversion,'-.',linewidth=diLnThk,color='Green')
                    #print('lonXYs*asConversion = '+repr(lonXYs*asConversion))
                    ## Draw Semi-Major axis
                    diMain.plot(semiMajorLocs[:,0]*asConversion,semiMajorLocs[:,1]*asConversion,'-',linewidth=diLnThk,color='Green')
                ## Draw larger star for primary star's location
                atot = asConversion*np.sqrt((semiMajorLocs[:,0][0]-semiMajorLocs[:,1][0])**2 + (semiMajorLocs[:,0][0]-semiMajorLocs[:,1][0])**2)
                #print('5*paramsDI[10] = '+str(5*paramsDI[10]))
                #print('atot/15.0 = '+str(atot/15.0))
                starWidth = atot/15.0
                #old width was 5*paramsDI[10]
                starPolygon = star(starWidth,0,0,color='yellow',N=6,thin=0.5)
                if parSetInt<1:
                    diMain.add_patch(starPolygon)
                ## plot predicted locations for the data points
                if True:
                    for i in range(0,len(predicted_rapa_model)):
                        diMain.plot(predicted_rapa_model[i]*asConversion,predicted_decsa_model[i]*asConversion,c='red',marker='.',markersize=diLnThk*11)#$$$$$$$$ Place for custimization
                        #print('plotted point ['+str(predictedDataDI[i,0]*asConversion)+', '+str(predictedDataDI[i,1]*asConversion)+']')
                ## Add DI data to plot
                (diMain,[xmin,xmax,ymin,ymax]) =  addDIdataToPlot(diMain,realDataDI,asConversion,errMult=diErrMult,thkns=diLnThk,pasa=pasa)#$$$$$$$$ Place for custimization
                ## set limits and other basics of plot looks
                
                ## update lims with custom values if provided
                if len(DIlims)>0:
                    #DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]
                    xLimsF=(DIlims[0][0][0],DIlims[0][0][1])
                    yLimsF=(DIlims[0][1][0],DIlims[0][1][1])
                    xLimsC=(DIlims[1][0][0],DIlims[1][0][1])
                    yLimsC=(DIlims[1][1][0],DIlims[1][1][1])
                else:
                    xLimsF = (np.min([xmin,np.min(real_rapa[:]*asConversion)]),np.max([xmax,np.max(real_rapa[:]*asConversion)]))
                    yLimsF = (np.min([ymin,np.min(real_decsa[:]*asConversion)]),np.max([ymax,np.max(real_decsa[:]*asConversion)]))
                    #print('xLimsF = '+repr(xLimsF))
                    #print('yLimsF = '+repr(yLimsF))
                    xLimsC = (xmin,xmax)
                    yLimsC = (ymin,ymax)
                ##fix FULL limits
                #force these to be square
                xR = xLimsF[1]-xLimsF[0]
                yR = yLimsF[1]-yLimsF[0]
                if xR>yR:
                    yLimsF = (yLimsF[0]-0.5*abs(xR-yR),yLimsF[1]+0.5*abs(xR-yR))
                elif xR<yR:
                    xLimsF = (xLimsF[0]-0.5*abs(xR-yR),xLimsF[1]+0.5*abs(xR-yR))
                #pad by 5%
                xLimsFull = (xLimsF[0]-(xLimsF[1]-xLimsF[0])*0.05,xLimsF[1]+(xLimsF[1]-xLimsF[0])*0.05)
                yLimsFull = (yLimsF[0]-(yLimsF[1]-yLimsF[0])*0.05,yLimsF[1]+(yLimsF[1]-yLimsF[0])*0.05)     
                #print('Full DI plot ranges: '+repr(xLimsFull[1]-xLimsFull[0])+' X '+repr(yLimsFull[1]-yLimsFull[0]))
                ##fix CROPPED limits
                #force these to be square
                xR = xLimsC[1]-xLimsC[0]
                yR = yLimsC[1]-yLimsC[0]
                if xR>yR:
                    yLimsC = (yLimsC[0]-0.5*abs(xR-yR),yLimsC[1]+0.5*abs(xR-yR))
                elif xR<yR:
                    xLimsC = (xLimsC[0]-0.5*abs(xR-yR),xLimsC[1]+0.5*abs(xR-yR))
                #pad by 5%
                xLimsCrop = (xLimsC[0]-(xLimsC[1]-xLimsC[0])*0.05,xLimsC[1]+(xLimsC[1]-xLimsC[0])*0.05)
                yLimsCrop = (yLimsC[0]-(yLimsC[1]-yLimsC[0])*0.05,yLimsC[1]+(yLimsC[1]-yLimsC[0])*0.05)     
                #print 'cropped DI plot ranges: '+repr(xLimsCrop[1]-xLimsCrop[0])+' X '+repr(yLimsCrop[1]-yLimsCrop[0])
                
        ## FLIP X-AXIS to match backawards Right Ascension definition
        a = diMain.axis()
        diMain.axis([a[1],a[0],a[2],a[3]])
        plt.minorticks_on()
        diMain.tick_params(axis='both',which='major',width=1,length=5,pad=10,direction='in',labelsize=30)
        diMain.tick_params(axis='both',which='minor',width=1,length=2,pad=10,direction='in')
        diMain.spines['right'].set_linewidth(1.0)
        diMain.spines['bottom'].set_linewidth(1.0)
        diMain.spines['top'].set_linewidth(1.0)
        diMain.spines['left'].set_linewidth(1.0)
        #[left,btm,width,height]
        limsList = [yLimsCrop[0],yLimsCrop[1],yLimsFull[0],yLimsFull[1]]
        if (-999>np.min(limsList))or(999<np.max(limsList)):
            diMain.set_position([0.19,0.129,0.76,0.76*10./9])
        else:
            diMain.set_position([0.17,0.129,0.76,0.76*10./9])
        xLabel = 'Relative RA  '+unitStr
        yLabel = 'Relative Dec  '+unitStr
        if latex:
            if asConversion==1.0:
                unitStr = ' $[^{\prime\prime}]$'
            else:
                unitStr = ' [mas]'
            xLabel = r'$\Delta\alpha$ '+unitStr
            yLabel = r'$\Delta\delta$ '+unitStr
        diMain.set_xlabel(xLabel, fontsize=35)
        diMain.set_ylabel(yLabel, fontsize=35)
        if len(legendStrs)>0:
            diMain.legend()
        ##
        ## Save full size fig, then cropped to file and maybe convert to pdf if fl_format=='eps'
        ##
        #plt.tight_layout()
        orientStr = 'landscape'
        if fl_format=='eps':
            orientStr = 'portrait'
        # save full size            
        diMain.axes.set_xlim((xLimsFull[1],xLimsFull[0]))
        diMain.axes.set_ylim(yLimsFull)
        plotFilenameFull = plotFnameBase+'-DI.'+fl_format
        
        if plotCustomInsets:
            #######################################################
            # Locations of the subplots.  Notice that I explicitly multiply by 9/10
            # rather than approximating 0.8444 as 0.84.
            #######################################################
            xloc = [0.25, 0.28, 0.22]
            yloc = [0.48, 0.69, 0.819]
            dx = [7, 7, 7]
            dy = [15, 7, 7]
            dxfull = np.abs(xLimsFull[1] - xLimsFull[0])/0.76
            dyfull = np.abs(yLimsFull[1] - yLimsFull[0])/0.76*9/10
            
            for i in range(len(xloc)):
                if i > 0:
                    xx = -1*predicted_rapa_model[i + 1]*asConversion
                    yy = predicted_decsa_model[i + 1]*asConversion
                else:
                    xx = np.mean(-1*predicted_rapa_model[:2]*asConversion)
                    yy = np.mean(predicted_decsa_model[:2]*asConversion)
                expand = 2.25
                # Actual x and y positions we're blowing up as insets.
                realxloc = 0.171 + (xx - dx[i]*expand/2. + np.amax(xLimsFull))/dxfull
                realyloc = 0.13 + (yy - dy[i]*expand/2. - np.amin(yLimsFull))/dyfull
                a = diFig.add_axes([realxloc, realyloc, expand*dx[i]/dxfull, expand*dy[i]/dyfull])
                # Draw clear boxes around these areas, expand:1 scale.
                a.patch.set_visible(False)
                plt.xticks([])
                plt.yticks([])
    
                # Now draw the opaque insets at 10:1 scale.  Draw the orbit.
                a = diFig.add_axes([xloc[i], yloc[i], 10*dx[i]/dxfull, 10*dy[i]/dyfull],alpha=0.5)
                a.patch.set_visible(True)
                plt.plot(-1*fit_rapa_model[:]*asConversion,fit_decsa_model[:]*asConversion,linewidth=diLnThk,color='Blue') 
    
                # Draw thicker dots in the insets
                for j in range(len(predicted_decsa_model)):
                    xval = -1*predicted_rapa_model[j]*asConversion
                    yval = predicted_decsa_model[j]*asConversion
                    dotsize = diLnThk
                    plt.plot(xval, yval,c='red',marker='.',markersize=10*dotsize)#$$$$$$$$ Place for custimization
    
                # Now draw the actual astrometric measurements.
                data_DI = realDataDI.copy()
                data_DI[:, 3] *= -1
                data_DI[:, 1::2] *= -1
                (a,[xmin,xmax,ymin,ymax]) =  addDIdataToPlot(a,data_DI,asConversion,errMult=1,thkns=diLnThk,pasa=pasa)
    
                # Set the limits and don't label the axes.
                plt.xlim(xx - 0.5*dx[i], xx + 0.5*dx[i])
                plt.ylim(yy - 0.5*dy[i], yy + 0.5*dy[i])
                plt.xticks([])
                plt.yticks([])

        if plotFilenameFull!='':
            plt.savefig(plotFilenameFull, dpi=800, orientation=orientStr)
            log.info("DI orbit plot (Full) saved to:\n"+plotFilenameFull)
        ##crop to limits of data and save
        diMain.axes.set_xlim((xLimsCrop[1],xLimsCrop[0]))
        diMain.axes.set_ylim(yLimsCrop)
        plotFilenameCrop = plotFnameBase+'-DI-cropped.'+fl_format
        if plotFilenameCrop!='':
            plt.savefig(plotFilenameCrop, dpi=800, orientation=orientStr)
            log.info("DI orbit plot (cropped) saved to:\n"+plotFilenameCrop)
        plt.close()
        if (fl_format=='eps')and True:
            log.debug('converting to PDF as well')
            try:
                    os.system("epstopdf "+plotFilenameFull)
                    os.system("epstopdf "+plotFilenameCrop)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")    
        ## log params used in DI plot
        log.info('\n'+"*"*50+"\nOrbital Elements used in DI plot:\n"+repr(paramsDI))
        log.info("\n with an omega value = "+str(paramsDI[9]+settingsDI["omegaFdi"])+'\n'+"*"*50+'\n')

    ################
    # Make RV plot #
    ################
    settingsRV = settings[0]
    ##get the real data
    realData = loadRealData(diFilename=settingsRV['di_dataFile'],rvFilename=settingsRV['rv_dataFile'],dataMode=settingsRV['data_mode'])
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in orbParams[0][0]:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    if settingsRV['data_mode']!='DI':     
        realDataRV = copy.deepcopy(realData)
        realDataRV = realDataRV[np.where(realDataRV[:,6]<1e6)[0],:]
        
        ## instantiate Model class for this set of settings
        Model = ExoSOFTmodel(settingsRV)
        paramsRV = copy.deepcopy(params)
        paramsRVraw = copy.deepcopy(Model.Params.stored_to_direct(paramsRV))
        real_epochs_rv = copy.deepcopy(Model.Data.epochs_rv)
        real_rv = copy.deepcopy(Model.Data.rv)
        #real_rv_err = copy.deepcopy(Model.Data.rv_err)
        # ExoSOFTmodel will calculate predicted/ExoSOFTmodel data for each epoch in Model.Data.epochs_di
        ## calculate the fit locations for the DI epochs to calculate 0-C
        _ = ln_posterior(paramsRVraw, Model)
        predicted_rv_model = copy.deepcopy(Model.Data.rv_model)
        
        ##Make ExoSOFTmodel data for 100~1000 points for plotting fit
        nPts = 500
        fakeEpochs = np.zeros((nPts),dtype=np.dtype('d'))
        Model.Data.rv = np.ones((nPts),dtype=np.dtype('d'))
        Model.Data.rv_err = np.ones((nPts),dtype=np.dtype('d'))
        Model.Data.rv_model = np.ones((nPts),dtype=np.dtype('d'))
        Model.Data.rv_inst_num = np.zeros((nPts),dtype=np.dtype('i'))
        #fakeEpochs[:] = paramsRV[6]-(days_per_year*paramsRV[7]/2.0)
        last_epoch = paramsRV[6]-(days_per_year*paramsRV[7]/2.0)
        for i in range(0,nPts):
            last_epoch += (days_per_year*paramsRV[7])/float(nPts)
            fakeEpochs[i] = last_epoch
        fakeEpochs[-1] = fakeEpochs[-2]#paramsRV[6]+(days_per_year*paramsRV[7]/2.0)
        Model.Data.epochs_rv = fakeEpochs
        _ = ln_posterior(paramsRVraw, Model)        
        fit_epochs = copy.deepcopy(Model.Data.epochs_rv)
        fit_rv_model = copy.deepcopy(Model.Data.rv_model)
        
        ##Need to subtract RV offsets from the RVs 
        ##The fakeRealData had all offsets set to zero, so realDataRV needs to be "zeroed" to match
        numOffsets = int(len(paramsRV)-13)
        if (numOffsets-1)!=np.max(realDataRV[:,7]):
            log.critical("Number of RV offsets in params does not match largest value in realData[:,7]")
            log.critical("# of offsets in params = "+str(numOffsets)+" != # max in realData = "+str(np.max(realDataRV[:,7])+1))
        else:
            log.debug("There was a matching number of RV offsets in realData and params, = "+str(numOffsets))
            #for i in range(0,realDataRV.shape[0]):
            #    if realDataRV[i,6]<1000000.0:
            #        print('before offset subtract realData = '+str(realDataRV[i,0])+', '+str(realDataRV[i,5])+", "+str(realDataRV[i,6])+", "+str(realDataRV[i,7]))
            zeroedRealDataRV = copy.deepcopy(realDataRV)
            for i in range(0,zeroedRealDataRV.shape[0]):
                #rvBefore = zeroedRealDataRV[i,5]
                zeroedRealDataRV[i,5]-=paramsRV[13+int(zeroedRealDataRV[i,7])]
                #print str(rvBefore)+' - '+str(paramsRV[13+int(zeroedRealDataRV[i,7])])+" = "+str(zeroedRealDataRV[i,5])
            
            ##convert epochs to phases for plotting
            phasesReal = epochsToPhases(copy.deepcopy(realDataRV[:,0]),paramsRV[6],paramsRV[7], halfOrbit=True)
            phasesFit = epochsToPhases(copy.deepcopy(fit_epochs[:]),paramsRV[6],paramsRV[7], halfOrbit=True)            
            
            ## determine if to plot [km/s] or [m/s]
            kmConversion = 1.0/1000.0
            unitStr = '[km/s]'
            if np.max(np.sqrt(zeroedRealDataRV[:,5]**2.0))<1500:
                kmConversion = 1.0
                unitStr = '[m/s]'
            ## start making figure for residual and fit plots
            figRV = plt.figure(3,figsize=(10,5))
            residualsPlot = figRV.add_subplot(212)
            fitPlot = figRV.add_subplot(211)
            xLabel = "Orbital Phase"
            fitYlabel = 'v '+unitStr
            residYlabel = 'O-C '
            if latex:
                residYlabel = r'${\rm O-C}$ '
                fitYlabel = r'$v$  ${\rm '+unitStr+'}$'
                xLabel = r"${\rm Orbital}$  ${\rm Phase}$"
            residualsPlot.axes.set_xlabel(xLabel,fontsize=30)
            residualsPlot.axes.set_ylabel(residYlabel,fontsize=20)
            fitPlot.xaxis.set_ticklabels([])#this is just a hack way of killing the tick labels
            fitPlot.axes.set_ylabel(fitYlabel,fontsize=30)
            
            ## real-ExoSOFTmodel=residual, then plot it
            residualData = copy.deepcopy(realDataRV)
            #residualData[:,5]-= predicted_rv_model[:]
            for i in range(0,len(predicted_rv_model)):
                epoch_i = residualData[i,0]
                ## Silly, but the new RVs are in a different order than those in 'realDataRv'.
                ## So, need to search for proper index in new array.
                for j in range(0,len(predicted_rv_model)):
                    if epoch_i == real_epochs_rv[j]:
                        break
                residualData[i,5]-= predicted_rv_model[j]
            
            ##plot fit epochsORphases,RVs,RVerrs
            fitPlot.plot(phasesFit,fit_rv_model[:]*kmConversion,c='Blue',linewidth=diLnThk*0.8,alpha=1.0)
            ##plot zero vel line
            residualsPlot.axhline(linewidth=diLnThk*0.8,c='Blue') #adds a x-axis origin line
            ## add real data to plots
            residualsPlot = addRVdataToPlot(residualsPlot,phasesReal,residualData[:,5]*kmConversion,residualData[:,6]*kmConversion,datasetInts=residualData[:,7],alf=0.1,markersize=15,plotErrorBars=True)
            fitPlot = addRVdataToPlot(fitPlot,phasesReal,zeroedRealDataRV[:,5]*kmConversion,zeroedRealDataRV[:,6]*kmConversion,datasetInts=residualData[:,7],alf=0.2,markersize=9,plotErrorBars=True)
            
            ## Find and set limits 
            xLims = (np.min([np.min(phasesFit),np.min(phasesReal)]),np.max([np.max(phasesFit),np.max(phasesReal)]))
            xLims = (xLims[0]-(xLims[1]-xLims[0])*.05,xLims[1]+(xLims[1]-xLims[0])*.05)
            fitYlims = (np.min([np.min(fit_rv_model[:]*kmConversion),np.min(zeroedRealDataRV[:,5]*kmConversion)]),np.max([np.max(fit_rv_model[:]*kmConversion),np.max(zeroedRealDataRV[:,5]*kmConversion)]))
            fitYrange = fitYlims[1]-fitYlims[0]
            fitYlims = (fitYlims[0]-fitYrange*.05,fitYlims[1]+fitYrange*.05)
            residYlims = (np.min(residualData[:,5]*kmConversion),np.max(residualData[:,5]*kmConversion))
            residYrange = residYlims[1]-residYlims[0]
            residYlims = (residYlims[0]-residYrange*.05,residYlims[1]+residYrange*.05)
            ## update lims with custom values if provided
            if len(RVlims)>0:
                #RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
                fitYlims = (RVlims[0][0],RVlims[0][1])
                residYlims = (RVlims[1][0],RVlims[1][1])
                xLims = (RVlims[2][0],RVlims[2][1])
            residualsPlot.axes.set_xlim(xLims)
            residualsPlot.axes.set_ylim(residYlims)
            #RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
            fitPlot.axes.set_xlim(xLims)
            fitPlot.axes.set_ylim(fitYlims)
            ## adjust plot locations to account for longer numbers
            #[left,btm,width,height]
            if (1000<abs(fitYlims[0])) or (1000<abs(fitYlims[1])):
                residualsPlot.set_position([0.17,0.17,0.81,0.23])
                fitPlot.set_position([0.17,0.39,0.81,0.57])
            else:
                residualsPlot.set_position([0.13,0.17,0.84,0.23])
                fitPlot.set_position([0.13,0.39,0.84,0.57])
            
            
            ##load resulting data to file for re-plotting by others
            #real [phases,JD,offset subtracted RV, residual, dataset#]
            outRVdataReal = []
            for i in range(0,len(phasesReal)):
                outRVdataReal.append([phasesReal[i],realDataRV[i,0],zeroedRealDataRV[i,5],residualData[i,5],realDataRV[i,7]])
            #print repr(outRVdataReal)
            #fit [phases,JD,RV]
            outRVdataFit = []
            for i in range(0,len(phasesFit)):
                outRVdataFit.append([phasesFit[i],fit_epochs[i],fit_rv_model[i]])
            fnameBase = os.path.join(os.path.dirname(plotDataDir),'RVplotData')
            np.savetxt(fnameBase+'-real.dat',outRVdataReal,header="[phases,JD,offset subtracted RV, residual, dataset#]")
            np.savetxt(fnameBase+'-fit.dat',outRVdataFit,header="[phases,JD,RV]")
            for i in range(0,int(np.max(residualData[:,7])+1)):
                ary = []
                chiSqr = 0
                print('')
                for j in range(0,len(phasesReal)):
                    if residualData[j,7]==float(i):
                        if residualData[j,5] not in ary:
                            ary.append(residualData[j,5])
                            chiSqr+=(residualData[j,5]**2.0)/(residualData[j,6]**2.0)
                            #if np.abs(residualData[j,5])>residualData[j,6]:
                                #print('dataset '+str(i)+', epoch '+str(realDataRV[j,0])+' had residual  '+str(np.abs(residualData[j,5])-residualData[j,6])+' greater than error')
                d = np.array(ary)
                #data = residualData[np.where(residualData[:,7]==float(i))[0],:]
                d = np.abs(d)
                #print('For RV dataset #'+str(i))
                #print('Max abs residual = '+str(np.max(d)))
                #print('Min abs residual = '+str(np.min(d)))
                #print('mean abs residual = '+str(np.mean(d)))
                #print('chiSqr for this set was = '+str(chiSqr))
                #print('chiSqr/numRVs for this set was = '+str(chiSqr/len(d)))
            
            ##clean up boarders, axis ticks and such 
            plt.minorticks_on()
            fitPlot.tick_params(axis='both',which='major',width=1,length=5,pad=8,direction='in',labelsize=25)
            fitPlot.tick_params(axis='both',which='minor',width=1,length=2,pad=8,direction='in')
            fitPlot.spines['right'].set_linewidth(1.0)
            fitPlot.spines['bottom'].set_linewidth(1.0)
            fitPlot.spines['top'].set_linewidth(1.0)
            fitPlot.spines['left'].set_linewidth(1.0)
            fitPlot.locator_params(axis='y',nbins=6) #fix number of y-axis label points
            plt.minorticks_on()
            #plot.axhline(linewidth=2.0) #adds a x-axis origin line
            #plot.axvline(linewidth=2.0) #adds a y-axis origin line
            plt.minorticks_on()
            residualsPlot.tick_params(axis='y',which='major',width=1,length=5,pad=8,direction='in',labelsize=15)
            residualsPlot.tick_params(axis='x',which='major',width=1,length=5,pad=8,direction='in',labelsize=25)
            residualsPlot.tick_params(axis='y',which='minor',width=2,length=4,pad=8,direction='in')
            residualsPlot.tick_params(axis='x',which='minor',width=2,bottom='on',length=4,pad=8,direction='in')
            residualsPlot.locator_params(axis='y',nbins=4) #fix number of y-axis label points
            plt.minorticks_on()
            residualsPlot.spines['right'].set_linewidth(1.0)
            residualsPlot.spines['bottom'].set_linewidth(1.0)
            residualsPlot.spines['top'].set_linewidth(1.0)
            residualsPlot.spines['left'].set_linewidth(1.0)
    
            ## save fig to file and maybe convert to pdf if fl_format=='eps'
            #plt.tight_layout()
            orientStr = 'landscape'
            if fl_format=='eps':
                orientStr = 'portrait'
            plotFilename = plotFnameBase+'-RV.'+fl_format
            if plotFilename!='':
                plt.savefig(plotFilename, dpi=300, orientation=orientStr)
                log.info("RV orbit plot saved to:\n"+plotFilename)
            plt.close()
            if (fl_format=='eps')and True:
                log.debug('converting to PDF as well')
                try:
                    os.system("epstopdf "+plotFilename)
                except:
                    log.warning("Seems epstopdf failed.  Check if it is installed properly.")
                ## log params used in RV plot
            log.info('\n'+"*"*50+"\nOrbital Elements used in RV plot:\n"+repr(paramsRV))
            log.info("\n with an omega value = "+str(paramsRV[9]+settingsRV["omegaFrv"])+'\n'+"*"*50+'\n')

def nodeEpochsCalc(paramsDI,omegaDIoffset):
    """
    Calculate the epochs for the Ascending and Descending nodes, might be in different orbital 
    periods and AN/DN might be the wrong order, but should work for plotting... I hope...
    """
    taAtNodes = [-1.0*paramsDI[9]+omegaDIoffset,180.0-1.0*paramsDI[9]+omegaDIoffset]
    nodeEpochs = []
    for ta in taAtNodes:
        if ta<0.0:
            ta =ta+360.0
        elif ta>360:
            ta =ta-360.0
        TA_s_rad = ta*(np.pi/180.0)
        top = np.sqrt(1.0-paramsDI[4])*np.sin(TA_s_rad/2.0)   
        btm = np.sqrt(1.0+paramsDI[4])*np.cos(TA_s_rad/2.0) 
        ATAN_rad = np.arctan2(top, btm)
        #NOTE: both math.atan2 and np.arctan2 tried with same results, both produce negatives rather than continuous 0-360 
        #thus, must correct for negative outputs
        if ATAN_rad<0:
            ATAN_rad = ATAN_rad+(2.0*np.pi)
        M_s_rad = ATAN_rad*2.0-paramsDI[4]*np.sin(ATAN_rad*2.0)
        delta_t = (M_s_rad*paramsDI[7]*days_per_year)/(2.0*np.pi)
        nodeEpochs.append(paramsDI[5]+delta_t)
    return nodeEpochs
 
def densConfInt(x, pdf, confidence_level):
    """copied directly from https://gist.github.com/adrn/3993992"""
    a = pdf[pdf > x].sum() - confidence_level
    #print str(confidence_level)+", "+str(pdf[pdf > x].sum())+" -> so far "+str(a)
    return a

def densityContourFunc(xdata, ydata, nbins, ax=None,ranges=None,bests=None):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        
    Copied directly, and heavily modified by Kyle after, 
    from https://gist.github.com/adrn/3993992    
    """
    import matplotlib.cm as cm
    from matplotlib.colors import from_levels_and_colors
    ## get data for 2D hist and a pdf from it
    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=nbins,range=ranges, normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins,1))
    H = ndimage.gaussian_filter(H, sigma=1)
    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    ## find contour levels and make color maps
    tiny_sigma= so.brentq(densConfInt, 0., 1., args=(pdf, 0.0001))
    #ptone_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.080))
    #oneQuarter_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.197))
    #oneHalf_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.383))
    #one_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.68))
    one_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, (1-np.exp(-0.5*(1.0**2)))))
    #two_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.95))
    two_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, (1-np.exp(-0.5*(2.0**2)))))
    #three_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.997))
    three_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, (1-np.exp(-0.5*(3.0**2)))))
    #four_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.99994))
    four_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, (1-np.exp(-0.5*(4.0**2)))))
    #five_sigma = so.brentq(densConfInt, 0., 1., args=(pdf, 0.999999426))
    levels3sigs = [three_sigma,two_sigma,one_sigma]
    #levels7lvls = [four_sigma,three_sigma, two_sigma, one_sigma, oneHalf_sigma, oneQuarter_sigma, ptone_sigma]
    c = []
    c2 = []
    for _ in range(len(levels3sigs)+1):
        c.append('k')
        c2.append('g')
    (black_cmap,n) = from_levels_and_colors(levels3sigs,c,extend='both')
    (solidColor_cmap,n) = from_levels_and_colors(levels3sigs,c2,extend='both')
    X, Y, Z = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1]), pdf.T
    if ax == None:
        # Lots of color map choices at (http://matplotlib.org/users/colormaps.html) and note that most have a reversed version by adding '_r' to the end.
        # this is because we flipped the levels, so the colors needed to also be flipped.
        #possible choices suitable here are: YlGnBu, Blues, PuBuGn, PuBu, afmhot_r, copper_r
        #main contour
        contour = plt.contourf(X, Y, np.log(Z+1e-10),levels=np.linspace(np.log(four_sigma),np.log(tiny_sigma),50),origin="lower",cmap=cm.afmhot_r)
        #add lines for 1-2-3sigmas
        contour = plt.contour(X, Y, Z, levels=levels3sigs, origin="lower", cmap=solidColor_cmap,linewidths=3,linestyles='dashed')
        #Plot lines for best values
        if bests!=None:
            contour = plt.plot([X.min(),X.max()], [bests[0],bests[0]],linewidth=3,color='blue')
            contour = plt.plot([bests[1],bests[1]],[Y.min(),Y.max()],linewidth=3,color='blue')
    else:
        contour = ax.contourf(X, Y, np.log(Z+1e-10),levels=np.linspace(np.log(four_sigma),np.log(tiny_sigma),50),origin="lower",cmap=cm.afmhot_r)
        #add lines for 1-2-3sigmas
        contour = ax.contour(X, Y, Z, levels=levels3sigs, origin="lower", cmap=solidColor_cmap,linewidths=3,linestyles='dashed')
        #Plot lines for best values
        if bests!=None:
            contour = ax.plot([X.min(),X.max()], [bests[1],bests[1]],linewidth=3,color='blue')
            contour = ax.plot([bests[0],bests[0]],[Y.min(),Y.max()],linewidth=3,color='blue')
        if True:
            #print('np.median(xdata) = '+repr(np.median(xdata)))
            #print('np.median(ydata) = '+repr(np.median(ydata)))
            contour = ax.plot([X.min(),X.max()], [np.median(ydata),np.median(ydata)],linewidth=3,color='red')
            contour = ax.plot([np.median(xdata),np.median(xdata)],[Y.min(),Y.max()],linewidth=3,color='red')
    return contour 

def densityPlotter2D(outputDataFilename,plotFilename,paramsToPlot=[],bestVals=None,ranges=None,rectanglePlot=True,smooth=True):
    """
    Will create a 2D density contour plot.
    Must pass in ONLY 2 params to plot.
    """    
    scatterTest = False
    latex=True
    plotFormat='eps'
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
    if len(paramsToPlot)==2:
        (head,data) = loadFits(outputDataFilename)
        if head!=False:  
            log.debug(' Inside densityPlotter')
            s= '\nCreating 2D density plot for file:\n'+outputDataFilename
            s=s+ '\nInput plotfilename:\n'+plotFilename
            log.info(s)
            ## check if the passed in value for plotFilename includes format extension
            if '.'+plotFormat not in plotFilename:
                plotFilename = plotFilename+"."+plotFormat
                log.debug('updating plotFilename to:\n'+plotFilename)
            else:
                plotFilename = plotFilename
            ## Get strings representing axes titles and plot filenames in latex and standard formats
            (_,paramStrs,_) = getParStrs(head,latex=latex,getALLpars=True)
            ## modify x labels to account for DI only situations where M1=Mtotal
            if np.var(data[:,1])==0:
                paramStrs[0] = r'$m_{total}$ [$M_{\odot}$]'
            ## check if a subset is to be plotted or the whole set
            ## remake lists of params to match subset.
            if len(paramsToPlot)!=0:
                paramStrsUse = []
                for par in paramsToPlot:
                    paramStrsUse.append(paramStrs[par])
                paramStrs = paramStrsUse
                
            ## convert m2 to Mjup if necessary
            in_m_jup = False
            if 1 in paramsToPlot:
                if np.max(data[:,1])<0.02:
                    data[:,1] = data[:,1]*(const.M_sun.value/const.M_jup.value)
                    bestVals[1] = bestVals[1]*(const.M_sun.value/const.M_jup.value)
                    in_m_jup = True
                    
            xdata = data[:,paramsToPlot[0]]
            ydata = data[:,paramsToPlot[1]]           
            
            nbins=50
            ## update lims with custom values if provided
            if ranges!=None:
                if len(ranges)>0:
                    #ranges=[[xMin,xMax],[yMin,yMax]]
                    xLims=ranges[0]
                    yLims=ranges[1]
            else:
                xLims = [np.min(xdata),np.max(xdata)]
                yLims = [np.min(ydata),np.max(ydata)]
            rangesOrig=[xLims,yLims]
            #force ranges to be square
            xR = xLims[1]-xLims[0]
            yR = yLims[1]-yLims[0]
            if xR>yR:
                yLims = [yLims[0]-0.5*abs(xR-yR),yLims[1]+0.5*abs(xR-yR)]
            elif xR<yR:
                xLims = [xLims[0]-0.5*abs(xR-yR),xLims[1]+0.5*abs(xR-yR)]
            rangesSqr=[xLims,yLims]
            sqBools=[True,False]
            for sqBool in sqBools:
                x = 10.0
                y = 10.0
                if rectanglePlot and (sqBool==False): 
                    xR = rangesOrig[0][1]-rangesOrig[0][0]
                    yR = rangesOrig[1][1]-rangesOrig[1][0]
                    if xR>yR:
                        y=(yR/xR)*y
                    elif xR<yR:
                        x=(xR/yR)*x
                    ########## Hack to get bottom axis label to show properly #$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    ########## when y range is small
                    if True:
                        y=1.05*y
                        hgt=0.8/1.05
                        btm=0.17+(0.8-hgt)
                #print('Figure size is ('+str(x)+', '+str(y)+")")
                fig = plt.figure(figsize=(x,y))
                subPlot = fig.add_subplot(111)
                xLabel = paramStrs[0]
                yLabel = paramStrs[1]
                ## Tweak plot axis and labels to look nice
                subPlot.locator_params(axis='x',nbins=7) # maximum number of x labels
                subPlot.locator_params(axis='y',nbins=7) # maximum number of y labels
                subPlot.tick_params(axis='x',which='major',width=0.5,length=3,pad=4,direction='in',labelsize=25)
                subPlot.tick_params(axis='y',which='major',width=0.5,length=3,pad=4,direction='in',labelsize=25)
                subPlot.spines['right'].set_linewidth(0.7)
                subPlot.spines['bottom'].set_linewidth(0.7)
                subPlot.spines['top'].set_linewidth(0.7)
                subPlot.spines['left'].set_linewidth(0.7)                        
                if rectanglePlot and (sqBool==False): 
                    #[left,btm,width,height]
                    subPlot.set_position([0.17,btm,0.80,hgt])
                else:
                    #[left,btm,width,height]
                    subPlot.set_position([0.17,0.14,0.80,0.80])
                # add axes labels
                fsize=34
                fsizeY = fsize
                fsizeX = fsize
                ## update font size for ecc label
                if xLabel in ['e','$e$']:
                    fsizeX =fsize+10
                elif yLabel in ['e','$e$']:
                    fsizeY =fsize+10
                ## update label for m2 if it is in mjup
                if in_m_jup:
                    if xLabel in ['m2','$m_2{\rm [M}_{\odot}{\rm ]}$']:
                        xLabel='m2 [Mjupiter]'
                        if latex:
                            xLabel=r'$m_2$ [$M_{J}$]'
                    elif yLabel in ['m2','$m_2{\rm [M}_{\odot}{\rm ]}$']:
                        yLabel='m2 [Mjupiter]'
                        if latex:
                            yLabel=r'$m_2$ [$M_{J}$]'
                if latex:
                    subPlot.axes.set_xlabel(xLabel,fontsize=fsizeX)
                    subPlot.axes.set_ylabel(yLabel,fontsize=fsizeY)
                else:
                    subPlot.axes.set_xlabel(xLabel,fontsize=fsizeX)
                    subPlot.axes.set_ylabel(yLabel,fontsize=fsizeY)
                ranges = rangesOrig
                if sqBool:
                    ranges = rangesSqr
                ## call densityContour func to fill up subplot with density/contour plot
                if scatterTest==False:
                    subPlot = densityContourFunc(xdata, ydata, nbins, ax=subPlot,ranges=ranges,bests=bestVals)
                ##$$$$$$$$$$$$$$$$$$$$$$$ sanity check with a 1/1000 scatter plot $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                if scatterTest:
                    subPlot.scatter(xdata[::1000],ydata[::1000],c='green',edgecolors='k')
                    #print('np.median(xdata) = '+repr(np.median(xdata)))
                    #print('np.median(ydata) = '+repr(np.median(ydata)))
                    subPlot.plot([ranges[0][0],ranges[0][1]], [np.median(ydata),np.median(ydata)],linewidth=3,color='red')
                    subPlot.plot([np.median(xdata),np.median(xdata)],[ranges[1][0],ranges[1][1]],linewidth=3,color='red')
                    subPlot.plot([ranges[0][0],ranges[0][1]], [bestVals[1],bestVals[1]],linewidth=3,color='blue')
                    subPlot.plot([bestVals[0],bestVals[0]],[ranges[1][0],ranges[1][1]],linewidth=3,color='blue')
                    subPlot.axes.set_xlim((ranges[0][0],ranges[0][1]))
                    subPlot.axes.set_ylim((ranges[1][0],ranges[1][1]))
                ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                ## Save file if requested.
                log.debug('\nStarting to save density contour figure:')
                if plotFilename!='':
                    plotnm = plotFilename
                    if rectanglePlot and (sqBool==False): 
                        plotnm = plotFilename[:-4]+'-rectangular.eps'
                    elif sqBool:
                        plotnm = plotFilename[:-4]+'-squareRanges.eps'
                    plt.savefig(plotnm,format=plotFormat)
                    s= 'density contour plot saved to: '+plotnm
                    log.info(s)
                plt.close()
                if (plotFormat=='eps') and True:
                    log.debug('converting to PDF as well')
                    try:
                        os.system("epstopdf "+plotnm)
                    except:
                        log.warning("Seems epstopdf failed.  Check if it is installed properly.")
    else:
        log.critical(repr(len(paramsToPlot))+" params requested to be plotted, yet only 2 is acceptable.")
    
def cornerPlotter(outputDataFilename,plotFilename,paramsToPlot=[],bestVals=[],smooth=True):
    """
    make a triangle/corner plot by using the corner package written by dfm:
    https://github.com/dfm/corner.py
    NOTE: the contours of the density plots in here are [ 0.1175031 ,  0.39346934,  0.67534753,  0.86466472]
    """
    from corner import corner
    
    latex=True
    plotFormat='eps'
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
        
    (head,data) = loadFits(outputDataFilename)
    
    if head!=False:  
        log.debug(' Inside tranglePlotter')
        s= '\nCreating corner plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        log.info(s)
        
        ## check if the passed in value for plotFilename includes fl_format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            log.debug('updating plotFilename to:\n'+plotFilename)
        else:
            plotFilename = plotFilename
                
        ## get parameter lists and filter accordingly
        (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=True)
        (_,paramStrs2,_) = getParStrs(head,latex=False,getALLpars=True)
        # modify x labels to account for DI only situations where M1=Mtotal
        if (np.var(data[:,1])==0)and (0 in paramList):
            paramStrs2[0] = 'm total [Msun]'
            paramStrs[0] = r'$m_{\rm total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'm-total'
        # check if a subset is to be plotted or the whole set
        # remake lists of params to match subset.
        dataUse = data
        if len(paramsToPlot)!=0:
            dataUse = data[:,paramsToPlot]
            paramStrs2Use = []
            paramStrsUse = []
            paramFileStrsUse = []
            paramListUse = []
            bestValsUse = []
            for par in paramsToPlot:
                paramStrs2Use.append(paramStrs2[par])
                paramStrsUse.append(paramStrs[par])
                paramFileStrsUse.append(paramFileStrs[par])
                paramListUse.append(par)
                bestValsUse.append(bestVals[par])
            paramStrs2 = paramStrs2Use
            paramStrs = paramStrsUse
            paramFileStrs = paramFileStrsUse 
            paramList = paramListUse
            bestVals = bestValsUse
            
        log.info("will try to make a triangle/corner plot for data of shape: "+repr(dataUse.shape))
        
        bests=None
        if len(bestVals)==len(paramsToPlot)!=0:
            bests = bestVals
        if latex:
            paramStrs = paramStrs2
        ##########################################################    
        ##call triangle plot function corner to make the figure
        ##########################################################   
        log.debug('About to call corner func')
        tic=timeit.default_timer()
        ## Create empty figure to be filled up with plots
        _ = corner(dataUse, bins=50, range=None, color="k",
                           smooth=smooth,labels=paramStrs,
                           truths=bests, truth_color="#4682b4",
                           verbose=False, fig=None,
                           max_n_ticks=5, top_ticks=False)
        log.debug("back from corner func")
        toc = timeit.default_timer()
        log.info("corner plotting took a total of "+timeStrMaker(toc-tic))
    
        #plt.tight_layout()
        ## Save file if requested.
        log.debug('\nStarting to save corner figure:')
        if plotFilename!='':
            plt.savefig(plotFilename,fl_format=plotFormat)
            log.info('Corner plot saved to: '+plotFilename)
        plt.close()
        toc2 = timeit.default_timer()
        log.info("Saving took a total of "+timeStrMaker(toc2-toc))
        if (plotFormat=='eps') and True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
    
def progressPlotter(outputDataFilename,plotFilename,paramToPlot,yLims=[],xLims = [],expectedVal=None,emcee_stage=False):
    """
    Plots progress of a single parameter's chain over one stage of simulation, AND the 
    reduced chi squared as a time series.
    
    Make sure to set emcee_stage=True if plotting results from emcee.
    """
    downSample = True # down sample the data to clean up plot
    cutDownBy = 1 # down sample by using ever #th data point
    scatter = True
    bigFig = True
    latex=True
    plotFormat='eps'
    matplotlib.rcParams['ps.useafm']= True
    matplotlib.rcParams['pdf.use14corefonts'] = True
    matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        matplotlib.rc('text', usetex=True)
        matplotlib.rcParams['text.latex.unicode']=True 
        matplotlib.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
    else:
        matplotlib.rc('font',family='serif')
        matplotlib.rc('text', usetex=False)
        
    lblSz = 20
    fntSz = 25
    lnWdth = 2
    if bigFig:
        lblSz = 50
        fntSz = 60
        lnWdth = 6
    (head,data) = loadFits(outputDataFilename)
    
    if head!=False:  
        log.debug(' Inside progressPlotter')
        s= '\nCreating progress plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        log.info(s)
        
        #Find best orbit params
        bestPars = findBestOrbit(outputDataFilename,bestToFile=False,by_ln_prob=emcee_stage)
        #print('back from findBestOrbit')
        
        ## check if the passed in value for plotFilename includes fl_format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            log.debug('updating plotFilename to:\n'+plotFilename)
        else:
            plotFilename = plotFilename
                
        (_,paramStrs,paramFileStrs) = getParStrs(head,latex=latex,getALLpars=True)
        (_,paramStrs2,_) = getParStrs(head,latex=False,getALLpars=True)
        nu =  head['NU']
        
        ## modify y labels to account for DI only situations where M1=Mtotal
        if np.var(data[:,1])==0:
            paramStrs2[0] = 'm total [Msun]'
            paramStrs[0] = '$m_{total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'm-total'
        ## check if m2 and if it should be in jupiter masses
        elif paramToPlot==1:
            if np.max(data[:,1])<0.02:
                data[:,1]=data[:,paramToPlot]*(const.M_sun.value/const.M_jup.value)
                bestPars[1] = bestPars[1]*(const.M_sun.value/const.M_jup.value)
                #valRange = np.max(data[:,1])-np.min(data[:,1])
                paramStrs2[1]='m2 [Mjupiter]'
                paramStrs[1]='$m_2$ [$M_{J}$]'

        ##make progress plots for parameter requested and reduced chi squared
        saveInt = int(head['SAVEINT'])
        #samples = range(0,len(data[:,11])-1)*saveInt
                
        ## Create empty figure to be filled up with plots
        if bigFig:
            fig = plt.figure(figsize=(16,10),dpi=200)
        else:
            fig = plt.figure(figsize=(8,5),dpi=100)  
        #print('made fig')
        #plot requested param  
        subPlot = fig.add_subplot(2,1,1)
        #print('len(data[:,paramToPlot]) = '+repr(len(data[:,paramToPlot])))
        samples = np.arange(0,len(data[:,paramToPlot])*saveInt,saveInt)
        #print('len(samples) = '+repr(len(samples)))
        if downSample:
            d = data[::cutDownBy,paramToPlot]
            samples = samples[::cutDownBy]
        else:
            d = data[:,paramToPlot]
        if scatter:
            subPlot.scatter(samples,d,c='green',edgecolors='green')
        else:
            subPlot.plot(samples,d,color='blue',linewidth=1)
        #print('##plotted data##')
        #print repr([samples[0],samples[-1]])+', '+repr([bestPars[paramToPlot],bestPars[paramToPlot]])
        
        # plot line for best fit
        if len(bestPars)>1:
            subPlot.plot([samples[0],samples[-1]],[bestPars[paramToPlot],bestPars[paramToPlot]],color='k',linewidth=lnWdth)
        # plot line for expected value if provided
        if expectedVal is not None:
            subPlot.plot([samples[0],samples[-1]],[expectedVal,expectedVal],color='blue',linewidth=lnWdth)
        #print('##plotted best##')
        if latex:
            subPlot.axes.set_ylabel(r''+paramStrs[paramToPlot],fontsize=fntSz)
        else:
            subPlot.axes.set_ylabel(paramStrs2[paramToPlot],fontsize=fntSz)           
        #place limits on y axis if provided
        if len(yLims)>0:
            subPlot.axes.set_ylim((yLims[0][0],yLims[0][1]))
        #place limits on x axis if provided
        if len(xLims)>0:
            subPlot.axes.set_xlim((xLims[0],xLims[1]))
        # change label sizes
        subPlot.tick_params(axis='both',labelsize=lblSz)
        #plot chi squareds
        subPlot = fig.add_subplot(2,1,2)
        #print('len(data[:,11]*(1.0/nu)) = '+repr(len(data[:,11]*(1.0/nu))))
        #samples = np.arange(0,len(data[:,11])*saveInt,saveInt)
        #print('len(samples) = '+repr(len(samples)))
        #subPlot.plot(samples,data[:,11]*(1.0/nu),color='k',linewidth=1)
        if downSample:
            d = data[::cutDownBy,11]*(1.0/nu)
        else:
            d = data[:,11]*(1.0/nu)
        if scatter:
            subPlot.scatter(samples,d,c='green',edgecolors='green')
        else:
            subPlot.plot(samples,d,color='k',linewidth=1)
        #print('made subplot')
        if latex:
            subPlot.axes.set_xlabel(r''+'$Sample$',fontsize=fntSz)
            #subPlot.axes.set_ylabel(r''+'$\chi^2_{\nu}$',fontsize=20)
            subPlot.axes.set_ylabel(r''+'$reduced$ $\chi^2$',fontsize=fntSz)
        else:
            subPlot.axes.set_xlabel('Sample',fontsize=fntSz)
            subPlot.axes.set_ylabel('reduced chi sqr',fontsize=fntSz)
        #print('labeled plot')   
        #place limits on y axis if provided
        if len(yLims)>0:
            subPlot.axes.set_ylim((yLims[1][0],yLims[1][1]))
        #place limits on x axis if provided
        if len(xLims)>0:
            subPlot.axes.set_xlim((xLims[0],xLims[1]))
        # change label sizes
        subPlot.tick_params(axis='both',labelsize=lblSz)
        plt.tight_layout()
        ## Save file if requested.
        log.debug('\nStarting to save param progress figure:')
        if plotFilename!='':
            plt.savefig(plotFilename,fl_format=plotFormat)
            s= 'progress plot saved to: '+plotFilename
            log.info(s)
        plt.close()
        if (plotFormat=='eps') and True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
#END OF FILE

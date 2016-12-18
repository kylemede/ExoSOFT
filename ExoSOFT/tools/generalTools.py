#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp  or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
#import numpy as np
#import timeit
import copy
import os
import shutil
import datetime
import glob
import numpy as np
from scipy import interpolate
import sys
from astropy.io import fits as pyfits
from astropy import constants as const
import warnings
import jdcal
import emcee
import KMlogger
from six.moves import range

## import from modules in ExoSOFT ##
#from . import constants as const
from .cytools import mean_corr_len
from .model import ExoSOFTmodel, ln_posterior
from .readWriteTools import loadFits, writeFits, rmFiles, dataReader

days_per_year = 365.2422
sec_per_year = 60*60*24*days_per_year

warnings.simplefilter("error")
log = KMlogger.getLogger('main.genTools',lvl=100,addFH=False)  
    
def mcmcEffPtsCalc(outputDataFilename):
    """
    Calculate average correlation length and the number of effective steps for each parameter 
    that was varying during the simulation.  The results are put into the log.
    
    Using Correlation Length based on Tegmark2004.  This is the number of steps until the 
    variance (aka correlation) is half that of the total chain.  We perform this in a "box car"
    style that starts calculating it again on the step just after the previously found correlation 
    length step.  This way it produces an average value that is more reliable.
    """
    log.info("Starting to calculate correlation lengths")
    (head,data) = loadFits(outputDataFilename)
    numSteps = data.shape[0]
    (paramList,paramStrs,_) = getParStrs(head,latex=False)
    completeStr=""
    try:
        completeStr+= '\n'+'-'*47+'\nThe mean correlation lengths of all params are:\n'+'(based on Tegmark et. al. 2004)\n' +'-'*47
        completeStr+='\nTotal # of steps stored '+str(numSteps)+'\nparam #, param name, mean correlation length'
        completeStr+= ' -> total # of steps/mean correlation length = number of effective points\n'
        for i in range(0,len(paramList)):
            log.debug( "*"*60+"\n"+'starting to mean calculate corr length for '+paramStrs[i]+' with CPP')
            
            if False:
                print("\n"+paramStrs[i]+":")
                meanCorrLength = np.cov(data[:,paramList[i]])
                print(" np output = "+str(meanCorrLength))
            if True:
                dataC = np.ascontiguousarray(data[:,paramList[i]],dtype=np.dtype('d'))
                #print("Calling cy mean_corr_len")
                meanCorrLength = mean_corr_len(dataC)
                #print(" cy output = "+str(meanCorrLength))
            
            currParamStr = str(paramList[i])+', '+paramStrs[i]+", "+str(meanCorrLength)
            currParamStr+=    ' -> '+str(numSteps)+'/'+str(meanCorrLength)+' = '+str(numSteps/meanCorrLength)+'\n'
            completeStr+=currParamStr
            log.debug(currParamStr)
        log.debug(completeStr)
    except:
        log.error("An error occurred when trying to calculate the correlation lengths, # effective steps")
    return completeStr

def autocorr(outputDataFilename,fast=True):
    """
    This directly calls the integrated_time function from the emcee package.
    
    https://github.com/dfm/emcee/blob/master/emcee/autocorr.py
    ' This estimate uses the iterative procedure described on page 16 of `Sokal's
    notes <http://www.stat.unc.edu/faculty/cji/Sokal.pdf>`_ to determine a
    reasonable window size.'
    """
    log.info("Starting to calculate autocorrelation with emcee.autocorr.integrated_time")
    (head,data) = loadFits(outputDataFilename)
    #numSteps = data.shape[0]
    (paramList,paramStrs,_) = getParStrs(head,latex=False)
    completeStr=""
    completeStr+= '\n'+'-'*63+'\nThe longest integrated autocorrelation times of all params are:\n'+"(using emcee.autocorr.integrated_time)\n"+'-'*63+'\n'
    
    all_passed = True
    #print("\n\n in autocorr \n\n")
    for i in range(0,len(paramList)):
        #print('\ni = '+repr(i))
        x = data[:,paramList[i]]
        #print("\nx = "+repr(x))
        ac = "too long"
        try:
            ac = emcee.autocorr.integrated_time(x,fast=fast)
        except:
            s = "\nError thrown while calculating autocorr for parameter, "+paramStrs[i]
            s+= "\nThis is most likely due to the chain not having converged yet, try a larger number of samples, or ignore this warning."
            log.debug(s)
            all_passed = False
        #print("\n\n ac: "+repr(ac)+'\n\n')
        currParamStr = str(paramList[i])+', '+paramStrs[i]+" = "+str(ac)+'\n'
        completeStr+=currParamStr
        log.debug(currParamStr)
    if all_passed==False:
        s = "Was unable to calculate autocorr for at least one of the parameters.  Most likely as the chains hadn't converged."
        log.info(s)
        completeStr+=s+'\n'
    #except:
    #    log.error("An error occurred when trying to calculate the correlation lengths, # effective steps")
    #print('done autocorr')
    return completeStr

def burnInCalc(mcmcFnames,combinedFname):
    """
    NOTE: ExoSOFT was designed to start the full MCMC chain from the last point of the 
        Sigma Tuning stage.  As this stage effectively acts as a form of burn-in period
        the burn-in value found from the pure MCMC tends to be very short.
 
    Calculate the burn in for a set of MCMC chains following the formulation of Tegmark.
     
    Burn-in is defined as the first point in a chain where the likelihood is greater than 
    the median value of all the chains.  Thus, there MUST be more than 1 chain to perform this calculation.
    """
    log.info("Starting to calculate burn-in.")
     
    #chiSquaredsALL = np.array([])
    burnInLengths = []
    # calculate median of combined data ary
    (head,data) = loadFits(combinedFname)
    #nu = float(head0['NU'])
    chiSqs = data[:,11]
    if type(chiSqs)!=np.ndarray:
        chiSqs = np.array(chiSqs)
    likelihoods = np.exp(-chiSqs/2.0)
    log.debug("likelihoodsALL min = "+repr(np.min(likelihoods)))
    log.debug("likelihoodsALL max = "+repr(np.max(likelihoods)))
    medainALL = np.median(likelihoods)         
    log.debug("medainALL = "+str(medainALL))
    s =21*'-'+'\nBurn-In lengths were:\n'+21*'-'
    s+='\nmedian value for all chains = '+str(medainALL)
    ## calculate location of medianALL in each chain
    for filename in mcmcFnames:
        if os.path.exists(filename):
            (head,data) = loadFits(filename)
            chiSqs = data[:,11]
            likelihoods = np.exp(-chiSqs/2.0)
            #medianChain = np.median(chiSquaredsChain)
            burnInLength = len(likelihoods)
            i=0
            while i<(len(likelihoods)-1):
                i+=1
                if likelihoods[i]>medainALL:
                    #print 'chiSqs[i] = '+str(chiSqs[i])
                    burnInLength = i+1
                    break
            burnInLengths.append(burnInLength)
            s2 = "\nfor chain #"+str(head['chainNum'])
            s2 += "\nTotal number of points in the chain = "+str(len(chiSqs))+"\n"
            s2 += "Burn-in length = "+str(burnInLength)+"\n"
            s+=s2
            log.debug(s2)
    log.debug(s)
    return (s,burnInLengths)

def burnInStripper(mcmcFnames,burnInLengths):
    """
    Strip the initial burn-in off each chain.
    """
    newFnames=[]
    for i in range(0,len(mcmcFnames)):
        filename = mcmcFnames[i]
        burnIn = burnInLengths[i]
        if os.path.exists(filename):
            (head,data) = loadFits(filename)
            ##strip burn-in and write to new fits     
            log.debug("Before stripping burn-in, file had "+str(len(data[:,0]))+" samples")        
            hdu = pyfits.PrimaryHDU(data[burnIn:,:])
            hdulist = pyfits.HDUList([hdu])
            newHead = hdulist[0].header
            log.debug("Before stripping burn-in, file had "+str(len(hdulist[0].data[:,0]))+" samples")
            for key in head:
                newHead[key]=(head[key],head.comments[key])
            n = os.path.splitext(os.path.basename(filename))
            newFname = os.path.join(os.path.dirname(mcmcFnames[0]),n[0]+'_BIstripped.fits')
            hdulist.writeto(newFname)
            log.info("output file written to:below\n"+newFname)
            hdulist.close()
            newFnames.append(newFname)
            log.debug("burn-in stripped file written to:\n"+newFname)
    return newFnames

def gelmanRubinCalc(mcmcFileList,nMCMCsamp=1,returnStrOnly=True):
    """
    Calculate Gelman-Rubin statistic for each varying param.
    Input MUST be the list of more than one MCMC chain files.
    """
    GRs=[]
    Ts = []
    grStr = '\n'+'-'*21+"\nGelman-Rubin Results:\n"+'-'*21+'\n'
    try:
        Lcfloat = float(nMCMCsamp)
        if os.path.exists(mcmcFileList[0]):
            log.info("Starting to calculate R&T")
            ###########################################################
            ## stage 1 ->  load up values for each param in each chain.
            ## allStg1vals = [chain#, param#, (mean,variance,Lc)]
            ## stage 2 ->  Use them to compare between chains 
            ##             then calc R and T.
            ###########################################################
            (head,data) = loadFits(mcmcFileList[0])
            (paramList,paramStrs,_) = getParStrs(head,latex=False)
            
            Nc = len(mcmcFileList)
            ##start stage 1
            allStg1vals=np.zeros((Nc,len(paramList),3))
            for i in range(0,len(mcmcFileList)):
                log.debug("Starting to calc chain #"+str(i)+' GR values')
                (head,data) = loadFits(mcmcFileList[i])
                allStg1vals[i,:,2]=data.shape[0]
                for j in range(0,len(paramList)):
                    log.debug("calculating stage 1 of GR for chain #"+str(i)+", param: "+paramStrs[j])
                    allStg1vals[i,j,0]=np.mean(data[:,paramList[j]])
                    allStg1vals[i,j,1]=np.var(data[:,paramList[j]])
            ##start stage 2         
            rHighest = 0
            tLowest = 1e9       
            rHighStr = ''
            tLowStr = ''
            for j in range(0,len(paramList)):
                try:
                    log.debug("Starting stage 2 for param: ("+str(paramList[j])+"/"+str(len(paramList))+"), "+paramStrs[j])
                    Ncfloat = float(Nc)
                    ##calc R
                    W = 0
                    for i in range(0,Nc):
                        #Lcfloat = float(allStg1vals[i,j,2])
                        W+=(Lcfloat/(Lcfloat-1.0))*allStg1vals[i,j,1]
                    W=W/Ncfloat
                    V = np.mean(allStg1vals[:,j,1])+(Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])
                    R=np.NaN
                    if W!=0:
                        R = np.sqrt(V/W)
                    GRs.append(R)
                    ##calc T
                    #B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0]*allStg1vals[:,j,2])
                    B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])*Lcfloat
                    #Uses the mean of Lc values
                    T = np.NAN
                    if B!=0:
                        #T = np.mean(allStg1vals[:,j,2])*Ncfloat*np.min([(V/B),1.0])
                        T = Lcfloat*Ncfloat*np.min([(V/B),1.0])
                    Ts.append(T)       
                    grStr+=paramStrs[j]+" had R = "+str(R)+", T = "+str(T)+'\n'
                    if T<tLowest:
                        tLowest=T
                        tLowStr="Lowest T = "+str(T)+' for '+paramStrs[j]+'\n'
                    if R>rHighest:
                        rHighest=R
                        rHighStr="Highest R = "+str(R)+' for '+paramStrs[j]+'\n'
                except:
                    log.critical("an error occured while performing stage 2 of GR calc for param "+paramStrs[j])
            grStr+='\nWorst R&T values were:\n'+rHighStr+tLowStr+'\n'
        else:
            log.critical("Gelman-Rubin stat can NOT be calculated as file does not exist!!:\n"+mcmcFileList[0])
    except:
        log.critical("Gelman-Rubin stat FAILED to be calculated for some reason")
    if returnStrOnly:
        return grStr
    else:
        return (GRs,Ts,grStr)
def jdToGcal(jd):
    "Convert standard Julian Date to a double representing its Gregorian date."
    (y,m,d,s) = jdcal.jd2gcal(2400000.5, jd-2400000.5)
    yrs = float('%.2f'%(y+(m/12.0)+(d/days_per_year)+(s/sec_per_year)))
    return yrs
    
def timeStrMaker(deltaT):
    """
    Convert a time in seconds into a nicer string.
    """
    timeStr = ''
    if type(deltaT) is not str:
        if deltaT>60:
            if deltaT>(60*60):
                hr = int(deltaT//(60*60))
                minutes = int((deltaT-hr*60*60)/60.0)
                timeStr = str(hr)+" hours and "+str(minutes)+" minutes"
            else:
                minutes = int(deltaT//(60))
                timeStr = str(minutes)+" minutes and "+str(int(deltaT-minutes*60))+" seconds"
        else:
            timeStr = str(int(deltaT))+' seconds'
    return timeStr

def dateStrMaker(now,numberSecondsLater,militaryTime=False):
    """
    Returns now+numberSecondsLater in a nice string.
    
    param now: datetime.datetime object
    """
    extra = datetime.timedelta(seconds=numberSecondsLater)
    laterDate = now+extra
    minute = ''
    if laterDate.minute<10:
        minute = '0'+str(laterDate.minute)
    else:
        minute = str(laterDate.minute)
    s =""
    if laterDate.day!=now.day:
        d = {1:'st',2:"nd",3:"rd",21:'st',22:"nd",23:'rd',31:'st'}
        for i in range(1,30):
            if i not in d:
                d[i]='th'
        s+=" on the "+str(laterDate.day)+d[laterDate.day]+","
    else:
        s+=" today,"
    if militaryTime:
        s+=" at about "+str(laterDate.hour)+":"+minute
    else:
        ampm = 'AM'
        hr = laterDate.hour
        if laterDate.hour>12:
            ampm='PM'
            hr = laterDate.hour-12
        s+=" at about "+str(hr)+":"+minute+" "+ampm
    return s

def getParStrs(head,latex=True,getALLpars=False):
    """
    Return matching paramList, paramStrs, paramFileStrs for provided header.
    
    latex=True will return latex formated strs for use in Python code.
    getALLpars=True will signal to get all params, not just the ones that were varying.
    """
    paramList = getParInts(head)    
    paramFileStrs = ['m1','m2','parallax','Omega','e','To', 'Tc','P','i','omega','a_total','chiSquared','K']
    paramStrs = ['m1 [Msun]','m2 [Msun]','Parallax [mas]','Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]','Probability','K [m/s]']
    #paramStrs = ['m1 [Msun]','m2 [Msun]','Parallax [mas]','Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]','chiSquared','K [m/s]']
    if latex:
        paramStrs = [r'$m_1{\rm [M}_{\odot}{\rm ]}$',r'$m_2{\rm [M}_{\odot}{\rm ]}$',r'$\varpi{\rm [mas]}$',r'$\Omega{\rm [deg]}$',r'$e$',r'$T_o{\rm  [JD]}$', r'$T_c{\rm  [JD]}$',r'$P{\rm  [Yrs]}$',r'$i{\rm  [deg]}$',r'$\omega{\rm  [deg]}$',r'$a_{{\rm total}} {\rm [AU]}$',r'$posterior_{prob}$',r'$K{\rm  [m/s]}$']
        #paramStrs = [r'$m_1{\rm [M}_{\odot}{\rm ]}$',r'$m_2{\rm [M}_{\odot}{\rm ]}$',r'$\varpi{\rm [mas]}$',r'$\Omega{\rm [deg]}$',r'$e$',r'$T_o{\rm  [JD]}$', r'$T_c{\rm  [JD]}$',r'$P{\rm  [Yrs]}$',r'$i{\rm  [deg]}$',r'$\omega{\rm  [deg]}$',r'$a_{{\rm total}} {\rm [AU]}$',r'$\chi^2$',r'$K{\rm  [m/s]}$']
    if head["nRVdsets"]>0:
        for dataset in range(1,head["nRVdsets"]+1):
            paramFileStrs.append('offset_'+str(dataset))
            if latex:
                paramStrs.append(r"$\gamma_{{\rm "+str(dataset)+"}}{\\rm [m/s]}$")
            else:
                paramStrs.append('offset '+str(dataset)+' [m/s]')        
    ## clean up lists if not returning ALL   
    ## else set paramList to contain its for ALL params
    if (len(paramFileStrs)>len(paramList))and(getALLpars==False):
            paramStrsUse = []
            paramFileStrsUse = []
            for par in paramList:
                paramStrsUse.append(paramStrs[par])
                paramFileStrsUse.append(paramFileStrs[par])
            paramStrs = paramStrsUse
            paramFileStrs = paramFileStrsUse 
    elif getALLpars:
        paramList=np.arange(0,len(paramStrs))
    return (paramList,paramStrs,paramFileStrs)
    
def cleanUp(settings,stageList,allFname):
    """
    Clean up final directory after simulation completes
    """
    #os.mkdir(settings['finalFolder'])
    
    ## write best orbit to a fits file for minimal customPost.py plotting
    bst = findBestOrbit(allFname, bestToFile=False, findAgain=False,by_ln_prob=stageList[-1]=='emcee')
    outFname = os.path.join(settings['finalFolder'],'bestFit.fits')
    writeFits(outFname,bst,settings,clob=False)
    
    delFiles = []
    fnames = glob.glob(os.path.join(settings['finalFolder'],"pklTemp-*"))
    for i in range(0,len(fnames)):
        delFiles.append(fnames[i])
    ##get chain data filenames to delete
    if settings["delChains"]:
        for stage in stageList:
            fnames = glob.glob(os.path.join(settings['finalFolder'],"outputData"+stage+"*.fits"))
            for i in range(0,len(fnames)):
                delFiles.append(fnames[i])
            if stage=='SA':
                fnames = glob.glob(os.path.join(settings['finalFolder'],"SAtempData*.fits"))
                for i in range(0,len(fnames)):
                    delFiles.append(fnames[i])
    ##get combined data filename to delete
    if settings["delCombined"]:
        delFiles.append(allFname)
    if (settings['rmBurn'])and(settings['nChains']>1):
        ##the burn-in was stripped in the final file, so kill the non-stripped version if it exists
        nm = os.path.join(os.path.dirname(allFname),'combinedMCMCdata.fits')
        if os.path.exists(nm):
            delFiles.append(nm)
            
    ##try to delete files
    rmFiles(delFiles)
    
def summaryFile(settings,stageList,finalFits,clStr,burnInStr,bestFit,grStr,effPtsStr,iacStr,allTime,postTime,durationStrings,MCmpo,SAmpo,STmpo,MCMCmpo,emcee_mpo):
    """
    Make a txt file that summarizes the results nicely.
    """
    summaryFname = os.path.join(settings['finalFolder'],'RESULTS.txt')
    if os.path.exists(summaryFname):
        f = open(summaryFname,'a')
    else:
        f = open(summaryFname,'w')
    head = loadFits(finalFits,noData=True)
    totalSamps = head['NSAMPLES']
    (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False,getALLpars=True)
    (paramListCleaned,paramStrsCleaned,paramFileStrsCleaned) = getParStrs(head,latex=False)
    t = datetime.date.today()
    f.write("Date ExoSOFT completed this run: "+t.strftime('%b %d, %Y')+'\n')
    f.write("\n"+"*"*80+"\noutRoot:  "+settings['outRoot']+"\n"+"*"*80+"\n")
    f.write('\n'+'-'*7+'\nBasics:\n'+'-'*7)
    f.write("\nComplete set of parameters that varied directly:\n"+'-'*48)
    f.write('\nTheir integers:\n'+repr(paramListCleaned))
    f.write('\nTheir name+units:\n'+repr(paramStrsCleaned))
    f.write('\nTheir file name postpends:\n'+repr(paramFileStrsCleaned))
    f.write("\n\nIntegers of those used in Astrometry/DI nu calc: ")
    try:
        f.write(repr(settings['DIvars']))
        f.write("\nIntegers of those used in RV nu calc: "+repr(settings['RVvars']))
    except:
        log.debug("no DIvars or RVvars it seems")
    try:
        ## try to make and write the more advanced summary strings to the file
        nusStr = "\nnu values were: [total,DI,RV] = ["+str(head['NU'])+", "+str(head['NUDI'])+", "+str(head['NURV'])+"]\n"
        f.write(nusStr)
        numEpochsStr = "Number of epochs in data [total, DI, RV, num RV datasets] = "
        numEpochsStr+="["+str(settings['n3Depoch'])+", "+str(settings['nDIepoch'])+", "+str(settings['nRVepoch'])+", "+str(settings['nRVdsets'])+"]\n"
        f.write(numEpochsStr)
        stgNsampStrDict = {"MC":"nSamples","SA":"nSAsamp","ST":"nSTsamp","MCMC":"nSamples","emcee":'nSamples'}
        try:
            totSampStored = int(effPtsStr.split("Total # of steps stored ")[1].split('\n')[0])
        except:
            log.debug('did not work yet')
            totSampStored=1
        numFilesStr = '\nTotal # of files that finished each stage were:\n'
        chiSquaredsStr = '\nBest Reduced Chi Squareds for each stage were:\n'
        flSzStr='\nFile sizes for each stage were:\n'
        for stage in stageList:
            fnames = np.sort(glob.glob(os.path.join(settings['finalFolder'],"outputData"+stage+"*.fits")))
            if (stage=="MCMC")and settings["rmBurn"]:
                fnames = np.sort(glob.glob(os.path.join(settings['finalFolder'],"outputData"+stage+"*BIstripped.fits")))
            numFilesStr+=stage+' = '+str(len(fnames))+", each with "+str(settings[stgNsampStrDict[stage]])+" samples\n"
            if len(fnames)>0:
                chiSquaredsStr+=stage+" = ["
                flSzStr+=stage+" = ["
                for fname in fnames: 
                    try:
                        bestFit2 = findBestOrbit(fname,bestToFile=False,findAgain=True,by_ln_prob=stage=='emcee')
                        
                        #### calc chi squared for these params
                        Model = ExoSOFTmodel(settings)
                        paramsLast = copy.deepcopy(Model.Params.stored_to_direct(bestFit2))
                        _ = ln_posterior(paramsLast, Model)
                        #Model.chi_squared_3d
                        #Model.chi_squared_di
                        #Model.chi_squared_rv
                        chiSquaredsStr+=str(Model.chi_squared_3d/float(head['NU']))+', '
                        #chiSquaredsStr+=str(bestFit2[11]/float(head['NU']))+', '
                        flSzStr+=fileSizeHR(fname)+', '
                    except:
                        log.error("A problem occurred while trying to find best fit of:\n"+fname)
                chiSquaredsStr = chiSquaredsStr[:-2]+']\n'
                flSzStr = flSzStr[:-2]+']\n'
        numFilesStr+="\n"+"*"*61
        numFilesStr+="\nThe final combined file was for a total of "+str(totalSamps)+" samples"
        numFilesStr+="\nThe 'saveInt' was "+str(settings['saveInt'])+" ie. every "+str(settings['saveInt']-1)+" was skipped. And burn-in stripped if requested."
        numFilesStr+="\n i.e. "+str(len(fnames))+"*"+str(settings[stgNsampStrDict[stage]])+"*(1/"+str(settings['saveInt'])+")-(burn-in) = "
        numFilesStr+=str(totSampStored)+" final stored samples.\n"+"*"*61+'\n'
        flSzStr+="Final fits data file was "+fileSizeHR(finalFits)+'\n'
        f.write(numFilesStr)
        f.write(flSzStr)
        f.write(chiSquaredsStr)
        bestStr = '\n'+'-'*21+"\nBest fit values were:\n"+'-'*21+'\n'
        ############################################
        ## calculate chi squareds for the best fit #
        ############################################
        Model = ExoSOFTmodel(settings)
        _params = copy.deepcopy(Model.Params.stored_to_direct(bestFit))
        _ = ln_posterior(_params, Model)
        reducedDI = Model.chi_squared_di/head['NUDI']
        reducedRV = Model.chi_squared_rv/head['NURV']
        reduced3D = Model.chi_squared_3d/head['NU']
        
        for i in range(len(bestFit)):
            if i==2:
                bestStr+=paramStrs[2]+" = "+str(bestFit[2])
                if (bestFit[2]!=0):
                    bestStr+=", OR  "+str(1.0/(bestFit[2]/1000.0))+'[PC]\n'
                else:
                    bestStr+='\n'
            elif i==1:
                if bestFit[1]<0.1:
                    mJupMult=(const.M_sun.value/const.M_jup.value)
                    bestStr+=paramStrs[1]+" = "+str(bestFit[1])+", OR "+str(bestFit[1]*mJupMult)+' in [Mjupiter]\n'
                else:
                    bestStr+=paramStrs[1]+" = "+str(bestFit[1])+'\n'
            elif i==0:
                if bestFit[1]==0:
                    bestStr+="m1 = m_total [Msun] = "+str(bestFit[0])+'\n'
                else:
                    bestStr+=paramStrs[0]+" = "+str(bestFit[0])+'\n'
            elif i in [5,6]:
                bestStr+=paramStrs[i]+" = "+str(bestFit[i])+", OR "+str(bestFit[i]-2400000.5)+' in [MJD]\n'
            else:
                bestStr+=paramStrs[i]+" = "+str(bestFit[i])+'\n'
        bestStr+='-'*36+"\nValues differing for the two bodies:\n"+'-'*36+'\n'
        bestStr+="(Masses, Parallax, e, To, P, i, chiSquared and RV offsets)"
        bestStr+=" all SAME as above."
        bestStr+="\nONLY Omega, omega and semi-major axis values differ.\n"
        omega1 = bestFit[3]+180.0
        if omega1>360.0:
            omega1-=360.0
        bestStr+="Omega_1 [deg] = "+str(omega1)+'\n'
        bestStr+="Omega_2 [deg] = "+str(bestFit[3])+'\n'
        omega1 = bestFit[9]+180.0
        if omega1>360.0:
            omega1-=360.0
        bestStr+="omega_1 [deg] = "+str(omega1)+'\n'
        bestStr+="omega_2 [deg] = "+str(bestFit[9])+'\n'
        if bestFit[1]!=0:
            a1 = bestFit[10]/(1.0+(bestFit[0]/bestFit[1]))
            bestStr+="a_1 [AU] = "+str(a1)+'\n'
            bestStr+="a_2 [AU] = "+str(bestFit[10]-a1)+'\n'
        bestStr+='\n'+'*'*90+'\nBEST REDUCED CHISQUAREDS: [total,DI,RV] = ['+str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]\n"+'*'*90
        f.write(bestStr)
    except:
        log.critical("A problem occured while trying to produce advanced summary strings.")
    f.write('\n'+clStr)
    f.write('\n'+burnInStr)
    f.write('\n'+grStr)
    ## Add note about GR validity if starting all MCMC chains at same position.
    if (settings['stages']=='MCMC') or (settings['strtMCMCatBest']==True):
        s ="\n NOTE: as the MCMC chains were all started at the same position, "+\
        "it questions the validity of the resulting Gelman-Rubin statistics.  "+\
        "\nThus, they must only be considered a loose indicator/estimator of convergence." 
        f.write(s)
    f.write(effPtsStr)
    f.write('\n'+iacStr)
    f.write('\n'+'-'*40+'\nAverage acceptance rates for each stage:\n'+'-'*40+'\n')
    if MCmpo!=None:
        f.write('For MC stage: '+str(np.mean(MCmpo.avgAcceptRates))+'\n')
        #(outFname,params,sigmas,bstChi,avgAcceptRates,acceptStrs) = MCmpo.getBest()
    if SAmpo!=None:
        f.write('For SA stage: '+str(np.mean(SAmpo.avgAcceptRates))+'\n')
        #(outFname,params,sigmas,bstChi,avgAcceptRates,acceptStrs) = SAmpo.getBest()
    if STmpo!=None:
        f.write('For ST stage: '+str(np.mean(STmpo.avgAcceptRates))+'\n')
        #(outFname,params,sigmas,bstChi,avgAcceptRates,acceptStrs) = STmpo.getBest()
    if MCMCmpo!=None:
        f.write('For MCMC stage: '+str(np.mean(MCMCmpo.avgAcceptRates))+'\n')
    if emcee_mpo!=None:
        f.write('For emcee stage: '+str(np.mean(emcee_mpo.avgAcceptRates))+'\n')
        #(outFname,params,sigmas,bstChi,avgAcceptRates,acceptStrs) = MCMCmpo.getBest()
    f.write('\n\n'+'-'*24+'\nDurations of each stage:\n'+'-'*24+'\n'+durationStrings)
    f.write('\nPost-Processing took: '+timeStrMaker(postTime)+'\n')
    f.write('Total simulation took: '+timeStrMaker(allTime)+'\n')
    f.write('\n\nEND OF RESULTS :-D')
    f.close()
    log.info("Summary file written to:\n"+summaryFname)  

def fileSizeHR(filename):
    """
    convert number of bytes to human readable form.
    code modified from post at:
    http://stackoverflow.com/questions/14996453/python-libraries-to-calculate-human-readable-filesize-from-bytes
    """
    nbytes = os.path.getsize(filename)
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def keplersThird(p=0,atot=0,mtot=0):
    """
    Kepler's Third rule.  
    Find the missing value provided you know 2 of them.
    p in [years]
    atot in [AU]
    mtot in [Msun]
    """
    if (atot==0)and(mtot!=0)and(p!=0):
        atot = (((p**2)*(sec_per_year**2)*const.G.value*const.M_sun.value*mtot)/(4.0*np.pi**2))**(1.0/3.0)
        atot = atot/const.au.value
    elif (atot!=0)and(mtot==0)and(p!=0):
        mtot = ((((atot*const.au.value)**3)*4.0*np.pi**2)/((p**2)*(sec_per_year**2)*const.G.value*const.M_sun.value))
    elif (atot!=0)and(mtot!=0)and(p==0):
        p = np.sqrt((((atot*const.au.value)**3)*4.0*(np.pi**2))/(const.G.value*mtot*(sec_per_year**2)*const.M_sun.value))
    else:
        log.critical('More than 1 parameter was zero, so I can not calc K3')
    
    return (p,atot,mtot)
    
def recheckFit3D(orbParams,settings,finalFits='',nus=[]):
    if finalFits!='':
        (head,data) = loadFits(finalFits)
        nu = head['NU']
        nuDI = head['NUDI']
        nuRV = head['NURV']
    elif len(nus)>1:
        nu = nus[0]
        nuDI = nus[1]
        nuRV = nus[2]
    else:
        log.error("nus and finalFits not defined, so cannont calc reduced chiSquareds")
        nu = 1.0
        nuDI = 1.0
        nuRV = 1.0
        
        
    ############################################
    ## calculate chi squareds for the best fit #
    ############################################
    Model = ExoSOFTmodel(settings)
    _params = copy.deepcopy(Model.Params.stored_to_direct(orbParams))
    _ = ln_posterior(_params, Model)
    raw3D = Model.chi_squared_3d
    reducedDI = Model.chi_squared_di/nuDI
    reducedRV = Model.chi_squared_rv/nuRV
    reduced3D = Model.chi_squared_3d/nu  
        
    ## OLD WAY below, delete it soon! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#     ##get the real data
#     realData = loadRealData(diFilename=settings['di_dataFile'],rvFilename=settings['rv_dataFile'],dataMode=settings['data_mode'])
#     ## Make Orbit cpp obj
#     Orbit = tools.cppTools.Orbit()
#     try:
#         pasa = settings["pasa"]
#     except:
#         pasa = False
#     Orbit.loadStaticVars(settings['omegaFdi'],settings['omegaFrv'],settings['lowEcc'],pasa)
#     Orbit.loadConstants(const.G.value,np.pi,const.M_sun.value, days_per_year,sec_per_year,const.au.value)
#     ## ensure orbParams are in required format for Orbit
#     params = []
#     for par in orbParams:
#         params.append(par)
#     params=np.array(params,dtype=np.dtype('d'),order='C')
#     Orbit.loadRealData(realData)
#     predictedData = np.ones((realData.shape[0],3),dtype=np.dtype('d'),order='C')
#     Orbit.calculate(predictedData,params)
#     ## Calculate chi squareds for 3D,DI,RV and update bestPars and bestSumStr if this is better than the best
#     (raw3D, reducedDI, reducedRV, reduced3D) = chiSquaredCalc3D(realData,predictedData,nuDI,nuRV,nu)

    print('(raw3D, reducedDI, reducedRV, reduced3D) = ',repr((raw3D, reducedDI, reducedRV, reduced3D)))
    
def predictLocation(orbParams,settings,epochs=[]):
    
    
    ############################################
    ## calculate chi squareds for the best fit #
    ############################################
    Model = ExoSOFTmodel(settings)
    
    ## make empty inputs for measured data
    nPts = len(epochs)
    predEpochs = np.array(epochs,dtype=np.dtype('d'))
    Model.Data.rapa = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.rapa_err = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.decsa = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.decsa_err = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.rapa_model = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.decsa_model = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.epochs_di = predEpochs

    Model.Data.rv = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.rv_err = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.rv_model = np.ones((nPts),dtype=np.dtype('d'))
    Model.Data.rv_inst_num = np.zeros((nPts),dtype=np.dtype('i'))
    Model.Data.epochs_rv = predEpochs
    
    ## call model to predict data for given epochs
    _params = copy.deepcopy(Model.Params.stored_to_direct(orbParams))
    _ = ln_posterior(_params, Model)
    
    ## resulting measurable astrometry and rv for those epochs and orbital elements
    fit_decsa_model = copy.deepcopy(Model.Data.decsa_model)
    fit_rapa_model = copy.deepcopy(Model.Data.rapa_model)
    fit_rv_model = copy.deepcopy(Model.Data.rv_model)
    
    ## OLD WAY, delete soon!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#     ##get the real data
#     realData = loadRealData(diFilename=settings['di_dataFile'],rvFilename=settings['rv_dataFile'],dataMode= settings['data_mode'])
#     ## Make Orbit cpp obj
#     Orbit = tools.cppTools.Orbit()
#     try:
#         pasa = settings["pasa"]
#     except:
#         pasa = False
#     Orbit.loadStaticVars(settings['omegaFdi'],settings['omegaFrv'],settings['lowEcc'],pasa)
#     Orbit.loadConstants(const.G.value,np.pi,const.M_sun.value, days_per_year,sec_per_year,const.au.value)
#     ## ensure orbParams are in required format for Orbit
#     params = []
#     for par in orbParams:
#         params.append(par)
#     params=np.array(params,dtype=np.dtype('d'),order='C')
#     fakeData = np.ones((len(epochs),7),dtype=np.dtype('d'),order='C')
#     fakeData[:,0]=epochs[:]
#     Orbit.loadRealData(fakeData)
#     predictedData = np.ones((len(epochs),3),dtype=np.dtype('d'),order='C')
#     print("fakeData are:\n"+repr(fakeData))
#     print("predicted epochs data before are:\n"+repr(predictedData))
#     Orbit.calculate(predictedData,params)
#     print("predicted epochs data are:\n"+repr(predictedData))

    print("Predicted astrometry and RV for the epochs: "+repr(predEpochs))
    print("Astrometry (Dec or SA, depending on 'pasa' setting): "+repr(fit_decsa_model))
    print("Astrometry (RA or PA, depending on 'pasa' setting): "+repr(fit_rapa_model))
    print("Radial Velocity: "+repr(fit_rv_model))

    
def chiSquaredCalc3D(realData,modelData,nuDI,nuRV,nu3D,pasa=False): 
    """
    Based on definition, chiSquared=sum((modelVal_i-dataVal_i)^2/(dataError_i^2)) over all values of 'i'.
    This function will do so for DI, RV and 3D sets of data and provide the reduced chi squared for each.
    The raw 3D value will also be returned.
    NOTES: realData is the standard 7 parameter format from loadRealData function, and modelData is 
           the standard 3 param format.
           
    returned (raw3D, reducedDI, reducedRV, reduced3D)
    """   
    ## convert E,N to SA,PA?
    if pasa:
        (PA,PA_error,SA,SA_error) = ENtoPASA(modelData[:,0], 0, modelData[:,1], 0)
        modelData[:,0] = PA
        modelData[:,1] = SA
    diffs = np.concatenate(((realData[:,1]-modelData[:,0]),(realData[:,3]-modelData[:,1]),(realData[:,5]-modelData[:,2])))
    errors = np.concatenate((realData[:,2],realData[:,4],realData[:,6]))
    raw3D = np.sum((diffs**2)/(errors**2))
    diffsDI = np.concatenate(((realData[:,1]-modelData[:,0]),(realData[:,3]-modelData[:,1])))
    errorsDI = np.concatenate((realData[:,2],realData[:,4]))
    diffsRV = (realData[:,5]-modelData[:,2])
    errorsRV = realData[:,6][np.where(diffsRV!=0)]
    rawDI = np.sum((diffsDI[np.where(diffsDI!=0)]**2)/(errorsDI[np.where(diffsDI!=0)]**2))
    rawRV = np.sum((diffsRV[np.where(diffsRV!=0)]**2)/(errorsRV**2))
    return (raw3D,rawDI/nuDI,rawRV/nuRV,raw3D/nu3D)

def copyToDB(settings):
    """
    Copy vital results files to Dropbox.
    """
    fnamesALL = []
    for extension in ['pdf','txt','log']:
        fnames = glob.glob(os.path.join(settings['finalFolder'],"*."+extension))
        log.debug('found files to copy to DB:\n'+repr(fnames))
        for name in fnames:
            fnamesALL.append(name)
    dbDir = os.path.join(settings['dbFolder'],settings['outRoot'])
    if os.path.exists(dbDir):
        log.critical("drobpox directory already exists, so just copying into it.  Please remove manually if you want.")
    else:
        os.mkdir(dbDir)
    log.debug('DB dir is:\n'+repr(dbDir))
    for f in fnamesALL:
        try:
            log.debug('trying to copy file:\n'+repr(os.path.basename(f)))
            shutil.copy(f,os.path.join(dbDir,os.path.basename(f)))
        except:
            log.error('failed to move file:\n'+f+'\nintto DB folder:\n'+dbDir)
    log.info("vital results files copied to DB folder:\n"+dbDir)
    
def getParInts(head):
    """
    convert string version of paramInts into a list again.
    """
    s = head['parInts'] 
    ints = s.split("[")[1].split("]")[0].split(',')
    parInts = []
    for i in ints:
        parInts.append(int(i))  
    return parInts        

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
        
def confLevelFinder(filename, colNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column 
    centered on the median value.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.    
    
    NOTE: This function calculates the confidence levels based on the median, while the
          more common standard is to centere it on the mean.  Doing that though requires a 
          time consuming loop over the data to total up that around the mean...
          Can't think of any faster way, so not doing it for now.
    """
    verboseInternal = False
    log.debug('Inside confLevelFinder')
    outStr=''
    if os.path.exists(filename):
        (dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = dataReader(filename, colNum)
        if len(dataAry>0) or (dataValueStart!=dataValueMid!=dataValueEnd):
            #print 'ln688:confLevelFinder'
            histAry = histMakeAndDump(chiSquareds,dataAry,weight=False,retHist=True)
            #print 'ln690:confLevelFinder'
            [conf68Vals,conf95Vals,range99,range100] = histConfLevels(histAry)
            #print 'ln692:confLevelFinder'
            
#             #Convert data array to a sorted numpy array
#             dataAry = np.sort(dataAry)
#             size = dataAry.size
#             mid=size//2
#             minVal = np.min(dataAry)
#             maxVal = np.max(dataAry)
#                 
#             minLoc68=mid-int(float(size)*0.683)//2
#             if minLoc68<0:
#                 minLoc68 = 0
#             maxLoc68 = mid+int(float(size)*0.683)//2
#             if maxLoc68>(size-1):
#                 maxLoc68 = size
#             minLoc95=mid-int(float(size)*0.958)//2
#             if minLoc95<0:
#                 minLoc95 = 0
#             maxLoc95= mid+int(float(size)*0.958)//2
#             if maxLoc95>(size-1):
#                 maxLoc95 = size
#             
#             conf68Vals = [dataAry[minLoc68],dataAry[maxLoc68]]
#             conf95Vals = [dataAry[minLoc95],dataAry[maxLoc95]]
#             conf68ValsRough=[]
#             conf95ValsRough=[]
#             
#             if ((len(conf68Vals)==0) or (len(conf95Vals)==0)):
#                 if (len(conf68Vals)==0):
#                     log.error('confLevelFinder: ERROR!!! No FINE 68.3% confidence levels were found')
#                     if (len(conf68ValsRough)==0):
#                         log.error('confLevelFinder: ERROR!!! No ROUGH 68% confidence levels were found, so returning [0,0]')
#                         conf68Vals = [0,0]
#                     else:
#                         conf68Vals = conf68ValsRough
#                         log.error("confLevelFinder: Had to use ROUGH 68% [68,69] as no FINE 68.3% was found. So, using range "+repr(conf68Vals))                
#                 if (len(conf95Vals)==0):
#                     log.error('confLevelFinder: ERROR!!! No FINE 95.4% confidence levels were found')
#                     if (len(conf95ValsRough)==0):
#                         log.error('confLevelFinder: ERROR!!! No ROUGH 95% confidence levels were found, so returning [0,0]')
#                         conf95Vals = [0,0]
#                     else:
#                         conf95Vals = conf95ValsRough
#                         log.error("confLevelFinder: Had to use ROUGH 95% [95,96] as no FINE 95.4% was found. So, using range "+repr(conf95Vals))
        else:
            ## There was no useful data, so return values indicating that
            dataAry=bestDataVal=dataMedian=dataValueStart
            conf68Vals = [dataValueStart,dataValueStart]
            conf95Vals = [dataValueStart,dataValueStart]
            chiSquareds = 0
            log.error("confLevelFinder: Entire column had a constant value of "+str(dataValueStart))
        #print 'ln741:confLevelFinder'
        mJupMult=(const.M_sun.value/const.M_jup.value)
        s = "\nFinal Range values:\nTOTAL "+repr([range100[0],range100[1]])
        s+= '\n95%   '+repr(conf95Vals)+"\n68%   "+repr(conf68Vals)+'\n'
        s+= "   median,      68.3% error above,   68.3% error below\n"
        s+= str(dataMedian)+',  +'+str(conf68Vals[1]-dataMedian)+',   '+str(conf68Vals[0]-dataMedian)
        s+= "\nAverage 68.3% error = +/- "+str((conf68Vals[1]-conf68Vals[0])/2.0)
        if False:
            s+= "\nWidth Estimates:\n((68.3% error range)/median)x100% = "+str(((conf68Vals[1]-conf68Vals[0])/abs(dataMedian))*100.0)+"%"
            s+= "\n((68.3% error range)/(95.4% error range))x100% = "+str(((conf68Vals[1]-conf68Vals[0])/(conf95Vals[1]-conf95Vals[0]))*100.0)+"%"
            s+= "\n((68.3% error range)/(Total range))x100% = "+str(((conf68Vals[1]-conf68Vals[0])/(range100[1]-range100[0]))*100.0)+"%"
        if (colNum==1) and (dataMedian<0.1):
            s+='\n'+"~"*55+'\nOR in units of Mjupiter:\n'
            s+=str(dataMedian*mJupMult)+',  +'+str(mJupMult*(conf68Vals[1]-dataMedian))+',   '+str(mJupMult*(conf68Vals[0]-dataMedian))
            s+="\nAverage 68.3% error = +/- "+str(mJupMult*((conf68Vals[1]-conf68Vals[0])/2.0))+'\n'
            s+="~"*55
        if colNum in [5,6]:
            s+='\n'+"~"*55+'\nOR in units of MJD:\n'
            s+=str(dataMedian-2400000.5)+',  +'+str(conf68Vals[1]-dataMedian)+',   '+str(conf68Vals[0]-dataMedian)
            s+='\n'+"~"*55
        outStr+=s
        s=s+75*'-'+'\n Leaving confLevelFinder \n'+75*'-'+'\n'
        log.debug('\n'+s)
        ## return requested form of results
        if (returnData and returnChiSquareds and (returnBestDataVal==False)):
            returnList =  ([conf68Vals,conf95Vals],dataAry, chiSquareds)
        elif (returnData and returnChiSquareds and returnBestDataVal):
            returnList =   ([conf68Vals,conf95Vals],dataAry, chiSquareds, bestDataVal,outStr) ##MODIFIED
        elif (returnData and (returnChiSquareds==False)and (returnBestDataVal==False)):
            returnList =   ([conf68Vals,conf95Vals],dataAry)
        elif (returnData and (returnChiSquareds==False) and returnBestDataVal):
            returnList =   ([conf68Vals,conf95Vals],dataAry, bestDataVal,outStr) ##MODIFIED
        elif ((returnData==False) and returnChiSquareds):
            returnList =   ([conf68Vals,conf95Vals], chiSquareds)
        elif ((returnData==False) and returnChiSquareds and returnBestDataVal):
            returnList =   ([conf68Vals,conf95Vals], chiSquareds, bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and returnBestDataVal):
            returnList = ([conf68Vals,conf95Vals], bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and (returnBestDataVal==False)):
            returnList =   [conf68Vals,conf95Vals]
        #print 'ln800:confLevelFinder'
        return returnList 
        
    else:
        log.critical( "confLevelFinder: ERROR!!!! file doesn't exist")            
         
def histConfLevels(histAry):
    """
    Calculates the locations on the x-axis where 68.5% and 95.4% of the data
    lie above a certain probability.  
    x-axis points need to be centered on the middle of the bins, rather than
    the edges.
    
    Returns the min and max of the data points for both percentages.
    [range68,range95,range100]
    """
    ##avoid float comparison errors for near identicle numbers by rounding all
    if True:
        for i in range(0,len(histAry[:,0])):
            histAry[i,0] = round(histAry[i,0],12)
    halfStepSize = round((histAry[1,0]-histAry[0,0])*0.5,12)
    res = abs(histAry[1,0]-histAry[0,0])/100.0
    xs=np.arange(histAry[0,0],histAry[-1,0],res)
    
    ## interpolate
    f = interpolate.interp1d(histAry[:,0], histAry[:,1],kind='cubic')
    if xs[-1]>np.max(histAry[:,0]):
        xs[-1]=np.max(histAry[:,0])
        #print 'new xs[-1] = '+str(xs[-1])
    ys=f(xs)
    #print 'ln831:histConfLevels'
    ysorted=np.sort(ys)
    #print 'ln832:histConfLevels'
    ysorted=ysorted[::-1]
    [got68,got95,got99,vls68,vls95,vls99]=[False,False,False,[],[],[]]
    #print 'ln833:histConfLevels'
    for i in range(0,len(xs)):
        vls = ys[np.where(ys>ysorted[i])]
        perc = np.sum(vls)/np.sum(ysorted)
        if (perc>0.683) and (got68==False):
            got68=True
            vls68=xs[np.where(ys>ysorted[i])]
        if (perc>0.954) and (got95==False):
            got95=True
            vls95=xs[np.where(ys>ysorted[i])]
        if (perc>0.99) and (got99==False):
            got99=True
            vls99=xs[np.where(ys>ysorted[i])]
    #print 'ln846:histConfLevels'
    range68=[np.min(vls68),np.max(vls68)]
    range95=[np.min(vls95),np.max(vls95)]
    range99=[np.min(vls99),np.max(vls99)]
    range100 = [np.min(xs),np.max(xs)]
    retAry = [range68,range95,range99,range100]
    #print 'ln852:histConfLevels'
    ##round to top or bottom value if confidence levels are next to them
    for i in range(0,4):
        retAry[i][0] = round(retAry[i][0],12)
        retAry[i][1] = round(retAry[i][1],12)
        if round(retAry[i][0],12)==round(np.min(xs),12):
            retAry[i][0]=histAry[0,0]-halfStepSize
        if round(retAry[i][1],12)==round(np.max(xs),12):
            #print 'rounding up '+str(retAry[i][1])+' to '+str(retAry[i][1]+halfStepSize)
            retAry[i][1]=halfStepSize+histAry[-1,0]
    #print 'ln862:histConfLevels'
    return retAry
             
def findBestOrbit(filename,bestToFile=True,findAgain=False, by_ln_prob=True):        
    """
    Find the orbital elements for the best fit in a ExoSOFT format fits file.
    """             
    bestFname = os.path.join(os.path.dirname(filename),'bestOrbitParams.txt')
    #print   bestFname
    gotIt = False
    if os.path.exists(bestFname)and(findAgain==False):
        try:
            orbBest = np.loadtxt(bestFname,delimiter=',')
            log.debug("Using previously found best orbit in file:\n"+bestFname)   
            gotIt=True
        except:
            log.error("Tried to load previously found best orbit from file, but failed, so will find it from data again.")
    if gotIt==False:
        log.debug("trying to find best orbit in file:\n"+filename)   
        (head,data) = loadFits(filename)
        if by_ln_prob:
            probBest = np.max(data[:,11])
            loc = np.where(data[:,11]==probBest)
        else:
            chiBest = np.min(data[:,11])
            loc = np.where(data[:,11]==chiBest)
        orbBest = data[loc[0][0],:]
        log.info("Best fit found to be:\n"+repr(orbBest))
        if bestToFile:
            f = open(bestFname,'w')
            f.write(nparyTolistStr(orbBest,brackets=False)+'\n')
            f.close()
            log.info("Best fit params written to :\n"+bestFname)
    return orbBest
                            
def unitlessSTD(ary):
    """
    Calculate the bias corrected standard deviation, then divide by the mean to make it unitless.
    """
    if len(ary)<1:
        log.error('no elements passed to unitlessSTD, so returning 1e6.')
        return 1e6
    if len(ary)==1:
        log.error('only 1 element passed into unitlessSTD, so returning 0.')
        return 0
    else:
        if type(ary)!=np.ndarray:
            if type(ary)==list:
                ary = np.array(ary)
            else:
                return 0.0
        bcstd = np.sqrt((1.0/(len(ary)-1.0))*np.sum(abs(ary-ary.mean())**2))
        return bcstd/np.mean(ary)

def nparyTolistStr(ary,brackets=True,dmtr=','):
    s=''
    if brackets:
        s+='['
    for val in ary:
        s+=str(val)+dmtr
    s=s[:-1]
    if brackets:
        s+=']'
    return s

def ENtoPASA(E, E_error, N, N_error):
    """
    Will calculate the Separation and Position Angles for a given East and North, including their errors.
    PA and error will be in [deg], with SA and error in ["]
    :returns: (PA,PA_error,SA,SA_error)
    """
    verbose = False
    PA = np.degrees(np.arctan2(E,N))
    #NOTE: both math.atan2 and np.arctan2 tried with same results, both produce negatives rather than continuous 0-360 
    #thus, must correct for negative outputs
    if type(PA)!=np.ndarray:
        if PA<0:
            PA = PA+360.0
    else:
        for i in range(0,len(PA)):
            if PA[i]<0:
                PA[i]=PA[i]+360.0
    
    SA = np.sqrt(E**2.0 + N**2.0)
    
    PA_error=0
    SA_error=0
    if (E_error==0)or(N_error==0):
        if False:
            print("either the E or N error value was zero, so setting the PA and SA return errors to zero!!")
    else:
        top = abs(E/N)*np.sqrt((E_error/E)**2.0 + (N_error/N)**2.0)
        btm = 1.0+(E/N)**2.0
        PA_error = abs(np.degrees(top/btm))

        top = SA*(abs(E*E_error)+abs(N*N_error))
        btm = E**2.0+N**2.0
        SA_error = abs(top/btm)
    if verbose:
        print(repr((E, E_error, N, N_error))+" -> "+repr((PA,PA_error,SA,SA_error)))   
    return (PA,PA_error,SA,SA_error)

def PASAtoEN(PA,PA_error,SA,SA_error):
    """
    Convert provided Position Angle and Separation Angle, and their errors, into 
    RA and DEC with errors.  These are the same equations for calculating 
    x and y in the Thiele-Innes orbit fitting.  Remember that x and y are 
    flipped in that fitting approach due to how Thiele defined the coord 
    system when deriving the equations used.
    
    NOTE: this can also be used to calculate x and y used in Thiele-Innes
          With East=RA=y and North=DEC=x.  
    
    :returns: (E, E_error, N, N_error)
    """
    verbose = False
    printForExcel = False
    N = SA*np.cos(np.radians(PA))
    E = SA*np.sin(np.radians(PA))
    
    E_error=0
    N_error=0
    if (SA_error==0)or(PA_error==0):
        if verbose:
            print("either the PA and SA error value was zero, so setting the E or N return errors to zero!!")
    else:
        tempA = (SA_error/SA)**2.0
        tempB = ((np.cos(np.radians(PA+PA_error))-np.cos(np.radians(PA))) / np.cos(np.radians(PA)))**2.0
        N_error = abs(N*np.sqrt(tempA+tempB))
        
        # Another way to calculate the error, but the one above is belived to be more correct 
        tempA2 = (SA_error*np.cos(np.radians(PA)))**2.0
        tempB2 = (SA*np.sin(np.radians(PA))*np.radians(PA_error))**2.0
        N_error2 = np.sqrt(tempA2+tempB2)
        
        tempC = (SA_error/SA)**2.0
        tempD = ((np.sin(np.radians(PA+PA_error))-np.sin(np.radians(PA))) / np.sin(np.radians(PA)))**2.0
        E_error = abs(E*np.sqrt(tempC+tempD))
        
        # Another way to calculate the error, but the one above is belived to be more correct 
        tempC2 = (SA_error*np.sin(np.radians(PA)))**2.0
        tempD2 = (SA*np.cos(np.radians(PA))*np.radians(PA_error))**2.0
        E_error2 = np.sqrt(tempC2+tempD2)
        
        if verbose:
            print('N_error2-N_error = '+str(N_error2-N_error))
            print('E_error2-E_error = '+str(E_error2-E_error))
            print('E_error2 = '+str(E_error2)+', E_error = '+str(E_error))
            print('N_error2 = '+str(N_error2)+', N_error = '+str(N_error)+"\n")
    if printForExcel:
        print('E = '+str(E))
        print('E error= '+str(E_error))
        print('N = '+str(N))
        print('N error = '+str(N_error))
        
        
    return (E, E_error, N, N_error)

def m2siniCalc(K,p,m1,e):
    """To calculate the commonly quoted m2sin(i) value assuming m2<<m1.
    follows:
    m2sin(i) = K*m1^(2/3)*sqrt(1-e^2)*[P/2piG]^(1/3)
    units:
    K [m/s]
    m1 [Msun]
    p [yrs]
    """
    pSecs = p*sec_per_year
    m1KG = m1*const.M_sun.value
    A = K*((pSecs/(2*np.pi*const.G.value))**(1.0/3.0))
    B = (m1KG**(2.0/3.0))*np.sqrt(1-e**2)
    C = (1.0/const.M_jup.value)
    m2sini=A*B*C
    return m2sini

def m2siniRangeCalc():
    """
    A basic custom tool to find the ranges of m2sin(i) values from RV only 
    runs.  Need to go into code below to update e, p, and K values along with 
    fixed m1.
    """
    m1 = 1.09 #Msun
    es = [0.1036,0.1556]
    ps = [10.501,13.096] #yrs
    ks = [705.328,721.005] #m/s
    mn = 1e6
    mx = 0
    for e in es:
        for p in ps:
            for k in ks:
                m2sini = m2siniCalc(k, p,m1 , e)
                if m2sini>mx:
                    mx = m2sini
                if m2sini<mn:
                    mn = m2sini
    print('min = '+str(mn)+", max = "+str(mx))             
    
def semiMajAmp(m1,m2,inc,ecc,p):
    """
    K = [(2*pi*G)/p]^(1/3) [m2*sin(i)/m2^(2/3)*sqrt(1-e^2)]
    units:
    K [m/s]
    m1 [Msun]
    m2 [Mj]
    p [yrs]
    inc [deg]
    """
    pSecs = p*sec_per_year
    m1KG = m1*const.M_sun.value
    A = ((2.0*np.pi*const.G.value)/pSecs)**(1.0/3.0)
    B = (m2*const.M_jup.value*np.sin(np.radians(inc)))
    C = m1KG**(2.0/3.0)*np.sqrt(1.0-ecc**2.0)
    print('Resulting K is '+repr(A*(B/C)))
    #return A*(B/C)
                    
    
    
    
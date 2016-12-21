#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp  or kylemede@gmail.com
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from astropy.io import fits as pyfits
import os
import re
import sys
import pickle
import shutil
import warnings
import KMlogger
from six.moves import range

warnings.simplefilter("error")
log = KMlogger.getLogger('main.rwTools',lvl=100,addFH=False)  


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
    
def copytree(src, dst):
    """
    Recursively copy a directory and its contents to another directory.
    
    WARNING: this is not advised for higher level folders as it can also copy subfolders 
    thus leading to a very large copy command if not careful.
    
    Code taken and simplified from:
    http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
    """
    skipStrs = [".git",".pyc",".py~",".sty",".so",".cxx",".o",".h~",".cc~"]
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            try:
                shutil.copytree(s, d)
                log.debug("Copying directory:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))
        else:
            try:
                #Check if filepath contains one of the skip 
                #strs and copy only if it is fine.
                fine = True
                for skipStr in skipStrs:
                    if skipStr in item:
                        fine=False
                if fine:
                    shutil.copy2(s, d)
                    log.debug("Copying file:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))      

def copyCodeFiles(src, dst,settingsFiles=None):
    """
    For copying the code and settings files used to the output directory.
    """
    #First copy the code directory/tree
    copytree(src, dst)
    #now copy the settings files
    if settingsFiles is not None:
        if (type(settingsFiles)!=list)and(type(settingsFiles)==str):
            settingsFiles = [settingsFiles]
        setdst = os.path.join(dst,'settingsFilesUsed')
        os.mkdir(setdst)
        for f in settingsFiles:
            s=f
            d=os.path.join(setdst, os.path.basename(f))
            try:
                shutil.copy2(s, d)
                log.debug("Copying:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))

def dataReader(filename, colNum=0):
    """
    Read in the data for a single column of data.
    """
    try:
        verboseInternal = False
        ## First get ranges of param and ChiSquared values
        log.debug('\nOpening and finding ranges for data in column # '+str(colNum))
        
        ## Check if file has useful data for that column#
        (head,data) = loadFits(filename)
        if head!=False:
            TotalSamples=data.shape[0]
            dataAry = data[:,colNum]
            chiSquareds = data[:,11]
            bestDataVal = dataAry[np.where(chiSquareds==np.min(chiSquareds))][0]          
            return (dataAry,chiSquareds,[bestDataVal,np.median(dataAry),dataAry[0],dataAry[len(dataAry)//2],dataAry[-1]])  
    except:
        log.critical("a problem occured while trying to load data file:\n"+filename)
        return False

def loadDIdata(filename):
    """
    Load the astrometry data into a numpy array.
    
    file format:
    title 
    column headers
    data
    .
    .
    .
    
    Data must be in the columns:
    obsDate[JD] x["] x_error["] y["] y_error["]
    OR if pasa key in settingsAdvanced ==True, then:
    obsDate[JD] PA[deg] PA_error[deg] SA["] SA_error["]
    """
    diData = []
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        for line in lines:
            #print "line was:'"+line+"'"
            #log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[3].replace('.','',1).replace('-','',1).isdigit():
                    diData.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])  
                    #print repr([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
        diData = np.array(diData)
    except:
        log.critical("a problem occured while trying to load DI data. \nPlease check it is formatted correctly.")
    return diData

def loadDIdata_cy(filename):
    """
    Load the astrometry data into a numpy array.
    
    file format:
    title 
    column headers
    data
    .
    .
    .
    
    Data must be in the columns:
    obsDate[JD] x["] x_error["] y["] y_error["]
    OR if pasa key in settingsAdvanced ==True, then:
    obsDate[JD] PA[deg] PA_error[deg] SA["] SA_error["]
    """
    epochs_di = []
    rapa = []
    rapa_err = []
    decsa = []
    decsa_err = []
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        for line in lines:
            #print "line was:'"+line+"'"
            #log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[3].replace('.','',1).replace('-','',1).isdigit():
                    epochs_di.append(float(line.split()[0]))
                    rapa.append(float(line.split()[1]))
                    rapa_err.append(float(line.split()[2]))
                    decsa.append(float(line.split()[3]))
                    decsa_err.append(float(line.split()[4]))
                    #print repr([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
        epochs_di = np.array(epochs_di,dtype=np.dtype('d'))
        rapa = np.array(rapa,dtype=np.dtype('d'))
        rapa_err = np.array(rapa_err,dtype=np.dtype('d'))
        decsa = np.array(decsa,dtype=np.dtype('d'))
        decsa_err = np.array(decsa_err,dtype=np.dtype('d'))
    except:
        log.critical("a problem occured while trying to load DI data. \nPlease check it is formatted correctly.")
    return (epochs_di, rapa, rapa_err, decsa, decsa_err)
    
def loadRVdata(filename):
    """
    Load the radial velocity data into a numpy array.  Provided jitter values will be added in quadrature with 
    the errors.
    
    file format:
    title 
    column headers
    data
    .
    .
    Empty line between data sets
    data
    .
    .
    
    Data must be in the columns:
    obsDate[JD] RV[m/s] RV_error[m/s] jitter[m/s] datasetNumber[int]
    NOTE: datasetNumber is optional, if not provided they will be automatically set to 0,1,2... following the order of the data in the file.
          If jitter is not provided, it will be assumed zero.
    """
    rvData=[]
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        rvData = []
        datasetNumLast = 0
        jitterLast = 0
        lastWasDataLine=False
        thisIsDataLine = False
        for line in lines:
            #print "line = "+line
            lastWasDataLine=thisIsDataLine
            thisIsDataLine=False
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[1].replace('.','',1).replace('-','',1).isdigit():
                    thisIsDataLine=True
                    curDataAry = [float(line.split()[0]),float(line.split()[1])]
                    #if jitter was provided on first line of data set
                    if len(line.split())>3:
                        try:
                            jitterLast = float(line.split()[3])
                        except:
                            log.error("could not convert 4th element of split into jitter.  4th element was: "+str(line.split()[3]))
                    curDataAry.append(np.sqrt(float(line.split()[2])**2+jitterLast**2))
                    #if datasetNum was provided on first line of data set
                    if len(line.split())>4:
                        try:
                            datasetNumLast = float(line.split()[4])
                        except:
                            log.error("could not convert 5th element of split into datasetNum.  5th element was: "+str(line.split()[4]))
                    curDataAry.append(datasetNumLast)
                    #print repr(curDataAry)
                    rvData.append(curDataAry)
            if lastWasDataLine and (thisIsDataLine==False):
                jitterLast = 0
                datasetNumLast+=1
                #print 'incrementing datasetNum'
        rvData = np.array(rvData)
    except:
        log.critical("a problem occured while trying to load RV data.  \nPlease check it is formatted correctly.")
    return rvData

def loadRVdata_cy(filename):
    """
    Load the radial velocity data into a numpy array.  Provided jitter values 
    will be added in quadrature with the errors.
    
    file format:
    title 
    column headers
    data
    .
    .
    Empty line between data sets
    data
    .
    .
    
    Data must be in the columns:
    obsDate[JD] RV[m/s] RV_error[m/s] jitter[m/s] datasetNumber[int]
    NOTE: datasetNumber is optional, if not provided they will be automatically set to 0,1,2... following the order of the data in the file.
          If jitter is not provided, it will be assumed zero.
    """
    epochs_rv = []
    rv = []
    rv_err = []
    rv_inst_num = []
    try:
        if filename[-4:]!='.dat':
            filename = filename+'.dat'
        fl = open(filename, 'r')
        lines = fl.readlines()
        fl.close()
        datasetNumLast = 0
        jitterLast = 0
        lastWasDataLine=False
        thisIsDataLine = False
        for line in lines:
            #print "line = "+line
            lastWasDataLine=thisIsDataLine
            thisIsDataLine=False
            if len(line.split())>2:
                if line.split()[0].replace('.','',1).isdigit() and line.split()[1].replace('.','',1).replace('-','',1).isdigit():
                    thisIsDataLine=True
                    epochs_rv.append(float(line.split()[0]))
                    rv.append(float(line.split()[1]))
                    #if jitter was provided on first line of data set
                    if len(line.split())>3:
                        try:
                            jitterLast = float(line.split()[3])
                        except:
                            log.error("could not convert 4th element of split into jitter.  4th element was: "+str(line.split()[3]))
                    rv_err.append(np.sqrt(float(line.split()[2])**2+jitterLast**2))
                    #if datasetNum was provided on first line of data set
                    if len(line.split())>4:
                        try:
                            datasetNumLast = float(line.split()[4])
                        except:
                            log.error("could not convert 5th element of split into datasetNum.  5th element was: "+str(line.split()[4]))
                    rv_inst_num.append(datasetNumLast)
            if lastWasDataLine and (thisIsDataLine==False):
                jitterLast = 0
                datasetNumLast+=1
                #print 'incrementing datasetNum'
        epochs_rv = np.array(epochs_rv,dtype=np.dtype('d'))
        rv = np.array(rv,dtype=np.dtype('d'))
        rv_err = np.array(rv_err,dtype=np.dtype('d'))
        rv_inst_num = np.array(rv_inst_num,dtype=np.dtype('i'))
    except:
        log.critical("a problem occured while trying to load RV data.  \nPlease check it is formatted correctly.")
    return (epochs_rv, rv, rv_err, rv_inst_num)
    
def loadRealData(diFilename='',rvFilename='',dataMode='3D'):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data,sorted into cronological order.
    Both filenames must be the full path to the files.
    """
    realData = False
    try:
        diEpochs = []
        rvEpochs = []
        if dataMode!='RV':
            #print 'using diFilename = '+diFilename        
            if os.path.exists(diFilename):
                diData = loadDIdata(diFilename)
                diEpochs = diData[:,0]
        if dataMode!='DI':
            #print 'using rvFilename = '+rvFilename
            if os.path.exists(rvFilename):
                rvData = loadRVdata(rvFilename)
                rvEpochs = rvData[:,0]
        #print 'rvData = '+repr(rvData)
        #for i in range(0,rvData.shape[0]):
        #    print 'ORIG rv data = '+str(rvData[i,0])+', '+str(rvData[i,1])+", "+str(rvData[i,2])+", "+str(rvData[i,3])
        #print repr(diData)
        ##load in epochs from both sets, sort and kill double entries
        epochsTemp = np.concatenate((diEpochs,rvEpochs))
        epochsTemp.sort()
        epochs = []
        for epoch in epochsTemp:
            if epoch not in epochs:
                epochs.append(epoch)
        epochs = np.array(epochs)
        realData = np.zeros((epochs.shape[0],8))
        #print 'ln171'
        ##set error values to 1e6 which signals not to calculate the predicted version in orbit.cc
        realData[:,2]=1e6
        realData[:,4]=1e6
        realData[:,6]=1e6
        realData[:,0]=epochs[:]
        #print 'ln177'
        #print repr(diEpochs)
        #print repr(diData[:,0])
        #print repr(epochs)
        for i in range(epochs.shape[0]):
            if len(diEpochs)>0:
                if epochs[i] in diData[:,0]:
                    pos = np.where(diData[:,0]==epochs[i])[0]
                    if len(pos)>1:
                        pos = [pos[0]]
                        log.critical('More than 1 set of DI data for epoch '+\
                                     str(epochs[i])+'.  Only using first!!')
                    realData[i,1:5]=diData[pos,1:]
            
            if len(rvEpochs)>0:
                if epochs[i] in rvData[:,0]:
                    pos = np.where(rvData[:,0]==epochs[i])[0]
                    if len(pos)>1:
                        pos = [pos[0]]
                        log.critical('More than 1 set of RV data for epoch '+\
                                     str(epochs[i])+'.  Only using first!!')
                    realData[i,5:]=rvData[pos,1:]
        #print 'dataMode'+dataMode+'->realData = '+repr(realData)
        #for i in range(0,realData.shape[0]):
        #    print 'realData = '+str(realData[i,0])+', '+str(realData[i,5])+", "+str(realData[i,6])+", "+str(realData[i,7])
    except:
        log.critical("An error occured while trying to load data!!")
    return realData
            
def load_settings(settings_in,advanced_settings_in,priors_in):
    """  
    This is to load the settings file from the path provided, then load in the 
    ExoSOFTpriors object.  It will first see if there is a priors.py in same 
    folder as the settings.py, else, this will load the default priors.py 
    packaged with ExoSOFT.
    NOTE: filenames must be exactly settings.py and priors.py.
    """
    
    if priors_in==None:
        ## load the default priors object
        from .priors import ExoSOFTpriors as priors_in
        #print(repr(ExoSOFTpriors))
        
    if advanced_settings_in==None:
        log.debug("no advanced settings provided by user, loading in defaults.")
        from .advanced_settings import advanced_settings_dict
        advanced_settings_in = advanced_settings_dict

    settings = settings_in
    
    # merge two settings dicts
    for key in advanced_settings_in:
        settings[key] = advanced_settings_in[key]
    
    # push priors into settings
    settings['ExoSOFTpriors'] = priors_in
            
    #######################################################
    ## determine argPeriOffsetRV and argPeriOffsetDI values
    #######################################################
    omegaFdi = 0
    ## ExoSOFT assumes the RV data was measured from the primary's spectra, and 
    ## the orbit of the companion is being fit, NOT the orbit of the primary.
    ## Thus, there is a 180deg shift forced to account for this.
    omegaFrv=180.0
    #now update due to fixed argPeriPlus values
    omegaFdi+=settings['omega_offset_di']
    omegaFrv+=settings['omega_offset_rv']
    settings['omegaFdi'] = (omegaFdi,"Total fixed val added to DI omega in model")
    settings['omegaFrv'] = (omegaFrv,"Total fixed val added to RV omega in model")
    log.debug("Setting fixed omega offsets to:\nomegaFdi = "+str(omegaFdi)+"\nomegaFrv = "+str(omegaFrv))
    
    #for key in settings:
    #    print(key+' = '+repr(settings[key]))
    #sys.exit('shirt')
    return settings

def loadFits(filename,noData=False):
    """
    Load in a fits file written by ExoSOFT.
    Return (header dict, data)
    """
    if os.path.exists(filename):
        f = pyfits.open(filename,'readonly')
        head = f[0].header
        if noData==False:
            data = f[0].data
        f.close()
    else:
        log.critical("fits file does not exist!!! filename =\n"+str(filename))
        head=data=False
    if noData:
        return head
    else:
        return (head,data)

def writeFits(baseFilename,data,settings,clob=False):
    """
    Data will be written to a fits file with a single PrimaryHDU,
    with the .header loaded up with the tuples from the settings 
    and .data = provided data.
    File will be stored in the 'finalFolder' directory from the settings.
    If data variable is a string, this function will assume it is a filename 
    of where the data is stored in a .npy file, and load it in.
    
    NOTE: lots of error msg checking that can be removed if further testing 
          does not show same issues that caused me to put them in.
    """
    outFname=''
    errMsg = ''
    notThere=False
    try:
        ##check if data is a .npy filename
        if type(data)==str:
            if os.path.exists(data):
                dataFname = data
                errMsg='about to try to load data from input filename'
                #print "input datafile was: "+dataFname
                if '.npy' in dataFname:
                    data = np.load(dataFname)
                elif ('.dat' in dataFname)or('.txt' in dataFname):
                    print('.dat or .txt file recognised and will try to use loadtxt')
                    try:
                        data = np.loadtxt(dataFname)
                    except:
                        print("about to try and load as normal file")
                        print(os.fstat(dataFname))
                        f = open(dataFname,'r')
                        for line in f.readlines:
                            print(line)
                    #print 'successfully loaded with loadtxt'
                    #print repr(data.shape())
                else:
                    log.critical("input data type was a filename, but type was not .npy .txt or .dat")
                os.remove(dataFname)
                log.debug("just removed data file from disk:\n"+dataFname)
            else:
                notThere=True
        if (len(data)>0)and(notThere==False):
            errMsg="got data and about to make outFname"
            if '.fits' not in baseFilename:
                baseFilename=baseFilename+'.fits'
            outFname = os.path.join(settings['finalFolder'],baseFilename)
            errMsg = 'Got outFname, about to make PrimaryHDU for data of type:\n'+repr(type(data))+'\noutFname was:\n'+outFname
            hdu = pyfits.PrimaryHDU(data)
            errMsg = 'Got outname and data into a PrimaryHDU, about to make an HDUList. outFname was:\n'+outFname
            hdulist = pyfits.HDUList([hdu])
            header = hdulist[0].header
            errMsg = 'HDUList created, about to try and load up header with keys from commentsDict'
            ##load up header with tuples from settings
            commentsDict = settings['commentsDict']
            #print repr(commentsDict)
            for key in settings:
                if key in commentsDict:  
                    good_key = True
                    #print repr(key)
                    #print key+' = '+repr(settings[key])
                    #try to store key as is, else store as a string          
                    try:    
                        #print 'trying straight up'
                        header[key]=settings[key]
                    except:
                        try:
                            #print 'trying with repr'
                            header[key]=repr(settings[key])
                            log.debug("Key '"+key+"' had a value that could not be stored as is, into the fits header, so it was converted to a string.")
                        except:
                            #print 'trying 3rd version'
                            good_key = False
                            if len(key)>8:
                                log.debug("Key '"+key+"' could not be set as a header key as it was too long! ie. greater than 8 characters.  So ignoring it.")
                            else:
                                log.debug("Key '"+key+"' could not be be set as a header key for some reason...")
                    if good_key:
                        com = commentsDict[key]
                        if len(com)>47:
                            s = "comment too long for pyfits headers:"+com
                            s+=". still writting to file, but cutting >47th characters off"
                            log.warning(s)
                            com = com[:46]
                        header.comments[key] = com
                        #print key+' = '+repr((header[key],header.comments[key]))
            errMsg = 'keys from commentsDict loaded into header, about to try to write hdu to disk at:\n'+outFname
            if os.path.exists(outFname):
                if clob:
                    log.debug("clob==True, but file exists, so deleting file:\n"+outFname)
                    rmFiles([outFname])
                else:
                    log.error("That file already exists on disk!!\n")
            if os.path.exists(outFname)==False:
                hdulist.writeto(outFname)
                log.debug("output file written to:below\n"+outFname)
                hdulist.close()
                errMsg=" all done, should be no errors after this."
            ## check resulting fits file header
            if False:
                f = pyfits.open(os.path.join(settings['finalFolder'],baseFilename),'readonly')
                head = f[0].header
                f.close()
                if False:
                    for key in head:
                        print(key+' = '+repr((header[key],header.comments[key])))
                        #print 'type(header[key] = '+repr(type(header[key]))
                print('\n\nEntire Header as a repr:\n'+repr(head))
        else:
            log.info("No data to write to file:\n"+baseFilename+'\n')
    except:
        m = "Could not write file to disk for some reason.  Last errMsg was:\n"
        log.error(m+errMsg+'\n')
    return outFname

def periodicDataDump(filename,d):
    """
    dump a ndarray to disk.  If first time, just dump it.
    Else, load current ary and cat d to it before dumping.
    """
    old=True
    if len(d)!=0:
        if os.path.exists(filename):
            if old:
                d0 = np.load(filename)
                np.save(filename,np.concatenate((d0,d)))
            else:
                with open(filename,'a') as outfile:
                    for i in range(0,d.shape[0]):
                        outstr = ''
                        for val in d[i]:
                            outstr+='%.14g  ' % (val)
                        outstr += '\n'
                        outfile.write(outstr)#nparyTolistStr(d[i],brackets=False,dmtr=' ')+'\n')
                        #outfile.write(re.sub("\n ","\n",re.sub("[\\[\\]]","",np.array2string(d,precision=16))))
        else:
            if old:
                np.save(filename,d)
            else:
                with open(filename,'w') as outfile:
                    for i in range(0,d.shape[0]):
                        outstr = ''
                        for val in d[i]:
                            outstr+='%.14g  ' % (val)
                        outstr += '\n'
                        outfile.write(outstr)#nparyTolistStr(d[i],brackets=False,dmtr=' '))
                        #outfile.write(re.sub("\n ","\n",re.sub("[\\[\\]]","",np.array2string(d))))
                        #print nparyTolistStr(d[i],brackets=False,dmtr=' ')
                    #raise IOError('\n\n'+s)
                


def pklIt(settings,dataObj, rootFnm):
    """
    Pickle something to a file and truncate file if it already exists.
    """
    pklFname = os.path.join(settings["pklDir"],rootFnm+".pkl")
    with open(pklFname,'w+') as f:
        pickle.dump(dataObj, f)

def unPklIt(settings,rootFnm):
    """
    To load back in pickled objects during post-processing of ExoSOFT.
    Useful for customPost work.
    """
    pklFname = os.path.join(settings["pklDir"],rootFnm+".pkl")
    #print pklFname
    with open(pklFname,'r') as f:
        o = pickle.load(f)
    return o

def reloadMpoROs(settings):
    pklDir = settings["pklDir"]
    pklDict = {'MCMC':None,'MC':None,'SA':None,'ST':None}
    for key in pklDict:
        if os.path.exists(os.path.join(pklDir,key+"mpoRO.pkl")):
            pklDict[key]=unPklIt(settings, key+'mpoRO')
    return [pklDict['MC'],pklDict['SA'],pklDict['ST'],pklDict['MCMC']]

def combineFits(filenames,outFname):
    """
    combine the data in multiple ExoSOFT fits files together.
    Used primarily for after multi-process runs.
    """
    nFiles = len(filenames)
    (head0,dataALL) = loadFits(filenames[0])
    for filename in filenames[1:]:
        (head,data) = loadFits(filename)
        dataALL = np.concatenate((dataALL,data))
    hdu = pyfits.PrimaryHDU(dataALL)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    for key in head0:
        if key=='NSAMPLES':
            ##load in total number of samples for this combined file
            header['NSAMPLES'] = (int(head0['NSAMPLES'])*nFiles,head0.comments['NSAMPLES'])
        else:
            header[key] = (head0[key],head0.comments[key])
    hdulist.writeto(outFname)
    hdulist.close()
    log.info("output file written to:below\n"+outFname)
        
def renameFits(filenameIn,filenameOut,killInput=True,overwrite=True):
    """
    Load input into a new fits HDU, write to output name, close, rm input file.
    """
    goodToGo = True
    if os.path.exists(filenameOut):
        if overwrite:
            rmFiles([filenameOut])
        else:
            goodToGo=False
            log.critical("File already exists, either set overwrite==True, or there be more careful.")
    if goodToGo:
        (head,data) = loadFits(filenameIn)
        hdu = pyfits.PrimaryHDU(data)
        hdulist = pyfits.HDUList([hdu])
        header = hdulist[0].header
        for key in head:
            header[key] = (head[key],head.comments[key])
        hdulist.writeto(filenameOut)
        hdulist.close()
        log.info("output file written to:below\n"+filenameOut)
        if killInput:
            rmFiles([filenameIn])
    
def rmFiles(files):
    if type(files)==str:
        files = [files]
    if type(files)==list:
        ##try to delete files
        for fname in files:
            try:
                if os.path.exists(fname):
                    log.debug('Deleting file: '+os.path.basename(fname))
                    os.remove(fname) 
            except:
                log.error('Failed to delete file: '+os.path.basename(fname))
    else:
        log.debug("files passed into rmFiles was not of list or str type")
    
def writeBestsFile(settings,pars,sigs,bstChiSqr,stage):

    filename = os.path.join(settings['finalFolder'],'best'+stage+'paramsAndSigs.txt')
    f = open(filename,'w')
    f.write("Best-fit between all "+stage+" chains had a reduce chi squared of "+str(bstChiSqr)+'\n')
    f.write("\nIts parameters were:\n")
    f.write(nparyTolistStr(pars)+'\n')
    f.write("\n\nIts sigmas were:\n")
    #double check clean up sigs of pars that were not varying
    paramInts = settings['paramInts']
    if stage=='MC':
        sigs_out = np.zeros(len(pars))
    else:
        try:
            sigs_out = [0]*len(sigs)
            for i in range(0,len(pars)):
                if i not in paramInts:
                    sigs_out[i] = 0
        except:
            sigs_out = sigs
            log.debug("Failed to zero non-varying sigmas, but no biggie. It's a benign bug.")
    f.write(nparyTolistStr(sigs_out)+'\n')
    f.close()
    log.info("Best fit params and sigmas from "+stage+" stage were written to :\n"+filename)    
    
    
#END OF FILE
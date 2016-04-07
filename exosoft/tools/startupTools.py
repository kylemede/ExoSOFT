import sys
import shutil
import os
import copy
import numpy as np
import exoSOFTlogger
import generalTools as genTools
import readWriteTools as rwTools
import constants as const
import warnings
warnings.simplefilter("error")

log = exoSOFTlogger.getLogger('main.suTools',lvl=100,addFH=False) 

def startup(argv,ExoSOFTdir,rePlot=False):
    """
    Perform the following vital start up steps (and return filled out and cleaned up settings):
    -Figure out important directories
    -Copy settings files to temp directory to combine them into the master settings dict
    -Get master settings dict
    -Make output folder
    -Remake the SWIG tools?
    -Copy all exosoft code into output dir for emergencies.
    -check data exists    
    -push all comments from tuples into a sub dictionary with key 'commentsDict'
    -Check range min and max values, including updates for 'lowecc' mode.
    -find which parameters will be varying.
    -Check if start startParams and startSigmas in dictionary make sense.
    NOTE: have ExoSOFTdir handled with setup.py??
    """    
    ## Pull in settings filename prepend from command line args, if provided
    inArg = ''
    if len(argv)>1:
        try:
            inArg = argv[1]
        except:
            print '\nWarning: Only able to handle a single argument on the cmd !!\n'    
    ## Load up the required specific directory paths in dict
    if '/' in inArg:
        #assume full path
        settFilePath = inArg
    else:
        examplesDir = os.path.join(ExoSOFTdir.split('exosoft')[0],'examples/')
        if "settings.py" in inArg:
            settFileNm = inArg
        else:
            settFileNm = inArg+"settings.py"
        settFilePath = os.path.join(examplesDir,settFileNm)

    if os.path.exists(settFilePath)==False:
        log.critical('Critical error occured while trying to load in the settings.  Quiting exosoft!!')
        #***************************************************************************************************
        s = "THAT SETTINGS FILE DOES NOT EXIST!!"
        s+= "\nARGUMENT PROVIDED TO exosoft WAS: "+inArg
        s+= "\nTHIS MUST BE A PREPEND TO A SETTINGS FILE IN THE STANDARD DIRECTORY, OR THE FULL PATH TO A SETTINGS FILE."
        s+="\n\n!!EXITING exosoft!!"
        sys.exit(s)
        #***************************************************************************************************
    else:
        settings = rwTools.loadSettings(ExoSOFTdir,settFilePath)
        log.setStreamLevel(settings['logLevel'])
        settings['ExoSOFTdir']=ExoSOFTdir
        settings['settingsDir']=os.path.dirname(settFilePath)
        settings['settFilePath']= settFilePath
        ## Make a directory (folder) to place all the files from this simulation run
        settings['finalFolder'] = os.path.join(settings['outDir'],settings['outRoot'])
        ##if not doing a re-post analysis with customPost.py
        if rePlot==False:
            if os.path.exists(settings['finalFolder']):
                if settings['logLevel']<50: ## Handle this with a 'clob' bool in dict??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    print '\n'+'$'*50
                    print 'WARNING!! the folder:\n"'+settings['finalFolder']+'"\nALREADY EXISTS!'
                    print 'You can overwrite the data in it, or exit this simulation.'
                    YN = raw_input('OVERWRITE current folder (y/n):')
                else:
                    YN = 'y'
                if (('y' in YN) or ('Y' in YN)):
                    shutil.rmtree(settings['finalFolder'])
                    os.mkdir(settings['finalFolder'])
                    dbDir = os.path.join(settings['dbFolder'],settings['outRoot'])
                    if os.path.exists(dbDir):
                        shutil.rmtree(dbDir)
                else: #elif (('n' in YN) or ('N' in YN)):
                    sys.exit()
                if settings['logLevel']<50:
                    print '$'*50+'\n'
            else:
                os.mkdir(settings['finalFolder'])
            if False:
                for key in settings:
                    print key+' = '+repr(settings[key])
            ## copy all of current code to output directory
            codeCopyDir = os.path.join(settings['finalFolder'],'codeUsed')
            os.mkdir(codeCopyDir)
            log.debug('Copying all files in the RESULTS folder over to output folder:\n '+codeCopyDir)
            setFiles = [settings['settFilePath'],settings['RVdataFile'],settings['DIdataFile']]
            genTools.copyCodeFiles(settings['ExoSOFTdir'], codeCopyDir,setFiles)
        #########################################################################################
        ## Check parameter range settings make sense for data provided and mode of operation.   #
        ## Then load up a list of the parameters to vary during simulation.                     #
        ## Note: Originally this was done in simulator startup, but thought better to move here.#
        #########################################################################################
        realData = rwTools.loadRealData(diFilename=settings['DIdataFile'],rvFilename=settings['RVdataFile'],dataMode=settings['dataMode'])
        if (type(realData)!=list)and(type(realData)!=np.ndarray):
            log.critical('Critical error occured while trying to load real data files.  Quiting exosoft!!')
            #***************************************************************************************************
            s = "THERE WAS A PROBLEM LOADING THE REAL DATA."
            s+= "IF 3D OR DI DATA MODE RQUESTED: MAKE SURE A DI FILENAME IN SETTINGS FILES EXISTS."
            s+= "IF 3D OR RV DATA MODE RQUESTED: MAKE SURE A RV FILENAME IN SETTINGS FILES EXISTS."
            s+= "IF THEY EXIST, MAKE SURE THEIR FORMATS MATCH THAT DESCRIBED IN THE readWriteTools.loadDIdata AND loadRVdata FUNCTIONS."
            s+="\n\n!!EXITING exosoft!!"
            sys.exit(s)
            #***************************************************************************************************
        ## push all comments from tuples into a sub dictionary to ensure all  
        ## values for requested keys are just the value with no comments.
        commentsDict = {}
        for key in settings:
            if type(settings[key])==tuple:
                    commentsDict[key] = settings[key][1]
                    settings[key] = settings[key][0]
        settings['commentsDict'] = commentsDict
        ##check there are matching number of RV datasets and provided min/max vals for offsets
        if np.min(realData[:,6])<1e6:
            numVmins=len(settings['vMINs'])
            if numVmins==0:
                numVmins=1
            if np.max(realData[:,7])!=(numVmins-1):
                log.error("THE NUMBER OF vMINs DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
                               "please check the vMINs/vMAXs arrays in the simple settings file\n"+\
                               "to make sure they have matching lengths to the number of RV datasets.")
        if  settings['TMAX']==settings['TMIN']==-1:
            ## set T range to [earliest Epoch-max period,earliest epoch]
            settings['TMAX']=np.min(realData[:,0])
            settings['TMIN']=np.min(realData[:,0])-settings['PMAX']*const.daysPerYear
        ##In DI mode can only find Mtotal, thus push all mass into M1 and kill M2
        if settings['dataMode']=='DI':
            settings['mass1MIN']=settings['mass1MIN']+settings['mass2MIN']
            settings['mass1MAX']=settings['mass1MAX']+settings['mass2MAX']
            settings['mass2MIN']=0
            settings['mass2MAX']=0
            log.debug("DI dataMode, so pushed all mass range vals into M1 and set ones for M2 to zero")
        if settings['OmegaMIN']==None:
            settings['OmegaMIN'] = 0.0
        if settings['OmegaMAX']==None:
            if settings['dataMode']=='3D':
                settings['OmegaMAX'] = 360.0
            else:
                settings['OmegaMAX'] = 180.0
        if settings['eMIN']==None:  
            settings['eMIN'] = 0.0
        if settings['eMAX']==None:
            settings['eMAX'] = 0.98
        if settings['incMIN']==None:
            settings['incMIN'] = 0.0
        if settings['incMAX']==None:
            settings['incMAX'] = 180.0  
        if settings['omegaMIN']==None:
            settings['omegaMIN'] = 0.0
        if settings['omegaMAX']==None:
            settings['omegaMAX'] = 360.0
        ##load up range min,max and sigma arrayS
        rangeMaxs = [settings['mass1MAX'],\
                     settings['mass2MAX'],\
                     settings['paraMAX'],\
                     settings['OmegaMAX'],\
                     settings['eMAX'],\
                     settings['TMAX'],\
                     settings['TMAX'],\
                     settings['PMAX'],\
                     settings['incMAX'],\
                     settings['omegaMAX'],\
                     0,\
                     0,\
                     settings['KMAX']]
        rangeMins = [settings['mass1MIN'],\
                     settings['mass2MIN'],\
                     settings['paraMIN'],\
                     settings['OmegaMIN'],\
                     settings['eMIN'],\
                     settings['TMIN'],\
                     settings['TMIN'],\
                     settings['PMIN'],\
                     settings['incMIN'],\
                     settings['omegaMIN'],\
                     0,\
                     0,\
                     settings['KMIN']]
        ##start with uniform sigma values
        sigSize = settings['strtSig']
        sigmas = [sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,0,0,sigSize]
        if len(settings['vMINs'])!=len(settings['vMAXs']):
            log.critical("THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!")
            #***************************************************************************************************
            s="THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!\n"
            s+="PLEASE CHECK THE ADVANCED SETTINGS FILES AND FIX THIS"
            s+="!\n\n!!EXITING exosoft!!"
            sys.exit(s)
            #***************************************************************************************************
        for i in range(0,len(settings['vMINs'])):
            sigmas.append(sigSize)
            rangeMins.append(settings['vMINs'][i])
            rangeMaxs.append(settings['vMAXs'][i])
        rangeMaxs = np.array(rangeMaxs)
        rangeMins = np.array(rangeMins)
        ##For lowEcc case, make Raw min/max vals for param drawing during MC mode
        rangeMaxsRaw = copy.deepcopy(rangeMaxs)
        rangeMinsRaw = copy.deepcopy(rangeMins)
        if settings['lowEcc']:
            ## run through the possible numbers for e and omega to find min/max for RAW versions
            ## Only relevent for MC and the begining jumps of SA,  
            ## Otherwise the values are converted back to omega and e and the original ranges are used for the check.
            fourMin=1e6
            fourMax=-1e6
            nineMin=1e6
            nineMax=-1e6
            for omeg in range(int(rangeMins[9]*10),int(rangeMaxs[9]*10),1):
                omega = float(omeg)/10.0
                for e in range(int(rangeMins[4]*100),int(rangeMaxs[4]*100),1):
                    ecc = float(e)/100.0
                    four = np.sqrt(ecc)*np.sin((np.pi/180.0)*omega)
                    nine = np.sqrt(ecc)*np.cos((np.pi/180.0)*omega)
                    if four>fourMax:
                        fourMax = four
                    if four<fourMin:
                        fourMin = four
                    if nine>nineMax:
                        nineMax = nine
                    if nine<nineMin:
                        nineMin = nine
            rangeMaxsRaw[9] = nineMax
            rangeMaxsRaw[4] = fourMax
            rangeMinsRaw[9] = nineMin
            rangeMinsRaw[4] = fourMin
        ## figure out which parameters are varying in this run.
        ## Don't vary atot or chiSquared ever, and take care of TcEqualT and Kdirect cases
        paramInts = []
        for i in range(0,len(rangeMins)):
            if (i!=10)and(i!=11):
                if (i>12):
                    if settings['dataMode']!='DI':
                        if rangeMaxs[i]!=0:
                            paramInts.append(i) 
                elif (i==8)or(i==12):
                    if (settings['dataMode']!='RV')or(settings['Kdirect']==False):
                        if (rangeMaxs[8]!=0)and(i==8):
                            paramInts.append(8)
                    elif settings['Kdirect']:
                        if (rangeMaxs[12]!=0)and(i==12):
                            paramInts.append(12)                                           
                elif i<4:
                    if (settings['dataMode']!='RV'):
                        if(rangeMaxs[i]!=0):
                            paramInts.append(i)
                    elif ([settings,'Kdirect']==False)and((i!=3)and(i!=2)):
                        if(rangeMaxs[i]!=0):
                            paramInts.append(i)
                elif rangeMaxs[i]!=0:
                    if (i==5)or(i==6):
                        if settings['TcStep']and(i!=5):
                                paramInts.append(i)
                        elif (settings['TcStep']==False)and(i!=6):
                                paramInts.append(i)
                    else:
                        paramInts.append(i)
            
        ## push all these important parameter related items into the dict for later use.
        settings['realData'] = realData
        settings['rangeMinsRaw'] = rangeMinsRaw
        settings['rangeMaxsRaw'] = rangeMaxsRaw
        settings['rangeMins'] = rangeMins
        settings['rangeMaxs'] = rangeMaxs
        settings['paramInts'] = np.array(paramInts)
        ## use modePrep to make sure all is ready for the stages requested
        settings = modePrep(settings,sigmas)
        
        return settings
        

def modePrep(settings,sigmas):
    """
    Check if start startParams and startSigmas in dictionary make sense.
    Arrays must exist, be the right length, and have non-zeros values for all varying parameters.     
    Mode logic:
    
    """
    startParams = settings['startParams']
    startSigmas = settings['startSigmas']
    paramInts = settings['paramInts']
    rangeMaxs = settings['rangeMaxs']
    autoMode = settings['autoMode']
    
    ##check if startParams in settings file are useful
    gotParams = False
    if (type(startParams)==list)or(type(startParams)==np.ndarray):
        if type(startParams)==list:
            startParams = np.array(startParams)
        if len(startParams)==len(rangeMaxs):
            i=0
            gotParams = True
            while (i<len(startParams))and(gotParams==True):
                if i in paramInts:
                    if startParams[i]==0:
                        gotParams=False
                i+=1
        else:
            gotParams=False
    if gotParams==False:
        log.info("Original startParams in settings files were not usable, so setting to False.")
        startParams = False
    #check if startSigmas in settings file are useful
    gotSigmas = False
    if (type(startSigmas)==list)or(type(startSigmas)==np.ndarray):
        if len(startSigmas)==len(rangeMaxs):
            i=0
            gotSigmas = True
            while (i<len(startSigmas))and(gotSigmas==True):
                if i in paramInts:
                    if startSigmas[i]==0:
                        gotSigmas=False
                i+=1
        else:
            gotSigmas=False       
    ## check out logic between provided parameters and requested mode/stages
    if gotParams==False:
        if autoMode==False:
            #No useful params proided, so check if manual settings make sense
            if (settings['stages'] in ['ST','MCMC']):
                log.critical('ST or MCMC mode requested, but not starting parameters provided.  Quiting exosoft!!')
                #***************************************************************************************************
                s="MUST PROVIDE USEFUL STARTPARAMS IN SIMPLE SETTINGS DICT FOR ST or MCMC MODE, ELSE NOTHING TO START CHAINS WITH"
                s+="\n\nRUN IN AUTO MODE, OR PERFORM A ROUND OF SA OR SAST TO GET USEFUL VALUES TO STARTING VALUES."
                s+="!\n\n!!EXITING exosoft!!"
                sys.exit(s)
                #***************************************************************************************************
        else:
            log.critical("Auto mode and no params provided, so run default stages: SASTMCMC.")
            settings['stages']='SASTMCMC'  
        if gotSigmas==False:
            startSigmas = sigmas
    elif gotSigmas==False:
        #Got params, but no useful sigmas in settings file, so check if there should be and update stages or exit.
        if autoMode==False:
            #check if manual settings make sense
            if (settings['stages'] in ['MCMC']):
                log.critical('MCMC mode requested, but not starting sigmas provided.  Quiting exosoft!!')
                #***************************************************************************************************
                s = "MUST PROVIDE USEFUL STARTSIGMAS IN SIMPLE SETTINGS DICT FOR MCMC MODE, ELSE NOTHING TO START CHAINS WITH."
                s+="\n\nRUN IN AUTO MODE OR PERFORM ST TO GET STARTING VALUES.\n\n!!EXITING exosoft!!"
                sys.exit(s)
                #***************************************************************************************************
        else:
            if type(startParams)!=np.ndarray:
                log.info("Auto mode and no params or sigmas provided, so run default stages: SASTMCMC.")
                settings['stages']='SASTMCMC'
            else:
                log.info("Auto mode and params, but no sigmas provided, so run default stages: STMCMC.")
                settings['stages']='STMCMC'
            log.info("Original startSigmas in settings files were not usable, so setting to default values.")
        startSigmas = sigmas
    else:
        #got params and sigmas.
        log.info("Both startParams and startSigmas passed examination :-D")
        if autoMode:
            #So we can skip all the initialization and go right to the juicy MCMC stage :-D
            settings['stages']='MCMC'
            log.info("Skipping initialization stages and going straight to MCMC.")
    log.info("Logic checks on mode/stages requested, and startParams/startSigmas provided, the resulting stages to run are: "+settings['stages'])
    # clean up sigs of pars that were not varying, and that ary is a ndarray.
    for i in range(0,len(startSigmas)):
        if i not in paramInts:
            startSigmas[i] = 0
    if type(startSigmas)!=np.ndarray:
        startSigmas = np.array(startSigmas)
        
    ##make list of stages to run
    stgLstDict = {'MC':['MC'],'SA':['SA'],'SAST':['SA','ST'],'ST':['ST'],'SASTMCMC':['SA','ST','MCMC'],'STMCMC':['ST','MCMC'],'MCMC':['MCMC']}
    stageList = stgLstDict[settings['stages']]
    ## take care of initialization settings if in autoMode
    if autoMode:
        uSTDdict = {'loose':0.1,'enough':0.05,'tight':0.02}
        settings['maxUstd'] = uSTDdict[settings['initCrit']]
        nSTsampDict = {'loose':10000,'enough':100000,'tight':500000}
        settings['nSTsamp']= nSTsampDict[settings['initCrit']]
        settings['commentsDict']['nSTsamp'] = "Num ST samples"
        
    settings['startParams'] = startParams
    settings['startSigmas'] = startSigmas
    settings['stageList'] = stageList
    
    return settings

#END OF FILE
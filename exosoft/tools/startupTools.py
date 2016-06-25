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
#from IPython.core.prompts import cwd_filt
warnings.simplefilter("error")

log = exoSOFTlogger.getLogger('main.suTools',lvl=100,addFH=False) 

def startup(argv,ExoSOFTdir,rePlot=False):
    """
    Ways to start ExoSOFT:
    If shebang at top of ExoSOFT.py matches your system and location of ExoSOFT/exosoft
    has been added to your PATH, then follow type 'a' below.  Else, use 'b' version.
    
    
    1. Provide full path to settings file.
        a: $ExoSOFT.py /path/you/want/settings.py
        b: cd into /../ExoSOFT/exosoft/ first, then: 
           $python ExoSOFT.py /path/you/want/settings.py
    
    2. From current directory where custom settings files exist.  
       Replace 'settingsfilename' with name of file you want to run.
        a: $ExoSOFT.py settingsfilename.py
        b: $python /../ExoSOFT/exosoft/ExoSOFT.py settingsfilename.py
        
    3. Basic.  User will be prompted to select example 
       ExoSOFT/examples/settings.py, or provide full path to specific settings 
       file.
        a:  $ExoSOFT.py
        b: cd into /../ExoSOFT/exosoft/ first, then: $python ExoSOFT.py
            For both 'a' and 'b':
            To use example: type 'y' and press enter.
            Else, type 'n' and press enter, then type 
            '/path/you/want/settings.py' and press enter.
        
    
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
    settFilePath = getSettFilePath(argv,ExoSOFTdir)
    # extra double check file exists
    if os.path.exists(settFilePath)==False:
        log.critical('Critical error occured while trying to load in the settings.  Quiting exosoft!!')
        #*********************************************************************
        s = "\nSETTINGS FILE NOT FOUND!"
        s+= "\nPLEASE DOUBLE CHECK RULES FOR THE 3 WAYS TO START ExoSOFT WITH "
        s+="INPUT ARGUMENT --help\n\n!!EXITING ExoSOFT!!"
        log.raisemsg(s)
        raise IOError('\n\n'+s)
        #*********************************************************************
    else:
        settings = rwTools.loadSettings(ExoSOFTdir,settFilePath)
        log.setStreamLevel(settings['logLevel'])
        settings['ExoSOFTdir']=ExoSOFTdir
        settings['settingsDir']=os.path.dirname(settFilePath)
        settings['settFilePath']= settFilePath
        ## check if outDir is defined, else to write outputs in cwd
        if settings['outDir']==None:
            print '*'*35+"\n* 'outDir' parameter was None."
            cwd = os.getenv('PWD')
            print '* your current working directory is: '+cwd
            yn = raw_input("* Do you want to write output files here? (y/n): ")
            if (('y' in yn) or ('Y' in yn)):
                settings['outDir'] = cwd
                print '*'*35
            else:
                #*********************************************************************
                s = "\nPlease change the value for the 'outDir' key in your\n"
                s+="settings file to where you want to write the outputs to."
                s+="\n\n!!EXITING ExoSOFT!!"
                log.raisemsg(s)
                raise IOError('\n\n'+s)
                #*********************************************************************
        ## Make a directory (folder) to place all the files from this simulation run
        settings['finalFolder'] = os.path.join(settings['outDir'],settings['outRoot'])
        pklDir = os.path.join(settings['finalFolder'],'pklDir')
        settings['pklDir'] = pklDir
        ##if not doing a re-post analysis with customPost.py
        if rePlot==False:
            if os.path.exists(settings['finalFolder']):
                if settings['logLevel']<50: ## Handle this with a 'clob' bool in dict??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    print '$'*50
                    print '$ WARNING!! the folder:\n$ "'+settings['finalFolder']+'"\n$ ALREADY EXISTS!'
                    print '$ You can overwrite the data in it, or exit this simulation.'
                    yn = raw_input('$ OVERWRITE current folder (y/n): ')
                    if settings['logLevel']<50:
                        print '$'*50+'\n'
                else:
                    yn = 'y'
                if (('y' in yn) or ('Y' in yn)):
                    shutil.rmtree(settings['finalFolder'])
                    os.mkdir(settings['finalFolder'])
                    try:
                        dbDir = os.path.join(settings['dbFolder'],settings['outRoot'])
                        if os.path.exists(dbDir):
                            shutil.rmtree(dbDir)
                    except:
                        log.debug('seems copy to dropbox related settings keys are not there...')
                else: #elif (('n' in YN) or ('N' in YN)):
                    log.raisemsg("")
                    raise IOError('\n\n'+s)
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
            ## make folder for later copying pickle files for recovery if error
            ## or customPost work needs them.
            os.mkdir(pklDir)
        ## push all comments from tuples into a sub dictionary to ensure all  
        ## values for requested keys are just the value with no comments.
        commentsDict = {}
        for key in settings:
            if type(settings[key])==tuple:
                    commentsDict[key] = settings[key][1]
                    settings[key] = settings[key][0]
        settings['commentsDict'] = commentsDict
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
            log.raisemsg(s)
            raise IOError('\n\n'+s)
            #***************************************************************************************************
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
        
        if len(settings['vMINs'])!=len(settings['vMAXs']):
            log.critical("THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!")
            #***************************************************************************************************
            s="THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!\n"
            s+="PLEASE CHECK THE ADVANCED SETTINGS FILES AND FIX THIS"
            s+="!\n\n!!EXITING exosoft!!"
            log.raisemsg(s)
            raise ValueError('\n\n'+s)
            #***************************************************************************************************
        for i in range(0,len(settings['vMINs'])):
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
        ##start with uniform sigma values
        sigmas = np.zeros(rangeMinsRaw.shape)
        for i in paramInts:
            sigmas[i]=(rangeMaxsRaw[i]-rangeMinsRaw[i])*settings['strtSig']
        # Maximum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
        settings['sigMax'] = 1.0
        settings['commentsDict']['sigMax']= 'Max ratio of params range,for step size.'
        ## push all these important parameter related items into the dict for later use.
        settings['realData'] = realData
        settings['rangeMinsRaw'] = rangeMinsRaw
        settings['rangeMaxsRaw'] = rangeMaxsRaw
        settings['rangeMins'] = rangeMins
        settings['rangeMaxs'] = rangeMaxs
        settings['paramInts'] = np.array(paramInts)
        ## map some prior keys to updated values
        incDict = {True:'sin',False:False,None:False,'sin':'sin','cos':'cos'}
        settings['incPrior']=incDict[settings['incPrior']]
        m1Dict = {True:'PDMF',False:False,None:False,'IMF':'IMF','PDMF':'PDMF'}
        settings['M1Prior']=m1Dict[settings['M1Prior']]
        m2Dict = {True:'CMF',False:False,None:False,'IMF':'IMF','PDMF':'PDMF','CMF':'CMF'}
        settings['M2Prior']=m2Dict[settings['M2Prior']]
        ## use modePrep to make sure all is ready for the stages requested
        settings = modePrep(settings,sigmas)
        
        ## extra to be extra careful, but will be able to kill this soon
        if settings['nChains']<settings['nMCMCcns']:
            settings['nChains'] = settings['nMCMCcns']
        
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
                log.raisemsg(s)
                raise IOError('\n\n'+s)
                #***************************************************************************************************
        else:
            if settings['stages']!='SASTMCMC':
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
                log.raisemsg(s)
                raise IOError('\n\n'+s)
                #***************************************************************************************************
        else:
            if type(startParams)!=np.ndarray:
                if settings['stages']!='SASTMCMC':
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

def getSettFilePath(argv,ExoSOFTdir):
    ## Pull in settings filename prepend from command line args, if provided
    inArg = ''
    if len(argv)>1:
        try:
            inArg = argv[1]
        except:
            print '\nWarning: Only able to handle a single argument on the cmd !!\n'    
    if inArg=='--help':
        s ="\nWays to start ExoSOFT\n"+'-'*50
        s+="\nIf shebang at top of ExoSOFT.py matches "
        s+="\nyour system and location of ExoSOFT/exosoft has been added to your"
        s+="\nPATH, then follow type 'a' below.  Else, use 'b' version.\n"+'-'*50
        s+="\n1. Provide full path to settings file.\n      a: $ExoSOFT.py "
        s+="/path/you/want/settings.py\n      b: cd into /../ExoSOFT/exosoft/ "
        s+="first, then: \n      $python ExoSOFT.py /path/you/want/settings.py"
        s+="\n\n2. From current directory where custom settings files exist.\n" 
        s+="   Replace 'settingsfilename' with name of file you want to run.\n"
        s+="      a: $ExoSOFT.py settingsfilename.py\n      b: $python "
        s+="/../ExoSOFT/exosoft/ExoSOFT.py settingsfilename.py\n\n"
        s+="3. Basic.  User will be prompted to select example ExoSOFT/"
        s+="examples/settings.py, \n   or provide full path to specific settings "
        s+="file.\n     a:  $ExoSOFT.py\n     b: cd into /../ExoSOFT/exosoft/ "
        s+="first, then: $python ExoSOFT.py\n        For both 'a' and 'b':\n  "
        s+=" To use example: type 'y' and press enter.\n   Else, type 'n' and "
        s+="press enter, then type '/path/you/want/settings.py' and press "
        s+="enter.\n"
        s= '\n'+'*'*50+'\n'+s+'\n'+'*'*50+'\n'
        log.raisemsg(s)
        raise IOError('\n\n'+s)
            
    else: 
        ## Load up the required specific directory paths in dict
        # Starting method #1, provide full path to settings file          
        if '/' in inArg:
            #assume full path
            settFilePath = inArg
            if os.path.exists(settFilePath)==False:
                log.critical('Critical error occured while trying to load in the settings.  Exiting ExoSOFT!!')
                #*****************************************************************
                s = "\nTHAT SETTINGS FILE DOES NOT EXIST!!"
                s+= "\nARGUMENT PROVIDED TO ExoSOFT WAS: '"+inArg+"'"
                s+= "\nFOR METHOD #1 THE FULL PATH TO SETTINGS FILE IS EXPECTED\n"
                s+= "ie. $python ExoSOFT.py /path/you/want/settings.py"
                s+="\n\n!!EXITING ExoSOFT!!"
                log.raisemsg(s)
                raise IOError('\n\n'+s)
                #*****************************************************************
            else:
                log.debug("Starting ExoSOFT by method #1")
        else:
            gotIt = False
            # Starting method #2, settings file in cwd.
            if inArg!='':
                gotIt=True
                if os.path.exists(inArg):
                    settFilePath=inArg
                    log.debug("Starting ExoSOFT by method #2")
                else:
                    log.critical('Critical error occured while trying to load in the settings.  Exiting ExoSOFT!!')
                    #********************************************************
                    s = "\nTHAT SETTINGS FILE DOES NOT EXIST!!"
                    s+= "\nARGUMENT PROVIDED TO ExoSOFT WAS: '"+inArg+"'"
                    s+= "\nNO SETTINGS FILE WITHT THAT NAME EXISTS IN CWD."
                    s+= "\nFOR METHOD #2 THE SETTINGS FILENAME IS EXPECTED.\n"
                    s+= "ie. $python ExoSOFT.py settingsfilename.py"
                    s+="\n\n!!EXITING ExoSOFT!!"
                    log.raisemsg(s)
                    raise IOError('\n\n'+s)
                    #*********************************************************
            elif os.path.exists('settings.py'):
                print '*'*35
                print '* A settings.py exists in your cwd.'
                yn = raw_input('* Use it? (y/n): ')
                if (('y' in yn) or ('Y' in yn)):
                    gotIt=True
                    settFilePath='settings.py'
                    print '*'*35
                    log.debug("Starting ExoSOFT by method #2")
            if gotIt==False:
                # Starting method #3, prompt user to use defaults or provide a 
                # settings file path.
                print '*'*50
                print '* No settings file provided.  '
                yn = raw_input('* Use the example one? (y/n): ')
                if (('y' in yn) or ('Y' in yn)):
                    examplesDir = os.path.join(ExoSOFTdir.split('exosoft')[0],'examples/')
                    settFilePath = os.path.join(examplesDir,"settings.py")
                    if os.path.exists(settFilePath)==False:
                        log.critical('Critical error occured while trying to load in the settings.  Exiting ExoSOFT!!')
                        s = "\nTHAT SETTINGS FILE DOES NOT EXIST!!"
                        s+= "\nARGUMENT PROVIDED TO ExoSOFT WAS: '"+inArg+"'"
                        s+= "\nMETHOD #3 IN EXAMPLE MODE BROKEN!!!\n"
                        s+= 'PLEASE CHECK THAT DIRECTORY TO SEE IF FILE IS THERE.'
                        s+="\n\n!!EXITING ExoSOFT!!"
                        log.raisemsg(s)
                        raise IOError('\n\n'+s)
                    else:
                        print '*'*50
                        log.debug("Starting ExoSOFT by method #3 with defaults")
                else:
                    #prompt user to provide settings file with path
                    print '* Please provide full path to settings file of your choice.'
                    print '* ex: /path/you/want/settings.py'
                    settFilePath = raw_input('* : ')
                    if os.path.exists(settFilePath)==False:
                        log.critical('Critical error occured while trying to load in the settings.  Exiting ExoSOFT!!')
                        s = "\nTHAT SETTINGS FILE DOES NOT EXIST!!"
                        s+= "\nARGUMENT PROVIDED TO ExoSOFT WAS: '"+inArg+"'"
                        s+= "\nFOR METHOD #3 WITH USER SPECIFIED FILE."
                        s+= "\nPLEASE PROVIDE FULL PATH WHEN PROMPTED\n"
                        s+= "ie. /path/you/want/settings.py\n"
                        s+= 'PLEASE DOUBLE CHECK PATH AND FILE NAME.'
                        s+="\n\n!!EXITING ExoSOFT!!"
                        log.raisemsg(s)
                        raise IOError('\n\n'+s)
                    else:
                        print '*'*50
                        log.debug("Starting ExoSOFT by method #3 with user provided file and path") 
        return settFilePath

#END OF FILE
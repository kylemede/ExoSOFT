#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp  or kylemede@gmail.com
#import sys
#from IPython.core.prompts import cwd_filt
from __future__ import absolute_import
from __future__ import print_function
import shutil
import os
import copy
import numpy as np
import KMlogger
import warnings
import multiprocessing
#from astropy import constants as const
from six.moves import range
from six.moves import input

## import from modules in ExoSOFT ##
#from . import constants as const
from .readWriteTools import load_settings, loadRealData

daysPerYear = 365.2422
#secPerYear = 60*60*24*daysPerYear

warnings.simplefilter("error")
log = KMlogger.getLogger('main.suTools',lvl=100,addFH=False) 

def startup(settings_in,advanced_settings_in,priors_in,rePlot=False):
    """
    Ways to start ExoSOFT:
    If shebang at top of ExoSOFT.py matches your system and location of ExoSOFT/ExoSOFT
    has been added to your PATH, then follow type 'a' below.  Else, use 'b' version.
    
    1. Provide full path to settings file.
        a: $ExoSOFT.py /path/you/want/settings.py
        b: cd into /../ExoSOFT/ExoSOFT/ first, then: 
           $python ExoSOFT.py /path/you/want/settings.py
    
    2. From current directory where custom settings files exist.  
       Replace 'settingsfilename' with name of file you want to run.
        a: $ExoSOFT.py settingsfilename.py
        b: $python /../ExoSOFT/ExoSOFT/ExoSOFT.py settingsfilename.py
        
    3. Basic.  User will be prompted to select example 
       ExoSOFT/examples/settings.py, or provide full path to specific settings 
       file.
        a:  $ExoSOFT.py
        b: cd into /../ExoSOFT/ExoSOFT/ first, then: $python ExoSOFT.py
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
    -Copy all ExoSOFT code into output dir for emergencies.
    -check data exists    
    -push all comments from tuples into a sub dictionary with key 'commentsDict'
    -Check range min and max values, including updates for 'lowecc' mode.
    -find which parameters will be varying.
    -Check if start startParams and startSigmas in dictionary make sense.
    NOTE: have ExoSOFTdir handled with setup.py??
    """    
    #settFilePath = getSettFilePath(sett_file_path)
    # extra double check file exists
    #if os.path.exists(settFilePath)==False:
    if settings_in==None:
        log.critical('Critical error occured while trying to load in the settings.  Quiting ExoSOFT!!')
        #*********************************************************************
        s = "\nSETTINGS FILE NOT FOUND!"
        s+= "\nPLEASE DOUBLE CHECK RULES FOR THE 3 WAYS TO START ExoSOFT WITH "
        s+="INPUT ARGUMENT --help\n\n!!EXITING ExoSOFT!!"
        log.raisemsg(s)
        raise IOError('\n\n'+s)
        #*********************************************************************
    else:
        #settings = loadSettings(ExoSOFTdir,settFilePath)
        #settings = load_settings(settFilePath)
        settings = load_settings(settings_in,advanced_settings_in,priors_in)
        log.setStreamLevel(settings['logLevel'])
        #settings['ExoSOFTdir']=ExoSOFTdir
        #settings['settingsDir']=os.path.dirname(settFilePath)
        #settings['settFilePath']= settFilePath
        ## check if outDir is defined, else to write outputs in cwd
        if settings['outDir']==None:
            print('*'*35+"\n* 'outDir' parameter was None.")
            cwd = os.getenv('PWD')
            print('* your current working directory is: '+cwd)
            yn = input("* Do you want to write output files here? (y/n): ")
            if (('y' in yn) or ('Y' in yn)):
                settings['outDir'] = cwd
                print('*'*35)
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
                    print('$'*50)
                    print('$ WARNING!! the folder:\n$ "'+settings['finalFolder']+'"\n$ ALREADY EXISTS!')
                    print('$ You can overwrite the data in it, or exit this simulation.')
                    yn = input('$ OVERWRITE current folder (y/n): ')
                    if settings['logLevel']<50:
                        print('$'*50+'\n')
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
                    raise IOError("")
            else:
                os.mkdir(settings['finalFolder'])
            if False:
                for key in settings:
                    print(key+' = '+repr(settings[key]))
#             ## copy all of current code to output directory
#             codeCopyDir = os.path.join(settings['finalFolder'],'codeUsed')
#             os.mkdir(codeCopyDir)
#             log.debug('Copying all files in the RESULTS folder over to output folder:\n '+codeCopyDir)
#             setFiles = [settings['settFilePath'],settings['rv_dataFile'],settings['di_dataFile']]
#             copyCodeFiles(settings['ExoSOFTdir'], codeCopyDir,setFiles)
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
        #push in important settings to comments dict needed for later fits file headers
        settings['commentsDict']['nSamples'] = 'number of samples'
        #########################################################################################
        ## Check parameter range settings make sense for data provided and mode of operation.   #
        ## Then load up a list of the parameters to vary during simulation.                     #
        ## Note: Originally this was done in simulator startup, but thought better to move here.#
        #########################################################################################
        realData = loadRealData(diFilename=settings['di_dataFile'],rvFilename=settings['rv_dataFile'],dataMode=settings['data_mode'])
        if (type(realData)!=list)and(type(realData)!=np.ndarray):
            log.critical('Critical error occured while trying to load real data files.  Quiting ExoSOFT!!')
            #***************************************************************************************************
            s = "THERE WAS A PROBLEM LOADING THE REAL DATA."
            s+= "IF 3D OR DI DATA MODE RQUESTED: MAKE SURE A DI FILENAME IN SETTINGS FILES EXISTS."
            s+= "IF 3D OR RV DATA MODE RQUESTED: MAKE SURE A RV FILENAME IN SETTINGS FILES EXISTS."
            s+= "IF THEY EXIST, MAKE SURE THEIR FORMATS MATCH THAT DESCRIBED IN THE readWriteTools.loadDIdata AND loadRVdata FUNCTIONS."
            s+="\n\n!!EXITING ExoSOFT!!"
            log.raisemsg(s)
            raise IOError('\n\n'+s)
            #***************************************************************************************************
        ##check there are matching number of RV datasets and provided min/max vals for offsets
        if np.min(realData[:,6])<1e6:
            numVmins=len(settings['offset_mins'])
            if numVmins==0:
                numVmins=1
            if np.max(realData[:,7])!=(numVmins-1):
                log.error("THE NUMBER OF offset_mins DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
                               "please check the offset_mins/offset_maxs arrays in the simple settings file\n"+\
                               "to make sure they have matching lengths to the number of RV datasets.")
        if  settings['t_max']==settings['t_min']==-1:
            ## set T range to [earliest Epoch-max period,earliest epoch]
            settings['t_max']=np.min(realData[:,0])
            settings['t_min']=np.min(realData[:,0])-settings['p_max']*daysPerYear
        ## check data_mode setting makes sense
        data_modes = ['3D','DI','RV']
        dm = settings['data_mode']
        if type(dm)!=str:
            s = "Setting for 'data_mode' is not valid!!.\nAvailable choices are '3D','DI' and 'RV'."
            s+="\nTrying a default of '3D', but might cause crash..."
            log.critical(s)
            dm = '3D'
        elif dm.upper() not in data_modes:
            s = "Setting for 'data_mode' is not valid!!.\nAvailable choices are '3D','DI' and 'RV'."
            s+="\nTrying a default of '3D', but might cause crash..."
            log.critical(s)
            dm = '3D'
        else:
            dm = dm.upper()
        settings['data_mode'] = dm
        
        ##In DI mode can only find Mtotal, thus push all mass into M1 and kill M2
        if settings['data_mode']=='DI':
            settings['m1_min']=settings['m1_min']+settings['m2_min']
            settings['m1_max']=settings['m1_max']+settings['m2_max']
            settings['m2_min']=0
            settings['m2_max']=0
            log.debug("DI dataMode, so pushed all mass range vals into M1 and set ones for M2 to zero")
        if settings['long_an_min']==None:
            settings['long_an_min'] = 0.0
        if settings['long_an_max']==None:
            if settings['data_mode']=='3D':
                settings['long_an_max'] = 360.0
            else:
                settings['long_an_max'] = 180.0
        if settings['ecc_min']==None:  
            settings['ecc_min'] = 0.0
        if settings['ecc_max']==None:
            settings['ecc_max'] = 0.98
        if settings['inc_min']==None:
            settings['inc_min'] = 0.0
        if settings['inc_max']==None:
            settings['inc_max'] = 180.0  
        if settings['arg_peri_min']==None:
            settings['arg_peri_min'] = 0.0
        if settings['arg_peri_max']==None:
            settings['arg_peri_max'] = 360.0
        ##load up range min,max and sigma arrayS
        rangeMaxs = [settings['m1_max'],\
                     settings['m2_max'],\
                     settings['para_max'],\
                     settings['long_an_max'],\
                     settings['ecc_max'],\
                     settings['t_max'],\
                     settings['p_max'],\
                     settings['inc_max'],\
                     settings['arg_peri_max']]
        rangeMins = [settings['m1_min'],\
                     settings['m2_min'],\
                     settings['para_min'],\
                     settings['long_an_min'],\
                     settings['ecc_min'],\
                     settings['t_min'],\
                     settings['p_min'],\
                     settings['inc_min'],\
                     settings['arg_peri_min']]   
        ## Check all range parameters have a suitable value     
        rangeMaxs_check = [  10,\
                             5,\
                             300,\
                             720,\
                             0.98,\
                             3000000,\
                             3000,\
                             180,\
                             720]
        rangeMins_check = [  0.00001,\
                             0.000000001,\
                             1,\
                             -720,\
                             0,\
                             2000000,\
                             0.00001,\
                             -180,\
                             -720] 
        for i in range(len(rangeMaxs)):
            ## force cast all range values to floats just to be sure
            rangeMaxs[i] = float(rangeMaxs[i])
            rangeMins[i] = float(rangeMins[i])
            rangeMins_check[i] = float(rangeMins_check[i])
            rangeMaxs_check[i] = float(rangeMaxs_check[i])
            if rangeMaxs[i]>rangeMaxs_check[i]:
                s = "Max parameter value was out of range.\n So, it was"
                s+= "changed from "+str(rangeMaxs[i])+" to "+str(rangeMaxs_check[i])
                rangeMaxs[i]=rangeMaxs_check[i]
                log.debug(s)
            if rangeMins[i]<rangeMins_check[i]:
                s = "Min parameter value was out of range.\n So, it was"
                s+= "changed from "+str(rangeMins[i])+" to "+str(rangeMins_check[i])
                rangeMins[i]=rangeMins_check[i]
                log.debug(s)
        ## load in RV instrument offsets and check their values are suitable
        if len(settings['offset_mins'])!=len(settings['offset_maxs']):
            log.critical("THE NUMBER OF offset_mins NOT EQUAL TO NUMBER OF offset_maxs!!!")
            #***************************************************************************************************
            s="THE NUMBER OF offset_mins NOT EQUAL TO NUMBER OF offset_maxs!!!\n"
            s+="PLEASE CHECK THE ADVANCED SETTINGS FILES AND FIX THIS"
            s+="!\n\n!!EXITING ExoSOFT!!"
            log.raisemsg(s)
            raise ValueError('\n\n'+s)
            #***************************************************************************************************
        settings['num_offsets'] = len(settings['offset_maxs'])
        for i in range(0,len(settings['offset_mins'])):
            v_min = settings['offset_mins'][i]
            v_max = settings['offset_maxs'][i]
            if v_min<-50000:
                s = "Min velocity offset parameter #"+str(i)+" value was out of"
                s+= " range.\nSo, it was changed from "+str(v_min)+" to "+str(-50000)
                log.debug(s)
                v_min = -50000
            if v_max>50000:
                s = "Max velocity offset parameter #"+str(i)+" value was out of"
                s+= " range.\nSo, it was changed from "+str(v_min)+" to "+str(50000)
                log.debug(s)
                v_max = 50000
            rangeMins.append(v_min)
            rangeMaxs.append(v_max)
        rangeMaxs = np.array(rangeMaxs)
        rangeMins = np.array(rangeMins)
        ##For low_ecc case, make Raw min/max vals for param drawing during MC mode
        rangeMaxsRaw = copy.deepcopy(rangeMaxs)
        rangeMinsRaw = copy.deepcopy(rangeMins)
    
        if settings['low_ecc']:
            rangeMaxsRaw[8] = rangeMaxs[4]
            rangeMaxsRaw[4] = rangeMaxs[4]
            rangeMinsRaw[8] = (-1.0*rangeMaxs[4])
            rangeMinsRaw[4] = (-1.0*rangeMaxs[4])
            #print(repr(rangeMinsRaw))
            ## Cutting out this loop below as it is running very slow on pip installed version...
            if False:
                ## run through the possible numbers for e and omega to find min/max for RAW versions
                ## Only relevent for MC and the begining jumps of SA,  
                ## Otherwise the values are converted back to omega and e and the original ranges are used for the check.
                fourMin=1e6
                fourMax=-1e6
                eightMin=1e6
                eightMax=-1e6
                for omeg in range(int(rangeMins[9]*10),int(rangeMaxs[9]*10),1):
                    omega = float(omeg)/10.0
                    for e in range(int(rangeMins[4]*100),int(rangeMaxs[4]*100),1):
                        ecc = float(e)/100.0
                        four = np.sqrt(ecc)*np.sin((np.pi/180.0)*omega)
                        eight = np.sqrt(ecc)*np.cos((np.pi/180.0)*omega)
                        if four>fourMax:
                            fourMax = four
                        if four<fourMin:
                            fourMin = four
                        if eight>eightMax:
                            eightMax = eight
                        if eight<eightMin:
                            eightMin = eight
                # check values are within [-e_max,e_max] and push them into RAW arrays
                if eightMax>rangeMaxs[4]:
                    eightMax = rangeMaxs[4]
                if fourMax>rangeMaxs[4]:
                    fourMax = rangeMaxs[4]
                if eightMin<(-1.0*rangeMaxs[4]):
                    eightMin = (-1.0*rangeMaxs[4])
                if fourMin<(-1.0*rangeMaxs[4]):
                    fourMin = (-1.0*rangeMaxs[4])
                rangeMaxsRaw[8] = eightMax
                rangeMaxsRaw[4] = fourMax
                rangeMinsRaw[8] = eightMin
                rangeMinsRaw[4] = fourMin
        ## figure out which parameters are varying in this run.
        ## basically, don't vary m1,m2,parallax if in RV mode
        paramInts = []
        for i in range(0,len(rangeMins)):
            if (i>8):
                # Only add Velocity offsets to varied params if NOT in DI mode
                if settings['data_mode']!='DI':
                    if rangeMaxs[i]!=0:
                        paramInts.append(i)                                         
            elif (i in [2,3]) and (settings['data_mode']!='RV'):
                # Don't vary parallax or Long_an if in RV mode
                if(rangeMaxs[i]!=0):
                    paramInts.append(i)
            elif rangeMaxs[i]!=0:
                paramInts.append(i)
        ##start with uniform sigma values
        sigmas = np.zeros(rangeMinsRaw.shape)
        for i in paramInts:
            sigmas[i]=(rangeMaxsRaw[i]-rangeMinsRaw[i])*settings['strtSig']
        # Maximum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
        settings['sigMax'] = 1.0
        #settings['commentsDict']['sigMax']= 'Max ratio of params range,for step size.'
        ## push all these important parameter related items into the dict for later use.
        settings['realData'] = realData
        settings['rangeMinsRaw'] = rangeMinsRaw
        settings['rangeMaxsRaw'] = rangeMaxsRaw
        settings['range_mins'] = rangeMins
        settings['range_maxs'] = rangeMaxs
        settings['paramInts'] = np.array(paramInts)
        ## map some prior keys to updated values
        eccDict = {True:'beta',False:'uniform',None:'uniform','2e':'2e','st08':'ST08','j08':'J08','rayexp':'RayExp','beta':'beta','uniform':'uniform'}
        ecc_prior = settings['ecc_prior']
        if type(ecc_prior)==str:
            ecc_prior = ecc_prior.upper()
        if ecc_prior not in eccDict:
            ecc_prior = True
        settings['ecc_prior'] = eccDict[ecc_prior]
        incDict = {True:'sin',False:False,None:False,'sin':'sin','cos':'cos'}
        inc_prior = settings['inc_prior']
        if type(inc_prior)==str:
            inc_prior = inc_prior.lower()
        if inc_prior not in incDict:
            inc_prior = True
        settings['inc_prior'] = incDict[inc_prior]
        m1Dict = {True:'PDMF',False:False,None:False,'IMF':'IMF','PDMF':'PDMF'}
        m1_prior = settings['m1_prior']
        if type(m1_prior)==str:
            m1_prior = m1_prior.upper()
        if m1_prior not in m1Dict:
            m1_prior = True
        settings['m1_prior']=m1Dict[m1_prior]
        m2Dict = {True:'CMF',False:False,None:False,'IMF':'IMF','PDMF':'PDMF','CMF':'CMF'}
        m2_prior = settings['m2_prior']
        if type(m2_prior)==str:
            m2_prior = m2_prior.upper()
        if m2_prior not in m2Dict:
            m2_prior = True
        settings['m2_prior']=m2Dict[m2_prior]
        ## for the basic True/False priors of period and parallax
        boolDict = {True:True,False:False,None:False}
        if settings['p_prior'] not in boolDict:
            settings['p_prior'] = True
        settings['p_prior'] = boolDict[settings['p_prior']]
        if settings['para_prior'] not in boolDict:
            settings['para_prior'] = True
        settings['para_prior'] = boolDict[settings['para_prior']]
        
        ### Check all other flag (bool) settings
        bool_setting_strs = ['pltDists','pltOrbit','delChains','delCombined','CalcBurn','rmBurn','calcCL','calcIAC','CalcGR','autoMode','strtMCMCatBest','pasa','vary_tc', 'tc_equal_to', 'Kdirect']
        bool_defaults = [True,      True,       True,       True,           True,       True,   True,   True,     True,   False,        False,            False,   False,      True,           False,]
        for i in range(len(bool_setting_strs)):
            if settings[bool_setting_strs[i]] not in boolDict:
                settings[bool_setting_strs[i]] = bool_defaults[i]
            settings[bool_setting_strs[i]] = boolDict[settings[bool_setting_strs[i]]]
                
        ## Check all other number settings
        num_setting_mins_dict = {'nSamples':1,\
                                 'n_wlkrs':1,\
                                 'logLevel':0,\
                                 'chiMAX':0.0001,\
                                 'chiMaxST':0.0001,\
                                 'cMaxMCMC':0.0001,\
                                 'saveInt':1,\
                                 'n_emcee_burn':0,\
                                 'thin_rate':0,\
                                 'nSumry':0,\
                                 'nGRcalc':1,\
                                 'nSAsamp':1,\
                                 'strtTemp':1,\
                                 'tempInt':1,\
                                 'maxUstd':0.01,\
                                 'nSTsamp':1,\
                                 'strtSig':0.001,\
                                 'sigInt':1,\
                                 'sigMin':0.0001,\
                                 'dmpInt':1,\
                                 'omega_offset_rv':-360,\
                                 'omega_offset_di':-360,\
                                 'KMAX':0,\
                                 'KMIN':0}
        num_setting_maxs_dict = {'nSamples':1e9,\
                                 'n_wlkrs':1e4,\
                                 'logLevel':100,\
                                 'chiMAX':1e5,\
                                 'chiMaxST':1e5,\
                                 'cMaxMCMC':1e5,\
                                 'saveInt':1e9,\
                                 'n_emcee_burn':1e9,\
                                 'thin_rate':1e9,\
                                 'nSumry':100,\
                                 'nGRcalc':1000,\
                                 'nSAsamp':1e8,\
                                 'strtTemp':1e5,\
                                 'tempInt':1e9,\
                                 'maxUstd':100,\
                                 'nSTsamp':1e7,\
                                 'strtSig':0.5,\
                                 'sigInt':1e3,\
                                 'sigMin':0.1,\
                                 'dmpInt':1e9,\
                                 'omega_offset_rv':360,\
                                 'omega_offset_di':360,\
                                 'KMAX':1e7,\
                                 'KMIN':1e7}
        ks = num_setting_maxs_dict.keys()
        for i in range(len(ks)):
            s = "Setting '"+ks[i]+"' value was out of range.\n"
            if settings[ks[i]]<num_setting_mins_dict[ks[i]]:
                s+= "So, it was changed from "+str(settings[ks[i]])+" to "+\
                str(num_setting_mins_dict[ks[i]])
                log.debug(s)
                settings[ks[i]] = num_setting_mins_dict[ks[i]]
            if settings[ks[i]]>num_setting_maxs_dict[ks[i]]:
                s+= "So, it was changed from "+str(settings[ks[i]])+" to "+\
                str(num_setting_maxs_dict[ks[i]])
                log.debug(s)
                settings[ks[i]] = num_setting_maxs_dict[ks[i]]
            
        
        ## use modePrep to make sure all is ready for the stages requested
        settings = modePrep(settings,sigmas)
        
        # load up # of cpus to use 
        ncpu = multiprocessing.cpu_count()
        if type(settings['nCPUs'])==int:
            if settings['nCPUs']<0:
                settings['nChains'] = ncpu-settings['nCPUs']
                settings['nMCMCcns'] = ncpu-settings['nCPUs']
            elif settings['nCPUs']>0:
                settings['nChains'] = settings['nCPUs']
                settings['nMCMCcns'] = settings['nCPUs']
            else:
                settings['nChains'] = ncpu
                settings['nMCMCcns'] = ncpu
        else:
            settings['nChains'] = ncpu
            settings['nMCMCcns'] = ncpu
        
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
    num_params_direct = len(settings['rangeMaxsRaw'])
    #print("len(settings['rangeMaxsRaw']) = "+repr())
    #rangeMaxs = settings['range_maxs']
    autoMode = settings['autoMode']
    num_params_stored = 13+settings['num_offsets']
    
    ##check if startParams in settings file are useful
    gotParams = False
    #print("startParams in settings originally: "+repr(startParams))
    #print('paramInts = '+repr(paramInts))
    if (type(startParams)==list)or(type(startParams)==np.ndarray):
        if type(startParams)==list:
            startParams = np.array(startParams)
        if len(startParams)==num_params_stored:
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
        if len(startSigmas)==num_params_direct:
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
                log.critical('ST or MCMC mode requested, but not starting parameters provided.  Quiting ExoSOFT!!')
                #***************************************************************************************************
                s="MUST PROVIDE USEFUL STARTPARAMS IN SIMPLE SETTINGS DICT FOR ST or MCMC MODE, ELSE NOTHING TO START CHAINS WITH"
                s+="\n\nRUN IN AUTO MODE, OR PERFORM A ROUND OF SA OR SAST TO GET USEFUL VALUES TO STARTING VALUES."
                s+="!\n\n!!EXITING ExoSOFT!!"
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
                log.critical('MCMC mode requested, but not starting sigmas provided.  Quiting ExoSOFT!!')
                #***************************************************************************************************
                s = "MUST PROVIDE USEFUL STARTSIGMAS IN SIMPLE SETTINGS DICT FOR MCMC MODE, ELSE NOTHING TO START CHAINS WITH."
                s+="\n\nRUN IN AUTO MODE OR PERFORM ST TO GET STARTING VALUES.\n\n!!EXITING ExoSOFT!!"
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
    stgLstDict = {'MC':['MC'],'SA':['SA'],'SAST':['SA','ST'],'ST':['ST'],\
                  'SASTMCMC':['SA','ST','MCMC'],'STMCMC':['ST','MCMC'],\
                  'MCMC':['MCMC'],'SAemcee':['SA','emcee'],'emcee':['emcee']}
    if settings['stages'] not in stgLstDict:
        s = "Setting for 'stages' "+str(settings['stages'])+"was not valid."
        s+="\n using the default of 'SASTMCMC'"
        s+="\n Options are:'MC','SA','ST','SAST','SASTMCMC,'MCMC','SAemcee','emcee'"
        settings['stages'] = 'SASTMCMC'
    stageList = stgLstDict[settings['stages']]
    ## take care of initialization settings if in autoMode
    if autoMode:
        uSTDdict = {'loose':0.1,'enough':0.05,'tight':0.02}
        if settings['initCrit'] not in uSTDdict:
            settings['initCrit'] = 'loose'
        settings['maxUstd'] = uSTDdict[settings['initCrit']]
        nSTsampDict = {'loose':10000,'enough':100000,'tight':500000}
        settings['nSTsamp']= nSTsampDict[settings['initCrit']]
        #settings['commentsDict']['nSTsamp'] = "Num ST samples"
    
    ## Check on 'accRates' setting values
    accRates = settings['accRates']
    if type(accRates)!=list:
        accRates = [0.25,0.35]
    else:
        if accRates[0]<0:
            accRates[0] = 0
        if accRates[0]>1:
            accRates[0] = 1
        if accRates[1]<0:
            accRates[1] = 0
        if accRates[1]>1:
            accRates[1] = 1
        
    
    if settings['thin_rate']==None:
        settings['thin_rate'] = 0
    if settings['thin_rate']<=settings['nSamples']/settings['n_wlkrs']:
        if 'emcee' in stageList:
            log.critical("thin_rate was less than nSamples/n_wlkers.  "+\
                         "Thus, thin_rate was set to 0.")
            settings['thin_rate'] = 0
    
    if settings['n_emcee_burn']==None:
        settings['n_emcee_burn'] = 0
    if settings['n_emcee_burn']>=settings['nSamples']/settings['n_wlkrs']:
        if 'emcee' in stageList:
            log.critical("n_emcee_burn was greater than nSamples/n_wlkers.  "+\
                             "Thus, n_emcee_burn was set to 0.")
        settings['n_emcee_burn'] = 0
        
    settings['startParams'] = startParams
    settings['startSigmas'] = startSigmas
    settings['stageList'] = stageList
    
    return settings

def getSettFilePath(sett_file_path):
    ## Pull in settings filename prepend from command line args, if provided
    inArg = ''
    if len(sett_file_path)>1:
        try:
            inArg =sett_file_path
        except:
            print('\nWarning: Only able to handle a single argument on the cmd !!\n')    
    if inArg=='--help':
        s ="\nWays to start ExoSOFT\n"+'-'*50
        s+="\nIf shebang at top of ExoSOFT.py matches "
        s+="\nyour system and location of ExoSOFT/ExoSOFT has been added to your"
        s+="\nPATH, then follow type 'a' below.  Else, use 'b' version.\n"+'-'*50
        s+="\n1. Provide full path to settings file.\n      a: $ExoSOFT.py "
        s+="/path/you/want/settings.py\n      b: cd into /../ExoSOFT/ExoSOFT/ "
        s+="first, then: \n      $python ExoSOFT.py /path/you/want/settings.py"
        s+="\n\n2. From current directory where custom settings files exist.\n" 
        s+="   Replace 'settingsfilename' with name of file you want to run.\n"
        s+="      a: $ExoSOFT.py settingsfilename.py\n      b: $python "
        s+="/../ExoSOFT/ExoSOFT/ExoSOFT.py settingsfilename.py\n\n"
        s+="3. Basic.  User will be prompted to select example ExoSOFT/"
        s+="examples/settings.py, \n   or provide full path to specific settings "
        s+="file.\n     a:  $ExoSOFT.py\n     b: cd into /../ExoSOFT/ExoSOFT/ "
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
                print('*'*35)
                print('* A settings.py exists in your cwd.')
                yn = input('* Use it? (y/n): ')
                if (('y' in yn) or ('Y' in yn)):
                    gotIt=True
                    settFilePath='settings.py'
                    print('*'*35)
                    log.debug("Starting ExoSOFT by method #2")
            if gotIt==False:
                # Starting method #3, prompt user to use defaults or provide a 
                # settings file path.
                print('*'*50)
                print('* No settings file provided.  ')
                yn = input('* Use the example one? (y/n): ')
                if (('y' in yn) or ('Y' in yn)):
                    #examplesDir = os.path.join(ExoSOFTdir.split('ExoSOFT')[0],'examples/')
                    #settFilePath = os.path.join(examplesDir,"settings.py")
                    if os.path.exists(settFilePath)==False:
                        log.critical('Critical error occured while trying to load in the settings.  Exiting ExoSOFT!!')
                        s = "PROVIDED/AUTO SETTINGS FILE WAS: "+settFilePath
                        s+= "\n\nTHAT SETTINGS FILE DOES NOT EXIST!!"
                        s+= "\nARGUMENT PROVIDED TO ExoSOFT WAS: '"+inArg+"'"
                        s+= "\nMETHOD #3 IN EXAMPLE MODE BROKEN!!!\n"
                        s+= 'PLEASE CHECK THAT DIRECTORY TO SEE IF FILE IS THERE.'
                        s+="\n\n!!EXITING ExoSOFT!!"
                        log.raisemsg(s)
                        raise IOError('\n\n'+s)
                    else:
                        print('*'*50)
                        log.debug("Starting ExoSOFT by method #3 with defaults")
                else:
                    #prompt user to provide settings file with path
                    print('* Please provide full path to settings file of your choice.')
                    print('* ex: /path/you/want/settings.py')
                    settFilePath = input('* : ')
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
                        print('*'*50)
                        log.debug("Starting ExoSOFT by method #3 with user provided file and path") 
        return settFilePath

#END OF FILE
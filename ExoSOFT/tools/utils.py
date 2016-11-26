#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp or kylemede@gmail.com
from __future__ import absolute_import
import numpy as np
import os
import shutil
import KMlogger
from six.moves import range

log = KMlogger.getLogger('main.tools',lvl=100,addFH=False)

def load_di_data(filename):
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


def load_rv_data(filename):
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
            #print "line = "+line  #$$$$$$$$$$$$$$$$$
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
                #print 'incrementing datasetNum'  #$$$$$$$$$$$$$$$$$
        epochs_rv = np.array(epochs_rv,dtype=np.dtype('d'))
        rv = np.array(rv,dtype=np.dtype('d'))
        rv_err = np.array(rv_err,dtype=np.dtype('d'))
        rv_inst_num = np.array(rv_inst_num,dtype=np.dtype('i'))
    except:
        log.critical("a problem occured while trying to load RV data.  \nPlease check it is formatted correctly.")
    return (epochs_rv, rv, rv_err, rv_inst_num)

def load_settings_dict(settings_filename):
    """
    A way to load in settings from a settings an ExoSOFT type settings file.
    
    Users can choose to copy and modify the one provided in the 'examples' 
    directory to match the system they are working on, or write their own 
    function like this to load in the required parameters/arrays... that 
    comprise all the inputs for ExoSOFTpriors, ExoSOFTdata and ExoSOFTparams.
    
    NOTE: THE DICTIONARY KEY WORDS USED HERE DO NOT WORK FOR OLD EXOSOFT SETTINGS FILES.  NEED TO CHANGE ALL OF THEM TO THIS FORMAT EVENTUALLY!!!!
    """   
    # First see if try to remove any old settings files.
    try:
        os.remove('./ExoSOFTmodel/temp/settings.py')
    except:
        pass
    cwd = os.getenv('PWD')
    tmp_dir = '.'
    #tmp_dir = os.path.dirname(os.path.abspath(settings_filename))
    tmp_settings_filename = os.path.join(tmp_dir,'__settings.py')
    tmp_init = os.path.join(tmp_dir,"__init__.py")
    kill_init = False
    if os.path.exists(tmp_init)==False:
        kill_init = True
        f = open(tmp_init,'w')
        f.close()
        
    # Copy settings file to temp dir
    shutil.copy(settings_filename,tmp_settings_filename)
    # cd into temp dir, import settings then cd back to pwd
    if tmp_dir!='.':
        os.chdir(tmp_dir)
    from __settings import settings as sd
    if kill_init:
        os.remove(tmp_init)
    os.remove(tmp_settings_filename)
    if tmp_dir!='.':
        os.chdir(cwd)
    
    ## Load in default ranges if any values in settings file were None.
    if sd['data_mode']=='DI':
        sd['m1_min']=sd['m1_min']+sd['m2_min']
        sd['m1_max']=sd['m1_max']+sd['m2_max']
        sd['m2_min']=0
        sd['m2_max']=0
        log.debug("DI dataMode, so pushed all mass range vals into M1 and set ones for M2 to zero")
    if sd['long_an_min']==None:
        sd['long_an_min'] = 0.0
    if sd['long_an_max']==None:
        if sd['data_mode']=='3D':
            sd['long_an_max'] = 360.0
        else:
            sd['long_an_max'] = 180.0
    if sd['ecc_min']==None:  
        sd['ecc_min'] = 0.0
    if sd['ecc_max']==None:
        sd['ecc_max'] = 0.98
    if sd['inc_min']==None:
        sd['inc_min'] = 0.0
    if sd['inc_max']==None:
        sd['inc_max'] = 180.0  
    if sd['arg_peri_min']==None:
        sd['arg_peri_min'] = 0.0
    if sd['arg_peri_max']==None:
        sd['arg_peri_max'] = 360.0
    
    ## load up modified versions of dictionary elements needed by ExoSOFTpriors, ExoSOFTdata and ExoSOFTparams.
    #[m1, m2, parallax, long_an, e, to/tc, p, inc, arg_peri]
    range_mins = [sd['m1_min'], sd['m2_min'], sd['para_min'], sd['long_an_min'], 
                  sd['ecc_min'], sd['t_min'], sd['p_min'], sd['inc_min'], 
                  sd['arg_peri_min']]
    for i in range(len(sd['offset_mins'])):
        range_mins.append(sd['offset_mins'][i])
    range_maxs = [sd['m1_max'], sd['m2_max'], sd['para_max'], sd['long_an_max'], 
                  sd['ecc_max'], sd['t_max'], sd['p_max'], sd['inc_max'], 
                  sd['arg_peri_max']]
    for i in range(len(sd['offset_maxs'])):
        range_maxs.append(sd['offset_maxs'][i])
    
    sd['range_maxs'] = range_maxs
    sd['range_mins'] = range_mins
    sd['num_offsets'] = len(sd['offset_mins'])
    
    # Move comments from any parameter in settings dict to sub dict            
    commentsDict = {}
    for key in sd:
        if type(sd[key])==tuple:
                commentsDict[key] = sd[key][1]
                sd[key] = sd[key][0]
    sd['commentsDict'] = commentsDict
            
    return sd
    
def make_starting_params(pars,n,scale=0.01):
    """
    Creates an array of starting guesses tightly centered on the input 
    parameters.
    """
    # pars: [m1,m2,parallax,long_an, e/sqrte_sinomega,to/tc,p,inc,arg_peri/sqrte_cosomega,v1,v2...]
    starting_params = []
    for i in range(n):
        p = []
        for j in range(len(pars)):
            s = scale
            if pars[j]>1e4:
                # make scale smaller for To/Tc
                s = scale*0.001
            p.append( np.random.normal(pars[j],abs(s*pars[j])) )
        starting_params.append(p)
    starting_params = np.array(starting_params,dtype=np.dtype('d'))
    return starting_params

#EOF
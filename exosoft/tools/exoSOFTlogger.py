#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import logging
import platform
import psutil
import sys
import os
#import traceback

log_dict={}
verbose = False # global variable like this is not recommended, just temporary till I do something better or kill it.

class ExoSOFTlogger(logging.getLoggerClass()):
    """
    This is the advanced logging object used throughout the exosoft.  
    It inherits from the standard 
    Python library 'logging' and provides added features.
    The default log level for the output file will be 1, ie ALL messages;
    while the default for the screen will be INFO, and can be changed easily 
    using the setStreamLevel(lvl) member function.
    """       
    
    def setStreamLevel(self,lvl=20):
        """
        Set/change the level for the stream handler for a logging object.
        Any file handlers will be left alone.
        All messages of a higher severity level than 'lvl' will be printed 
        to the screen.
        
        Args:    
            lvl (int): The severity level of messages printed to the screen with 
                    the stream handler, default = 20.
        
        +---------------------+
        |    Standard Levels  |
        +---------------+-----+
        |    Name       |Level|
        +===============+=====+
        |CRITICAL       |  50 | 
        +---------------+-----+
        |ERROR          |  40 | 
        +---------------+-----+
        |WARNING        |  30 |
        +---------------+-----+
        |INFO           |  20 | 
        +---------------+-----+
        |DEBUG          |  10 |
        +---------------+-----+
        |NOTSET         |  0  |
        +---------------+-----+
        """
        # Kill off the old handlers and reset them with the setHandlers func
        for i in range(0,len(self.handlers)):
            h = self.handlers[i]
            if isinstance(h,logging.StreamHandler):
                self.removeHandler(h)
                break
        addStreamHandler(self,lvl)
        
def getLogger(name='generalLoggerName',dir='',lvl=20,addFH=True,addSH=True,):
    """This will either return the logging object already
    instantiated, or instantiate a new one and return it.  
    **Use this function to both create and return any logger** to avoid 
    accidentally adding additional handlers by using the setUpLogger function 
    instead.
    
    Args:
        name (str): The name for the logging object and 
                    name.log will be the output file written to disk.
        lvl (int): The severity level of messages printed to the screen with 
                    the stream handler, default = 20.
        addFH (boolean): Add a file handler to this logger?  Default severity 
                         level for it will be 1, and it will be named following
                         name+'.log'.  Default = True.
        addSH (boolean): Add a stream handler to this logger? Severity set with 
                        the lvl argument.  Default = True.        
    Returns:
        log (CharisLogger object): A CharisLogger object that was either 
                                  freshly instantiated or determined to 
                                  already exist, then returned.
    """
    log = False
    try:
        log = log_dict[name]
        if verbose:
            print repr(log_dict)
            print 'found a log by the name already exists so returning it'
    except:
        if verbose:
            print 'No logger object found so creating one with the name '+name
        log = setUpLogger(name,dir,lvl,addFH,addSH)
    return log

def setUpLogger(name='generalLoggerName',dir='',lvl=20,addFH=True,addSH=True):
    """ This function is utilized by getLogger to set up a new logging object.
    It will have the default name 'generalLoggerName' and stream handler level
    of 20 unless redefined in the function call.  
    NOTE:
    If a file handler is added, it will have the lowest severity level by 
    default (Currently no need for changing this setting, so it will stay 
    this way for now).  Remember that any messages will be passed up to any 
    parent loggers, so children do not always need their own file handler.
    
    Args:
        name (str): The name for the logging object and 
                    name.log will be the output file written to disk.
        lvl (int): The severity level of messages printed to the screen with 
                    the stream handler, default = 20.
        addFH (boolean): Add a file handler to this logger?  Default severity 
                         level for it will be 1, and it will be named following
                         name+'.log'.  Default = True.
        addSH (boolean): Add a stream handler to this logger? Severity set with 
                        the lvl argument.  Default = True.
    Returns:
        log (SMODTlogger object): A SMODTlogger object that was freshly 
                                   instantiated.
    """
    logging.setLoggerClass(ExoSOFTlogger)
    log = logging.getLogger(name)
    log_dict[name]=log
    log.setLevel(1)
    # add the requested handlers to the log
    if addFH:
        addFileHandler(log,dir,lvl=1)
    # make a stream handler
    if addSH:
        addStreamHandler(log,lvl)
    return log    
        
    
def addFileHandler(log, dir='',lvl=1):
    """
    This function will add a file handler to a log with the provided level.
    
    Args:
        log (CharisLogger object): A ExoSOFTlogger object that was freshly 
                                   instantiated.
        lvl (int): The severity level of messages printed to the file with 
                    the file handler, default = 1.
    """
    if verbose:
        print 'Setting FileHandler level to '+str(lvl)
    fh = logging.FileHandler(os.path.join(dir,log.name+'.log'))
    fh.setLevel(lvl)
    frmtString = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    fFrmt = logging.Formatter(frmtString)
    fh.setFormatter(fFrmt)
    log.addHandler(fh)

def addStreamHandler(log,lvl=20):
    """
    This function will add a stream handler to a log with the provided level.
    
    Args:
        log (CharisLogger object): A ExoSOFTlogger object that was freshly 
                                   instantiated.
        lvl (int): The severity level of messages printed to the screen with 
                    the stream handler, default = 20.
    """
    if verbose:
        print 'Setting StreamHandler level to '+str(lvl)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(lvl)
    sFrmt = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    sh.setFormatter(sFrmt)
    log.addHandler(sh)   
    
def logSystemInfo(log):
    """ A function to be called just after a logging object is instantiated 
    for exosoft to load the log up with info about the computer it is 
    being ran on and the software version.  This function utilizes the 
    psutil and platform libraries, so they must be install for it to work.  
    For clarity of the log, it is suggested to perform immediately after 
    instantiation to put it at the top of the log file.
    
    Args:
        log (Python logging object): logging object to have the system's 
                                    info summarized in.
                                    
    The messages this prints to the log will look like:
    
    | System Information Summary:
    | OS type = Linux
    | OS Version = 3.9.10-100.fc17.x86_64
    | Machine UserName = xxxxxx.astron.s.u-tokyo.ac.jp
    | Machine Processor Type = x86_64
    | Number of cores = 8
    | Total RAM [GB] = 23.5403785706, % used = 15.9
    | Python Version = '2.7.3'

    """
    infoStr = ''
    infoStr+="\n"+"="*11+' System Information Summary '+'='*11
    infoStr+="\n"+'OS type = '+platform.uname()[0]
    infoStr+="\n"+'OS Version = '+platform.uname()[2]
    infoStr+="\n"+'Machine UserName = '+platform.uname()[1]
    infoStr+="\n"+'Machine Processor Type = '+platform.processor()
    infoStr+="\n"+'Number of cores = '+str(psutil.cpu_count())
    totMem = int(round(psutil.virtual_memory()[0]/1073741824.0))
    percentMem = int(round(psutil.virtual_memory()[2]))
    infoStr+="\n"+'Total RAM = '+str(totMem)+'[GB], with ~ '+str(percentMem)+"% already in use at simulation start"
    infoStr+="\n"+'Python Version = '+repr(platform.python_version())
    infoStr+="\n"+'='*50
    log.debug(infoStr)
        
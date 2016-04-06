

### run make for swig if requested
#if settingsDict['remake']:
#    cwd = os.getcwd()
#    log.debug("-"*45+" Starting to remake CPP/SWIG tools "+45*"-")
#    os.chdir(os.path.join(settingsDict['ExoSOFTdir'],'tools/cppTools/'))
#    os.system('make clean')
#    os.system('make')
#    os.chdir(cwd)
#    log.debug("-"*45+" Done re-making CPP/SWIG tools "+45*"-")
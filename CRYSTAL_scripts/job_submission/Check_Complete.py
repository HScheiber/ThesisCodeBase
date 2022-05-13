#!/global/software/python-3.5.0/ThesisCodeBase/python3
import os
import sys
import re
import datetime
    
MainName = sys.argv[1]
logfile = MainName + '.out'

# Check if file log exist
if os.path.isfile(logfile):
    logfound = True
else:
    logfound = False

# Open the output file if found
Error_Found = False
if logfound:
    f = open(logfile,'r')
    logtext = f.read()
    f.close()

    # Check for Abnormal End
    Tempfolder_re = re.compile(r'SCF abnormal end') 
    SCF_Error_matches = re.search(Tempfolder_re,logtext)
    if SCF_Error_matches is not None:
        Error_Found = True
        f = open('ERROR.FOUND', 'a+')
        f.write(str(datetime.datetime.now()) + ': SCF Abnormal End. I/O Error?\n')
        f.close()
    
    # Check for ERROR in log file
    Errorcheck_re = re.compile(r'ERROR')
    error_matches = re.search(Errorcheck_re,logtext)
    if error_matches is not None:
        Error_Found = True
        f = open('ERROR.FOUND', 'a+')
        f.write(str(datetime.datetime.now()) + ': Unspecified Error Detected, check log file.\n')
        f.close()
else:
    Error_Found = True
    f = open('ERROR.FOUND', 'a+')
    f.write(str(datetime.datetime.now()) + ': CRYSTAL Log File Not Found.\n')
    f.close()

if not Error_Found:
    f2 = open('RUN.COMPLETE', 'w+')
    f2.close()

sys.exit(0)
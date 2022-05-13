#!/global/software/python-3.5.0/ThesisCodeBase/python3
import os
import sys
import re
import shutil
import datetime
    
MainName = sys.argv[1]
Restart_WF = bool(sys.argv[2])
logfile = MainName + '.out'
Job_input_file = MainName + '.d12'

# Check for Complete run file (indicates run completed)
if os.path.isfile('RUN.COMPLETE'):
    Restart = '0';
    f = open('Restart.Response', 'w')
    f.write(Restart)
    f.close()
    sys.exit("Previous Run Complete. Do not restart job.")
else:
    Restart = '1';
    f = open('Restart.Response', 'w')
    f.write(Restart)
    f.close()

# Check if file log exist
if not os.path.isfile(logfile):
    logfound = False
    TempFound = False
else:
    logfound = True

# Open the output file if found
if logfound:
    f = open(logfile,'r')
    logtext = f.read()
    f.close()

    # Get temp folder from log file
    tempfolder_re = re.compile(r'(temporary directory )(.+?)(\n)') 
    folder_matches = re.search(tempfolder_re,logtext)

    # Missing folder
    if folder_matches is None:
        TempFound = False
    else:
        TempFound = True
    
    # Check for ERROR in log file
    errorcheck_re = re.compile(r'ERROR')
    error_matches = re.search(errorcheck_re,logtext)
    if error_matches is not None:
        f = open('ERROR.FOUND', 'a+')
        f.write(str(datetime.datetime.now()) + ': Unspecified Error Detected, check log file.\n')
        f.close()


if TempFound:
    Prev_calc_folder = folder_matches.group(2)
    if os.path.isdir(Prev_calc_folder):
        Prev_calc_files = os.listdir(Prev_calc_folder)

        # Copy Restart files over
        f20_copied = False
        for i in range(0, len(Prev_calc_files)):
            Current_file = Prev_calc_files[i]
            src = os.path.join(Prev_calc_folder,Current_file)
            extension = os.path.splitext(Current_file)[1]
            if extension == '.DAT':
                dst_dir= os.curdir
                shutil.copy(src,dst_dir)
            elif Current_file == 'fort.79' and os.path.getsize(src) > 0:
                dst = MainName + '.f20'
                shutil.copy(src, dst)
                f20_copied = True
            elif Current_file == 'fort.9' and os.path.getsize(src) > 0 and not f20_copied:
                src = os.path.join(Prev_calc_folder,Current_file)
                dst = MainName + '.f20'
                shutil.copy(src, dst)
                f20_copied = True
            elif Current_file == 'fort.20' and os.path.getsize(src) > 0 and not f20_copied:
                src = os.path.join(Prev_calc_folder,Current_file)
                dst = MainName + '.f20'
                shutil.copy(src, dst)
                f20_copied = True
            elif Current_file == 'fort.13' and os.path.getsize(src) > 0:
                # Copy and rename it to .f13 file type
                src = os.path.join(Prev_calc_folder,Current_file)
                dst = MainName + '.f13'
                shutil.copy(src, dst)
        # Delete old folder and files
        shutil.rmtree(Prev_calc_folder)

# Check if f20 file exists in folder
Current_files = os.listdir(os.curdir)
f20_exist = False
for i in range(0,len(Current_files)):
    FileName = Current_files[i]
    if FileName == (MainName + '.f20') and os.path.getsize(FileName) > 0:
        f20_exist = True
        break
    elif FileName == (MainName + '.f9') and os.path.getsize(FileName) > 0:
        f20_exist = True
        break

 # Get input file text
f3 = open(Job_input_file,'r')
Inputtxt = f3.read()

# Add restart of wavefunction using GUESSP if the fort.20 file was found and Restart_WF active
checkreWF = re.search(r'GUESSP',Inputtxt)
if f20_exist and (checkreWF is None) and Restart_WF:
    Inputtxt = Inputtxt.replace('99 0\nEND', '99 0\nEND\nGUESSP')
    
elif not f20_exist or not Restart_WF: # If no fort.20 file or restartWF disabled, dont use GUESSP
    Inputtxt = Inputtxt.replace('GUESSP\n', '')

# Figure out type of job restart: QHA, frequency, or optimization   
optgeom_re = re.compile(r'OPTGEOM') 
optgeom_matches = re.search(optgeom_re,Inputtxt)

freqcalc_re = re.compile(r'FREQCALC')
freqcalc_matches = re.search(freqcalc_re,Inputtxt)

QHA_re = re.compile(r'QHA')
QHA_matches = re.search(QHA_re,Inputtxt)

# Search for optimization
if optgeom_matches is not None:
    # If found, look for OPTINFO.DAT file (required for restart)
    Current_files = os.listdir(os.curdir)
    OPTINFO_Found = False
    for i in range(0,len(Current_files)):
        FileName = Current_files[i]
        if FileName == 'OPTINFO.DAT':
            OPTINFO_Found = True
            break
        elif FileName == (MainName + '.optinfo'):
            os.rename(MainName + '.optinfo', 'OPTINFO.DAT')
            OPTINFO_Found = True
            break
    
    # If optinfo.dat file is found, add in RESTART to CRYSTAL input in optgeom section if not present
    checkrestart_re = re.compile(r'OPTGEOM' + r'\n' + r'RESTART')
    checkrestart = re.search(checkrestart_re,Inputtxt)
    
    if OPTINFO_Found and checkrestart is None:
        Inputtxt_out = Inputtxt.replace('OPTGEOM', 'OPTGEOM' + '\n' + 'RESTART')
    else:
        Inputtxt_out = Inputtxt
# Search for QHA analysis
elif QHA_matches is not None:
    # If found, look for EOSINFO.DAT
    Current_files = os.listdir(os.curdir)
    EOSINFO_Found = False
    
    for i in range(0,len(Current_files)):
        FileName = Current_files[i]
        if FileName == 'EOSINFO.DAT':
            EOSINFO_Found = True
            break
        elif FileName == (MainName + '.eosinfo'):
            os.rename(MainName + '.eosinfo', 'EOSINFO.DAT')
            EOSINFO_Found = True
            break
    
    checkrestart_re = re.compile(r'QHA' + r'\n' + r'RESTART')
    checkrestart = re.search(checkrestart_re,Inputtxt)
    if EOSINFO_Found and checkrestart is None:
        Inputtxt_out = Inputtxt.replace('QHA','QHA' + '\n' + 'RESTART')
    else:
        Inputtxt_out = Inputtxt
# Search for vibrational analysis
elif freqcalc_matches is not None:
    # If found, look for FREQINFO.DAT
    Current_files = os.listdir(os.curdir)
    FREQINFO_Found = False
    
    for i in range(0,len(Current_files)):
        FileName = Current_files[i]
        if FileName == 'FREQINFO.DAT':
            FREQINFO_Found = True
            break
        elif FileName == (MainName + '.freqinfo'):
            os.rename(MainName + '.freqinfo', 'FREQINFO.DAT')
            FREQINFO_Found = True
            break
    
    checkrestart_re = re.compile(r'FREQCALC' + r'\n' + r'RESTART')
    checkrestart = re.search(checkrestart_re,Inputtxt)
    if FREQINFO_Found and checkrestart is None:
        Inputtxt_out = Inputtxt.replace('FREQCALC','FREQCALC' + '\n' + 'RESTART')
    else:
        Inputtxt_out = Inputtxt
else:
    Inputtxt_out = Inputtxt

# Re-Save input file
Inputtxt_out = re.sub(r'(RESTART\n)+', r'RESTART\n', Inputtxt_out)
f5 = open(MainName + '.d12', 'w');
f5.write(Inputtxt_out)
f5.close()
sys.exit()

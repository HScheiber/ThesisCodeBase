function Restart_CRYSTAL(MainName)

logfile = [MainName '.out'];

% Open the output file
fid = fopen(logfile,'rt');
logtext = fread(fid);
fclose(fid);
logtext = char(logtext.');

if isempty(logtext)
    error(['Log text from ' logfile.name ' is empty.'])
end

% Get temp folder
folder_matches = regexp(logtext,'(temporary directory )(.+?)(\n)','tokens','ONCE');

if isempty(folder_matches)
    error(['Unable to find temporary folder in ' logfile.name])
end

Prev_calc_folder = folder_matches{2};

Prev_calc_files = dir(Prev_calc_folder);

% If no old directory, run completed
if isempty(Prev_calc_files)
    % Exit bash script without running crystal
    Restart = false;
    fid = fopen('matlab_output.txt', 'w');
    fprintf(fid, '%d', Restart);
    fclose(fid);
    return
end

% Otherwise Copy files over
f20_copied = false;
for i = 1:length(Prev_calc_files)
    Current_file = Prev_calc_files(i).name;
    if ~isempty(regexp(Current_file,'\.DAT','ONCE'))
        copyfile([Prev_calc_folder filesep Current_file],pwd)
    elseif ~isempty(regexp(Current_file,'fort\.20','ONCE'))
        % Copy and rename it to .f20 file type
        copyfile([Prev_calc_folder filesep Current_file],pwd)
        movefile(Current_file,[MainName '.f20'])
        f20_copied = true;
    elseif ~isempty(regexp(Current_file,'fort\.13','ONCE'))
        % Copy and rename it to .f13 file type
        copyfile([Prev_calc_folder filesep Current_file],pwd)
        movefile(Current_file,[MainName '.f13'])
    end
end

% Get input file
Job_input_file = [MainName '.d12'];
fid = fopen(Job_input_file,'rt');
Inputtxt = fread(fid);
fclose(fid);
Inputtxt = char(Inputtxt.');

% Add restart of wavefunction using GUESSP if the fort.20 file was found
if f20_copied && isempty(regexp(Inputtxt,'GUESSP','ONCE'))
	Inputtxt = strrep(Inputtxt,['99 0' newline 'END'],['99 0' newline 'END' newline 'GUESSP']);
elseif ~f20_copied % If no fort.20 file found, dont use GUESSP
    Inputtxt = strrep(Inputtxt,['GUESSP' newline],'');
end

% Figure out type of job restart: frequency or optimization
% Search for optimization
if ~isempty(regexp(logtext,'OPTGEOM','ONCE'))
    
    % If found, look for OPTINFO.DAT file (required for restart)
    Current_files = dir(pwd);
    OPTINFO_Found = false;
    for i = 1:size(Current_files,1)
        FileName = Current_files(i).name;
        if strcmp(FileName,'OPTINFO.DAT')
            OPTINFO_Found = true;
            break
        elseif strcmpi(FileName,[MainName '.optinfo'])
            copyfile([MainName '.optinfo'],'OPTINFO.DAT')
            OPTINFO_Found = true;
            break
        end
    end
    
    % If optinfo.dat file is found, add in RESTART to CRYSTAL input in
    % optgeom section unless already added
    if OPTINFO_Found && isempty(regexp(Inputtxt,['OPTGEOM' newline 'RESTART'],'ONCE'))
        Inputtxt_out = strrep(Inputtxt,'OPTGEOM',['OPTGEOM' newline 'RESTART']);
    else
        Inputtxt_out = Inputtxt;
    end
    
% Search for vibrational analysis
elseif isempty(regexp(logtext,'FREQCALC','ONCE'))
    % If found, look for FREQINFO.DAT, fort.13, fort.20, fort.28, fort.80 files (required for restart)
    Current_files = dir(pwd);
    FREQINFO_Found = false;
    for i = 1:size(Current_files,1)
        FileName = Current_files(i).name;
        if strcmp(FileName,'FREQINFO.DAT')
            FREQINFO_Found = true;
            break
        elseif strcmpi(FileName,[MainName '.freqinfo'])
            copyfile([MainName '.freqinfo'],'FREQINFO.DAT')
            FREQINFO_Found = true;
            break
        end
    end
    
    % All necessary files for restart are present: set up restart if not
    % already set up
    if FREQINFO_Found && isempty(regexp(Inputtxt,['FREQCALC' newline 'RESTART'],'ONCE'))
        Inputtxt_out = strrep(Inputtxt,'FREQCALC',['FREQCALC' newline 'RESTART']);
    else
        Inputtxt_out = Inputtxt;
    end
else
    Inputtxt_out = Inputtxt;
end

% Signal to restart job
Restart = true;
fid = fopen('matlab_output.txt', 'w');
fprintf(fid, '%d', Restart);
fclose(fid);

% Save input file
delete([MainName '.d12'])
fid = fopen([MainName '.d12'], 'w');
fprintf(fid, '%s', Inputtxt_out);
fclose(fid);

end
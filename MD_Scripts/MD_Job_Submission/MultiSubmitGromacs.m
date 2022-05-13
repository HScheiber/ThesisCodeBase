% Script to submit N Jobs that link together, given an input script.
% Inputs: MultiSubmitGromacs(N,[Directory, Filename, CPTFilename])

% N is the number of links to create.

% Directory is the directory which you would like your job to run in

% Filename is the name of the batch submission script

% CPTfilename is the name of the cpt/checkpoint file (optional) which was given as a 
% checkpoint from a previous job

% If not using optional inputs, the following are defaults:

% Directory will be the present working directory (i.e. the directory which
% the script is run in)

% Filename will be found in the directory as the most recent file with
% extenstion of .sh or .subm

% CPTFilename will be found in the directory as the most recent file with
% extension of .cpt. If no such file exists, then mdrun will start without
% a checkpoint (i.e. from initial inputs such as gro, mdp, top, etc.)

% Must have a template submission script present or directed to by
% Filename. The script copies this template submission and links N jobs
% with the same settings.

function MultiSubmitGromacs(N,varargin)

% newline does not exist as a built in function for older versions of matlab
% if ~exist('newline','builtin')
%     newline = char(10); %#ok<CHARTEN>
% end

if nargin > 1
    Directory = varargin{1};
else
    Directory = pwd;
end

% Get the filename of the submission template script (*.sh or *.subm)
if nargin > 2
    SubmFilename = fullfile(Directory,varargin{2});
else
    % Find all files in directory
    Files = dir(Directory);

    % Find only files with .sh or .subm extention
    Index1 = arrayfun(@(a) regexp(a,'\.sh$|\.subm$'),{Files.name});

    Index2 = ~cellfun(@isempty,Index1);

    Files_Of_Interest = Files(Index2);
    sFOI = size(Files_Of_Interest,1);
    
    if sFOI > 1
        % If multiple submission files exist, choose most recent
        Times = zeros(1,sFOI);
        for y = 1:sFOI
            Times(y) = datenum(Files_Of_Interest(y).date,...
                'dd-mmm-yyyy HH:MM:SS');
        end
        [~,Idx] = max(Times);
        Files_Of_Interest = Files_Of_Interest(Idx);
    elseif sFOI == 0
        error('Unable to find submission file. Use .sh or .subm file extention.');
    end
    
    SubmFilename = fullfile(Directory,Files_Of_Interest.name);
end
[SubmFilename_Dir,SubmFilename_Name,SubmFilename_Ext] = fileparts(SubmFilename);

% Get checkpoint files
if nargin > 3
    CPTFilename = fullfile(Directory,varargin{3});
    if isempty(CPTFilename)
        CPT_Present = false;
        CPTFilename = fullfile(Directory,'state.cpt');
    else
        CPT_Present = true;
    end
else
    % Find only files with .cpt extention
    Index1 = arrayfun(@(a) regexp(a,'\.cpt$'),{Files.name});

    Index2 = ~cellfun(@isempty,Index1);

    Files_Of_Interest = Files(Index2);
    sFOI = size(Files_Of_Interest,1);
    
    if sFOI > 1
        % If multiple cpt files exist, choose most recent
        Times = zeros(1,sFOI);
        for y = 1:sFOI
            Times(y) = datenum(Files_Of_Interest(y).date,...
                'dd-mmm-yyyy HH:MM:SS');
        end
        [~,Idx] = max(Times);
        Files_Of_Interest = Files_Of_Interest(Idx);
    end
    
    if sFOI == 0
        CPT_Present = false;
        CPTFilename = fullfile(Directory,[SubmFilename_Name '.cpt']);
    else
        CPT_Present = true;
        CPTFilename = fullfile(Directory,Files_Of_Interest.name);
    end
end

% Determine scheduler type
svr = strtrim(getenv('HOSTNAME'));
if isempty(svr)
    svr = 'HOSTNAME NOT FOUND';
    Server = 'HOSTNAME NOT FOUND';
else
    Server = svr(1:3);
end

% SLURM Scheduler
if strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'gra') ...
        || strcmpi(Server,'bel')
    Slurm = true;
    disp('SLURM scheduler detected.')
% PBS Scheduler
elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') ...
        || strcmpi(Server,'sea') || strcmpi(Server,'pod')
    Slurm = false;
    disp('PBS scheduler detected.')
else
    disp(['This script does not know the scheduler type of the current host: ' svr])
    disp('Please choose from the following options:')
    prompt = ['Input 0 for PBS scheduler, or 1 for SLURM scheduler. [0]:' newline];
    str = '';
    while ~strcmpi(str,'0') && ~strcmpi(str,'1')
        str = input(prompt,'s');
        if isempty(str)
            str = 0;
            break
        elseif ~strcmpi(str,'0') && ~strcmpi(str,'1')
            disp('Invalid answer.')
        end
    end
    Slurm = strcmpi(str,'1');
    if Slurm
        disp('SLURM scheduler selected.')
    else
        disp('PBS scheduler selected.')
    end
end

% Read submission script
FileData = fileread(SubmFilename);

if Slurm
    % Grab Job Name and put placeholder for job name
    JobMainName = regexp(FileData,'#SBATCH --job-name=(.+?)\n','tokens','ONCE');
    if isempty(JobMainName)
        JobMainName = 'job';
    else
        JobMainName = strtrim(JobMainName{1});
        FileData = regexprep(FileData,'(#SBATCH --job-name=)(.+?)(\n)','$1##JOBNAME##$3');
    end
    
    % Grab error and standard output
    Stdo = regexp(FileData,'#SBATCH --output=(.+?)\n','tokens','ONCE');
    if isempty(Stdo)
        Stdo_Dir = Directory;
        Stdo_Name = SubmFilename_Name;
        Stdo_Ext = '.stdo';
    else
        Stdo = strtrim(Stdo{1});
        FileData = regexprep(FileData,'(#SBATCH --output=)(.+?)(\n)','$1##STDO##$3');
        [Stdo_Dir,Stdo_Name,Stdo_Ext] = fileparts(Stdo);
    end
    Stde = regexp(FileData,'#SBATCH --error=(.+?)\n','tokens','ONCE');
    if isempty(Stde)
        Stde_Dir = Directory;
        Stde_Name = SubmFilename_Name;
        Stde_Ext = '.stde';
    else
        Stde = strtrim(Stde{1});
        FileData = regexprep(FileData,'(#SBATCH --error=)(.+?)(\n)','$1##STDE##$3');
        [Stde_Dir,Stde_Name,Stde_Ext] = fileparts(Stde);
    end
else
    % Grab Job Name and put placeholder for job name
    JobMainName = regexp(FileData,'#PBS -N (.+?)\n','tokens','ONCE');
    if isempty(JobMainName)
        % default job name if none given
        JobMainName = 'job';
    else
        JobMainName = strtrim(JobMainName{1});
        FileData = regexprep(FileData,'(#PBS -N )(.+?)(\n)','$1##JOBNAME##$3');
    end
    
    % Grab error and standard output
    Stdo = regexp(FileData,'#PBS -o (.+?)\n','tokens','ONCE');
    if isempty(Stdo)
        Stdo_Dir = Directory;
        Stdo_Name = SubmFilename_Name;
        Stdo_Ext = '.stdo';
    else
        Stdo = strtrim(Stdo{1});
        FileData = regexprep(FileData,'(#PBS -o )(.+?)(\n)','$1##STDO##$3');
        [Stdo_Dir,Stdo_Name,Stdo_Ext] = fileparts(Stdo);
    end
    Stde = regexp(FileData,'#PBS -e (.+?)\n','tokens','ONCE');
    if isempty(Stde)
        Stde_Dir = Directory;
        Stde_Name = SubmFilename_Name;
        Stde_Ext = '.stde';
    else
        Stde = strtrim(Stde{1});
        FileData = regexprep(FileData,'(#PBS -e )(.+?)(\n)','$1##STDE##$3');
        [Stde_Dir,Stde_Name,Stde_Ext] = fileparts(Stde);
    end
end

% File parts of checkpoint filename
[CPT_Dir,CPT_Name,CPT_Ext] = fileparts(CPTFilename);

% Grab the gromacs mdrun command
gmdrun = regexp(FileData,'[^\n]*mdrun[^\n]*','match','once');

if isempty(gmdrun)
    error('GROMACS mdrun command not found in the batch script.')
end

% Check if there is a max time yet
if isempty(regexp(gmdrun,'-maxh','match','once'))
    if Slurm
        maxhour = regexp(FileData,'#SBATCH --time=([0-9|:]+)','tokens','once');
    else
        maxhour = regexp(FileData,'walltime=([0-9|:]+)','tokens','once');
    end
    maxhour = hours(duration(maxhour{1}));
    gmdrun_edit = [gmdrun ' -maxh ' num2str(maxhour)];
else
    gmdrun_edit = gmdrun;
end

% Check if there is already a cpt and/or cpo flag, remove them if present
if ~isempty(regexp(gmdrun_edit,' -cpi ','match','once'))
    cptline = regexp(gmdrun,'-cpi [^\s]+ *','match','once');
    gmdrun_edit = strrep(gmdrun_edit,cptline,'');
end

if ~isempty(regexp(gmdrun_edit,' -cpo ','match','once'))
    cpoline = regexp(gmdrun,'-cpo [^\s]+ *','match','once');
    gmdrun_edit = strrep(gmdrun_edit,cpoline,'');
end

% Check for the name of the final output gro file
if ~isempty(regexp(gmdrun_edit,' -c ','match','once'))
    outGeom = regexp(gmdrun,'-c [^\s]+ *','match','once');
else
    outGeom = 'confout.gro'; % gromacs default
end

for idx = 1:N
    Index = ['-' num2str(idx,'%03.0f')];
    
    JobName = [JobMainName Index];
    Stdo_i = fullfile(Stdo_Dir,[Stdo_Name Index Stdo_Ext]);
    Stde_i = fullfile(Stde_Dir,[Stde_Name Index Stde_Ext]);
    Filename_i = fullfile(SubmFilename_Dir,[SubmFilename_Name Index SubmFilename_Ext]);    
    cpt_output_i =  fullfile(CPT_Dir,[CPT_Name Index CPT_Ext]);
    
    % Add in cpt output name
    gmdrun_edit_i = [gmdrun_edit ' -cpo ' cpt_output_i];
    
    % Replace Strings
    JobData = strrep(FileData,'##JOBNAME##',JobName);
    JobData = strrep(JobData,'##STDO##',Stdo_i);
    JobData = strrep(JobData,'##STDE##',Stde_i);
    
    if idx == 1
        if CPT_Present
            gmdrun_edit_i = [gmdrun_edit_i ' -cpi ' CPTFilename]; %#ok<*AGROW>
        end
        if Slurm
            qsub_cmd = ['sbatch ' Filename_i];
        else
            qsub_cmd = ['qsub ' Filename_i];
        end
    else
        gmdrun_edit_i = [gmdrun_edit_i ' -cpi ' cpt_output_prev];
        if Slurm
            qsub_cmd = ['sbatch --depend=afterany:' PrevJobID ' ' Filename_i ];
        else
            qsub_cmd = ['qsub -W depend=afterany:' PrevJobID ' ' Filename_i ];
        end
        
        % Add in BASH condition to only run if final output file is not present
        % AND a cpt file from the previous run is detected
        Cond1 = ['if [[ ! -f "' outGeom '" && -f "' cpt_output_prev '" ]]; then'];
        Cond2 = 'fi';
        gmdrun_edit_i = [Cond1 newline '  ' gmdrun_edit_i newline Cond2];
    end
    
    % Replace Jobdata gmx mdrun command with edited one containing
    % dependancy
    JobData = strrep(JobData,gmdrun,gmdrun_edit_i);
    
    % Save the submission script to file
    fid = fopen(Filename_i,'wt');
    fwrite(fid,regexprep(JobData,'\r',''));
    fclose(fid);
    
    disp(['Submitting Job ' JobName])
    [errcode,output] = system(qsub_cmd);

    if errcode ~= 0
        disp(output);
        error(['Error submitting Job: ' newline qsub_cmd]);
    end
    PrevJobID = regexp(output,'[0-9]+','match','ONCE');
    cpt_output_prev = cpt_output_i;
end

disp('All Jobs Submitted Sucessfully!');
end
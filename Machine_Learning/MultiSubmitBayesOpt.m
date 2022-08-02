% Script to submit N Jobs that link together, given an input script.
% Inputs: MultiSubmitBayesOpt(N,[Directory, Filename])

% N is the number of links to create.

% Directory is the directory which you would like your job to run in

% Filename is the name of the batch submission script

% If not using optional inputs, the following are defaults:

% Directory will be the present working directory (i.e. the directory which
% the script is run in)

% Filename will be found in the directory as the most recent file with
% extenstion of .sh or .subm

% Must have a template submission script present or directed to by
% Filename. The script copies this template submission and links N jobs
% with the same settings.

function MultiSubmitBayesOpt(N,varargin)

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
        || strcmpi(Server,'bel') || strcmpi(Server,'nar')
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

% Grab the matlab command
command = regexp(FileData,'[^\n]*matlab( -nodisplay)* -r[^\n]*','match','once');

if isempty(command)
    error('matlab run command command not found in the batch script, must start with "matlab -r".')
end

for idx = 1:N
    Index = ['-' num2str(idx,'%03.0f')];
    
    JobName = [JobMainName Index];
    Stdo_i = fullfile(Stdo_Dir,[Stdo_Name Index Stdo_Ext]);
    Stde_i = fullfile(Stde_Dir,[Stde_Name Index Stde_Ext]);
    Filename_i = fullfile(SubmFilename_Dir,[SubmFilename_Name Index SubmFilename_Ext]);    
    
    % Replace Strings
    JobData = strrep(FileData,'##JOBNAME##',JobName);
    JobData = strrep(JobData,'##STDO##',Stdo_i);
    JobData = strrep(JobData,'##STDE##',Stde_i);
    
    if idx == 1
        if Slurm
            qsub_cmd = ['sbatch ' Filename_i];
        else
            qsub_cmd = ['qsub ' Filename_i];
        end
        command_i = command;
    else
        if Slurm
            qsub_cmd = ['sbatch --depend=afterany:' PrevJobID ' ' Filename_i ];
        else
            qsub_cmd = ['qsub -W depend=afterany:' PrevJobID ' ' Filename_i ];
        end
        
        % Add in BASH condition to only run if final output file is not present
        % AND a cpt file from the previous run is detected
        Cond1 = 'if ! ls *fullopt.mat 1> /dev/null 2>&1  && ls *opt.mat 1> /dev/null 2>&1; then';
        Cond2 = 'fi';
        command_i = [Cond1 newline '  ' command newline Cond2];
    end
    
    % Replace Jobdata gmx mdrun command with edited one containing
    % dependancy
    JobData = strrep(JobData,command,command_i);
    
    % Save the submission script to file
    fid = fopen(Filename_i,'wt');
    fwrite(fid,regexprep(JobData,'\r',''));
    fclose(fid);
    
    disp(['Submitting Job ' JobName])
    [errcode,output] = system(qsub_cmd);


    if errcode ~= 0 && ~Slurm
        disp(output);
        error(['Error submitting Job: ' newline qsub_cmd]);
    end
    
    % Attempt to recover from failed job submissions on SLURM servers
    while errcode ~= 0
        disp(output);
        disp(['Error submitting Job: ' newline qsub_cmd]);
        disp('Retrying.');
        
        errcode2 = 1;
        while errcode2 ~= 0
            pause(5);
            
            % Check if the job was actually submit or not
            check_cmd = 'squeue -u $USER -o "%10i %60j"';
            [errcode2,output] = system(check_cmd);
        end
        
        % Check if the previous job is in the queue
        regcheck = ['([0-9]+) +' JobName];
        jobid = regexp(output,regcheck,'once','tokens');
        
        if isempty(jobid) % Job did not submit successfully
            % Retry submitting
            disp(['Submitting Job ' JobName])
            [errcode,output] = system(qsub_cmd);
            
        else % Job did submit successfully (but failed to output properly)
            output = jobid{1};
            errcode = 0;
        end
    end

    PrevJobID = regexp(output,'[0-9]+','match','ONCE');
    pause(0.5);
end

disp('All Jobs Submitted Sucessfully!');
end
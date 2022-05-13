%%%%% Melting_Point_Batch_Submit_Sockeye %%%%%%
% NOTE THIS DOES NOT WORK YET - Requires MATLAB version 2021 and gromacs,
% which is not really compatible with sockeye
%% Global Calculation settings
clear;
Calculation.N_Calc = 2; % Number of chained calculations
Calculation.Hours = 3; % Max time for each job (hours)
Calculation.Mins = 0; % Max time for job (minutes)
Calculation.Nodes = 1; % Minimum number of cores to request for calculation.
Calculation.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Calculation.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Calculation.SinglePrecision = false; % choose true for single precision mode, false for double
Calculation.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Calculation.openMP = false; % for Sockeye only

% Load the batch script template and remove unnecessary fields
[Batch_Template,~,~,mdrun] = MD_Batch_Template(Calculation);
Batch_Template = strrep(Batch_Template,['##PREMIN##' newline],'');
Batch_Template = strrep(Batch_Template,['##CLEANUP##' newline],'');
Batch_Template = strrep(Batch_Template,mdrun,'##MDRUN##');
Batch_Template = regexprep(Batch_Template,' -maxh [0-9]+','','once');
Batch_Template = strrep(Batch_Template,['##EXT1##' newline],'');
Batch_Template = strrep(Batch_Template,['##EXT2##' newline],'');


% Set up calculations below
skip_calculations = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs

% Some shared calculation settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; %[5 3 2] Domain decomposition
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Shared_Settings.MeltFreezeThreshold = 0.15; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Shared_Settings.Optimizer = 'MPSearcher';
Shared_Settings.lb = 500; % K, lower bound on MP search
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.InitialMeshSize = 10; % Initial step size
Shared_Settings.BracketThreshold = 2; % [K] Sets the target bracket for the melting point
Shared_Settings.MinStepSize = 0.5; % [K] Sets the minimum step size for MPsearcher algorithm
Shared_Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Shared_Settings.Liquid_Fraction = 0.51;
Shared_Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Shared_Settings.MaxCheckTime = 1000; % ps. Max time for melting/freezing runs

% Initial calculation index
idx=0;

%% Set12: Repeated Measurements [6x6x7 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
Salt = 'NaCl';
Theory = 'JC';
T0 = 1073.2 + 216.8;
XY_Size = 6; % nm
Z_Size = 7; % nm
Reps = 1;
for kdx = 1:length(Reps)
    Rep = num2str(Reps(kdx));
    
    rng(Reps(kdx))
    T0_i = T0 + (rand-1/2)*2*20; % add a random number between -20 and +20 to the initial temperature
    
    idx = idx+1;
    Settings_array(idx) = Shared_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = 'Rocksalt'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = ['Set12_Rep_' Rep]; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Manual_Box = true; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
    Settings_array(idx).Manual_Box_a = XY_Size; % Box length a in nm
    Settings_array(idx).Manual_Box_b = XY_Size; % Box length b in nm
    Settings_array(idx).Manual_Box_c = Z_Size; % Box length c in nm
    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
    Settings_array(idx).MeltFreezeThreshold = 0.25;
end

%% Check for already running jobs
if check_running && isunix
    [~,~,~,slurm,~,~,~,~,~,~] = find_home;
    if slurm
    
        er = 1;
        jdx = 1;
        while er ~= 0
            [er,out] = system('squeue -u $USER -o %100j | tail -n +2');
            jdx = jdx+1;
            if jdx > 10
                disp('Failed to check current running jobs.')
                out = '';
                break
            end
        end

        % Grab the currently running jobs
        out = strtrim(out);
        c_jobs = strtrim(split(out,newline));

        % Remove the job number tags
        c_jobs_unique = unique(regexprep(c_jobs,'-([0-9]){3}',''));
    else % PBS without squeue
        er = 1;
        jdx = 1;
        while er ~= 0
            [er,out] = system('qstat -u $USER');
        end
        
        if isempty(out)
            % No jobs running
            c_jobs_unique = {};
        else
            er = 1;
            jdx = 1;
            while er ~= 0
                [er,out] = system('qstatf () { qstat -f $*| tr ''\n'' ''&'' | sed -e ''s/&\t//g'' | tr ''&'' ''\n''; }; qstatfi() { grep -e "$*" | sed -e "s/^\(\S*\).*,$*=\([^,]*\).*/\1 \2/g" ;}; qstatf | qstatfi PBS_O_WORKDIR | grep $USER');
                jdx = jdx+1;
                if jdx > 10
                    disp('Failed to check current running jobs.')
                    out = '';
                    break
                end
            end
            c_jobs = strtrim(splitlines(strtrim(out)));
            
            % Remove the outer directory and find unique jobs
            c_jobs_unique = strrep(unique(c_jobs),fullfile(Shared_Settings.project,Shared_Settings.Project_Directory_Name,filesep),'');
            
            % Replace file separator with '_'
            c_jobs_unique = strrep(c_jobs_unique,filesep,'_');
        end
    end
else
    c_jobs_unique = {''};
end

%% Prepare batch scripts
[~,~,~,Slurm] = find_home;
calc_cmd = 'matlab -nodisplay -r "Find_Melting_Point(##INPUTFILE##)" >> ##LOGFILE##';

% Loop through all models and build their submission file
for idx = 1:length(Settings_array)
    if any(idx == skip_calculations)
        continue
    end
    Settings = Settings_array(idx);
    
    [WorkDir,JobName,Full_Model_Name] = GetMDWorkdir(Settings);
    Batch_Template_idx = strrep(Batch_Template,'##DIRECTORY##',WorkDir);
    TaskName = [Settings.Salt '_' JobName];
    
    % Generate and move to the submission directory
    if ~exist(WorkDir, 'dir')
       mkdir(WorkDir)
    end
    cd(WorkDir);
    
    % Check if job is already complete
    ResultsFile = fullfile(WorkDir,[JobName '_MPResults.mat']);
    if check_complete && isfile(ResultsFile)
        disp([TaskName ': Job already completed. Skipping Job Submission.'])
        continue
    end
    
    % Check if Model is already running
    if check_running && ismember(TaskName,c_jobs_unique) 
        disp([TaskName ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    % Save input file
    inp_file = fullfile(WorkDir,[JobName '.inp']);
    save(inp_file,'Settings','-mat')
    
    % Create input command
    calc_cmd_idx = strrep(calc_cmd,'##INPUTFILE##',['''' JobName '.inp' '''']);
    calc_cmd_idx = strrep(calc_cmd_idx,'##LOGFILE##',[JobName '.MPlog']);
    
    % Set up job links
    for jdx = 1:Calculation.N_Calc
        
        if jdx == 1
            calc_cmd_idx_jdx = calc_cmd_idx;
        else
            EXT1 = ['if [[ ! -f "' ResultsFile '" ]]; then'];
            EXT2 = 'fi';
            calc_cmd_idx_jdx = [EXT1 newline calc_cmd_idx newline EXT2];
        end
        Batch_Template_idx_idx = strrep(Batch_Template_idx,'##MDRUN##',calc_cmd_idx_jdx);
        
        % Copy over batch template and modify
        TaskName_jdx = [TaskName '-' num2str(jdx,'%03.f')];
        Batch_Template_idx_jdx = strrep(Batch_Template_idx_idx,'##TASKNAME##',TaskName_jdx);
        Batch_Template_idx_jdx = strrep(Batch_Template_idx_jdx,'##ERROR##',TaskName_jdx);

        % Save job submission script
        subm_file_jdx = fullfile(WorkDir,[JobName '-' num2str(jdx,'%03.f') '.subm']);
        fid = fopen(subm_file_jdx,'wt');
        fwrite(fid,regexprep(Batch_Template_idx_jdx,{'\r'}',{''}));
        fclose(fid);

        % Submit job
        if jdx == 1
            if Slurm
                qsub_cmd = ['sbatch ' subm_file_jdx];
            else
                qsub_cmd = ['qsub ' subm_file_jdx];
            end
        else
            if Slurm
                qsub_cmd = ['sbatch --depend=afterany:' PrevJobID ' ' subm_file_jdx ];
            else
                qsub_cmd = ['qsub -W depend=afterany:' PrevJobID ' ' subm_file_jdx ];
            end
        end
        
        disp(['Submitting Job ' TaskName_jdx])
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
            regcheck = ['([0-9]+) +' TaskName_jdx];
            jobid = regexp(output,regcheck,'once','tokens');

            if isempty(jobid) % Job did not submit successfully
                % Retry submitting
                disp(['Submitting Job ' TaskName_jdx])
                [errcode,output] = system(qsub_cmd);

            else % Job did submit successfully (but failed to output properly)
                output = jobid{1};
                errcode = 0;
            end
        end
        
        PrevJobID = regexp(output,'[0-9]+','match','ONCE');
        pause(0.5);
    end
end
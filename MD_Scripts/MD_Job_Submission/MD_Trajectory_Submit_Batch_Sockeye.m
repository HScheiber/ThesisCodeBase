%%%%% MD_Trajectory_Submit_Batch_Sockeye %%%%%%

%% Global Calculation settings
clear;
% Load Settings object and set some shared calculation settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.BatchMode = true; % Sets up batch job when true, or runs immediately when false
Shared_Settings.JobSettings.N_Calc = 1; % Number of jobs to link together.
Shared_Settings.JobSettings.Hours = 3; % Max time for each job (hours)
Shared_Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Shared_Settings.JobSettings.Nodes = 1; % Minimum number of cores to request for calculation.
Shared_Settings.JobSettings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.JobSettings.MPI_Ranks = 32; % MPI ranks PER NODE!!!
Shared_Settings.JobSettings.OMP_Threads = 1; % OMP threads per MPI rank
Shared_Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Shared_Settings.MinMDP.nsteps_min = 1000;
Shared_Settings.TimePerFrame = 20; % post-processing time per frame check in ps

% Set up calculations below
skip_calculations = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs. NOTE: Does not work on sockeye

% Initial calculation index
idx=0;

%% Test simulation 1
idx = idx+1;
T0 = 1081; % Initial temperature
Settings_array(idx) = Shared_Settings;
Settings_array(idx).Theory = 'TF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = 'NaCl';
Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'TestSmall'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = []; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = []; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 5; % ns
Settings_array(idx).PreEquilibration = 10; % ps

% Thermostat options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Barostat options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = true; % Use this to raise or lower the barostat compressibility from the default experimental isothermal compressibility

% cutoff settings
Settings_array(idx).MDP.RVDW_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.9; % nm
Settings_array(idx).MDP.Disp_Correction = false;


%% Check for already running jobs
if check_running && isunix
    [~,~,~,slurm,~,~,~,~,~,~] = find_home;
    if slurm
        er = 1;
        jdx = 1;
        while er ~= 0
            [er,out] = system('squeue -u $USER -o %100Z | tail -n +2');
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

        % Remove the outer directory and find unique jobs
        c_jobs_unique = strrep(unique(c_jobs),fullfile(Shared_Settings.project,Shared_Settings.Project_Directory_Name,filesep),'');

        % Replace file separator with '_'
        c_jobs_unique = strrep(c_jobs_unique,filesep,'_');
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
    c_jobs_unique = {};
end

%% Prepare batch scripts
% Loop through all models and build their submission file
for idx = 1:length(Settings_array)
    if any(idx == skip_calculations)
        continue
    end
    Settings = Settings_array(idx);
    
    [WorkDir,JobName,~] = GetMDWorkdir(Settings);
    TaskName = [Settings.Salt '_' JobName];
    
    % Check if job is already complete (including postprocessing)
    OutConfFile = fullfile(WorkDir,[JobName '_OutConf.' Settings.CoordType]);
    OutConfFile2 = strrep(OutConfFile,'scratch','project');
    PostProcessFlagFile = fullfile(WorkDir,'POSTPROCESS_COMPLETE');
    PostProcessFlagFile2 = strrep(PostProcessFlagFile,'scratch','project');
    if check_complete && ( isfile(OutConfFile) || isfile(OutConfFile2) )
        Settings.JobSettings.N_Calc = 1;
        if (Settings.RunPostProcessor && (isfile(PostProcessFlagFile) || isfile(PostProcessFlagFile2)) ) || ~Settings.RunPostProcessor
            disp([TaskName ': Job already completed. Skipping Job Submission.'])
            continue
        end
    end
    
    % Check if Model is already running
    if check_running && ismember(TaskName,c_jobs_unique) 
        disp([TaskName ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    % Finish job setup
    Setup_LiX_Simulation(Settings)
end
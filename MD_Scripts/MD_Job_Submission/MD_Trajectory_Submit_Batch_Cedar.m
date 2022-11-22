%%%%% MD_Trajectory_Submit_Batch_Cedar %%%%%%

%% Global Calculation settings
clear;

% Set up calculations below
skip_calculations = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs. NOTE: Does not work on sockeye

% Load Settings object and set some shared calculation settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.N_Calc = 4; % Number of jobs to link together.
Shared_Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.BatchMode = true; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Hours = 3; % Max time for each job (hours)
Shared_Settings.Mins = 0; % Max time for job (minutes)
Shared_Settings.Nodes = 1; % Minimum number of cores to request for calculation.
Shared_Settings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.MPI_Ranks = -1;
Shared_Settings.OMP_Threads = 1;
Shared_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Shared_Settings.MinMDP.nsteps_min = 1000;
Shared_Settings.TimePerFrame = 1; % post-processing time per frame check in ps
Shared_Settings.Output_Coords = 1000; % output coords every 1 ps


% Initial calculation index
idx=0;

Experiment = Load_Experimental_Data;

% Example production runs for paper
%% LiI/CLJ polarized Alexandria model at 600 K starting in wurtzite crystal structure, 
% NPT, anisotropic barostat. Anneal structure slowly to 0 kelvin.
Salt = 'LiI';
Structure = 'Wurtzite';
Theory = 'JC';
vdW_Type = 'LJ_12-6';
T0 = 600;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'Anneal'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

% load the model
Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
pset = Initialize_MD_Settings;
pset.Salt = Salt;
[JC_MX,JC_MM,JC_XX] = JC_Potential_Parameters(pset);
Settings_array(idx).S.S.MM = Settings_array(idx).S.S.MM/JC_MM.sigma;
Settings_array(idx).S.S.XX = Settings_array(idx).S.S.XX/JC_XX.sigma;
Settings_array(idx).S.S.MX = Settings_array(idx).S.S.MX/JC_MX.sigma;
Settings_array(idx).S.E.MM = Settings_array(idx).S.E.MM/JC_MM.epsilon;
Settings_array(idx).S.E.XX = Settings_array(idx).S.E.XX/JC_XX.epsilon;
Settings_array(idx).S.E.MX = Settings_array(idx).S.E.MX/JC_MX.epsilon;
Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

Settings_array(idx).PreEquilibration = 0; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = [0  2000 2010]; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = [T0 0    0]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 2.01; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'anisotropic';
Settings_array(idx).Target_P = [1 1 1 0 0 0]; % Bar
Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = ones(1,6); % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = false;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

%% LiI/WBK polarized Alexandria model at 909 K starting in liquid crystal structure, Beredensen barostat (tc = 1 ps)
Salt = 'LiI';
Structure = 'Liquid';
Theory = 'BF';
vdW_Type = 'WBK';
T0 = 909;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'Berendsen_1ps'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

% load the model
Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

Settings_array(idx).PreEquilibration = 15; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).MDP.Trajectory_Time = 1.0; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = false;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

%% LiI/WBK polarized Alexandria model at 909 K starting in liquid crystal structure, Beredensen barostat (tc = 10 ps)
% NPT, isotropic barostat. Anneal structure slowly to 0 kelvin.
Salt = 'LiI';
Structure = 'Liquid';
Theory = 'BF';
vdW_Type = 'WBK';
T0 = 909;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'Berendsen_10ps'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

% load the model
Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

Settings_array(idx).PreEquilibration = 15; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).MDP.Trajectory_Time = 1.0; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = false;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

%% LiI/WBK polarized Alexandria model at 909 K starting in liquid crystal structure, PR barostat (tc = 1 ps)
% NPT, isotropic barostat. Anneal structure slowly to 0 kelvin.
Salt = 'LiI';
Structure = 'Liquid';
Theory = 'BF';
vdW_Type = 'WBK';
T0 = 909;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'PR_1ps'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

% load the model
Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

Settings_array(idx).PreEquilibration = 100; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).MDP.Trajectory_Time = 1.0; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = false;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

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
        Settings.N_Calc = 1;
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
    Setup_LiX_Simulation(Settings);
end
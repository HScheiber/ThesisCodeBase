%%%%% MD_Trajectory_Submit_Batch_Cedar %%%%%%

%% Global Calculation settings
clear;

% Set up calculations below
skip_calculations = [5:100];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs. NOTE: Does not work on sockeye

% Load Settings object and set some shared calculation settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.JobSettings.N_Calc = 8; % Number of jobs to link together.
Shared_Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.BatchMode = true; % Sets up batch job when true, or runs immediately when false
Shared_Settings.JobSettings.Hours = 3; % Max time for each job (hours)
Shared_Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Shared_Settings.JobSettings.Nodes = 1; % Minimum number of cores to request for calculation.
Shared_Settings.JobSettings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.MPI_Ranks = -1;
Shared_Settings.OMP_Threads = 1;
Shared_Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.JobSettings.npme = []; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = []; % Domain decomposition
Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Shared_Settings.MinMDP.nsteps_min = 2000;
Shared_Settings.TimePerFrame = 20; % post-processing time per frame check in ps

% Initial calculation index
idx=0;

Experiment = Load_Experimental_Data;

%% Example production runs for paper
% Run 1: NaCl nucleation with JC model
Salt = 'NaCl';
Structure = 'Liquid';
Theory = 'JC';
T0 = 1200;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).Output_Coords = 1000; % output coords every 1 ps
Settings_array(idx).JobSettings.N_Calc = 10; % Number of chained calculations
Settings_array(idx).JobSettings.Hours = 3; % Max time for each job (hours)
Settings_array(idx).JobSettings.MPI_Ranks = 4; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings_array(idx).JobSettings.OMP_Threads = 8; % Set the number of OMP threads per MPI rank

Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'ExampleNuc'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

Settings_array(idx).PreEquilibration = 100; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = [0   20000   21000]; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = [T0  T0-1000 T0-1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 21; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = true;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Run 2: LiI nucleation with TF model
Salt = 'LiI';
Structure = 'Liquid';
Theory = 'TF';
T0 = 1400;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).Output_Coords = 1000; % output coords every 1 ps
Settings_array(idx).JobSettings.N_Calc = 14; % Number of chained calculations
Settings_array(idx).JobSettings.Hours = 3; % Max time for each job (hours)
Settings_array(idx).JobSettings.MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings_array(idx).JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank

Settings_array(idx).MDP.RVDW_Cutoff = 1.4; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.4; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.4; % nm
Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'ExampleNuc'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

Settings_array(idx).PreEquilibration = 100; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = [0   20000   21000]; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = [T0  T0-1000 T0-1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 21; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = true;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Run 3: LiBr nucleation with TF model
Salt = 'LiBr';
Structure = 'Liquid';
Theory = 'TF';
T0 = 1400;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).Output_Coords = 1000; % output coords every 1 ps
Settings_array(idx).JobSettings.N_Calc = 14; % Number of chained calculations
Settings_array(idx).JobSettings.Hours = 3; % Max time for each job (hours)
Settings_array(idx).JobSettings.MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings_array(idx).JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank

Settings_array(idx).MDP.RVDW_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.9; % nm

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'ExampleNuc'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

Settings_array(idx).PreEquilibration = 100; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = [0   20000   21000]; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = [T0  T0-1000 T0-1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 21; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = true;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Run 4: LiBr nucleation with TF-like model
Salt = 'CsCl';
Structure = 'Liquid';
Theory = 'TF';
T0 = 1400;

idx = idx+1;
Settings_array(idx) = Shared_Settings;

Settings_array(idx).Output_Coords = 1000; % output coords every 1 ps
Settings_array(idx).JobSettings.N_Calc = 14; % Number of chained calculations
Settings_array(idx).JobSettings.Hours = 3; % Max time for each job (hours)
Settings_array(idx).JobSettings.MPI_Ranks = 32; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings_array(idx).JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank

Settings_array(idx).MDP.RVDW_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RCoulomb_Cutoff = 1.9; % nm
Settings_array(idx).MDP.RList_Cutoff = 1.9; % nm

Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings_array(idx).JobID = 'ExampleNuc'; % An ID that is tacked onto the folder name of all current jobs
Settings_array(idx).N_atoms = 20000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings_array(idx).c_over_a = 1;

Settings_array(idx).PreEquilibration = 100; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Settings_array(idx).Annealing_Times = [0   20000   21000]; % [ps] A list with the number of annealing reference/control points used
Settings_array(idx).Annealing_Temps = [T0  T0-1000 T0-1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Settings_array(idx).MDP.Trajectory_Time = 21; % ns

% Barostat Options
Settings_array(idx).Isotropy = 'isotropic';
Settings_array(idx).Target_P = 1; % Bar
Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings_array(idx).UseMoltenCompressibility = true;

% Thermostat Options
Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities

%% Generating reference structure for liquid and solid near melting point
Salt = 'NaCl';
Structures = {'Rocksalt' 'Liquid'};
Theories = {'JC' 'TF'};
Temps = {[1073.2 1289] [1073.2 1081]}; % At both experimental MP and model MP
for jdx = 1:length(Structures)
    Structure = Structures{jdx};
    
    for kdx = 1:length(Theories)
        Theory = Theories{kdx};
        Temp_vec = Temps{kdx};
        
        for tdx = 1:length(Temp_vec)
            T0 = Temp_vec(tdx);
            
            %% Thermal expansion for solid starting at 0 K and heating to 2200 K at 100 K / ns
            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = ['StrucRefT' num2str(T0,'%.0f')]; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).N_atoms = 8000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
            Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
            Settings_array(idx).T0 = T0; % K, Initial temperature
            Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
            Settings_array(idx).MDP.Trajectory_Time = 5; % ns
            Settings_array(idx).Thermal_Solid = true;
            Settings_array(idx).Output_Coords = 1000; % output coords every 1 ps
        end
    end
end

%% [Not Yet assigned] Liquid near the melting point, extracting out inherent structures
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Models = {{'JC' ''} {'TF' ''} {'JC' 'EQ'} {'BH' 'FE'}};
Nums_JC = {'2' '2' '1' '5' ''};
Nums_BH = {'5' '6' '7' '5' ''};

for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    Num_JC = Nums_JC{sidx};
    Num_BH = Nums_BH{sidx};
    
    for kdx = 1:length(Models)

        Theory = Models{kdx}{1};
        Model = Models{kdx}{2};
        
        MP = Load_Model_MP(Salt,Theory,Model);
        if ~isempty(Model) && strcmp(Theory,'JC') && isempty(Num_JC)
            continue
        elseif ~isempty(Model) && strcmp(Theory,'BH') && isempty(Num_BH)
            continue
        elseif ~isempty(Model) && strcmp(Theory,'JC')
            Model = [Model Num_JC];
        elseif ~isempty(Model) && strcmp(Theory,'BH')
            Model = [Model Num_BH];
        end
    
        %% Generating inherent structures near the melting point
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = Model; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = 'InhStrMP'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = MP; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = MP; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = MP; % K, Initial temperature
        Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
        Settings_array(idx).Annealing_Times = [0  999.999 1000 31000]; % [ps] A list with the number of annealing reference/control points used
        Settings_array(idx).Annealing_Temps = [MP MP      0    2000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
        Settings_array(idx).MDP.Trajectory_Time = 31; % ns
        Settings_array(idx).Cutoff_Buffer = 1.25;
        if strcmp(Salt,'LiF')
            Settings_array(idx).JobSettings.N_Calc = Settings_array(idx).JobSettings.N_Calc + 4;
        end
    end
end

%% [Not Yet assigned] Inherent structures for CsCl
Salts = {'CsCl'};
Structures = {'Liquid'}; %   'Rocksalt' 'Wurtzite' 'BetaBeO' 'NiAs'  'Sphalerite' 'FiveFive' 'AntiNiAs' 'CsCl' 
Models = {'JC' 'TF'};
MP = 918.2;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    for kdx = 1:length(Models)

        Theory = Models{kdx};
        
        %% Generating inherent structures near the melting point
        idx = idx+1;
        Settings_array(idx) = Shared_Settings;
        Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
        Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
        Settings_array(idx).JobID = 'IntStrMP'; % An ID that is tacked onto the folder name of all current jobs
        Settings_array(idx).N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
        Settings_array(idx).Target_T = MP; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
        Settings_array(idx).MDP.Initial_T = MP; % Initial termpature at which to generate velocities
        Settings_array(idx).T0 = MP; % K, Initial temperature
        Settings_array(idx).Annealing = 'single'; % Options: 'no' 'single' 'periodic'
        Settings_array(idx).Annealing_Times = [0  999.999 1000 31000]; % [ps] A list with the number of annealing reference/control points used
        Settings_array(idx).Annealing_Temps = [MP MP      0    2000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
        Settings_array(idx).MDP.Trajectory_Time = 31; % ns
        Settings_array(idx).Cutoff_Buffer = 1.25;
        if strcmp(Salt,'LiF')
            Settings_array(idx).JobSettings.N_Calc = Settings_array(idx).JobSettings.N_Calc + 4;
        end
    end
end

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
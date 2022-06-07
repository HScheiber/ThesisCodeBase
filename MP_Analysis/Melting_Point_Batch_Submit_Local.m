%%%%% Melting_Point_Batch_Submit_Cedar %%%%%%

%% Global Calculation settings
clear;

% Set up calculations below
skip_calculations = [];

% Load shared resource and mdrun settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 4; % Set the number of OMP threads per MPI rank

% Shared calculation parameters
Shared_Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Shared_Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Shared_Settings.Optimizer = 'MPSearcher';
Shared_Settings.lb = 0; % K, lower bound on MP search
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.BracketThreshold = 1; % [K] Sets the target bracket for the melting point
Shared_Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Shared_Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Shared_Settings.Liquid_Fraction = 0.50;
Shared_Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Shared_Settings.MaxWarn = 2;

% Initial calculation index
idx=0;

%% Prod2 Redos: Melting points of all alkali halides with "default" parameters + Rc = 1.4 + dispersion correction ON
% Structure: rocksalt
% N = 11000 atoms
% Cutoff 1.4 nm with dispersion correction added
% The compressibility is the same in all spatial directions
% Additional 1ps relaxation added to start of simulation

% LiF JC FiveFive: Tm = 1336.8
% LiBr TF Rocksalt: Tm = 802.9
% NaCl TF Rocksalt: Tm = 1081.4
% CsCl TF Rocksalt: Tm = 1179.9

Experiment = Load_Experimental_Data;
JobStuff = {{'LiF' 'JC' 'FiveFive' 1336.8} ...
            {'LiBr' 'TF' 'Rocksalt' 802.9} ...
            {'NaCl' 'TF' 'Rocksalt' 1081.4} ...
            {'CsCl' 'TF' 'Rocksalt' 1179.9}};
        
for jdx = 1:length(JobStuff)
    Salt = JobStuff{jdx}{1};
    Theory = JobStuff{jdx}{2};
    Structure = JobStuff{jdx}{3};
    T0 = JobStuff{jdx}{4};
    
    Theory_Settings = Shared_Settings;
    switch Theory
        case 'TF'
            Theory_Settings.JobSettings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
            Theory_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
            Theory_Settings.MDP.RVDW_Cutoff = 1.9; % nm
            Theory_Settings.MDP.RCoulomb_Cutoff = 1.9; % nm
            Theory_Settings.MDP.RList_Cutoff = 1.9; % nm
            Theory_Settings.Cutoff_Buffer = 1.05; % Not using verlet-buffer-tolerance

        otherwise
            Theory_Settings.JobSettings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
            Theory_Settings.JobSettings.OMP_Threads = 4; % Set the number of OMP threads per MPI rank
            Theory_Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
            Theory_Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
            Theory_Settings.MDP.RVDW_Cutoff = 1.4; % nm
            Theory_Settings.MDP.RCoulomb_Cutoff = 1.4; % nm
            Theory_Settings.MDP.RList_Cutoff = 1.4; % nm
    end

    idx = idx+1;
    Settings_array(idx) = Theory_Settings;
    Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).Structure = Structure; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Settings_array(idx).InitialMeshSize = 5; % K
    Settings_array(idx).Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings_array(idx).JobID = 'Prod2'; % An ID that is tacked onto the folder name of all current jobs
    Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0; % K, Initial temperature
    if contained_in_cell(Structure,{'Rocksalt' 'Sphalerite' 'CsCl'})
        Settings_array(idx).Isotropy = 'semiisotropic';
        Settings_array(idx).Target_P = [1 1]; % Bar
    else
        Settings_array(idx).Isotropy = 'anisotropic';
        Settings_array(idx).Target_P = [1 1 1 1 1 1]; % Bar
    end
    Settings_array(idx).c_over_a = 2;
    Settings_array(idx).N_atoms = 11000;

    Settings_array(idx).BracketThreshold = 1; % K
    Settings_array(idx).MinStepSize = 0.25;
    Settings_array(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs

    Settings_array(idx).Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).Equilibrate_Liquid = 5; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
    Settings_array(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.

    % Barostat Options
    Settings_array(idx).Isotropy = 'semiisotropic';
    Settings_array(idx).Target_P = [1 1]; % Bar
    Settings_array(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings_array(idx).ScaleCompressibility = 1;

    % Thermostat Options
    Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings_array(idx).Target_T = T0_i; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings_array(idx).MDP.Initial_T = T0_i; % Initial termpature at which to generate velocities
    Settings_array(idx).T0 = T0_i; % K, Initial temperature
end


% Loop through all models and build their submission file
for idx = 1:length(Settings_array)
    if any(idx == skip_calculations)
        continue
    end
    Settings = Settings_array(idx);
    Find_Melting_Point(Settings);
end
%% Global Calculation settings
clear;
[home,project,computer,slurm] = find_home;

%% Load Settings object and set some shared calculation settings
Settings = Initialize_MD_Settings;
Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Settings.Nodes = 1; % Minimum number of cores to request for calculation.
Settings.Cores = 8; % Minimum number of cores to request for calculation. Set to -1 for entire node
Settings.MPI_Ranks = 8;
Settings.OMP_Threads = 1;
Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Settings.SinglePrecision = false; % choose true for single precision mode, false for double
Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Settings.MinMDP.nsteps_min = 1000;
Settings.Output_Coords = 1000; % output coords every 1 ps
Settings.MinMDP.Parallel_Min = false;

%% NaCl/WBK polarized vs unpolarized Alexandria model at Tm starting in rocksalt structure: NVT
Settings.Theory = 'BF'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings.Salt = 'NaCl';
Settings.Structure = 'Rocksalt';
a_RS = 5.857263144; % Angstoms
T0 = 1075.168;
Settings.RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};

Settings.MDP.RVDW_Cutoff = 1.0; % nm
Settings.MDP.RCoulomb_Cutoff = 1.1; % nm
Settings.MDP.RList_Cutoff = 1.1; % nm
Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

Settings.Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings.JobID = 'ChkNVT'; % An ID that is tacked onto the folder name of all current jobs
Settings.N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings.c_over_a = 1;
Settings.Geometry = Default_Crystal(Settings,'Center_Coordinates',true);
Settings.Geometry.a = a_RS;
Settings.Geometry.b = a_RS;
Settings.Geometry.c = a_RS;
Settings.Skip_Minimization = true;

% load the model
Settings = Alexandria_Potential_Parameters(Settings,'vdW_Type','WBK');
Settings.GaussianCharge = true; % Turn on Gaussian distributed charges when true
Settings.Polarization = false; % Turn on polarizible Drude model when true

Settings.PreEquilibration = 25; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings.Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings.MDP.Trajectory_Time = 1.0; % ns

% Barostat Options
Settings.Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

% Thermostat Options
Settings.Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings.Nsttcouple = Get_nstcouple(Settings.Time_Constant_T,Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings.Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings.MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Postprocessing
Settings.SavePredictionsImage = true;
Settings.ML_TimeLength = 0;
Settings.ML_TimeStep = 0;
Settings.TimePerFrame = 1; % ps
Settings.SaveFeatures = false; % Save structure fraction vs time image when true for each temperature check
Settings.SavePredictions = false; % Save structure fraction vs time image when true for each temperature check
Settings.Qlm_Average = true;

Setup_LiX_Simulation(Settings)
MD_Postprocessor(Settings)
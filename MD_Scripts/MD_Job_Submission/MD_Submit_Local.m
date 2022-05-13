clear;
% Load Settings object and set some shared calculation settings
T0 = 1000; % Initial temperature

Shared_Settings = Initialize_MD_Settings;
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.JobSettings.MPI_Ranks = 2;
Shared_Settings.JobSettings.OMP_Threads = 4;
Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder

Shared_Settings.Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Shared_Settings.Salt = 'LiBr'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Shared_Settings.JobID = 'Thermals'; % An ID that is tacked onto the folder name of all current jobs
Shared_Settings.N_atoms = 10000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Shared_Settings.c_over_a = 1;
Shared_Settings.RunPostProcessor = true;

% Barostat options
Shared_Settings.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Shared_Settings.Isotropy = 'isotropic';
Shared_Settings.Thermal_Solid = true;
Shared_Settings.PreEquilibration = 10;
Shared_Settings.Target_P = 1.0; % Target pressure in bar
Shared_Settings.Time_Constant_P = 1.0; % [ps] time constant for coupling P. Should be 20*Nstpcouple*timestep
Shared_Settings.Nstpcouple = Get_nstcouple(Shared_Settings.Time_Constant_P,Shared_Settings.MDP.dt);

% Cutoff options
Shared_Settings.MDP.RVDW_Cutoff = 1.4; % nm
Shared_Settings.MDP.RCoulomb_Cutoff = 1.4; % nm
Shared_Settings.MDP.RList_Cutoff = 1.4; % nm
Shared_Settings.MDP.Disp_Correction = true;

% Thermostat Options
Shared_Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Shared_Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Shared_Settings.Nsttcouple = Get_nstcouple(Shared_Settings.Time_Constant_T,Shared_Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 
Shared_Settings.Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Shared_Settings.MDP.Initial_T = T0; % Initial termpature at which to generate velocities
Shared_Settings.T0 = T0; % K, Initial temperature


Shared_Settings.Annealing = 'single'; % Options: 'no' 'single' 'periodic'
Shared_Settings.Annealing_Times = [0    1000 2000 3000 4000]; % [ps] A list with the number of annealing reference/control points used
Shared_Settings.Annealing_Temps = [1000 1033 1033 1073 1073]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
Shared_Settings.MDP.Trajectory_Time = 4; % ns

% Postprocessing
Shared_Settings.SavePredictionsImage = true;
Shared_Settings.ML_TimeLength = 21;
Shared_Settings.ML_TimeStep = 5;
Shared_Settings.TimePerFrame = 1;
Settings.Output_Coords
            'SaveTrajectory',Shared_Settings.SaveTrajectory,...
            'SaveFeatures', Shared_Settings.SaveFeatures,...
            'SavePredictions',Shared_Settings.SavePredictions));


for idx = 1:length()
Shared_Settings.Structure = 'Liquid'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Shared_Settings.UseMoltenCompressibility = false; % Use this to raise or lower the barostat compressibility from the default experimental isothermal compressibility

Setup_LiX_Simulation(Shared_Settings)
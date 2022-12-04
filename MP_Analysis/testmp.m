Settings = Initialize_MD_Settings;
Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Settings.Liquid_Fraction = 0.50; % Only meaninful when Liquid_Interface = true. Sets the approximate fraction of the total number of atoms that will initialize as Liquid
Settings.MP_Liquid_Test_Time = 50; % ps. simulation time to sample the liquid for MSD calculation involvedin MP
Settings.MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MP_Equilibrate_Liquid = 20; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings.Cores = 8;
Settings.MPI_Ranks = 8;
Settings.OMP_Threads = 1;
setenv('OMP_NUM_THREADS',num2str(Settings.OMP_Threads))

Settings.Skip_Minimization = false;
Settings.Theory = 'BF'; % Input model(s) to use: JC, JC3P, JC4P, JCSD (for NaCl), TF, BH
Settings.Salt = 'LiI'; % Input model(s) to use: JC, JC3P, JC4P, JCSD (for NaCl), TF, BH
%[JCMX,JCMM,JCXX] = JC_Potential_Parameters(Settings);
Settings = Alexandria_Potential_Parameters(Settings,'vdW_Type','WBK');
% Settings.S.E.MM = Settings.S.E.MM/JCMM.epsilon;
% Settings.S.E.XX = Settings.S.E.XX/JCXX.epsilon;
% Settings.S.E.MX = Settings.S.E.MX/JCMX.epsilon;
% Settings.S.S.MM = Settings.S.S.MM/JCMM.sigma;
% Settings.S.S.XX = Settings.S.S.XX/JCXX.sigma;
% Settings.S.S.MX = Settings.S.S.MX/JCMX.sigma;
Settings.GaussianCharge = true;
Settings.Polarization = true;
Settings.niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
Settings.emtol_polarization = 1e-2; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence
Settings.Structure = 'Rocksalt'; % One of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'Liquid' 'Previous'
Settings.RefStructure = 'Rocksalt'; % Reference structure used for determination of melting or freezing (usually liquid)
Settings.Model = 'Alexandria'; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings.JobID = 'Test'; % An ID that is tacked onto the folder name of all current jobs
Settings.N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings.MDP.Disp_Correction = true;
Settings.MDP.vdw_modifier = 'Potential-shift';
Settings.MinMDP.Disp_Correction = true;
Settings.Find_Min_Params = true;
Settings.Delete_Equil = false;
Settings.initial_opt_type = true;
Settings.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Settings.ScaleInitialLiqDensity = 0.8; % Sets the scale of liquid density
Settings.CheckAmorphousHalide = false; 
Settings.AmorphousDiffThreshold = 1e-6; % [cm^2/s] this sets the diffusion threshold for amorphous vs liquid

% Barostat options
Settings.Isotropy = 'semiisotropic';
Settings.Target_P = [1 1]; % Bar
Settings.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings.Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings.Nstpcouple = Get_nstcouple(Settings.Time_Constant_P,Settings.MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 

% Thermostat Options
T0 = 900; % K %BestKnownMP(Settings);
Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings.Nsttcouple = Get_nstcouple(Settings.Time_Constant_T,Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 
Settings.Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings.MDP.Initial_T = T0; % Initial termpature at which to generate velocities

% Choose interface type
Settings.GenCluster = false; % Switch from FLAT INTERFACE to SPHERICAL INTERFACE
Settings.Cluster_Frac = 0.1; % Cluster fraction, must be between (0,1). Keep this fraction LOW to prevent cluster border issues (does not apply to manual radius jobs)

Settings.CheckTime = 25; % ps. Time between checking for melting/freezing
Settings.MaxCheckTime = 2000; % ps. Max time for melting/freezing runs
Settings.T0 = T0; % K, Initial temperature
Settings.MeltFreezeThreshold = 0.25; % Required CHANGE in fraction or number of atoms either frozen or melted to flag the end of the simulation
Settings.SaveTrajectory = true; % Save processed xyz trajectory file for each temperature check when true
Settings.SavePredictionsImage = true; % Save structure fraction vs time image when true for each temperature check
Settings.c_over_a = 2;
Settings.MaxTDiff = 0.01; % Max change in temperature before generating new initial conditions

% Debugging
Settings.Output_Coords = 1000; % Number of steps between outputting coordinates
Settings.Output_Velocity = 0; % Number of steps between outputting velocities
Settings.Expand_LP = true; % when true, allows expansion of lattice parameters of the solid.
Settings.MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
Settings.MDP.Ewald_rtol = 1e-05; %Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
Settings.MDP.Fourier_Spacing = 0.12; % used 0.1 in minimization. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
Settings.MaxWarn = 2;

% optimizer settings
Settings.Optimizer = 'MPSearcher';
Settings.lb = 0; % K, lower bound on MP search
Settings.ub = 2200; % K, upper bound on MP search
Settings.InitialMeshSize = 25; % Initial step size
Settings.MeshSizeMultiplier = 1;
Settings.BracketThreshold = 10; % [K] Sets the target bracket for the melting point
Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Settings.SlopeThreshold = 1e10; % A second requirement: the change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]

Settings.MDP.RVDW_Cutoff = 1.00; % nm
Settings.MDP.RCoulomb_Cutoff = 1.10; % nm
Settings.MDP.RList_Cutoff = 1.10; % nm
Settings.Cutoff_Buffer = 1.20;


Settings.MinMDP.RVDW_Cutoff = 1.90; % nm
Settings.MinMDP.RCoulomb_Cutoff = 1.90; % nm
Settings.MinMDP.RList_Cutoff = 1.90; % nm
Settings.MinMDP.VerletBT = -1;

Find_Melting_Point(Settings);
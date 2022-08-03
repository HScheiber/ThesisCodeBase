%%%%% Bayesian_Optimize_Local_Submit %%%%%%
%% INFO %%
%% Primary loss options
% Model.Loss_Options.regularization = 'L2'; % Set the regularization scheme
% Model.Loss_Options.Rocksalt.LE = 1;
% Model.Loss_Options.Wurtzite.RLE = 1;
% Model.Loss_Options.NiAs.RLE = 1;
% Model.Loss_Options.Sphalerite.RLE = 1;
% Model.Loss_Options.FiveFive.RLE = 1;

% Model.Loss_Options.Rocksalt.a = 1/10;
% Model.Loss_Options.Wurtzite.a = 1/20;
% Model.Loss_Options.Wurtzite.c = 1/20;
% Model.Loss_Options.NiAs.a = 1/20;
% Model.Loss_Options.NiAs.c = 1/20;
% Model.Loss_Options.Sphalerite.a = 1/10;
% Model.Loss_Options.FiveFive.a = 1/20;
% Model.Loss_Options.FiveFive.c = 1/20;

%% Incorporating energy gaps into loss function
% Value: The target gap between reference and current structure
% Gap.Value < 0 -> "reference structure" is favoured
% Gap.Value > 0 -> "current structure" is favoured
% Gap.Value = 0 -> structures are equal in energy
% Model.Loss_Options.(Structure).Gap.Value = 0;

% Weight: The weighting for this gap in the overall loss function: 
% larger values means more weight. Do not use negative weights!
% Default weight is zero (excluded from loss function)
% Model.Loss_Options.(Structure).Gap.Weight = 0;

% Type: Comparison function. Pick one of: lt | gt | eq | ge | le | ne
% Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
% Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
% Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
% Model.Loss_Options.(Structure).Gap.Type = @lt;

% Ref: The reference structure for the gap. The difference between the
% current structure and the reference structure is considered.
% Model.Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

% Defaults:
% Model.Loss_Options.(Structure).Gap.Value = 0;
% Model.Loss_Options.(Structure).Gap.Weight = 0;
% Model.Loss_Options.(Structure).Gap.Type = @lt;
% Model.Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

%% Set up Model
clear;

%% Global Calculation settings
clear;
% Shared calculation parameters
Model = Initialize_LiX_BO_Settings;
Model.Project_Directory_Name = 'Model_Building';
Model.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Model.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Model.JobSettings.N_Calc = 4; % Number of chained calculations
Model.JobSettings.Hours = 3; % Max time for each job (hours)
Model.JobSettings.Mins = 0; % Max time for job (minutes)
Model.JobSettings.Nodes = 0; % Minimum number of cores to request for calculation.
Model.JobSettings.Cores = 12; % Minimum number of cores to request for calculation. Set to -1 for entire node
Model.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Model.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Model.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Model.MaxWarn = 2;

% MP / Finite T Settings
Model.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Model.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Model.Optimizer = 'MPSearcher';
Model.lb = 0; % K, lower bound on MP search
Model.ub = 2200; % K, upper bound on MP search
Model.BracketThreshold = 5; % [K] Sets the target bracket for the melting point
Model.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Model.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Model.Liquid_Fraction = 0.50;
Model.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Model.Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
Model.Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Model.Delete_Equil = true; % switch to delete temporary calculation folders for finite T calcs

% More Finite T setings
Model.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Model.MDP.RVDW_Cutoff = 1.00; % nm
Model.MDP.RCoulomb_Cutoff = 1.1; % nm
Model.MDP.RList_Cutoff = 1.1; % nm
Model.Cutoff_Buffer = 1.20;
Model.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Model.c_over_a = 2;
Model.N_atoms = 2000;
Model.BracketThreshold = 5; % K
Model.MinStepSize = 0.25;
Model.MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
Model.MeltFreezeThreshold = 0.25;
Model.Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Model.Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Model.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Model.InitialMeshSize = 20;
Model.MeshSizeMultiplier = 5;
Model.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Model.ScaleInitialLiqDensity = 0.8;

% Barostat Options
Model.Isotropy = 'semiisotropic';
Model.Target_P = [1 1]; % Bar
Model.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Model.Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Model.Nstpcouple = Get_nstcouple(Model.Time_Constant_P,Model.MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Model.ScaleCompressibility = 1;

% Thermostat Options
Model.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Model.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Model.Nsttcouple = Get_nstcouple(Model.Time_Constant_T,Model.MDP.dt); %[ps] The frequency for coupling the temperature. 

% Ewald options
Model.MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
Model.MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
Model.MDP.Fourier_Spacing = 0.12;
Model.MDP.VerletBT = -1;

% Set up models below
Exp = Load_Experimental_Data;

Model.Parallel_Bayesopt = false;
Model.Parallel_Struct_Min = false;
Model.Parallel_LiX_Minimizer = false;
Model.MinMDP.Parallel_Min = false;
Model.JobSettings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Model.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Model.JobSettings.npme = 2; % Number of rank assigned to PME
Model.JobSettings.dd = [1 2 3]; % Domain decomposition
Salt = 'LiF';
Theory = 'BH';
Rep = '3';

% Set initial MP temperature
Model.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Model.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
Model.T0 = Exp.(Salt).mp; % K, Initial temperature

%% Model BH: IB3
Model.Salt = Salt;
Model.Theory = Theory;
Model.Trial_ID = ['IB' Rep];
Model.final_opt_type = 'patternsearch';
Model.switch_final_opt = false;
Model.Loss_Convergence = 1e-6;
Model.Param_Convergence = 1e-3;

% Loss
Model.Loss_Options.Rocksalt.LE = 1;
Model.Loss_Options.Rocksalt.a = 1;
Model.Loss_Options.Wurtzite.RLE = 1;

Model.Structures = Auto_Structure_Selection(Model.Loss_Options);
Model.SigmaEpsilon = true;
Model.Fix_Charge = true;
Model.Additivity = true;
Model.MinExpWallHeight = 100; % kJ/mol

%% Submit jobs
workdir = pwd;
% Loop through all models and build their submission file

Calc_Name = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
Calc_Dir = [workdir filesep Calc_Name];
Diary_Loc = [Calc_Dir filesep Calc_Name '.log'];

% Move to new directory
if ~isfolder(Calc_Dir)
    mkdir(Calc_Dir)
end
cd(Calc_Dir);

% Turn diary on and submit job
diary(Diary_Loc);
Bayesian_Optimize_LiX_Parameters(Model)
diary('off')
close all % closes figures

cd(workdir);
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
Settings = Initialize_LiX_BO_Settings;
Settings.Project_Directory_Name = 'Model_Building';
Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Settings.N_Calc = 4; % Number of chained calculations
Settings.Hours = 3; % Max time for each job (hours)
Settings.Mins = 0; % Max time for job (minutes)
Settings.Nodes = 0; % Minimum number of cores to request for calculation.
Settings.Cores = 12; % Minimum number of cores to request for calculation. Set to -1 for entire node
Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Settings.SinglePrecision = false; % choose true for single precision mode, false for double
Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Settings.MaxWarn = 2;

% MP / Finite T Settings
Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Settings.Optimizer = 'MPSearcher';
Settings.lb = 0; % K, lower bound on MP search
Settings.ub = 2200; % K, upper bound on MP search
Settings.BracketThreshold = 5; % [K] Sets the target bracket for the melting point
Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Settings.Liquid_Fraction = 0.50;
Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Settings.Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
Settings.Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Settings.Delete_Equil = true; % switch to delete temporary calculation folders for finite T calcs

% More Finite T setings
Settings.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Settings.MDP.RVDW_Cutoff = 1.00; % nm
Settings.MDP.RCoulomb_Cutoff = 1.1; % nm
Settings.MDP.RList_Cutoff = 1.1; % nm
Settings.Cutoff_Buffer = 1.20;
Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings.c_over_a = 2;
Settings.N_atoms = 2000;
Settings.BracketThreshold = 5; % K
Settings.MinStepSize = 0.25;
Settings.MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
Settings.MeltFreezeThreshold = 0.25;
Settings.MP_Liquid_Test_Time = 100; % ps. Time used for calculation of liquid MSD in melting point calculations.
Settings.MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MP_Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings.InitialMeshSize = 20;
Settings.MeshSizeMultiplier = 5;
Settings.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Settings.ScaleInitialLiqDensity = 0.8;

% Barostat Options
Settings.Isotropy = 'semiisotropic';
Settings.Target_P = [1 1]; % Bar
Settings.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Settings.Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Settings.Nstpcouple = Get_nstcouple(Settings.Time_Constant_P,Settings.MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Settings.ScaleCompressibility = 1;

% Thermostat Options
Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Settings.Nsttcouple = Get_nstcouple(Settings.Time_Constant_T,Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 

% Ewald options
Settings.MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
Settings.MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
Settings.MDP.Fourier_Spacing = 0.12;
Settings.MDP.VerletBT = -1;

% Set up models below
Exp = Load_Experimental_Data;

Settings.Parallel_Bayesopt = false;
Settings.Parallel_Struct_Min = false;
Settings.Parallel_LiX_Minimizer = false;
Settings.MinMDP.Parallel_Min = false;
Settings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Settings.npme = 2; % Number of rank assigned to PME
Settings.dd = [1 2 3]; % Domain decomposition
Salt = 'LiF';
Theory = 'BH';
Rep = '3';

% Set initial MP temperature
Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
Settings.T0 = Exp.(Salt).mp; % K, Initial temperature

%% Model BH: IB3
Settings.Salt = Salt;
Settings.Theory = Theory;
Settings.Trial_ID = ['IB' Rep];
Settings.final_opt_type = 'patternsearch';
Settings.switch_final_opt = false;
Settings.Loss_Convergence = 1e-6;
Settings.Param_Convergence = 1e-3;

% Loss
Settings.Loss_Options.Rocksalt.LE = 1;
Settings.Loss_Options.Rocksalt.a = 1;
Settings.Loss_Options.Wurtzite.RLE = 1;

Settings.Structures = Auto_Structure_Selection(Settings.Loss_Options);
Settings.SigmaEpsilon = true;
Settings.Fix_Charge = true;
Settings.Additivity = true;
Settings.MinExpWallHeight = 100; % kJ/mol

%% Submit jobs
workdir = pwd;
% Loop through all models and build their submission file

Calc_Name = [Settings.Salt '_' Settings.Theory '_Model_' Settings.Trial_ID];
Calc_Dir = [workdir filesep Calc_Name];
Diary_Loc = [Calc_Dir filesep Calc_Name '.log'];

% Move to new directory
if ~isfolder(Calc_Dir)
    mkdir(Calc_Dir)
end
cd(Calc_Dir);

% Turn diary on and submit job
diary(Diary_Loc);
Bayesian_Optimize_LiX_Parameters(Settings)
diary('off')
close all % closes figures

cd(workdir);
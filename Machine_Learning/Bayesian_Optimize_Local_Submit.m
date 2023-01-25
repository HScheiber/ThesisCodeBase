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

%% Set up fixed settings
clear;
close all

% Shared calculation parameters
Settings = Initialize_LiX_BO_Settings;
Settings.Project_Directory_Name = 'Model_Building';
Settings.Verbose = true;
Settings.SigmaEpsilon = true;

% Job settings
Settings.Cores = 8; % Minimum number of cores to request for calculation. Set to -1 for entire node
Settings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Settings.npme = []; % Number of rank assigned to PME
Settings.dd = []; % Domain decomposition

% Bayesian Optimization Stopping Criteria
Settings.Max_Bayesian_Iterations = 400;
Settings.Max_Secondary_Iterations = 50;

% Local optimization settings convergence settings
Settings.Loss_Convergence = 1e-6;
Settings.Param_Convergence = 1e-3;
Settings.switch_final_opt = false;
Settings.Max_Local_Iterations = 10;
Settings.final_opt_type = 'none';

% Parallel Settings
Settings.Parallel_Bayesopt = false;
Settings.Parallel_Struct_Min = true;
Settings.Parallel_LiX_Minimizer = false;

% Bayesian optimization constraints
Settings.MinSkipLoss = 2; % Minimum loss value required before skipping further computation
Settings.BadFcnLossPenalty = 1000; % Penalty to give bad potentials
Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Settings.MaxRepWellDepth = 0; % [kJ/mol] This is the maximum allowed depth of a well between like-like interactions before a loss penalty is applied
Settings.MaxAttWellDepth = -1000; % [kJ/mol] This is the maximum allowed depth of a well between MX interactions before a loss penalty is applied
Settings.MinModelVolume = 10; % [A^3/molecule] minimum allowed volume per molecule of the model solid before finite T calculations are skipped
Settings.MaxModelVolume = 1024; % [A^3/molecule] maximum allowed volume per molecule of the model solid before finite T calculations are skipped
Settings.MinMDP.E_Unphys = -2000; % [kJ/mol] Unphysical energy cutoff
Settings.EnforceRR = true; % enforce radius ratio < 1 (anion is larger)
Settings.MaxMXWellR = 8; % [A] maximum allowed distance for well minima.
Settings.MinMXWellR = 0.5; % [A] minimum allowable distance for well minima.

% Non-MP Finite T calculation settings
Settings.Liquid_Test_Time = 200; % ps. simulation time to sample the liquid for enthalpy / MSD calculations
Settings.Liquid_Equilibrate_Time = 25; % ps. time spent relaxing the liquid for enthalpy / MSD calculations
Settings.Solid_Test_Time = 50; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Settings.CheckAmorphousLiquid = false; % When this is on, calculations that are deemed amorphous solid receive extra loss penalty

% MP-specific calculation settings
Settings.c_over_a = 2;
Settings.MaxCheckTime = 5000; % ps. Max time for MP simulation points
Settings.BracketThreshold = 10; % [K] Sets the target bracket for the melting point
Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Settings.Liquid_Fraction = 0.50;
Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Settings.MP_Liquid_Test_Time = 100; % ps. Simulation time to sample the liquid for MSD calculation
Settings.MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MP_Equilibrate_Liquid = 20; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Settings.InitialMeshSize = 100; % Initial step size for MP calcs
Settings.MeshSizeMultiplier = 2;
Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Settings.Optimizer = 'MPSearcher';
Settings.lb = 0; % K, lower bound on MP search
Settings.ub = 2200; % K, upper bound on MP search

% Shared Finite T calculation Settings
Settings.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Settings.ScaleInitialLiqDensity = 0.8; % Sets the scale of liquid density
Settings.CheckAmorphousHalide = false; 
Settings.AmorphousDiffThreshold = 1e-6; % [cm^2/s] this sets the diffusion threshold for amorphous vs liquid
Settings.Delete_Equil = false; % switch to delete temporary calculation folders for finite T calcs
Settings.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Settings.MDP.RVDW_Cutoff = 1.00; % nm
Settings.MDP.RCoulomb_Cutoff = 1.1; % nm
Settings.MDP.RList_Cutoff = 1.1; % nm
Settings.Cutoff_Buffer = 1.20;
Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
Settings.N_atoms = 2000;

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

% Other stuff
Settings.Initial_N_Multiplier = 20; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Settings.Acquisition_Function = 'expected-improvement-plus';
Settings.ExplorationRatio = 2;
Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Settings.GaussianCharge = false;
Settings.Polarization = false;


%% Test Model Particular parameter
Settings.ShowPlots = false;
Settings.GPActiveSetSize = 1000; % Also applies to final optimization
Settings.Salt = 'LiI';
Settings.Theory = 'BF';
Settings.InnerRange = false;
Settings.Trial_ID = 'QX1';
Settings.UseCoupledConstraint = false;
Settings.Initialize_From_Model = {};

%Settings = Alexandria_Potential_Parameters(Settings,'Coulomb_Only',true); % Loads Gaussian charge parameters

% Loss function
Settings.Loss_Options.Rocksalt.LE   = 1;
Settings.Loss_Options.Rocksalt.a    = 1;
Settings.Loss_Options.Wurtzite.RLE  = 1;
Settings.Loss_Options.FiveFive.RLE  = 0;
Settings.Loss_Options.CsCl.RLE      = 0;
Settings.Loss_Options.Fusion_Enthalpy  = 0; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
Settings.Loss_Options.Liquid_DM_MP = 0; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
Settings.Loss_Options.MP_Volume_Change = 0; % Fitting the experimental change in volume due to melting at the experimental MP
Settings.Loss_Options.Liquid_MP_Volume = 0; % Fitting the experimental volume per formula unit at the experimental MP
Settings.Loss_Options.Solid_MP_Volume  = 0; % Fitting the experimental volume of the experimental solid structure at the experimental MP

% Add gaps
Settings.Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
Settings.Loss_Options.Wurtzite.Gap.Weight = 1000;
Settings.Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
Settings.Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

Settings.Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
Settings.Loss_Options.FiveFive.Gap.Weight = 1000;
Settings.Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
Settings.Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

Settings.Loss_Options.CsCl.Gap.Value = 0; % Negative value:
Settings.Loss_Options.CsCl.Gap.Weight = 1000;
Settings.Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
Settings.Loss_Options.CsCl.Gap.Ref = 'Rocksalt';

% Other loss options options
Settings.Fix_Charge = true;
Settings.Additivity = true;
Settings.Comb_rule = 'Kong'; % One of: 'Lorentz-Berthelot', 'Kong', 'Hogervorst[-wbk]', 'GROMACS'. Only applies to BH and BF models

% Auto structure selection
Settings.Structures = Auto_Structure_Selection(Settings);

%% Start job
[home,project,computer,~,~,~,~,~,~,~,scratch] = find_home;
Calc_Name = [Settings.Salt '_' Settings.Theory '_Model_' Settings.Trial_ID];
Model_Name_abrv = [Settings.Theory '_Model_' Settings.Trial_ID];

if contains(pwd,Calc_Name)
    Calc_Dir = pwd;
else
    Calc_Dir = [pwd filesep Calc_Name];
    if ~isfolder(Calc_Dir)
        mkdir(Calc_Dir)
    end
    cd(Calc_Dir);
end

Settings.Diary_Loc = [Calc_Dir filesep Calc_Name '.log'];
Settings.scratch_dir = fullfile(scratch,Settings.Project_Directory_Name,...
    Settings.Salt,Model_Name_abrv);

% Turn diary on and submit job
diary(Settings.Diary_Loc);
setenv('OMP_NUM_THREADS',num2str(Settings.OMP_Threads))
Bayesian_Optimize_LiX_Parameters(Settings)
diary('off')
close all % closes figures
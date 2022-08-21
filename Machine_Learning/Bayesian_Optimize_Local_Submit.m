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

% Shared calculation parameters
Shared_Settings = Initialize_LiX_BO_Settings;
Shared_Settings.Max_Bayesian_Iterations = 400;
Shared_Settings.Max_Secondary_Iterations = 100;
Shared_Settings.MaxFunEvals = 100; % Only applies to the 'fminsearchbnd' method
Shared_Settings.Loss_Convergence = 1e-6;
Shared_Settings.Param_Convergence = 1e-3;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.final_opt_type = 'fminsearchbnd'; % One of 'none', 'patternsearch', 'fminsearch', 'fminsearchbnd', or fmincon (uses gradients!)
Shared_Settings.switch_final_opt = false;
Shared_Settings.JobSettings.Cores = 8;
Shared_Settings.JobSettings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = []; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 2 4]; % Domain decomposition
Shared_Settings.Cutoff_Buffer = 1.2; % This affects Structure_Minimization as well as other aspects of code
Shared_Settings.MaxWarn = 2;
Shared_Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Shared_Settings.MaxRepWellDepth = 0; % [kJ/mol] This is the maximum allowed depth of a well between like-like interactions before a loss penalty is applied


% MP / Finite T Settings
Shared_Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Shared_Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Shared_Settings.Optimizer = 'MPSearcher';
Shared_Settings.lb = 0; % K, lower bound on MP search
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.BracketThreshold = 5; % [K] Sets the target bracket for the melting point
Shared_Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Shared_Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Shared_Settings.Liquid_Fraction = 0.50;
Shared_Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Shared_Settings.Liquid_Test_Time = 100; % ps. simulation time to sample the liquid for enthalpy / MSD calculations
Shared_Settings.Liquid_Equilibrate_Time = 25; % ps. time spent relaxing the liquid for enthalpy / MSD calculations
Shared_Settings.Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Shared_Settings.Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
Shared_Settings.Delete_Equil = false; % switch to delete temporary calculation folders for finite T calcs

% More Finite T setings
Shared_Settings.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Shared_Settings.MDP.RVDW_Cutoff = 1.00; % nm
Shared_Settings.MDP.RCoulomb_Cutoff = 1.1; % nm
Shared_Settings.MDP.RList_Cutoff = 1.1; % nm
Shared_Settings.Cutoff_Buffer = 1.20;
Shared_Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Shared_Settings.c_over_a = 2;
Shared_Settings.N_atoms = 2000;
Shared_Settings.BracketThreshold = 5; % K
Shared_Settings.MinStepSize = 0.25;
Shared_Settings.MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
Shared_Settings.MeltFreezeThreshold = 0.25;
Shared_Settings.Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Shared_Settings.Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Shared_Settings.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Shared_Settings.InitialMeshSize = 20;
Shared_Settings.MeshSizeMultiplier = 5;
Shared_Settings.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Shared_Settings.ScaleInitialLiqDensity = 0.8;

% Barostat Options
Shared_Settings.Isotropy = 'semiisotropic';
Shared_Settings.Target_P = [1 1]; % Bar
Shared_Settings.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
Shared_Settings.Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
Shared_Settings.Nstpcouple = Get_nstcouple(Shared_Settings.Time_Constant_P,Shared_Settings.MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
Shared_Settings.ScaleCompressibility = 1;

% Thermostat Options
Shared_Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
Shared_Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
Shared_Settings.Nsttcouple = Get_nstcouple(Shared_Settings.Time_Constant_T,Shared_Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 

% Ewald options
Shared_Settings.MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
Shared_Settings.MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
Shared_Settings.MDP.Fourier_Spacing = 0.12;
Shared_Settings.MDP.VerletBT = -1;

idx = 0;

%% BH Model (MP Test) XZ
Salts = {'LiBr'};
Theories = {'BH'};
Replicates = 1;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model BH: XZ
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['XX' Rep];
            
%             % T=0 Loss
%             Models(idx).Loss_Options.Rocksalt.LE = 1;
%             Models(idx).Loss_Options.Rocksalt.a = 1;
%             Models(idx).Loss_Options.Wurtzite.RLE = 1;
            
            % Finite T loss
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Models(idx).Loss_Options.MP  = 1; % Fitting the experimental MP, using the experimental structure as the solid
            
            % Aux options
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
        end
    end
end

%% Submit jobs
workdir = pwd;
% Loop through all models and build their submission file
for idx = 1:length(Models)
    Model = Models(idx);
    Calc_Name = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
    Calc_Dir = [workdir filesep Calc_Name];
    Model.Diary_Loc = [Calc_Dir filesep Calc_Name '.log'];
    
    % Move to new directory
    if ~isfolder(Calc_Dir)
        mkdir(Calc_Dir)
    end
    cd(Calc_Dir);
    
    % Turn diary on and submit job
    diary(Model.Diary_Loc);
    Bayesian_Optimize_LiX_Parameters(Model)
    diary('off')
    close all % closes figures
end
cd(workdir);
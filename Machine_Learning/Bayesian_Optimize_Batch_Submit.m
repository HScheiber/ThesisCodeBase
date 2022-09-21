%%%%% Bayesian_Optimize_Batch_Submit %%%%%%
%% INFO %%
%% Primary loss options
% Models(idx).Loss_Options.regularization = 'L2'; % Set the regularization scheme
% Models(idx).Loss_Options.Rocksalt.LE = 1;
% Models(idx).Loss_Options.Wurtzite.RLE = 1;
% Models(idx).Loss_Options.NiAs.RLE = 1;
% Models(idx).Loss_Options.Sphalerite.RLE = 1;
% Models(idx).Loss_Options.FiveFive.RLE = 1;

% Models(idx).Loss_Options.Rocksalt.a = 1/10;
% Models(idx).Loss_Options.Wurtzite.a = 1/20;
% Models(idx).Loss_Options.Wurtzite.c = 1/20;
% Models(idx).Loss_Options.NiAs.a = 1/20;
% Models(idx).Loss_Options.NiAs.c = 1/20;
% Models(idx).Loss_Options.Sphalerite.a = 1/10;
% Models(idx).Loss_Options.FiveFive.a = 1/20;
% Models(idx).Loss_Options.FiveFive.c = 1/20;

%% Incorporating energy gaps into loss function
% Value: The target gap between reference and current structure
% Gap.Value < 0 -> "reference structure" is favoured
% Gap.Value > 0 -> "current structure" is favoured
% Gap.Value = 0 -> structures are equal in energy
% Models(idx).Loss_Options.(Structure).Gap.Value = 0;

% Weight: The weighting for this gap in the overall loss function: 
% larger values means more weight. Do not use negative weights!
% Default weight is zero (excluded from loss function)
% Models(idx).Loss_Options.(Structure).Gap.Weight = 0;

% Type: Comparison function. Pick one of: lt | gt | eq | ge | le | ne
% Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
% Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
% Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
% Models(idx).Loss_Options.(Structure).Gap.Type = @lt;

% Ref: The reference structure for the gap. The difference between the
% current structure and the reference structure is considered.
% Models(idx).Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

% Defaults:
% Models(idx).Loss_Options.(Structure).Gap.Value = 0;
% Models(idx).Loss_Options.(Structure).Gap.Weight = 0;
% Models(idx).Loss_Options.(Structure).Gap.Type = @lt;
% Models(idx).Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

%% Additional_GAdjust options: 'MM' 'MX' 'XX'
% Introduces additional inverted Gaussians with 3 additional parameters each. 
% Leave empty to exclude
% e.g.: Settings.Additional_GAdjust = {'MM'};
%
% Additional_GAdjust_Ranges: Must be a 1D cell array of the same length as Additional_GAdjust
% Each cell Is a 3 (row) x 2 (column) matrix with:
% First row as the range of the Gaussian depth in kJ/mol (generally should be negative)
% Second row as the range of the Gaussian center in nm (must be strictly positive)
% Third row is the Gaussian standard deviation in nm (must be strictly positive and non-zero)
% e.g. Settings.Additional_GAdjust_Ranges = {[-10 0; 0.3 0.6; 0.01 0.1]};

% Models(idx).Fix_Charge = false;
% Models(idx).Additivity = false;
% Models(idx).Additional_MM_Disp = false;
% Models(idx).Additional_GAdjust = {};
% Models(idx).Additional_GAdjust_Ranges = {};

%% Incorporating dispersion "rational" damping
% Models(idx).C6Damp.MM = 1; % Place rational damping on MM dispersion interaction
% Models(idx).C6Damp.MX = 1; % Place rational damping on MX dispersion interaction
% Models(idx).C6Damp.XX = 1; % Place rational damping on XX dispersion interaction

%% Global Calculation settings
clear;
[home,project,computer,Slurm,~,~,~,~,~,~,scratch] = find_home;

% Some options
skip_models = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs
continue_completed = false; % If a job is already complete, but you wish to continue, this will rename the previous *fullopt.mat file and restart. Must be used with check_complete = false

% Shared calculation parameters
Shared_Settings = Initialize_LiX_BO_Settings;
Shared_Settings.Project_Directory_Name = 'Model_Building';
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.Verbose = false;
Shared_Settings.SigmaEpsilon = true;
Shared_Settings.InnerRange = true; % Sets domain of BH/TF

% Shared job settings
Shared_Settings.JobSettings.N_Calc = 5; % Number of chained calculations
Shared_Settings.JobSettings.Hours = 3; % Max time for each job (hours)
Shared_Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Shared_Settings.JobSettings.Nodes = 0; % Minimum number of cores to request for calculation.
Shared_Settings.JobSettings.Cores = 12; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.

% Local optimization settings convergence settings
Shared_Settings.Loss_Convergence = 1e-6;
Shared_Settings.Param_Convergence = 1e-3;
Shared_Settings.switch_final_opt = false;
Shared_Settings.final_opt_type = 'fminsearchbnd';

% Bayesian optimization constraints
Shared_Settings.MinSkipLoss = 2; % Minimum loss value required before skipping further computation
Shared_Settings.BadFcnLossPenalty = 1000; % Penalty to give bad potentials
Shared_Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Shared_Settings.MaxRepWellDepth = 0; % [kJ/mol] This is the maximum allowed depth of a well between like-like interactions before a loss penalty is applied
Shared_Settings.MaxAttWellDepth = -1500; % [kJ/mol] This is the maximum allowed depth of a well between MX interactions before a loss penalty is applied
Shared_Settings.MinModelVolume = 10; % [A^3/molecule] minimum allowed volume per molecule of the model solid before finite T calculations are skipped
Shared_Settings.MaxModelVolume = 1024; % [A^3/molecule] maximum allowed volume per molecule of the model solid before finite T calculations are skipped
Shared_Settings.MinMDP.E_Unphys = -2000; % [kJ/mol] Unphysical energy cutoff
Shared_Settings.MaxMXWellR = 8; % [A] maximum allowed distance for well minima.
Shared_Settings.MinMXWellR = 0.5; % [A] minimum allowable distance for well minima.

% Non-MP Finite T calculation settings
Shared_Settings.Liquid_Test_Time = 200; % ps. simulation time to sample the liquid for enthalpy / MSD calculations
Shared_Settings.Liquid_Equilibrate_Time = 25; % ps. time spent relaxing the liquid for enthalpy / MSD calculations
Shared_Settings.Solid_Test_Time = 50; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Shared_Settings.CheckAmorphousLiquid = false; % When this is on, calculations that are deemed amorphous solid receive extra loss penalty

% MP-specific calculation settings
Shared_Settings.c_over_a = 2;
Shared_Settings.MaxCheckTime = 5000; % ps. Max time for MP simulation points
Shared_Settings.BracketThreshold = 10; % [K] Sets the target bracket for the melting point
Shared_Settings.MinStepSize = 0.25; % [K] Sets the minimum step size for MPsearcher algorithm
Shared_Settings.SlopeThreshold = 1e10; % The change in the % fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [% Structure Fraction/ps]
Shared_Settings.Liquid_Fraction = 0.50;
Shared_Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Shared_Settings.Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Shared_Settings.Equilibrate_Liquid = 20; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Shared_Settings.PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
Shared_Settings.InitialMeshSize = 100; % Initial step size for MP calcs
Shared_Settings.MeshSizeMultiplier = 2;
Shared_Settings.Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Shared_Settings.MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
Shared_Settings.Optimizer = 'MPSearcher';
Shared_Settings.lb = 0; % K, lower bound on MP search
Shared_Settings.ub = 2200; % K, upper bound on MP search

% Shared Finite T calculation Settings
Shared_Settings.QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
Shared_Settings.ScaleInitialLiqDensity = 0.8; % Sets the scale of liquid density
Shared_Settings.CheckAmorphousHalide = false; 
Shared_Settings.AmorphousDiffThreshold = 1e-6; % [cm^2/s] this sets the diffusion threshold for amorphous vs liquid
Shared_Settings.Delete_Equil = false; % switch to delete temporary calculation folders for finite T calcs
Shared_Settings.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Shared_Settings.MDP.RVDW_Cutoff = 1.00; % nm
Shared_Settings.MDP.RCoulomb_Cutoff = 1.1; % nm
Shared_Settings.MDP.RList_Cutoff = 1.1; % nm
Shared_Settings.Cutoff_Buffer = 1.20;
Shared_Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
Shared_Settings.N_atoms = 2000;

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

% Experimental Data
Exp = Load_Experimental_Data;

%% Initial Model index
idx=0;
switch lower(computer)
    case 'cedar'
        %% Shared_Settings
        Shared_Settings.Max_Bayesian_Iterations = 300;
        Shared_Settings.Max_Secondary_Iterations = 200;
        Shared_Settings.Max_Local_Iterations = 100;
        Shared_Settings.Parallel_Bayesopt = false;
        Shared_Settings.Parallel_Struct_Min = true;
        Shared_Settings.Parallel_LiX_Minimizer = false;
        Shared_Settings.UseCoupledConstraint = true;
        Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Shared_Settings.InnerRange = true; % Sets domain of BH
        Shared_Settings.EnforceRR = true;
        
        %% BH [Gamma<6], BD, [Gamma<6] and Mie Models: OA, OB, OC on LiX with proper radius ratios
        Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
        Theories = {'BH' 'BD' 'Mie'};
        Replicates = 1:5;
        for tidx = 1:length(Theories)
            Theory = Theories{tidx};
            for sidx = 1:length(Salts)
                Salt = Salts{sidx};
                
                % Set initial MP temperature
                Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
                Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature
                
                for ridx = 1:length(Replicates)
                    Rep = num2str(Replicates(ridx));
                    
                    %% Model OA
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['OA' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model OB
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['OB' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 2;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 0.2; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                end
            end
        end
        
    case 'narval'
        %% Shared_Settings
        Shared_Settings.Max_Bayesian_Iterations = 300;
        Shared_Settings.Max_Secondary_Iterations = 200;
        Shared_Settings.Max_Local_Iterations = 200;
        Shared_Settings.Parallel_Bayesopt = false;
        Shared_Settings.Parallel_Struct_Min = true;
        Shared_Settings.Parallel_LiX_Minimizer = false;
        Shared_Settings.UseCoupledConstraint = true;
        Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Shared_Settings.InnerRange = true; % Sets domain of BH
        Shared_Settings.EnforceRR = false;
        Shared_Settings.final_opt_type = 'patternsearch';
        
        %% BH Models: MA, MB, MC, MF on LiX
        Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
        Theories = {'BH'};
        Replicates = 1:5;
        for tidx = 1:length(Theories)
            Theory = Theories{tidx};
            for sidx = 1:length(Salts)
                Salt = Salts{sidx};
                
                % Set initial MP temperature
                Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
                Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature
                
                for ridx = 1:length(Replicates)
                    Rep = num2str(Replicates(ridx));
                    
%                     %% Model MA
%                     idx = idx+1;
%                     Models(idx) = Shared_Settings;
%                     Models(idx).Salt = Salt;
%                     Models(idx).Theory = Theory;
%                     Models(idx).Trial_ID = ['MA' Rep];
%                     
%                     % Loss function
%                     Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
%                     Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
%                     Models(idx).Loss_Options.Rocksalt.LE  = 1;
%                     Models(idx).Loss_Options.Rocksalt.a  = 1;
%                     
%                     
%                     Models(idx).Structures = Auto_Structure_Selection(Models(idx));
%                     Models(idx).Fix_Charge = true;
%                     Models(idx).Additivity = true;
%                     
%                     %% Model MB
%                     idx = idx+1;
%                     Models(idx) = Shared_Settings;
%                     Models(idx).Salt = Salt;
%                     Models(idx).Theory = Theory;
%                     Models(idx).Trial_ID = ['MB' Rep];
%                     
%                     % Loss function
%                     Models(idx).Loss_Options.Rocksalt.LE  = 1;
%                     Models(idx).Loss_Options.Rocksalt.a  = 1;
%                     Models(idx).Loss_Options.Wurtzite.RLE  = 1;
%                     Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
%                     Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
%                     
%                     Models(idx).Structures = Auto_Structure_Selection(Models(idx));
%                     Models(idx).Fix_Charge = true;
%                     Models(idx).Additivity = true;
                    
                    
                    %% Model MC
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['MC' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
                    Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
                    Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model MF
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['MF' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.FiveFive.RLE  = 1;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
                    Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
                    Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                end
            end
        end
        
        %% Shared_Settings
        Shared_Settings.Max_Bayesian_Iterations = 300;
        Shared_Settings.Max_Secondary_Iterations = 200;
        Shared_Settings.Max_Local_Iterations = 100;
        Shared_Settings.Parallel_Bayesopt = false;
        Shared_Settings.Parallel_Struct_Min = true;
        Shared_Settings.Parallel_LiX_Minimizer = false;
        Shared_Settings.UseCoupledConstraint = false;
        Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Shared_Settings.InnerRange = false; % Sets domain of BH
        Shared_Settings.EnforceRR = true;
        
        %% BD[Gamma>7], BE[Gamma>7] Models: OA, OB on LiX
        Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
        Theories = {'BD' 'BE'};
        Replicates = 1:5;
        for tidx = 1:length(Theories)
            Theory = Theories{tidx};
            for sidx = 1:length(Salts)
                Salt = Salts{sidx};
                
                % Set initial MP temperature
                Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
                Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature
                
                for ridx = 1:length(Replicates)
                    Rep = num2str(Replicates(ridx));
                    
                    %% Model OA
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['OA' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model OB
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['OB' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 2;
                    Models(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
                    Models(idx).Loss_Options.Liquid_DM_MP = 0.2; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                end
            end
        end
        
    case 'graham'
        %% Shared_Settings
        Shared_Settings.Max_Bayesian_Iterations = 600;
        Shared_Settings.Max_Secondary_Iterations = 200;
        Shared_Settings.Max_Local_Iterations = 1000;
        Shared_Settings.switch_final_opt = true;
        Shared_Settings.Parallel_Bayesopt = true;
        Shared_Settings.Parallel_Struct_Min = false;
        Shared_Settings.Parallel_LiX_Minimizer = false;
        Shared_Settings.UseCoupledConstraint = true;
        Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Shared_Settings.InnerRange = true; % Sets domain of BH/TF
        Shared_Settings.EnforceRR = false;
        
        %% BH Models: NA NB NC ND NE NF NG (inner range - targets on T=0 crystals only)
        Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
        Theories = {'BH'};
        Replicates = 1:5;
        for tidx = 1:length(Theories)
            Theory = Theories{tidx};
            for sidx = 1:length(Salts)
                Salt = Salts{sidx};
                
                % Set initial MP temperature
                Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
                Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature

                for ridx = 1:length(Replicates)
                    Rep = num2str(Replicates(ridx));

                    %% Model NA
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NA' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model NB
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NB' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model NC
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NC' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = false;
                    Models(idx).Additivity = true;
                    
                    %% Model ND
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['ND' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = false;
                    
                    
                    %% Model NE
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NE' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.NiAs.RLE  = 1;
                    Models(idx).Loss_Options.AntiNiAs.RLE  = 1;
                    Models(idx).Loss_Options.BetaBeO.RLE  = 1;
                    Models(idx).Loss_Options.Sphalerite.RLE  = 1;
                    Models(idx).Loss_Options.FiveFive.RLE  = 1;
                    Models(idx).Loss_Options.CsCl.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model NF
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NF' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.NiAs.RLE  = 1;
                    Models(idx).Loss_Options.AntiNiAs.RLE  = 1;
                    Models(idx).Loss_Options.BetaBeO.RLE  = 1;
                    Models(idx).Loss_Options.Sphalerite.RLE  = 1;
                    Models(idx).Loss_Options.FiveFive.RLE  = 1;
                    Models(idx).Loss_Options.CsCl.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = false;
                    Models(idx).Additivity = true;
                    
                    %% Model NG
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NG' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    Models(idx).Loss_Options.NiAs.RLE  = 1;
                    Models(idx).Loss_Options.AntiNiAs.RLE  = 1;
                    Models(idx).Loss_Options.BetaBeO.RLE  = 1;
                    Models(idx).Loss_Options.Sphalerite.RLE  = 1;
                    Models(idx).Loss_Options.FiveFive.RLE  = 1;
                    Models(idx).Loss_Options.CsCl.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = false;
                end
            end
        end
        
    otherwise % Place jobs here for later assignment
        %% Shared_Settings
        Shared_Settings.Max_Bayesian_Iterations = 600;
        Shared_Settings.Max_Secondary_Iterations = 200;
        Shared_Settings.Max_Local_Iterations = 1000;
        Shared_Settings.switch_final_opt = true;
        Shared_Settings.Parallel_Bayesopt = true;
        Shared_Settings.Parallel_Struct_Min = false;
        Shared_Settings.Parallel_LiX_Minimizer = false;
        Shared_Settings.UseCoupledConstraint = false;
        Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
        Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
        Shared_Settings.InnerRange = false; % Sets domain of BH/TF
        Shared_Settings.EnforceRR = true;
        
        %% BH Models: NA NB NC ND (outer range - targets on T=0 crystals only, RR enforced)
        Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
        Theories = {'BD' 'BE' 'Mie'};
        Replicates = 1:5;
        for tidx = 1:length(Theories)
            Theory = Theories{tidx};
            for sidx = 1:length(Salts)
                Salt = Salts{sidx};
                
                % Set initial MP temperature
                Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
                Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature

                for ridx = 1:length(Replicates)
                    Rep = num2str(Replicates(ridx));

                    %% Model NA
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NA' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model NB
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NB' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = true;
                    
                    %% Model NC
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['NC' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = false;
                    Models(idx).Additivity = true;
                    
                    %% Model ND
                    idx = idx+1;
                    Models(idx) = Shared_Settings;
                    Models(idx).Salt = Salt;
                    Models(idx).Theory = Theory;
                    Models(idx).Trial_ID = ['ND' Rep];
                    
                    % Loss function
                    Models(idx).Loss_Options.Rocksalt.LE  = 1;
                    Models(idx).Loss_Options.Rocksalt.a  = 1;
                    Models(idx).Loss_Options.Wurtzite.RLE  = 1;
                    
                    Models(idx).Structures = Auto_Structure_Selection(Models(idx));
                    Models(idx).Fix_Charge = true;
                    Models(idx).Additivity = false;
                    
                end
            end
        end
        
end

%% Check for already running jobs
if check_running && Slurm
    er = 1;
    jdx = 1;
    while er ~= 0
        [er,out] = system('squeue -u $USER -o %100j | tail -n +2');
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
    
    % Remove the job number tags
    c_jobs_unique = unique(regexprep(c_jobs,'-([0-9]){3}',''));
else
    c_jobs_unique = {''};
end

%% Prepare batch scripts
% Find the scheduler and get the batch script template file

calc_cmd = 'matlab -nodisplay -r "Bayesian_Optimize_LiX_Parameters(''##TASKNAME##.inp'')"';
clean_cmd =  'matlab -r "cleanup_BO_log(''##TASKNAME##.log'')"';

% Loop through all models and build their submission file
for idx = 1:length(Models)
    if any(idx == skip_models)
        continue
    end
    Model = Models(idx);
    
    Batch_Template = MD_Batch_Template(Model.JobSettings);
    Batch_Template = strrep(Batch_Template,['##PREMIN##' newline],'');
    Batch_Template = strrep(Batch_Template,'##MDRUN##',calc_cmd);
    Batch_Template = strrep(Batch_Template,'##CLEANUP##',clean_cmd);
    Batch_Template = strrep(Batch_Template,['##EXT1##' newline],'');
    Batch_Template = strrep(Batch_Template,['##EXT2##' newline],'');
    
    
    Model_Name = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
    Model_Name_abrv = [Model.Theory '_Model_' Model.Trial_ID];
    Batch_Template = strrep(Batch_Template,'##TASKNAME##',Model_Name);
    Batch_Template = strrep(Batch_Template,'##ERROR##',Model_Name);
    
    % Generate and move to the submission directory
    submit_dir = fullfile(project,Model.Project_Directory_Name,Model.Salt,Model_Name_abrv);
    Model.scratch_dir = fullfile(scratch,Model.Project_Directory_Name,Model.Salt,Model_Name_abrv);
    Model.OuterDir = submit_dir;
    Batch_Template = strrep(Batch_Template,'##DIRECTORY##',submit_dir);
    Model.Diary_Loc = fullfile(submit_dir,[Model_Name '.log']);
    if ~exist(submit_dir, 'dir')
       mkdir(submit_dir)
    end
    cd(submit_dir);
    
    % Check if job is already complete
    if check_complete
        obs = dir(['*Model_' Model.Trial_ID '*fullopt.mat']);
        if ~isempty(obs)
            disp([Model_Name ': Job already completed. Skipping Job Submission.'])
            continue
        end
    end
    
    % Check if Model is already running
    if check_running && ismember(Model_Name,c_jobs_unique) 
        disp([Model_Name ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    if continue_completed
        obs = dir(['*Model_' Model.Trial_ID '*fullopt.mat']);
        if ~isempty(obs)
            disp([Model_Name ': Job already completed. Continuing completed Job.'])
            src = fullfile(obs.folder,obs.name);
            dest = fullfile(obs.folder,strrep(obs.name,'fullopt','oldopt'));
            movefile(src,dest);
        end
        obs = dir(['*Model_' Model.Trial_ID '*final_point.mat']); 
        if ~isempty(obs)
            src = fullfile(obs.folder,obs.name);
            dest = fullfile(obs.folder,strrep(obs.name,'final_point','old_finp'));
            movefile(src,dest);
        end
        
        BPT_Folder = fullfile(submit_dir,'BestPoint_Thermal');
        BPT_New = fullfile(submit_dir,'OldPoint_Thermal');
        if isfolder(BPT_Folder)
            movefile(BPT_Folder,BPT_New)
        end
    end
    
    % Save input file
    inp_file = fullfile(submit_dir,[Model_Name '.inp']);
    save(inp_file,'Model','-mat')
    
    % Save job submission script
    subm_file = fullfile(submit_dir,[Model_Name '.subm']);
    fid = fopen(subm_file,'wt');
    fwrite(fid,regexprep(Batch_Template,{'\r'}',{''}));
    fclose(fid);
    
    % Submit job
    MultiSubmitBayesOpt(Model.JobSettings.N_Calc,submit_dir,[Model_Name '.subm'])
    delete(subm_file); % This is used as a template and can be removed
end

% Return
cd(project);

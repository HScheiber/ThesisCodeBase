%%%%% Bayesian_Optimize_Batch_Submit_Cedar %%%%%%
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
% Shared calculation parameters
Shared_Settings = Initialize_LiX_BO_Settings;
Shared_Settings.Project_Directory_Name = 'Model_Building';
Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Submit_Jobs = false; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.JobSettings.N_Calc = 4; % Number of chained calculations
Shared_Settings.JobSettings.Hours = 3; % Max time for each job (hours)
Shared_Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Shared_Settings.JobSettings.Nodes = 0; % Minimum number of cores to request for calculation.
Shared_Settings.JobSettings.Cores = 12; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.MaxWarn = 2;

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
Shared_Settings.Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
Shared_Settings.Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Shared_Settings.Delete_Equil = true; % switch to delete temporary calculation folders for finite T calcs

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

% Set up models below
Exp = Load_Experimental_Data;
skip_models = [1:60 61:240]; % 1:60 61:240 240:1000
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs
continue_completed = false; % If a job is already complete, but you wish to continue, this will rename the previous *fullopt.mat file and restart. Must be used with check_complete = false

% Initial Model index
idx=0;

%% JC Models EV EW, EX
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.JobSettings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 6; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 0; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 1 2]; % Domain decomposition
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    % Set initial MP temperature
    Shared_Settings.Target_T = Exp.(Salt).mp; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Shared_Settings.MDP.Initial_T = Exp.(Salt).mp; % Initial termpature at which to generate velocities
    Shared_Settings.T0 = Exp.(Salt).mp; % K, Initial temperature
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EV
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EV' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EW
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EW' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EW
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EW' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1/3;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
    end
end

%% BH Models IA, IB, IC, ID
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 2; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: IA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IA' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: IB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IB' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: IC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IC' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: ID
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['ID' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
        end
    end
end

%% BH Models IE, IF, IG, IH, II
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 2; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: IE
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IE' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: IF
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IF' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Wurtzite.a = 2/3;
            Models(idx).Loss_Options.Wurtzite.c = 1/3;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: IG
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IG' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: IH
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IH' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: II
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['II' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
        end
    end
end

%% BH Models JA, JB, JC
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.JobSettings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 2; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 2 5]; % Domain decomposition
Shared_Settings.final_opt_type = 'none'; % One of 'none', 'patternsearch', 'fminsearch', 'fminsearchbnd', or fmincon (uses gradients!)
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

            %% Model BH: JA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JA' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: JB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JB' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            
            %% Model BH: JC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JC' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
        end
    end
end

%% JC Models JA, JB, JC
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.JobSettings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.JobSettings.OMP_Threads = 6; % Set the number of OMP threads per MPI rank
Shared_Settings.JobSettings.npme = 0; % Number of rank assigned to PME
Shared_Settings.JobSettings.dd = [1 1 2]; % Domain decomposition
Shared_Settings.final_opt_type = 'none'; % One of 'none', 'patternsearch', 'fminsearch', 'fminsearchbnd', or fmincon (uses gradients!)
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'JC'};
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

            %% Model BH: JA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JA' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            %% Model BH: JB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JB' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
            
            %% Model BH: JC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JC' Rep];
            Models(idx).final_opt_type = 'patternsearch';
            Models(idx).switch_final_opt = false;
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            Models(idx).MinExpWallHeight = 100; % kJ/mol
            
        end
    end
end

%% Check for already running jobs
if check_running && isunix
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
[home,project,~,Slurm] = find_home;
calc_cmd = 'matlab -nodisplay -r "Bayesian_Optimize_LiX_Parameters(''##TASKNAME##.inp'')" >> ##TASKNAME##.log';
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
    submit_dir = fullfile(project,'Model_Building',Model.Salt,Model_Name_abrv);
    Batch_Template = strrep(Batch_Template,'##DIRECTORY##',submit_dir);
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
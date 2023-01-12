function Settings = Initialize_LiX_BO_Settings
% Script to generate default LiX bayesian optimization calculation settings

%% Calculation parameters
Settings = Initialize_MD_Settings;
Settings.Trial_ID = ''; % Optional ID to add
Settings.Salt = 'LiI';
Settings.Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
Settings.Theory = 'JC'; % JC, TF
Settings.Continue_Calculation = true; % Continue the calculation from where it left off
Settings.Fix_Charge = true; % When true, fixes coulombic charges at 1
Settings.Additivity = true; % Only applies to JC model. When true, use combining rules
Settings.Comb_rule = 'Lorentz-Berthelot'; % (case insensitive). One of: 'Lorentz-Berthelot', 'Kong', 'Hogervorst', 'GROMACS'. Only applies to BH and BF models
Settings.Additional_MM_Disp = false; % When Additivity is activated with the JC model, tacks on an additional non-additive metal-metal interaction
Settings.SigmaEpsilon = false; % For all models, recasts the search space in terms of sigma/epsilon rather than dispersion/repulsion
Settings.InnerRange = false; % When SigmaEpislon is on, use this to set the TF/BH inner gamma range vs outer gamma range
Settings.Additional_Function = Init_Additional_Function; % Adding additional function options
Settings.Fix_Alpha = false; % Active for TF and BH models only when NOT using SigmaEpsilon form. When active, the value of exponential repulsion parameter, which is related to compressibility, is fixed
Settings.Fix_C8 = false; % Active for TF model only when NOT using SigmaEpsilon form. When active, the value of C8 is fixed by the value of C6, and is no longer a free parameter
Settings.Fix_Mie_n = true; % When using the Mie model, this fixes the value of the Born exponent, when false the value is an optimizable parameter
Settings.Diary_Loc = '';

% This feature allows new calculations to start from the outputs of previous completed calculations
% Note that the previous calculations must have the same targets, or at least cover all of the current targets
Settings.Initialize_From_Model = {};
Settings.Initialize_From_Model_Subsample = Inf;

% Charge
Settings.Q_Range = [0.95 1.05]; % Default range for charge scaling. Only meaningful when Fix_Charge = false
Settings.Q_value = 1.0; % Default value for the charge scale. Only meaningful when Fix_Charge = true

% Additional_GAdjust options: 'MM' 'MX' 'XX'
% Introduces additional inverted Gaussians with 3 additional parameters each. 
% Leave empty to exclude
Settings.Additional_GAdjust = {};
% Additional_GAdjust_Ranges: Must be a 1D cell array of the same length as Additional_GAdjust
% Each cell Is a 3 (row) x 2 (column) matrix with:
% First row as the range of the Gaussian depth in kJ/mol (generally should be negative)
% Second row as the range of the Gaussian center in nm (must be strictly positive)
% Third row is the Gaussian standard deviation in nm (must be strictly positive and non-zero)
Settings.Additional_GAdjust_Ranges = {};

% Dispersion damping
% C6Damping:
% 0 = no (default) damping. This is default of JC/TF model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)
Settings.Extra_Properties = true; % calculate extra properties for the final computation

% Loss and constraint options
Settings.Loss_Options = init_loss_options; % default loss options
Settings.UseCoupledConstraint = false; % Turns on coupled constraints
Settings.CheckBadFcn = true; % Switch to turn on or off the checking of pathological functions, adding a loss penalty for such functions
Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Settings.MaxRepWellDepth = 0; % [kJ/mol] This is the maximum allowed depth of a well between like-like interactions before a loss penalty is applied
Settings.MaxAttWellDepth = -1000; % [kJ/mol] This is the maximum allowed depth of a well between MX interactions before a loss penalty is applied
Settings.MinModelVolume = 10; % [A^3/molecule] minimum allowed volume per molecule of the model solid before finite T calculations are skipped
Settings.MaxModelVolume = 2000; % [A^3/molecule] maximum allowed volume per molecule of the model solid before finite T calculations are skipped
Settings.MaxMXWellR = 10; % [A] maximum allowed distance for well minima.
Settings.MinMXWellR = 0.5; % [A] minimum allowable distance for well minima.
Settings.BadFcnLossPenalty = 1000; % Penalty to give when (1) the function shape is deemed pathological or (2) the optimized crystal volume is too small or large
Settings.MinSkipLoss = 2; % Minimum loss value required before skipping further computation
Settings.EnforceRR = true; % Enforce radius ratios

% Intial optimization settings
Settings.GPActiveSetSize = 1000; % Also applies to final optimization
Settings.initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Settings.Initial_N_Multiplier = 10; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Settings.Max_Bayesian_Iterations = 800;
Settings.Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Settings.MinSurrogatePoints = 20; % Minimum number of random sample points to create at the start of a surrogate creation phase. Used only by surrogateopt
Settings.MinSampleDistance = 1e-1; % Controls how often surrogateopt resets but setting the minimum sample distance between points along any dimension.
Settings.ExplorationRatio = 2; % Exploration ratio used in 'plus' Acquisition Functions, or kappa used in lower-confidence-bound
Settings.KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
Settings.ShowPlots = true; % Switch to turn on/off BO plots
% Available Kernels:
% 'exponential', 'squaredexponential', 'matern32', 'matern52', 
% 'rationalquadratic', 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Acquisition Functions
% 'expected-improvement'
% 'expected-improvement-plus' with ExplorationRatio determining exploration/exploitation trade off
% 'lower-confidence-bound' with ExplorationRatio = kappa
% 'probability-of-improvement'


% (optional) secondary bayesian optimization settings
Settings.second_opt_type = 'bayesopt'; % One of 'bayesopt' or 'none'
Settings.Max_Secondary_Iterations = 200;
Settings.Secondary_Acquisition_Function = 'lower-confidence-bound'; % The acquisition function used in the secondary bayesian optimization
Settings.ExplorationRatio_Secondary = 1; % Exploration ratio used in 'plus' Acquisition Functions, or kappa used in lower-confidence-bound
Settings.SecondaryKernelFunction = 'ardmatern52'; % The covariance function used in the secondary bayesian optimization

% Final optimization settings
Settings.Max_Local_Iterations = 1000; % max iterations of the final optimizations
Settings.MaxFunEvals = 100; % Only applies to the 'fminsearchbnd' method
Settings.Loss_Convergence = 1e-6; % options for the final local optimization
Settings.Param_Convergence = 1e-3; % Note: with patternsearch only one of Loss_Convergence or Param_Convergence must be satisfied, but with fminsearch both must be satisfied simultaneously
Settings.final_opt_type = 'fminsearchbnd'; % One of 'none', 'patternsearch', 'fminsearch', 'fminsearchbnd', or fmincon (uses gradients!)
Settings.switch_final_opt = true; % When true, switch from 'fminsearchbnd' to 'patternsearch' after Settings.MaxFunEvals is reached

% Parallel settings (only one can be true at a time)
Settings.Parallel_LiX_Minimizer = false; % Run the parallel version of the LiX_Minimizer subroutine when true
Settings.Parallel_Struct_Min = false; % Run the parallel version of the LiX_Minimizer subroutine when true (this is generally less efficient)
Settings.Parallel_Bayesopt = true; % Run the parallel version of Bayesian Optimization when true

end


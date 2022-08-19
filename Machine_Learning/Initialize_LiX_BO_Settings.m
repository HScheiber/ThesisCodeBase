function Settings = Initialize_LiX_BO_Settings
% Script to generate default LiX bayesian optimization calculation settings

%% Calculation parameters
Settings = Initialize_MD_Settings;
Settings.Trial_ID = ''; % Optional ID to add
Settings.Salt = 'LiI';
Settings.Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
Settings.Theory = 'JC'; % JC, TF
Settings.Restart_Calculation = true; % restart the calculation from where it left off
Settings.Fix_Charge = false; % When true, fixes coulombic charges at 1
Settings.Additivity = false; % Only applies to JC model. When true, use combining rules
Settings.Additional_MM_Disp = false; % When Additivity is activated with the JC model, tacks on an additional non-additive metal-metal interaction
Settings.SigmaEpsilon = false; % For the JC model, recasts the search space in terms of sigma/epsilon rather than dispersion/repulsion
Settings.Additional_Function = Init_Additional_Function; % Adding additional function options
Settings.Fix_Alpha = false; % Active for TF and BH models only when NOT using SigmaEpsilon form. When active, the value of exponential repulsion parameter, which is related to compressibility, is fixed to its default
Settings.Fix_C8 = false; % Active for TF model only when NOT using SigmaEpsilon form. When active, the value of C8 is fixed by the value of C6, and is no longer a free parameter
Settings.Diary_Loc = '';

% Parameter Ranges
Settings.SDMM_Range = [0 1000]; % Default range for MM dispersion scaling in JC model.
Settings.SDXX_Range = [0 40]; % Default range for XX dispersion scaling in JC/BH model.
Settings.SDMX_Range = [0 40]; % Default range for MX dispersion scaling in JC/BH model.
Settings.SDMM2_Range = [0 1000]; % Default range for MM dispersion scaling in JC model. Only used when Additional_MM_Disp = true

Settings.SD6MM_Range = [0 1000]; % Default range for C6 MM dispersion scaling in TF model.
Settings.SD6XX_Range = [0 10]; % Default range for C6 XX dispersion scaling in TF model.
Settings.SD6MX_Range = [0 100]; % Default range for C6 MX dispersion scaling in TF model.

Settings.SD8MM_Range = [0 1000]; % Default range for C8 MM dispersion scaling in TF model.
Settings.SD8XX_Range = [0 20]; % Default range for C8 XX dispersion scaling in TF model.
Settings.SD8MX_Range = [0 200]; % Default range for C8 MX dispersion scaling in TF model.

Settings.SRTFMM_Range = [0.1 500]; % Default range for MM repulsive prefactor scaling in TF model.
Settings.SRTFXX_Range = [0.1 20]; % Default range for XX repulsive prefactor scaling in TF model.
Settings.SRTFMX_Range = [0.1 20]; % Default range for MX repulsive prefactor scaling in TF model.

Settings.SAMM_Range = [0.1 10]; % Default range for MM repulsive exponent scaling in TF model.
Settings.SAXX_Range = [0.1 5]; % Default range for XX repulsive exponent scaling in TF model.
Settings.SAMX_Range = [0.1 5]; % Default range for MX repulsive prefactor scaling in TF model.

Settings.SRMM_Range = [0.3 10]; % Default range for MM repulsive prefactor scaling in JC/BH model.
Settings.SRXX_Range = [0.3 5]; % Default range for XX repulsive prefactor scaling in JC/BH model.
Settings.SRMX_Range = [0.3 5]; % Default range for MX repulsive prefactor scaling in JC/BH model.

% Exp-C6 or Exp-C6-C8: "sigma-epsilon" version parameters
Settings.Sr0MM_Range = [0.1 0.4];  % [nm] Default range for MM r0 parameter
Settings.Sr0XX_Range = [0.1 0.5]; % [nm] Default range for XX r0 parameter
Settings.Sr0MX_Range = [0.1 0.5]; % [nm] Default range for MX r0 parameter

Settings.SepsMM_Range = [0 1000]; % [kJ/mol] Default range for MM epsilon parameter
Settings.SepsXX_Range = [0 1000]; % [kJ/mol] Default range for XX epsilon parameter
Settings.SepsMX_Range = [0 1000]; % [kJ/mol] Default range for MX epsilon parameter

Settings.SgamMM_Range = [1 48/7]; % [Unitless] Default range for MM gamma parameter in TF models
Settings.SgamXX_Range = [1 48/7]; % [Unitless] Default range for XX gamma parameter in TF models
Settings.SgamMX_Range = [1 48/7]; % [Unitless] Default range for MX gamma parameter in TF models
Settings.SgamMM_RangeBH = [6 50]; % [Unitless] Default range for MM gamma parameter in BH models
Settings.SgamXX_RangeBH = [6 50]; % [Unitless] Default range for XX gamma parameter in BH models
Settings.SgamMX_RangeBH = [6 50]; % [Unitless] Default range for MX gamma parameter in BH models

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

% Loss options
Settings.Loss_Options = init_loss_options; % default loss options
Settings.CheckBadFcn = true; % Switch to turn on or off the checking of pathological functions, adding a loss penalty for such functions
Settings.MinExpWallHeight = 100; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Settings.MaxRepWellDepth = 0; % [kJ/mol] This is the maximum allowed depth of a well between like-like interactions before a loss penalty is applied
Settings.MaxModelMismatch = 1; % [Fraction] This is the maximum allowed mismatch for the experimental vs model lattice energy before a loss penalty is applied and finite T calculations are skipped
Settings.BadFcnLossPenalty = 1000; % Penalty to give when function shape is deemed pathological

% Intial optimization settings
Settings.initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
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

% Finite T settings
Settings.Liquid_Test_Time = 50; % ps
Settings.Solid_Test_Time = 30; % ps

end


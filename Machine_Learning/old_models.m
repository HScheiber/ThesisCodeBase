
%% Model 1: LiI BF2 (targeting DFT/0.44 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF2';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';


Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [0 5.0000e-11];

%% Model 2: LiI BF3 (targeting 1 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF3';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -1;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 3: LiI BF4 (targeting 2 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF4';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -2;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 4: LiI BF5 (targeting 3 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF5';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -3;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 5: LiI BF6 (targeting 4 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF6';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -4;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 6: LiI BF7 (targeting 5 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF7';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -5;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 7: LiI BF8 (targeting 10 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BF8';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

% Loss function
Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -10;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.2900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;
Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-12 6e-11];

%% Model 8: LiI BG1
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG1';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';


Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-21 3e-20];

%% Model 9: LiI BG2
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG2';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 10: LiI BG3: (targetting 1 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG3';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -1;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 11: LiI BG4: (targetting 2 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG4';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -2;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 12: LiI BG5: (targetting 3 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG5';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -3;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 13: LiI BG6: (targetting 4 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG6';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -4;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 14: LiI BG7: (targetting 5 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG7';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -5;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 15: LiI BG8: (targetting 10 kJ/mol gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BG8';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 1;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -10;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Wurtzite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Rocksalt.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Rocksalt.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 4e-42];

%% Model 9: LiI BH1
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH1';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-21 3e-20];

%% Model 10: LiI BH2
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH2';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 11: LiI BH3: (targetting 1 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH3';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -1;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 12: LiI BH4: (targetting 2 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH4';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -2;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 13: LiI BH5: (targetting 3 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH5';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -3;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 14: LiI BH6: (targetting 4 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH6';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -4;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 15: LiI BH7: (targetting 5 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH7';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -5;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 16: LiI BH8: (targetting 10 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BH8';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -10;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 45: LiI BI1
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI1';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.LE = 10;
Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1900; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 50;
Models(idx).Additional_Function.MM.Range = [1e-21 3e-20];

%% Model 46: LiI BI2
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI2';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.LE = 10;
Models(idx).Loss_Options.NiAs.RLE = 10000;
Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;

Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 47: LiI BI3: (targetting 1 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI3';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -1;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 48: LiI BI4: (targetting 2 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI4';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -2;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 49: LiI BI5: (targetting 3 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI5';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -3;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 50: LiI BI6: (targetting 4 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI6';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -4;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 51: LiI BI7: (targetting 5 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI7';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -5;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% Model 52: LiI BI8: (targetting 10 kJ/mol RS/NiAs gap)
idx = idx+1;
Models(idx) = Initialize_LiX_BO_Settings;
Models(idx).Salt = 'LiI';
Models(idx).Theory = 'JC';
Models(idx).Trial_ID = 'BI8';

% Initial optimization
Models(idx).initial_opt_type = 'bayesopt'; % One of 'bayesopt' or 'surrogateopt'
Models(idx).Max_Bayesian_Iterations = 1000;
Models(idx).Acquisition_Function = 'expected-improvement-plus'; % The acquisition function used in bayesian optimization
Models(idx).KernelFunction = 'ardmatern52'; % The covariance function used in the primary bayesian optimization
% 'ardexponential', 'ardsquaredexponential', 
% 'ardmatern32', 'ardmatern52', 'ardrationalquadratic'

% Secondary optimization
Models(idx).second_opt_type = 'none';
Models(idx).final_opt_type = 'fminsearchbnd';

Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Experimental_LP = false; % Targets experimental rocksalt lattice energy rather than DFT
Models(idx).Loss_Options.Rocksalt.LE = 10;
Models(idx).Loss_Options.Rocksalt.a = 1/10;

Models(idx).Loss_Options.NiAs.a = 1/20;
Models(idx).Loss_Options.NiAs.c = 1/20;

% Targetting gaps
Models(idx).Loss_Options.Rocksalt.Gap.Value = -10;
Models(idx).Loss_Options.Rocksalt.Gap.Weight = 10000;
Models(idx).Loss_Options.Rocksalt.Gap.Ref = 'NiAs';
Models(idx).Loss_Options.Rocksalt.Gap.Type = @eq;

Models(idx).Loss_Options.Wurtzite.Gap.Value = 0;
Models(idx).Loss_Options.Wurtzite.Gap.Weight = 10000;
Models(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';
Models(idx).Loss_Options.Sphalerite.Gap = Models(idx).Loss_Options.Wurtzite.Gap;
Models(idx).Loss_Options.FiveFive.Gap = Models(idx).Loss_Options.Wurtzite.Gap;

Models(idx).Fix_Charge = true;
Models(idx).Additivity = true;
Models(idx).SDMM_Range = [0 20];

Models(idx).C6Damp.input_rvdw = true;
Models(idx).C6Damp.rvdw.Li = 0.1850; % nm
Models(idx).C6Damp.N.MM = 1;
Models(idx).C6Damp.MM = 1;

Models(idx).Additional_Function.MM.N = 100;
Models(idx).Additional_Function.MM.Range = [5e-44 3e-42];

%% LiBr JC Models C3a - C3e
Salts = {'LiBr'};
Replicates = {'a' 'b' 'c' 'd' 'e'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model 1
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['C3' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10000;
        Models(idx).Loss_Options.Rocksalt.a = 1/10;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 100];
    end
end

%% LiBr JC Models C5 testing Acquisition Functions
%'expected-improvement'
%'expected-improvement-plus' with ExplorationRatio 0.5, 1, 2
%'lower-confidence-bound' with ExplorationRatio (kappa) 1, 2, 3
%'probability-of-improvement'

Salts = {'LiBr'};
Replicates = {'a' 'b' 'c'};
for ridx = 1:length(Replicates)
    Rep = Replicates{ridx};

    %% Model 1
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_EI' Rep];
    Models(idx).Acquisition_Function = 'expected-improvement';
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 2
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_EIP0.5' Rep];
    Models(idx).Acquisition_Function = 'expected-improvement-plus';
    Models(idx).ExplorationRatio = 0.5;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 3
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_EIP1' Rep];
    Models(idx).Acquisition_Function = 'expected-improvement-plus';
    Models(idx).ExplorationRatio = 1;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 4
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_EIP2' Rep];
    Models(idx).Acquisition_Function = 'expected-improvement-plus';
    Models(idx).ExplorationRatio = 2;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 5
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_LCB1' Rep];
    Models(idx).Acquisition_Function = 'lower-confidence-bound';
    Models(idx).ExplorationRatio = 1;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 6
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_LCB2' Rep];
    Models(idx).Acquisition_Function = 'lower-confidence-bound';
    Models(idx).ExplorationRatio = 2;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 7
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_LCB3' Rep];
    Models(idx).Acquisition_Function = 'lower-confidence-bound';
    Models(idx).ExplorationRatio = 3;
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
    %% Model 8
    idx = idx+1;
    Models(idx) = Initialize_LiX_BO_Settings;
    Models(idx).Salt = 'LiBr';
    Models(idx).Theory = 'JC';
    Models(idx).Trial_ID = ['C5_POI' Rep];
    Models(idx).Acquisition_Function = 'probability-of-improvement';
    Models(idx).final_opt_type = 'fminsearchbnd';

    Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
    Models(idx).Loss_Options.Rocksalt.LE = 1;
    Models(idx).Loss_Options.Wurtzite.RLE = 10000;
    Models(idx).Loss_Options.Rocksalt.a = 1/10;
    Models(idx).Loss_Options.Wurtzite.a = 1/20;
    Models(idx).Loss_Options.Wurtzite.c = 1/20;

    Models(idx).Fix_Charge = false;
    Models(idx).Additivity = true;
    Models(idx).SDMM_Range = [0 100];
    
end

%% JC Models CA, CB, CC
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = {'a' 'b' 'c' 'd' 'e'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model 1
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CA' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
        
        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10000;
        Models(idx).Loss_Options.Wurtzite.RLE_freedom = 2; % kJ/mol
        Models(idx).Loss_Options.Rocksalt.a = 1/10;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model 2
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CB' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 10000;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1/10;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model 3
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CC' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 10000/10;
        Models(idx).Loss_Options.Wurtzite.a = 10000/20;
        Models(idx).Loss_Options.Wurtzite.c = 10000/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
    end
end

%% JC Models CD, CE
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = {'a' 'b' 'c' 'd' 'e'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model 1
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CD' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 10000;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model 2
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CE' Rep];
        Models(idx).second_opt_type = 'none';
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 10000;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 10000;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
    end
end

%% JC Models CF
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = {'a' 'b' 'c' 'd' 'e'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model 1
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CF' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10000;
        Models(idx).Loss_Options.Rocksalt.a = 1/10;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CG, CH, CI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CG' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fmincon';
        Models(idx).Structures = {'Rocksalt'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
        
        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 100;
        Models(idx).Loss_Options.Rocksalt.a = 100;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CH' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fmincon';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'
        
        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CI
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CI' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fmincon';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CJ, CK, CL
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CJ' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;
        

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CK' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CI
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CL' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;        

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CM, CN, CO
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CM
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CM' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CN
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CN' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 20;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 20;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 20;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.RLE = 20;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 20;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CO
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CO' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'FiveFive'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 20;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 20;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 20;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;        

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CP, CQ, CR
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CP
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CP' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1/5;
        Models(idx).Loss_Options.Wurtzite.c = 1/5;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model CQ
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CQ' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1/5;
        Models(idx).Loss_Options.Wurtzite.c = 1/5;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 100];
        
        %% Model CR
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CR' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1/5;
        Models(idx).Loss_Options.Wurtzite.c = 1/5;   

        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CS, CT
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CS
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CS' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 100];
        
        %% Model CT
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CT' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models DA, DB, DC, DD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model DB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model DC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        %% Model DD
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DD' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 1000];
    end
end

%% JC Models DE, DF, DG, DH
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DE
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DE' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model DF
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DF' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        %% Model DG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DG' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.V = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.V = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.V = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.V = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.V = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.V = 1;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.V = 1;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model DH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DH' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models CT
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CT
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['CT' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        
        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
    end
end

%% JC Models DA, DB, DC, DD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DD
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DD' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 1000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
    end
end

%% JC Models DE, DF, DG, DH
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DE
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DE' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DF
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DF' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DG' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.V = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.V = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.V = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.V = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.V = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.V = 1;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.V = 1;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['DH' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
    end
end

%% JC Models EA, EB, EC, ED, EE
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:10;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: ED
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ED' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        %% Model JC: EE
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EE' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
    end
end

%% JC Models EF, EG, EH, EI, EJ, EK, EL, EM
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EF
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EF' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EG' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        %% Model JC: EH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EH' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EI
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EI' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EJ
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EJ' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
        %% Model JC: EK
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EK' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EL
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EL' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EM
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EM' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
    end
end

%% JC Models EN, EO, EP
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:10;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EN
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EN' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EO
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EO' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EP
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EP' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];
        
    end
end

%% JC Models EQ, ER, ES, ET, EU
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EQ
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EQ' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: ER
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ER' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: ES
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ES' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 300];
        
        
        %% Model JC: ET
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ET' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EU
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EU' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        
    end
end

%% JC Models EQ, ER, ES, ET, EU
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: ET
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ET' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EU
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EU' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        
    end
end

%% JC Models EV EW, EX
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EV
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
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
        Models(idx) = Initialize_LiX_BO_Settings;
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
        Models(idx) = Initialize_LiX_BO_Settings;
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

%% TF Models CA, CB, CC, CD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = {'a' 'b' 'c' 'd' 'e'};
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10000;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 10000;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1/10;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 10000;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CD
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CD' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 10000;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 10000;
        Models(idx).Loss_Options.Wurtzite.a = 1/20;
        Models(idx).Loss_Options.Wurtzite.c = 1/20;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
    end
end

%% TF Models CG, CH, CI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CG' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 100;
        Models(idx).Loss_Options.Rocksalt.a = 100;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CH' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CI
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CI' Rep];
        Models(idx).second_opt_type = 'bayesopt';
        Models(idx).final_opt_type = 'fminsearchbnd';

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
    end
end

%% TF Models CJ, CK, CL
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CJ' Rep];
        Models(idx).final_opt_type = 'patternsearch'; % fminsearchbnd
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;
        

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CK' Rep];
        Models(idx).final_opt_type = 'patternsearch'; % fminsearchbnd
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CI
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CL' Rep];
        Models(idx).final_opt_type = 'patternsearch'; % fminsearchbnd
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;        

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
    end
end

%% TF Models CM, CN, CO
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model CM
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CM' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CN
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CN' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'


        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 20;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 20;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 20;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;
        Models(idx).Loss_Options.NiAs.RLE = 20;
        Models(idx).Loss_Options.NiAs.a = 1;
        Models(idx).Loss_Options.NiAs.c = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 20;
        Models(idx).Loss_Options.AntiNiAs.a = 1;
        Models(idx).Loss_Options.AntiNiAs.c = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        %% Model CO
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['CO' Rep];
        Models(idx).final_opt_type = 'patternsearch';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'FiveFive'}; % 'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 20;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        
        Models(idx).Loss_Options.Wurtzite.RLE = 20;
        Models(idx).Loss_Options.Wurtzite.a = 1;
        Models(idx).Loss_Options.Wurtzite.c = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 20;
        Models(idx).Loss_Options.FiveFive.a = 1;
        Models(idx).Loss_Options.FiveFive.c = 1;        

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
    end
end

%% TF Models DA, DB, DC, DD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Fix_Charge = false;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = false;
        Models(idx).SD6MM_Range = [0 2000];
        Models(idx).SD8MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DD
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DD' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        Models(idx).Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        Models(idx).SD8MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
    end
end

%% TF Models DE, DF, DG, DH
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = arrayfun(@num2str,1:5,'UniformOutput',false);
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = Replicates{ridx};
        
        %% Model DE
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DE' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        Models(idx).final_opt_type = 'patternsearch';
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DF
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DF' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = false;
        Models(idx).SD6MM_Range = [0 2000];
        Models(idx).SD8MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DG
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DG' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.V = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 10;
        Models(idx).Loss_Options.Wurtzite.V = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.V = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 10;
        Models(idx).Loss_Options.Sphalerite.V = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 10;
        Models(idx).Loss_Options.FiveFive.V = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 10;
        Models(idx).Loss_Options.AntiNiAs.V = 1;
        Models(idx).Loss_Options.BetaBeO.RLE = 10;
        Models(idx).Loss_Options.BetaBeO.V = 1;
        Models(idx).Loss_Options.CsCl.RLE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
        
        %% Model DH
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['DH' Rep];

        Models(idx).Loss_Options.Experimental_LE = true; % Targets experimental rocksalt lattice energy rather than DFT
        Models(idx).Loss_Options.Experimental_LP = true; % Targets experimental rocksalt lattice param rather than DFT
        
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 10;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 10;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 10;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 10;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 10;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 10;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 2000];
        
        % Used in older models
        Models(idx).CRDamp.MX.r_d = 0.15; % This is the value of the sigmoid's midpoint in nm. Set to a negative value to disable close range damping
        Models(idx).CRDamp.MX.b = 100; % sigmoid "steepness" for damping
        Models(idx).CRDamp.MM.r_d = 0.30; % LiCl = 0.21 , LiI = 0.24
        Models(idx).CRDamp.MM.b  = 75; % 75
        Models(idx).CRDamp.XX.r_d = 0.20; 
        Models(idx).CRDamp.XX.b  = 100;
    end
end

%% TF Models EA, EB, EC, ED, EE, EF
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:10;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model TF: EA
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['EA' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        Models(idx).Loss_Options.BetaBeO.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
        
        %% Model TF: EB
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['EB' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        Models(idx).Loss_Options.BetaBeO.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
        
        %% Model TF: EC
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['EC' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
        
        %% Model TF: ED
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['ED' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
        
        %% Model TF: EE
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['EE' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
        
        %% Model TF: EF
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'TF';
        Models(idx).Trial_ID = ['EF' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Max_Bayesian_Iterations = 1000;
            Models(idx).Max_Secondary_Iterations = 1000;
        else
            Models(idx).Max_Bayesian_Iterations = 800;
            Models(idx).Max_Secondary_Iterations = 200;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.RLE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.RLE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.RLE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.RLE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;
        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Fix_Alpha = true;
        Models(idx).Fix_C8 = true;
        Models(idx).SD6MM_Range = [0 3000];
    end
end

%% TF & BH Models FA, FB, FC, FD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:10;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: FA
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FA' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: FB
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FB' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: FC
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FC' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: FD
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FD' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% TF & BH Models FE, FF, FG, FH, FI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:10;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: FE
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FE' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: FF
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FF' Rep];
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
            Models(idx).Loss_Options.Wurtzite.a = 2/3;
            Models(idx).Loss_Options.Wurtzite.c = 1/3;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: FG
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FG' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: FH
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FH' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: FI
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['FI' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% TF & BH Models GA, GB, GC, GD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:5;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: GA
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GA' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GB
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GB' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GC
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GC' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: GD
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GD' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% TF & BH Models GE, GF, GG, GH, GI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:5;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: GE
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GE' Rep];
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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GF
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GF' Rep];
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
            Models(idx).Loss_Options.Wurtzite.a = 2/3;
            Models(idx).Loss_Options.Wurtzite.c = 1/3;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GG
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GG' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GH
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GH' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: GI
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GI' Rep];
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
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% BH Models JA, JB, JC
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 2; % Number of rank assigned to PME
Shared_Settings.dd = [1 2 5]; % Domain decomposition
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

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JB' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JC' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

        end
    end
end

%% JC Models JA, JB, JC
Shared_Settings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 6; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 0; % Number of rank assigned to PME
Shared_Settings.dd = [1 1 2]; % Domain decomposition
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

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JB' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JC' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

        end
    end
end

%% JC Models EA, EB, ED, EE
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};

    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));

        %% Model JC: EA
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EA' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EB
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EB' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: ED
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ED' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];

        %% Model JC: EE
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EE' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

    end
end

%% JC Models EG, EH, EJ, EK, EM
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};

    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));

        %% Model JC: EG
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EG' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];

        %% Model JC: EH
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EH' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EJ
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EJ' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];

        %% Model JC: EK
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EK' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EM
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EM' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];

    end
end

%% JC Models EN, EP
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};

    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));

        %% Model JC: EN
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EN' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: EP
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EP' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.LE = 1;
        Models(idx).Loss_Options.Wurtzite.a = 2/3;
        Models(idx).Loss_Options.Wurtzite.c = 1/3;
        Models(idx).Loss_Options.NiAs.LE = 1;
        Models(idx).Loss_Options.NiAs.a = 2/3;
        Models(idx).Loss_Options.NiAs.c = 1/3;
        Models(idx).Loss_Options.Sphalerite.LE = 1;
        Models(idx).Loss_Options.Sphalerite.a = 1;
        Models(idx).Loss_Options.FiveFive.LE = 1;
        Models(idx).Loss_Options.FiveFive.a = 2/3;
        Models(idx).Loss_Options.FiveFive.c = 1/3;
        Models(idx).Loss_Options.AntiNiAs.LE = 1;
        Models(idx).Loss_Options.AntiNiAs.a = 2/3;
        Models(idx).Loss_Options.AntiNiAs.c = 1/3;
        Models(idx).Loss_Options.BetaBeO.LE = 1;
        Models(idx).Loss_Options.BetaBeO.a = 2/3;
        Models(idx).Loss_Options.BetaBeO.c = 1/3;
        Models(idx).Loss_Options.CsCl.LE = 1;
        Models(idx).Loss_Options.CsCl.a = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 1000];

    end
end

%% JC Models EQ, ER, ES, ET
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};

    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));

        %% Model JC: EQ
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EQ' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: ER
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ER' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];

        %% Model JC: ES
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ES' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = false;
        Models(idx).Additivity = false;
        Models(idx).SDMM_Range = [0 300];

        %% Model JC: ET
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['ET' Rep];

        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.AntiNiAs.RLE = 1;

        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
    end
end

%% JC Models EV, EW, EX
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};

    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));

        %% Model JC: EV
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EV' Rep];

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

        %% Model JC: EX
        idx = idx+1;
        Models(idx) = Shared_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EX' Rep];

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

%% BH Models JD, JE, JF, JG
Shared_Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 2; % Number of rank assigned to PME
Shared_Settings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: JD
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JD' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JE
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JE' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JF
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JF' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

            %% Model BH: JG
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JG' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;

        end
    end
end

%% JC Models JD, JE, JF
Shared_Settings.MinExpWallHeight = 300; % [kJ/mol] in TF and BH models, this is the minimum allowed heighted of the repulsive wall before a loss penalty is applied
Shared_Settings.ub = 2200; % K, upper bound on MP search
Shared_Settings.MPI_Ranks = 2; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 6; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 0; % Number of rank assigned to PME
Shared_Settings.dd = [1 1 2]; % Domain decomposition
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

            %% Model BH: JD
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JD' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JE
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JE' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JF
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JF' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

        end
    end
end

%% BH Models IA, IB, IC, ID
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

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: IB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IB' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;


            %% Model BH: IC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IC' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;

            %% Model BH: ID
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['ID' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

        end
    end
end

%% BH Models IE, IF, IG, IH, II
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

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: IF
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IF' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Wurtzite.a = 2/3;
            Models(idx).Loss_Options.Wurtzite.c = 1/3;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: IG
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IG' Rep];

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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: IH
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['IH' Rep];

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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;

            %% Model BH: II
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['II' Rep];

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

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
        end
    end
end

%% BH Models JH, JI, JJ, JK, [not JL], JM
Shared_Settings.N_Calc = 10; % Number of chained calculations
Shared_Settings.Hours = 12; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 2; % Number of rank assigned to PME
Shared_Settings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: JH
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JH' Rep];

            % Loss
            Models(idx).Loss_Options.MP  = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JI
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JI' Rep];

            % Loss
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JJ
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JJ' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model BH: JK
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JK' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

            %% Model BH: JM
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JM' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;

        end
    end
end

%% JC Models JH, JI, JJ, JK, [not JL], JM
Shared_Settings.N_Calc = 10; % Number of chained calculations
Shared_Settings.Hours = 12; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 1; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 12; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
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

            %% Model JC: JH
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JH' Rep];

            % Loss
            Models(idx).Loss_Options.MP  = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model JC: JI
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JI' Rep];

            % Loss
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model JC: JJ
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JJ' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model JC: JK
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JK' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

            %% Model JC: JM
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JM' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.MP  = 2;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
        end
    end
end

%% BH Models JL
Shared_Settings.N_Calc = 5; % Number of chained calculations
Shared_Settings.Hours = 3; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 2; % Number of rank assigned to PME
Shared_Settings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: JL
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JL' Rep];

            % Loss
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

        end
    end
end

%% JC Models JL
Shared_Settings.N_Calc = 5; % Number of chained calculations
Shared_Settings.Hours = 3; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 1; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 12; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
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

            %% Model JC: JL
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JL' Rep];

            % Loss
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
        end
    end
end

%% BH Models JN, JO, JP
Shared_Settings.N_Calc = 5; % Number of chained calculations
Shared_Settings.Hours = 3; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = 2; % Number of rank assigned to PME
Shared_Settings.dd = [1 2 5]; % Domain decomposition
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

            %% Model BH: JN
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JN' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            %% Model BH: JO
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JO' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            %% Model BH: JP
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JP' Rep];

            Models(idx).CheckAmorphousHalide = true; % check to make sure halide ion is properly mobile in the liquid

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
        end
    end
end

%% JC Models JN, JO, JP
Shared_Settings.N_Calc = 5; % Number of chained calculations
Shared_Settings.Hours = 3; % Max time for each job (hours)
Shared_Settings.MPI_Ranks = 1; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 12; % Set the number of OMP threads per MPI rank
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
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

            %% Model JC: JN
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JN' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            %% Model JC: JO
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JO' Rep];

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            %% Model JC: JP
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['JP' Rep];

            Models(idx).CheckAmorphousHalide = true; % check to make sure halide ion is properly mobile in the liquid

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).SigmaEpsilon = false;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
        end
    end
end




%% New Models

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 200;
Shared_Settings.Max_Secondary_Iterations = 100;
Shared_Settings.Max_Local_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = false; % Sets domain of BH

%% NaCl - JC/BH Models: KA, KB, KC, KD, KE on NaCl
Salts = {'NaCl'}; % 'LiF' 'LiCl' 'LiBr' 'LiI' 
Theories = {'JC' 'BH'};
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

            %% Model KA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KA' Rep];

            % Loss function
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model KB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KB' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model KC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KC' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model KD
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KD' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

            %% Model KE
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KE' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 200;
Shared_Settings.Max_Secondary_Iterations = 100;
Shared_Settings.Max_Local_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = false; % Sets domain of BH

%% LiX - JC/BH Models: KA, KE
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'JC' 'BH'};
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

            %% Model KA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KA' Rep];

            % Loss function
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model KE
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KE' Rep];

            % Loss function
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 200;
Shared_Settings.Max_Secondary_Iterations = 100;
Shared_Settings.Max_Local_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH

%% NaCl - TF Model: KB and KF
Salts = {'NaCl'}; % 'LiF' 'LiCl' 'LiBr' 'LiI' 
Theories = {'TF'};
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

            %% Model KB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KB' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model KF
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['KF' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 2;
            Models(idx).Loss_Options.Rocksalt.a  = 2;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 300;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH

%% BH Models: MA, MB, MC on NaCl
Salts = {'NaCl'}; % 'LiF' 'LiCl' 'LiBr' 'LiI' 
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

            %% Model MA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['MA' Rep];

            % Loss function
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;


            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model MB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['MB' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;
            Models(idx).Loss_Options.Wurtzite.RLE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

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
        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 600;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 1000;
Shared_Settings.switch_final_opt = true;
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = false; % Sets domain of BH

%% JC/BH Models: LA, LB, LC, LD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' };
Theories = {'JC' 'BH'};
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

            %% Model LA
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['LA' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model LB
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['LB' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;
            Models(idx).Loss_Options.Wurtzite.RLE  = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model LC
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['LC' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Rocksalt.a  = 1;
            Models(idx).Loss_Options.Wurtzite.RLE  = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;

            %% Model LD
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['LD' Rep];

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

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 300;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH

%% BH Models: MD, ME on LiX
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

            %% Model MD
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['MD' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.LE  = 1;
            Models(idx).Loss_Options.Wurtzite.RLE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Models(idx).Structures = Auto_Structure_Selection(Models(idx));
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;

            %% Model ME
            idx = idx+1;
            Models(idx) = Shared_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['ME' Rep];

            % Loss function
            Models(idx).Loss_Options.Rocksalt.a  = 1;
            Models(idx).Loss_Options.Wurtzite.RLE  = 1;
            Models(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

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
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH
Shared_Settings.EnforceRR = true;

%% BH [Gamma<6], BD, [Gamma<6] and Mie Models: OA, OB on LiX with proper radius ratios
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
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['OA' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model OB
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['OB' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 2;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.2; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;
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
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = false; % Sets domain of BH
Shared_Settings.EnforceRR = true;

%% BH [Gamma>7], BD, [Gamma>7] and Mie Models: OC on LiX with proper radius ratios
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

            %% Model OC
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['OC' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.CsCl.RLE  = 1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 3; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;
        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 400;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH
Shared_Settings.EnforceRR = false;
Shared_Settings.final_opt_type = 'none';

%% BH [Gamma<6] and JC Models: MG
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'BH' 'JC'};
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

            %% Model MG
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MG' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 300;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 200;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH
Shared_Settings.EnforceRR = false;
Shared_Settings.final_opt_type = 'patternsearch';

%% BH [Gamma<6] Models: MA, MB, MC, MF on LiX
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
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MC' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model MF
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MF' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

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
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
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
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['OA' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model OB
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['OB' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 2;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 2; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.2; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;
        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 400;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH
Shared_Settings.EnforceRR = false;
Shared_Settings.final_opt_type = 'none';

%% Mie Models: MG
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'Mie'};
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

            %% Model MG
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MG' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 600;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 1000;
Shared_Settings.switch_final_opt = true;
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = true;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
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
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NA' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model NB
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NB' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model NC
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NC' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = false;
            Settings_Array(idx).Additivity = true;

            %% Model ND
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['ND' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = false;

            %% Model NE
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NE' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.NiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.AntiNiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.BetaBeO.RLE  = 1;
            Settings_Array(idx).Loss_Options.Sphalerite.RLE  = 1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 1;
            Settings_Array(idx).Loss_Options.CsCl.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model NF
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NF' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.NiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.AntiNiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.BetaBeO.RLE  = 1;
            Settings_Array(idx).Loss_Options.Sphalerite.RLE  = 1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 1;
            Settings_Array(idx).Loss_Options.CsCl.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = false;
            Settings_Array(idx).Additivity = true;

            %% Model NG
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NG' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;
            Settings_Array(idx).Loss_Options.NiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.AntiNiAs.RLE  = 1;
            Settings_Array(idx).Loss_Options.BetaBeO.RLE  = 1;
            Settings_Array(idx).Loss_Options.Sphalerite.RLE  = 1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 1;
            Settings_Array(idx).Loss_Options.CsCl.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = false;
        end
    end
end

%% Shared_Settings
Shared_Settings.Max_Bayesian_Iterations = 600;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Max_Local_Iterations = 1000;
Shared_Settings.switch_final_opt = true;
Shared_Settings.Parallel_Bayesopt = true;
Shared_Settings.Parallel_Struct_Min = false;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
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
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NA' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model NB
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['NB' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

            %% Model ND
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['ND' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE  = 1;
            Settings_Array(idx).Loss_Options.Rocksalt.a  = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 1;

            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = false;

        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.InnerRange = true; % Sets domain of BH/TF
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';

%% BF Models: MG
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'BF'};
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

            %% Model MI
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MG' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 1; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement';
Shared_Settings.ExplorationRatio = 1;
Shared_Settings.Max_Bayesian_Iterations = 50;
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.InnerRange = true; % Sets domain of BH
Shared_Settings.EnforceRR = true;
Shared_Settings.second_opt_type = 'none';
Shared_Settings.final_opt_type = 'none';
Shared_Settings.Liquid_Test_Time = 300; % ps. simulation time to sample the liquid for enthalpy / MSD calculations
Shared_Settings.Liquid_Equilibrate_Time = 25; % ps. time spent relaxing the liquid for enthalpy / MSD calculations
Shared_Settings.Solid_Test_Time = 300; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
Shared_Settings.Initialize_From_Model = {'MG'};

%% BH [Gamma<6] and JC Models: MH
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'BH' 'JC'};
Replicates = 1;
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

            %% Model MH
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MH' Rep];

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 10;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 10; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';

%% Mie Models: MI (n = 10)
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'Mie'};
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

            %% Model MI
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MI' Rep];
            Settings_Array(idx).S.n.MM = 10;
            Settings_Array(idx).S.n.XX = 10;
            Settings_Array(idx).S.n.MX = 10;

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Mie Models: MJ (n = 8)
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'Mie'};
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

            %% Model MJ
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MJ' Rep];
            Settings_Array(idx).S.n.MM = 8;
            Settings_Array(idx).S.n.XX = 8;
            Settings_Array(idx).S.n.MX = 8;

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 2;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end


%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.InnerRange = false; % Sets domain of BH/TF
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';
Shared_Settings.GaussianCharge = false;
Shared_Settings.Polarization = false;
Shared_Settings.Initialize_From_Model = {};

%% BF [Gamma>10] without gaussian charge Models: MH
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'BF'};
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

            %% Model MH
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MH' Rep];
            Settings_Array(idx) = Alexandria_Potential_Parameters(Settings_Array(idx),...
                'Coulomb_Only',true);

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 10;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 10; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.InnerRange = false; % Sets domain of BH/TF
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';
Shared_Settings.GaussianCharge = true;
Shared_Settings.Polarization = false;
Shared_Settings.Initialize_From_Model = {};

%% BF [Gamma>10] WITH gaussian charge Models: MI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; 
Theories = {'BF'};
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

            %% Model MH
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MI' Rep];
            Settings_Array(idx) = Alexandria_Potential_Parameters(Settings_Array(idx),...
                'Coulomb_Only',true);

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 10;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 10; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end


%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';
Shared_Settings.GaussianCharge = true;
Shared_Settings.Polarization = false;
Shared_Settings.Initialize_From_Model = {};

%% JC WITH gaussian charge Models: MI
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

            %% Model MI
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MI' Rep];
            Settings_Array(idx) = Alexandria_Potential_Parameters(Settings_Array(idx),...
                'Coulomb_Only',true);

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 10;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 10; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end

%% Shared_Settings
Shared_Settings.Initial_N_Multiplier = 40; % Multiply the number of input dimensions by this number to obtain the number of initial random points
Shared_Settings.Acquisition_Function = 'expected-improvement-plus';
Shared_Settings.ExplorationRatio = 2;
Shared_Settings.Max_Bayesian_Iterations = 800;
Shared_Settings.Max_Secondary_Iterations = 200;
Shared_Settings.Secondary_Acquisition_Function = 'expected-improvement'; % The acquisition function used in the secondary bayesian optimization
Shared_Settings.Parallel_Bayesopt = false;
Shared_Settings.Parallel_Struct_Min = true;
Shared_Settings.Parallel_LiX_Minimizer = false;
Shared_Settings.UseCoupledConstraint = false;
Shared_Settings.InnerRange = true; % Sets domain of BH/TF
Shared_Settings.MPI_Ranks = 12; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Shared_Settings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Shared_Settings.EnforceRR = true;
Shared_Settings.final_opt_type = 'none';
Shared_Settings.GaussianCharge = true;
Shared_Settings.Polarization = false;
Shared_Settings.Initialize_From_Model = {};

%% BH [Gamma<6] WITH gaussian charge Models: MI
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

            %% Model MI
            idx = idx+1;
            Settings_Array(idx) = Shared_Settings;
            Settings_Array(idx).Salt = Salt;
            Settings_Array(idx).Theory = Theory;
            Settings_Array(idx).Trial_ID = ['MI' Rep];
            Settings_Array(idx) = Alexandria_Potential_Parameters(Settings_Array(idx),...
                'Coulomb_Only',true);

            % Loss function
            Settings_Array(idx).Loss_Options.Rocksalt.LE   = 10;
            Settings_Array(idx).Loss_Options.Rocksalt.a    = 1;
            Settings_Array(idx).Loss_Options.Wurtzite.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.FiveFive.RLE  = 0.1;
            Settings_Array(idx).Loss_Options.CsCl.RLE      = 0.1;
            Settings_Array(idx).Loss_Options.Fusion_Enthalpy  = 10; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_DM_MP = 0.1; % Fitting the experimental metal ion diffusion constant of the molten salt at the experimental MP
            Settings_Array(idx).Loss_Options.MP_Volume_Change = 10; % Fitting the experimental change in volume due to melting at the experimental MP
            Settings_Array(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Settings_Array(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP

            % Add gaps
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.Wurtzite.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.FiveFive.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.FiveFive.Gap.Ref = 'Rocksalt';

            Settings_Array(idx).Loss_Options.CsCl.Gap.Value = 0; % Negative value:
            Settings_Array(idx).Loss_Options.CsCl.Gap.Weight = 1000;
            Settings_Array(idx).Loss_Options.CsCl.Gap.Type = @lt; % pick one of: lt | gt | eq | ge | le | ne
            Settings_Array(idx).Loss_Options.CsCl.Gap.Ref = 'Rocksalt';


            Settings_Array(idx).Structures = Auto_Structure_Selection(Settings_Array(idx));
            Settings_Array(idx).Fix_Charge = true;
            Settings_Array(idx).Additivity = true;

        end
    end
end
        
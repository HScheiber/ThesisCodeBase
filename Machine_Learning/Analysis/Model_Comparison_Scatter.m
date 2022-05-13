% Script for plotting an error comparison for multiple models, to show how
% each performs.
%#ok<*UNRCH>
%% Analysis parameters: Comparison 1 C1 - C6

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'C1' 'C4' 'C6' 'C3' 'C5' 'C2'};
% Model_Labels = {'DFT Target' ...
%     'Free Charge, No Additivity.' ...
%     'Fix Charge, No Additivity.' ...
%     'Free Charge, Partial Additivity' ...
%     'Fix Charge, Partial Additivity' ...
%     'Free Charge, Additivity.' ...
%     'Fix Charge, Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models C1 - C6}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_C1-C6_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};

%% Analysis Parameters: Comparison D - G

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'D' 'E' 'F' 'G'};
% Model_Labels = {'DFT Target' ...
%     'LE:[RS] \& RLE:[WZ,NA,SP,5-5]' ...
%     'RLE:[WZ,NA,SP,5-5]' ...
%     'LE:[RS] \& RLE:[WZ,NA]' ...
%     'RLE:[WZ,NA]'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models D - G}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_D-G_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison H - K
% Salt = 'LiI';
% Theory = 'JC';
% Models = {'H' 'I' 'J' 'K'};
% Model_Labels = {'DFT Target' ...
%     'RLE:[WZ]' ...
%     'LE:[RS] \& RLE:[WZ]' ...
%     'RLE:[NA]' ...
%     'LE:[RS] \& RLE:[NA]'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models H - K}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_H-K_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison H - K not H J
% Salt = 'LiI';
% Theory = 'JC';
% Models = {'I' 'K'};
% Model_Labels = {'DFT Target' ...
%     'LE:[RS] \& RLE:[WZ]' ...
%     'LE:[RS] \& RLE:[NA]'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models I and K}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_I-K_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {};
% Plot_ca = {};
% Plot_loss = {};


%% Analysis Parameters: Comparison L - Q2

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'M' 'N' 'O' 'P' 'Q' 'Q2' };
% Model_Labels = {'DFT Target' ...
%     'Fix Charge, Additivity.*' ...
%     'Free Charge, Additivity.*' ...
%     'Fix Charge, Additivity' ...
%     'Free Charge, Additivity' ...
%     'Fix Charge, No Additivity.' ...
%     'Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models M - Q2}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_M-Q_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};

%% Analysis Parameters: Comparison R - U

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'R' 'R2' 'S' 'T'};
% Model_Labels = {'DFT Target' ...
%     'Fix Charge, Additivity.' ...
%     'Fix Charge, Additivity.*' ...
%     'Fix Charge, No Additivity.' ...
%     'Free Charge, Additivity' ...
%     'Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models R - U}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_R-U_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};


%% Analysis Parameters: Comparison V - Y

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'V' 'W' 'X' 'Y'};
% Model_Labels = {'DFT Target' ...
%     'Fix Charge, Additivity.' ...
%     'Fix Charge, No Additivity.' ...
%     'Free Charge, Additivity' ...
%     'Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models V - Y}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_V-Y_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison Z1 - Z4

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'Z1' 'Z2' 'Z3' 'Z4'};
% Model_Labels = {'DFT Target' ...
%     'Fix Charge, Partial Additivity.' ...
%     'Free Charge, Partial Additivity.' ...
%     'Fix Charge, No Additivity' ...
%     'Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models Z1 - Z4}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_Z1-Z4_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AA1 - AA1b
% 
% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AA1' 'AA1b'};
% Model_Labels = {'DFT Target' ...
%     'G[0.53 nm], Free Charge, No Additivity' ...
%     'G[0.35 nm], Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models AA1 - AA1b}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_AA1-AA1b_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AA1 - AB1

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AA1' 'AA1b' 'AA2' 'AA3' 'AA4' 'AA5' 'AB1'};
% Model_Labels = {'DFT Target' ...
%     'G[0.53 nm], Free Charge, No Additivity' ...
%     'G[0.35 nm], Free Charge, No Additivity' ...
%     'G[0.31 nm], Fix Charge, No Additivity' ...
%     'G[0.33 nm], Free Charge, Additivity' ...
%     'G[0.34 nm], Fix Charge, Additivity' ...
%     'G[0.43 nm], Fix Charge, Additivity' ...
%     'G[0.35 nm], Fix Charge, Additivity*'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models AA1 - AB1}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = 'LiI_JC_AA1-AB1_Comparison.png';
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};


%% Analysis Parameters: Comparison AC1 - AC3

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AC1' 'AC2' 'AC3'};
% Model_Labels = {'DFT Target' ...
%     'RS and NiAs Target' ...
%     'RS, NiAs, and WZ Target' ...
%     'RS, NiAs, WZ, SP, and 5-5 Target'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = ['LiI_JC_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AD1 - AF3

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AD1' 'AD2' 'AD3' 'AE1' 'AE2' 'AE3' 'AF1' 'AF2' 'AF3'};
% Model_Labels = {'DFT Target' ...
%     'Gap = 0 kJ/mol, Fix Charge, Additivity' ...
%     'Gap = 0 kJ/mol, Free Charge, Additivity' ...
%     'Gap = 0 kJ/mol, Free Charge, No Additivity' ...
%     'Gap = 2 kJ/mol, Fix Charge, Additivity' ...
%     'Gap = 2 kJ/mol, Free Charge, Additivity' ...
%     'Gap = 2 kJ/mol, Free Charge, No Additivity' ...
%     'Gap = 5 kJ/mol, Fix Charge, Additivity' ...
%     'Gap = 5 kJ/mol, Free Charge, Additivity' ...
%     'Gap = 5 kJ/mol, Free Charge, No Additivity' ...
%     };
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = ['LiI_JC_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AG1 - AG2

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AG1' 'AG2'};
% Model_Labels = {'DFT Target' ...
%     'Gaussian Perturbation' ...
%     'Large Li-Li dispersion'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = ['LiI_JC_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AH1 - AJ6

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AH1' 'AH2' 'AH3' 'AI1' 'AI2' 'AI3' 'AJ1' 'AJ2' 'AJ3' 'AJ4' 'AJ5' 'AJ6'};
% Model_Labels = {'DFT Target' ...
%     'RS Target, Free Charge, No Additivity' ...
%     'RS Target, Free Charge, Additivity' ...
%     'RS Target, Fix Charge, Additivity' ...
%     'WZ Target, Free Charge, No Additivity' ...
%     'WZ Target, Free Charge, Additivity' ...
%     'WZ Target, Fix Charge, Additivity' ...
%     'NiAs Target, Free Charge, No Additivity*' ...
%     'NiAs Target, Free Charge, Additivity*' ...
%     'NiAs Target, Fix Charge, Additivity*' ...
%     'NiAs Target, Free Charge, No Additivity' ...
%     'NiAs Target, Free Charge, Additivity' ...
%     'NiAs Target, Fix Charge, Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = ['LiI_JC_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison AH1 - AJ6

% Salt = 'LiI';
% Theory = 'JC';
% Models = {'AR1'};
% Model_Labels = {'DFT Target' ...
%     'NiAs Target, Free Charge, No Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = ['LiI_JC_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison LiBr JC Model C5a-e: reproducibility test

% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C5a' 'C5b' 'C5c' 'C5d' 'C5e'};
% Model_Labels = {'DFT Target' 'C5a' 'C5b' 'C5c' 'C5d' 'C5e'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models ' Models{1} ' - ' Models{end} '}' '. Target = RS/WZ properties'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};


%% Analysis Parameters: Comparison LiBr JC Model C5a-e
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C5matern32' 'C5e' 'C5exponential' 'C5squaredexponential' 'C5rationalquadratic'};
% Model_Labels = {'DFT Target' 'Matern-3/2' 'Matern-5/2' 'Exponential' 'Squared Exponential' 'Rational Quadratic'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models. Target = RS/WZ properties}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};

%% Analysis Parameters: Comparison LiI JC Model BF1
% Salt = 'LiI';
% Theory = 'JC';
% Models = {'BF1'};
% Model_Labels = {'DFT Target' 'BF1'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models. Target = RS/WZ/NiAs properties}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {};

%% Analysis Parameters: Comparison LiBr JC C6
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5' 'C6b1' 'C6b2' 'C6b3' 'C6b4' 'C6b5'};
% Model_Labels = {'DFT Target' 'Matern52' 'Matern52' 'Matern52' 'Matern52' 'Matern52' ...
%     'Matern32' 'Matern32' 'Matern32' 'Matern32' 'Matern32'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models. Target = RS/WZ properties}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% 
% N_Models = length(Models);
% Colours = [0 0 0; cbrewer('seq','Reds',N_Models/2); cbrewer('seq','Blues',N_Models/2)];
% full_opt = false; % Set True to use the best results from the final local optimization, false for the bayesian optimization result

%% Analysis Parameters: Comparison LiBr JC C3 vs C6 set
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C3a' 'C3b' 'C3c' 'C3d' 'C3e' 'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5'};
% Model_Labels = {'DFT Target' 'No Additivity' 'No Additivity' 'No Additivity' 'No Additivity' 'No Additivity' ...
%     'Additivity' 'Additivity' 'Additivity' 'Additivity' 'Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
% Title_txt = '';
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% 
% N_Models = length(Models);
% Colours = [0 0 0; cbrewer('seq','Reds',N_Models/2); cbrewer('seq','Blues',N_Models/2)];
% full_opt = false; % Set True to use the best results from the final local optimization, false for the bayesian optimization result

%% Analysis Parameters: Comparison LiBr JC C3 vs C6 set
Salt = 'LiBr';
Theory = 'JC';
Models = {'C5a' 'C5b' 'C5c' 'C5d' 'C5e' 'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5'};
Model_Labels = {'DFT Target' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' ...
    'Fixed Charge' 'Fixed Charge' 'Fixed Charge' 'Fixed Charge' 'Fixed Charge'};
%Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
Title_txt = '';
Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];

Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
Plot_ca = {};
Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};

N_Models = length(Models);
Colours = [0 0 0; cbrewer('seq','Reds',N_Models/2); cbrewer('seq','Blues',N_Models/2)];
full_opt = true; % Set True to use the best results from the final local optimization, false for the bayesian optimization result



%% Analysis Parameters: Comparison LiBr JC C3 vs C6 set
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C3a' 'C3b' 'C3c' 'C3d' 'C3e' 'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5'};
% Model_Labels = {'DFT Target' 'No Add. A' 'No Add. B' 'No Add. C' 'No Add. D' 'No Add. E' ...
%     'Add. A' 'Add. B' 'Add. C' 'Add. D' 'Add. E'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% 
% N_Models = length(Models);
% Colours = [0 0 0; cbrewer('seq','Reds',N_Models/2); cbrewer('seq','Blues',N_Models/2)];
% full_opt = true; % Set True to use the best results from the final local optimization, false for the bayesian optimization result


%% Analysis Parameters: Comparison of acquisition functions with LiBr JC
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C5_EIa' 'C5_EIb' 'C5_EIc' ...
%           'C5_EIP0.5a' 'C5_EIP0.5b' 'C5_EIP0.5c' ...
%           'C5_EIP1a' 'C5_EIP1b' 'C5_EIP1c' ...
%           'C5_EIP2a' 'C5_EIP2b' 'C5_EIP2c' ...
%           'C5_LCB1a' 'C5_LCB1b' 'C5_LCB1c' ...
%           'C5_LCB2a' 'C5_LCB2b' 'C5_LCB2c' ...
%           'C5_LCB3a' 'C5_LCB3b' 'C5_LCB3c' ...
%           'C5_POIa' 'C5_POIb' 'C5_POIc'};
% Model_Labels = {'DFT Target' ...
%     'Exp. Improvement' 'Exp. Improvement' 'Exp. Improvement' ...
%     'Exp. Improv. Plus(0.5)' 'Exp. Improv. Plus(0.5)' 'Exp. Improv. Plus(0.5)' ...
%     'Exp. Improv. Plus(1)' 'Exp. Improv. Plus(1)' 'Exp. Improv. Plus(1)' ...
%     'Exp. Improv. Plus(2)' 'Exp. Improv. Plus(2)' 'Exp. Improv. Plus(2)' ...
%     'LCB ($\kappa = 1$)' 'LCB ($\kappa = 1$)' 'LCB ($\kappa = 1$)' ...
%     'LCB ($\kappa = 2$)' 'LCB ($\kappa = 2$)' 'LCB ($\kappa = 2$)' ...
%     'LCB ($\kappa = 3$)' 'LCB ($\kappa = 3$)' 'LCB ($\kappa = 3$)' ...
%     'Prob. of Improv.' 'Prob. of Improv.' 'Prob. of Improv.'};
% %Title_txt = ['\bf{Acqusition Function Comparison for ' Salt ' Models, Matern-5/2 Kernel}'];
% Title_txt = '';
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% 
% N_Models = length(Models);
% Colours = [0 0 0; repelem(cbrewer('qual','Set1',N_Models/3),3,1)];
% 
% full_opt = false; % Set True to use the best results from the final local optimization, false for the bayesian optimization result


%% Analysis Parameters: Comparison of kernel functions with LiBr JC
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'C5matern32' 'C5e' 'C5exponential' 'C5squaredexponential' 'C5rationalquadratic'};
% Model_Labels = {'DFT Target' 'Matern-3/2' 'Matern-5/2' 'Exponential' 'Squared Exponential' 'Rational Quadratic'};
% %Title_txt = ['\bf{Acqusition Function Comparison for ' Salt ' Models, Matern-5/2 Kernel}'];
% Title_txt = '';
% Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive'};
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Comparison.png'];
% 
% Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
% Plot_ca = {};
% Plot_loss = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
% 
% N_Models = length(Models);
% Colours = [0 0 0; cbrewer('qual','Set1',N_Models)];
% 
% full_opt = false; % Set True to use the best results from the final local optimization, false for the bayesian optimization result



%% Other parameters
savefile = false; % switch to save the final plots to file
show_as_percent_error = true; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)

% Calculates error in lattice energy with respect to experimental lattice energies when available
Target_Experimental_Energies = true; 
% Calculates error in lattice parameter with respect to experimental lattice parameters when available
Target_Experimental_latpar = true; 

fs = 22; % font size
%% Pre-allocate data arrays
N_Structures = length(Structures);

Energies = nan(N_Models,N_Structures);
Rel_Energies = nan(N_Models,N_Structures);
LatPars_a = nan(N_Models,N_Structures);
LatPars_c = nan(N_Models,N_Structures);
LatPars_ca = nan(N_Models,N_Structures);
Total_loss = nan(N_Models,1);

% Plot color scheme
%Colours = [0 0 0; cbrewer('qual','Set3',N_Models)];
%Colours = [0 0 0; cbrewer('qual','Paired',N_Models)];


% Load target data
DFT = Load_Best_DFT_Data;

if Target_Experimental_Energies || Target_Experimental_latpar
    EXP = Load_Experimental_Data;
    
    if Target_Experimental_Energies
        Correction = EXP.LiI.Rocksalt.E - DFT.(Salt).Rocksalt.Energy;
        for idx = 1:length(Structures)
            DFT.(Salt).(Structures{idx}).Energy = DFT.(Salt).(Structures{idx}).Energy + Correction;
        end
    end
    
    if Target_Experimental_latpar
        DFT.(Salt).Rocksalt.a = EXP.(Salt).Rocksalt.a;
        DFT.(Salt).Rocksalt.b = EXP.(Salt).Rocksalt.b;
        DFT.(Salt).Rocksalt.c = EXP.(Salt).Rocksalt.c;
        if ~isnan(EXP.(Salt).Wurtzite.a)
            DFT.(Salt).Wurtzite.a = EXP.(Salt).Wurtzite.a;
            DFT.(Salt).Wurtzite.b = EXP.(Salt).Wurtzite.b;
            DFT.(Salt).Wurtzite.c = EXP.(Salt).Wurtzite.c;
        end
    end
    clear('EXP')
end

% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

% Loop over chosen models
for idx = 1:N_Models
    Model = Models{idx};
    
    if full_opt
        % Find the fully optimized model
        sres = dir([ML_results_dir filesep '*' Salt '*' Theory '*' 'Model_' Model '_*fullopt.mat']);

        if length(sres) > 1
            warning(['Multiple results found for: ' Salt ', ' Theory ', Model ' Model '. Using first result.']);
            sres = sres(1);
        elseif isempty(sres)
            error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
        end

        model_filename = fullfile(sres.folder,sres.name);


        % Load data
        data = load(model_filename);
        try
            Minimization_Data = data.Minimization_Data;
            Loss = data.loss;
            clear('data')
        catch
            warning(['No saved crystal data for ' Salt ' ' Theory ' Model ' Model])
            continue
        end
    else
        % Find the bayesian optimized model
        sres = dir([ML_results_dir filesep '*' Salt '*' Theory '*' 'Model_' Model '_*bayesopt.mat']);

        if length(sres) > 1
            warning(['Multiple results found for: ' Salt ', ' Theory ', Model ' Model '. Using first result.']);
            sres = sres(1);
        elseif isempty(sres)
            error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
        end

        model_filename = fullfile(sres.folder,sres.name);
        
        % Load data
        data = load(model_filename);
        try
            best_idx = data.results.IndexOfMinimumTrace(end);
            Minimization_Data = data.results.UserDataTrace{best_idx};
            Loss = data.results.MinObjective;
            clear('data')
        catch
            warning(['No saved crystal data for ' Salt ' ' Theory ' Model ' Model])
            continue
        end
    end
    
    % Grab available structures from data
    Min_Dat_Structs = cell(1,length(Minimization_Data));
    for jdx = 1:length(Minimization_Data)
        Min_Dat_Structs{jdx} = Minimization_Data{jdx}.Structure;
    end
    
    % Loop through structures
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        % Absolute energies
        if contained_in_cell(Structure,Plot_Energies)
            if show_as_percent_error
                Energies(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                    DFT.(Salt).(Structure).Energy)*100/ ...
                    abs(DFT.(Salt).(Structure).Energy);
            else
                Energies(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                    DFT.(Salt).(Structure).Energy; 
            end
        end
        
        % Relative Energies
        if contained_in_cell(Structure,Plot_Rel_Energies)
            Rel_Energies(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                                    Minimization_Data{strcmpi(Min_Dat_Structs,'Rocksalt')}.E;
        end
        
        % lattice param a
        if contained_in_cell(Structure,Plot_a)
            if show_as_percent_error
                LatPars_a(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a - DFT.(Salt).(Structure).a)*100/ ...
                    DFT.(Salt).(Structure).a;
            else
                LatPars_a(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a - DFT.(Salt).(Structure).a;
            end
        end
        
        % lattice param c
        if contained_in_cell(Structure,Plot_c)
            if show_as_percent_error
                LatPars_c(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c - DFT.(Salt).(Structure).c)*100/ ...
                    DFT.(Salt).(Structure).c;
            else
                LatPars_c(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c - DFT.(Salt).(Structure).c;
            end
        end
        
        % lattice param c/a
        if contained_in_cell(Structure,Plot_ca)
            Target_ca = DFT.(Salt).(Structure).c / DFT.(Salt).(Structure).a;
            
            Model_ca  = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c / ...
                         Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a;
            
            if show_as_percent_error
                LatPars_ca(idx,jdx) = (Model_ca - Target_ca)*100 / Target_ca;
            else
                LatPars_ca(idx,jdx) = Model_ca - Target_ca;
            end
        end
    end
    
    % total loss
    if contained_in_cell(Structure,Plot_loss)
        Total_loss(idx) = Loss;
    end
end

% Add in DFT results
Rel_Energies_DFT = nan(1,N_Structures);
for jdx = 1:length(Structures)
    Structure = Structures{jdx};
    
    if contained_in_cell(Structure,Plot_Rel_Energies)
        Rel_Energies_DFT(jdx) = DFT.(Salt).(Structure).Energy - DFT.(Salt).Rocksalt.Energy;
    end
end

N_compares = 0;
if ~isempty(Plot_Energies)
    N_compares = N_compares+1;
end
if ~isempty(Plot_Rel_Energies)
    N_compares = N_compares+1;
end
if ~isempty(Plot_a)
    N_compares = N_compares+1;
end
if ~isempty(Plot_c)
    N_compares = N_compares+1;
end
if ~isempty(Plot_ca)
    N_compares = N_compares+1;
end
if ~isempty(Plot_loss)
    N_compares = N_compares+1;
end

Energies = [nan(1,N_Structures); Energies];
Rel_Energies = [Rel_Energies_DFT; Rel_Energies];
LatPars_a = [nan(1,N_Structures); LatPars_a];
LatPars_c = [nan(1,N_Structures); LatPars_c];
LatPars_ca = [nan(1,N_Structures); LatPars_ca];

all_y = {Total_loss, Energies, Rel_Energies, LatPars_a, LatPars_c, LatPars_ca};
bar_x = 1:N_Structures;
if full_opt
    t1 = 'Final Optimized Objective';
else
    t1 = 'Bayesian Optimized Objective';
end
Titles = {t1 ...
    'Error in $E_L$' ...
    '$E_{L}$[Structure] - $E_{L}$[Rocksalt]' ...
    'Error in $a$' ...
    'Error in $c$' ...
    'Error in $c/a$'};

if show_as_percent_error
    ylabs = {'log$(f(x) + 1)$' ...
        '[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
else
    ylabs = {'log$(f(x) + 1)$' ...
        '[kJ mol$^{-1}$]' ...
        '[kJ mol$^{-1}$]' ...
        '[\AA]' ...
        '[\AA]' ...
        '$(c/a)$'};
end


% Create figure and axis
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
t = tiledlayout(figh,N_compares,1,'TileSpacing','tight');

for idx = 1:length(all_y)
    bar_y = all_y{idx};
    if all(isnan(bar_y(2:end,:)),'all')
        continue
    end
    cTitle = Titles{idx};
    ylab = ylabs{idx};
    
    axobj = nexttile;
    if size(bar_y,2) == 1
        p = bar(axobj,bar_y,'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        p.CData = Colours(2:end,:);
        
        xticks(axobj,1:N_Models)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Models+0.5])
        xticklabels(axobj,[])
        %set(axobj, 'YScale', 'log')
    else
        p = bar(axobj,bar_x,bar_y','FaceColor','flat','Visible','on','BarWidth',1);
        % Apply colours and alpha
        for jdx=1:length(p)
            p(jdx).CData = Colours(jdx,:);
            p(jdx).LineWidth = 2;
        end
        p(1).EdgeColor = 'r';
        
        xticks(axobj,1:N_Structures)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        xticklabels(axobj,[])
        ylim(axobj,'padded') 
    end
    
    title(axobj, cTitle,'Fontsize',fs,'Interpreter','latex')
    grid(axobj,'minor');
    axobj.YGrid = 'on';
    axobj.XGrid = 'off';
    ylabel(axobj, ylab,'Fontsize',fs-1,'Interpreter','latex')
end


title(t, Title_txt,'Fontsize',fs+5,'Interpreter','latex')
%xlabel(t, 'Crystal Structure','Fontsize',fs,'Interpreter','latex')
xticklabels(axobj,strrep(Structures,'FiveFive','5-5'))


% legh =legend(axobj,Model_Labels,'Location','SouthOutside','Orientation','Horizontal',...
%     'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', 4);
legh = legend(axobj,p([1 6 11]),Model_Labels{[1 6 11]},'Location','EastOutside','Orientation','Horizontal',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', 1);

legh.Layout.Tile = 'East';
title(legh,'\bf{\underline{Acqusition Functions}}','Fontsize',fs,'Interpreter','latex')

% if show_as_percent_error
%     ylabel(t, 'Relative Error [\%]','Fontsize',fs,'Interpreter','latex')
% end

if savefile
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
end
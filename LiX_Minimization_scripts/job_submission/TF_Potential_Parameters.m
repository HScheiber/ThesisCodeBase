% Output units: energies are kJ/mol and distances are nm
function [OutputMX,OutputMM,OutputXX,S] = TF_Potential_Parameters(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'TF_Potential_Parameters';
addOptional(p,'PlotType','full')
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

%% Conversion factors
nm_per_Ang = 0.1; % nm per Angstrom
kj_per_erg = 1e-10; % kJ per erg
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_cm = 1e+7; % nm per cm
NA = 6.0221409e23; % Molecules per mole

%% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: MX = +-   MM = ++     XX = --
C.LiF.MM  = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.MX  = (0.8*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.XX  = (14.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiCl.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.MX = (2.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.XX = (111*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiBr.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.XX = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiI.MM = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.MX = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.XX = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaF.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.MX  = (4.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.XX  = (16.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaCl.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.MX  = (11.2*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.XX  = (116*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaBr.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.MX  = (14.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.XX  = (196*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaI.MM  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.MX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.XX  = (392*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KF.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.MX  = (19.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.XX  = (18.6*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KCl.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.MX  = (48*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.XX  = (124.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KBr.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.MX  = (60*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.XX  = (206*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KI.MM  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.MX  = (82*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.XX  = (403*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbF.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.MX  = (31*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.XX  = (18.9*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbCl.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.MX  = (79*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.XX  = (130*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbBr.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.MX  = (99*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.XX  = (215*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbI.MM  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.MX  = (135*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.XX  = (428*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsF.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.MX  = (52*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.XX  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsCl.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsCl.MX  = (129*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsCl.XX  = (129*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsBr.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsBr.MX  = (163*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsBr.XX  = (214*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsI.MM  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsI.MX  = (224*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsI.XX  = (424*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

%% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: MX = +-   MM = ++     XX = --
D.LiF.MM  = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.MX   = (0.6*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.XX    = (17*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiCl.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.MX = (2.4*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.XX   = (223*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiBr.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.XX = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiI.MM = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.MX = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.XX = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaF.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.MX  = (3.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.XX  = (20*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaCl.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.MX  = (13.9*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.XX  = (233*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaBr.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.MX  = (19*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.XX  = (450*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaI.MM  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.MX  = (31*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.XX  = (1100*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KF.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.MX  = (21*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.XX  = (22*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KCl.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.MX  = (73*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.XX  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KBr.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.MX  = (99*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.XX  = (470*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KI.MM  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.MX  = (156*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.XX  = (1130*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbF.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.MX  = (40*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbCl.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.MX  = (134*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.XX  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbBr.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.MX  = (180*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.XX  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbI.MM  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.MX  = (280*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.XX  = (1200*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsF.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.MX  = (78*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.XX  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsCl.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsCl.MX  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsCl.XX  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsBr.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsBr.MX  = (340*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsBr.XX  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsI.MM  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsI.MX  = (520*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsI.XX  = (1190*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

%% TF Repulsive Size Parameter sigma (AKA r+/-): M = +   X = -
% Metals
sigma.Li = 0.816*nm_per_Ang; % nm
sigma.Na = 1.170*nm_per_Ang; % nm
sigma.K  = 1.463*nm_per_Ang; % nm
sigma.Rb = 1.587*nm_per_Ang; % nm
sigma.Cs = 1.720*nm_per_Ang; % nm

% Halides
sigma.F  = 1.179*nm_per_Ang; % nm
sigma.Cl = 1.585*nm_per_Ang; % nm
sigma.Br = 1.716*nm_per_Ang; % nm
sigma.I  = 1.907*nm_per_Ang; % nm

%% TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
% Metals
valence.Li = 2;
valence.Na = 8;
valence.K = 8;
valence.Rb = 8;
valence.Cs = 8;

% Halides
valence.F = 8;
valence.Cl = 8;
valence.Br = 8;
valence.I = 8;

%% TF Hardness Parameter Rho
rho.LiF = 0.299*nm_per_Ang; % nm
rho.LiCl = 0.342*nm_per_Ang; % nm
rho.LiBr = 0.353*nm_per_Ang; % nm
rho.LiI = 0.430*nm_per_Ang; % nm

rho.NaF = 0.330*nm_per_Ang; % nm
rho.NaCl = 0.317*nm_per_Ang; % nm
rho.NaBr = 0.340*nm_per_Ang; % nm
rho.NaI = 0.386*nm_per_Ang; % nm

rho.KF = 0.338*nm_per_Ang; % nm
rho.KCl = 0.337*nm_per_Ang; % nm
rho.KBr = 0.335*nm_per_Ang; % nm
rho.KI = 0.355*nm_per_Ang; % nm

rho.RbF = 0.328*nm_per_Ang; % nm
rho.RbCl = 0.318*nm_per_Ang; % nm
rho.RbBr = 0.335*nm_per_Ang; % nm
rho.RbI = 0.337*nm_per_Ang; % nm

rho.CsF = 0.282*nm_per_Ang; % nm
rho.CsCl = 0.272*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsBr = 0.289*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsI = 0.291*nm_per_Ang; % nm THIS IS ESTIMATED

%% TF Parameter: q (charge)
q.Li =  Settings.S.Q; % atomic
q.Na =  Settings.S.Q; % atomic
q.K  =  Settings.S.Q; % atomic
q.Rb =  Settings.S.Q; % atomic
q.Cs =  Settings.S.Q; % atomic

q.F  = -Settings.S.Q; % atomic
q.Cl = -Settings.S.Q; % atomic
q.Br = -Settings.S.Q; % atomic
q.I  = -Settings.S.Q; % atomic

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2/valence.(Metal); % Unitless
beta.MX = 1 + 1/valence.(Metal) - 1/valence.(Halide); % Unitless
beta.XX = 1 - 2/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All*Settings.S.A.MM/rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All*Settings.S.A.MX/rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All*Settings.S.A.XX/rho.(Settings.Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All*Settings.S.R.MM*beta.MM*b*exp(2*sigma.(Metal)/rho.(Settings.Salt));
B.MX = Settings.S.R.All*Settings.S.R.MX*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Settings.Salt));
B.XX = Settings.S.R.All*Settings.S.R.XX*beta.XX*b*exp(2*sigma.(Halide)/rho.(Settings.Salt));

%% Scale Dispersion
C.(Settings.Salt).MM = Settings.S.D6D.All*Settings.S.D6D.MM*Settings.S.D.All*Settings.S.D.MM.*C.(Settings.Salt).MM;
D.(Settings.Salt).MM = Settings.S.D8D.All*Settings.S.D8D.MM*Settings.S.D.All*Settings.S.D.MM.*D.(Settings.Salt).MM;

C.(Settings.Salt).XX = Settings.S.D6D.All*Settings.S.D6D.XX*Settings.S.D.All*Settings.S.D.XX.*C.(Settings.Salt).XX;
D.(Settings.Salt).XX = Settings.S.D8D.All*Settings.S.D8D.XX*Settings.S.D.All*Settings.S.D.XX.*D.(Settings.Salt).XX;

C.(Settings.Salt).MX = Settings.S.D6D.All*Settings.S.D6D.MX*Settings.S.D.All*Settings.S.D.MX.*C.(Settings.Salt).MX;
D.(Settings.Salt).MX = Settings.S.D8D.All*Settings.S.D8D.MX*Settings.S.D.All*Settings.S.D.MX.*D.(Settings.Salt).MX;

%% Outputs
OutputMM.B = B.MM; % Repulsion parameter
OutputMM.C = C.(Settings.Salt).MM; % C6 dispersion constant
OutputMM.D = D.(Settings.Salt).MM; % C8 dispersion constant
OutputMM.alpha = alpha.MM; % Exponential steepness parameter

OutputXX.B = B.XX; % Repulsion parameter
OutputXX.C = C.(Settings.Salt).XX; % C6 dispersion constant
OutputXX.D = D.(Settings.Salt).XX; % C8 dispersion constant
OutputXX.alpha = alpha.XX; % Exponential steepness parameter

OutputMX.B = B.MX; % Repulsion parameter
OutputMX.C = C.(Settings.Salt).MX; % C6 dispersion constant
OutputMX.D = D.(Settings.Salt).MX; % C8 dispersion constant
OutputMX.alpha = alpha.MX; % Exponential steepness parameter

% Convert to sigma-epsilon form
OutputMM.sigma = sqrt((4/3)*(D.(Settings.Salt).MM/C.(Settings.Salt).MM));
OutputXX.sigma = sqrt((4/3)*(D.(Settings.Salt).XX/C.(Settings.Salt).XX));
OutputMX.sigma = sqrt((4/3)*(D.(Settings.Salt).MX/C.(Settings.Salt).MX));

OutputMM.gamma = alpha.MM*OutputMM.sigma;
OutputXX.gamma = alpha.XX*OutputXX.sigma;
OutputMX.gamma = alpha.MX*OutputMX.sigma;

OutputMM.epsilon = B.MM*exp(-OutputMM.gamma)*(1 - (7/48)*OutputMM.gamma);
OutputXX.epsilon = B.XX*exp(-OutputXX.gamma)*(1 - (7/48)*OutputXX.gamma);
OutputMX.epsilon = B.MX*exp(-OutputMX.gamma)*(1 - (7/48)*OutputMX.gamma);


S = Settings.S;
% Alpha (nm^-1)
S.A.MM = alpha.MM;
S.A.XX = alpha.XX;
S.A.MX = alpha.MX;

% Repulsive prefactor param % kJ/mol
S.R.MM = B.MM;
S.R.XX = B.XX;
S.R.MX = B.MX;

% C6 param
S.D6D.MM = C.(Settings.Salt).MM;
S.D6D.XX = C.(Settings.Salt).XX;
S.D6D.MX = C.(Settings.Salt).MX;

% C8 param
S.D8D.MM = D.(Settings.Salt).MM;
S.D8D.XX = D.(Settings.Salt).XX;
S.D8D.MX = D.(Settings.Salt).MX;

if Plotswitch
    figure;
    % Options
    lw=2;
    fs=25;
    
    %% Conversion Factors / Constants
    m_per_nm = 1e+9; % nm per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    J_per_kJ = 1000;
    epsilon_0 = (8.854187817620e-12)*J_per_kJ/(m_per_nm*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in C^-2 nm kJ mol^-1
    
    q.(Metal) = Settings.S.Q; % charge of M in C
    q.(Halide) = -Settings.S.Q; % charge of X in C
    
    r = Startpoint:Settings.Table_StepSize:Settings.Table_Length; % r vector
    
    U.MX.f = 1./r; % Electrostatics function f(r)
    U.MX.g = - C.(Settings.Salt).MX./(r.^6) - D.(Settings.Salt).MX./(r.^8); % Dispersion g(r)
    U.MX.h = B.MX*exp(-alpha.MX.*r); % Short range repulsion h(r)
    
    U.MX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MX.dg = - 6*C.(Settings.Salt).MX./(r.^7) - 8*D.(Settings.Salt).MX./(r.^9); % Dispersion -g(r)
    U.MX.dh = alpha.MX*B.MX*exp(-alpha.MX.*r); % Short range repulsion -h(r)
    
    U.XX.f = 1./r; % Electrostatics function f(r)
    U.XX.g = - C.(Settings.Salt).XX./(r.^6) - D.(Settings.Salt).XX./(r.^8); % Dispersion g(r)
    U.XX.h = B.XX*exp(-alpha.XX.*r); % Short range repulsion h(r)
    
    U.XX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.XX.dg = - 6*C.(Settings.Salt).XX./(r.^7) - 8*D.(Settings.Salt).XX./(r.^9); % Dispersion -g(r)
    U.XX.dh = alpha.XX*B.XX*exp(-alpha.XX.*r); % Short range repulsion -h(r)
    
    U.MM.f = 1./r; % Electrostatics function f(r)
    U.MM.g = - C.(Settings.Salt).MM./(r.^6) - D.(Settings.Salt).MM./(r.^8); % Dispersion g(r)
    U.MM.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion h(r)
    
    U.MM.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MM.dg = - 6*C.(Settings.Salt).MM./(r.^7) - 8*D.(Settings.Salt).MM./(r.^9); % Dispersion -g(r)
    U.MM.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion -h(r)
    
    h = cell(1,3);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'lj'
            h{1} = plot(r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g + U.MM.h,'Color','g','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g + U.XX.h,'Color','b','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Repulsive Potential';
        case 'coulomb'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Coulomb Potential';
    end
    
    title(['Plot of ' ttxt ' for ' Settings.Salt ' TF Model'],...
       'Interpreter','latex','fontsize',fs)
    
    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('Separation [\AA]','fontsize',fs,'Interpreter','latex');
    ylabel('Potential Energy [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');
    
    ylim(yl);
    xlim([Startpoint Settings.Table_Length]);
    
    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{[Metal '$^{+}$' ' - ' Halide '$^{-}$'] [Metal '$^{+}$' ' - ' Metal '$^{+}$'] [Halide '$^{-}$' ' - ' Halide '$^{-}$']});
    leg1.Interpreter = 'latex';
end

end
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND
% kJ/mol
% DS = dispersion scaling
% ES = Epsilon scaling (for JC only)
function Output = TF_Pair_Potential(Startpoint,Endpoint,Spacing,Salt,Scaling_Params,Potential_Type)

S_D = Scaling_Params(1);
S_R = Scaling_Params(2);
S_MMD = Scaling_Params(5);
S_XXD = Scaling_Params(6);
S_MXD = Scaling_Params(7);
S_C = Scaling_Params(8); % Scaling charges
S_D6 = Scaling_Params(9); % Scaling R6
S_D8 = Scaling_Params(10); % Scaling R8

S_MMD6 = S_D*S_D6*S_MMD;
S_XXD6 = S_D*S_D6*S_XXD;
S_MXD6 = S_D*S_D6*S_MXD;

S_MMD8 = S_D*S_D8*S_MMD;
S_XXD8 = S_D*S_D8*S_XXD;
S_MXD8 = S_D*S_D8*S_MXD;

S_MMR = S_R;
S_XXR = S_R;
S_MXR = S_R;

%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
Ang_per_cm = 1e+8; % Angstroms per cm
Ang_per_m = 1e+10; % Angstroms per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Angstrom^-1
k_0 = S_C/(4*pi*epsilon_0); % Coulomb constant in kJ Angstrom C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_Ang = 0.1; % nm per Angstrom

%% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: PM = +-   PP = ++     MM = --
C.LiF.MM  = (0.073*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiF.MX  = (0.8*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiF.XX  = (14.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.LiCl.MM = (0.073*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiCl.MX = (2.0*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiCl.XX   = (111*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.LiBr.MM = (0.073*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiBr.MX = (2.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiBr.XX = (185*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.LiI.MM = (0.073*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiI.MX = (3.3*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.LiI.XX = (378*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.NaF.MM  = (1.68*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaF.MX  = (4.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaF.XX  = (16.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.NaCl.MM  = (1.68*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaCl.MX  = (11.2*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaCl.XX  = (116*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.NaBr.MM  = (1.68*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaBr.MX  = (14.0*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaBr.XX  = (196*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.NaI.MM  = (1.68*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaI.MX  = (19.1*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.NaI.XX  = (392*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.KF.MM  = (24.3*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KF.MX  = (19.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KF.XX  = (18.6*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.KCl.MM  = (24.3*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KCl.MX  = (48*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KCl.XX  = (124.5*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.KBr.MM  = (24.3*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KBr.MX  = (60*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KBr.XX  = (206*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.KI.MM  = (24.3*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KI.MX  = (82*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.KI.XX  = (403*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.RbF.MM  = (59.4*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbF.MX  = (31*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbF.XX  = (18.9*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.RbCl.MM  = (59.4*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbCl.MX  = (79*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbCl.XX  = (130*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.RbBr.MM  = (59.4*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbBr.MX  = (99*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbBr.XX  = (215*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.RbI.MM  = (59.4*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbI.MX  = (135*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.RbI.XX  = (428*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

C.CsF.MM  = (152*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.CsF.MX  = (52*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6
C.CsF.XX  = (19.1*C_unit)*kj_per_erg*(Ang_per_cm^6)*NA; % kJ/mol Angstrom^6

%% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: PM = +-   PP = ++     MM = --
D.LiF.MM  = (0.03*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiF.MX   = (0.6*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiF.XX    = (17*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.LiCl.MM = (0.03*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiCl.MX = (2.4*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiCl.XX   = (223*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.LiBr.MM = (0.03*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiBr.MX = (3.3*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiBr.XX = (423*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.LiI.MM = (0.03*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiI.MX = (5.3*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.LiI.XX = (1060*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.NaF.MM  = (0.8*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaF.MX  = (3.8*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaF.XX  = (20*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.NaCl.MM  = (0.8*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaCl.MX  = (13.9*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaCl.XX  = (233*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.NaBr.MM  = (0.8*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaBr.MX  = (19*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaBr.XX  = (450*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.NaI.MM  = (0.8*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaI.MX  = (31*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.NaI.XX  = (1100*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.KF.MM  = (24*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KF.MX  = (21*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KF.XX  = (22*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.KCl.MM  = (24*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KCl.MX  = (73*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KCl.XX  = (250*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.KBr.MM  = (24*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KBr.MX  = (99*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KBr.XX  = (470*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.KI.MM  = (24*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KI.MX  = (156*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.KI.XX  = (1130*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.RbF.MM  = (82*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbF.MX  = (40*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbF.XX  = (23*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.RbCl.MM  = (82*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbCl.MX  = (134*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbCl.XX  = (260*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.RbBr.MM  = (82*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbBr.MX  = (180*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbBr.XX  = (490*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.RbI.MM  = (82*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbI.MX  = (280*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.RbI.XX  = (1200*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

D.CsF.MM  = (278*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.CsF.MX  = (78*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8
D.CsF.XX  = (23*D_unit)*kj_per_erg*(Ang_per_cm^8)*NA; % kJ/mol Angstrom^8

%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
% Metals
sigma.Li = 0.816; % Angstrom
sigma.Na = 1.170; % Angstrom
sigma.K  = 1.463; % Angstrom
sigma.Rb = 1.587; % Angstrom
sigma.Cs = 1.720; % Angstrom

% Halides
sigma.F  = 1.179; % Angstrom
sigma.Cl = 1.585; % Angstrom
sigma.Br = 1.716; % Angstrom
sigma.I  = 1.907; % Angstrom

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
rho.LiF = 0.299; % Angstrom
rho.LiCl = 0.342; % Angstrom
rho.LiBr = 0.353; % Angstrom
rho.LiI = 0.430; % Angstrom

rho.NaF = 0.330; % Angstrom
rho.NaCl = 0.317; % Angstrom
rho.NaBr = 0.340; % Angstrom
rho.NaI = 0.386; % Angstrom

rho.KF = 0.338; % Angstrom
rho.KCl = 0.337; % Angstrom
rho.KBr = 0.335; % Angstrom
rho.KI = 0.355; % Angstrom

rho.RbF = 0.328; % Angstrom
rho.RbCl = 0.318; % Angstrom
rho.RbBr = 0.335; % Angstrom
rho.RbI = 0.337; % Angstrom

rho.CsF = 0.282; % Angstrom

%% TF Parameter: q (charge)
q.Li =  1; % atomic
q.Na =  1; % atomic
q.K  =  1; % atomic
q.Rb =  1; % atomic
q.Cs =  1; % atomic

q.F  = -1; % atomic
q.Cl = -1; % atomic
q.Br = -1; % atomic
q.I  = -1; % atomic

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Generate range (r) in Angstrom
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Calculate Pauling Coefficients beta: PM = +-   PP = ++     MM = --
beta.MM = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.MX = 1 + q.(Metal)/valence.(Metal) + q.(Halide)/valence.(Halide); % Unitless
beta.XX = 1 + 2*q.(Halide)/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: PM = +-   PP = ++     MM = --
alpha.MM = 1/rho.(Salt); % Angstrom^-1
alpha.MX = 1/rho.(Salt); % Angstrom^-1
alpha.XX = 1/rho.(Salt); % Angstrom^-1

%% Calculate TF Repulsive Scaling Parameter B: PM = +-   PP = ++     MM = --
B.MM = S_MMR.*beta.MM*b*exp(2*sigma.(Metal)/rho.(Salt));
B.MX = S_MXR.*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Salt));
B.XX = S_XXR.*beta.XX*b*exp(2*sigma.(Halide)/rho.(Salt));

%% Build TF PES for LiF with dispersion scaling
if strcmpi(Potential_Type,'full')
    Output.MX = k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
        + B.MX*exp(-alpha.MX.*r) ...
        - S_MXD6.*C.(Salt).MX./(r.^6) ...
        - S_MXD8.*D.(Salt).MX./(r.^8);

    Output.MM = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
        + B.MM*exp(-alpha.MM.*r) ...
        - S_MMD6.*C.(Salt).MM./(r.^6) ...
        - S_MMD8.*D.(Salt).MM./(r.^8);

    Output.XX = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
        + B.XX*exp(-alpha.XX.*r) ...
        - S_XXD6.*C.(Salt).XX./(r.^6)...
        - S_XXD8.*D.(Salt).XX./(r.^8);

elseif strcmpi(Potential_Type,'vdw')    
    Output.MX = B.MX*exp(-alpha.MX.*r) ...
        - S_MXD6.*C.(Salt).MX./(r.^6) ...
        - S_MXD8.*D.(Salt).MX./(r.^8);

    Output.MM = B.MM*exp(-alpha.MM.*r) ...
        - S_MMD6.*C.(Salt).MM./(r.^6) ...
        - S_MMD8.*D.(Salt).MM./(r.^8);

    Output.XX = B.XX*exp(-alpha.XX.*r) ...
        - S_XXD6.*C.(Salt).XX./(r.^6)...
        - S_XXD8.*D.(Salt).XX./(r.^8);
elseif strcmpi(Potential_Type,'dispersion')
    Output.MX = - S_MXD6.*C.(Salt).MX./(r.^6) ...
        - S_MXD8.*D.(Salt).MX./(r.^8);

    Output.MM = - S_MMD6.*C.(Salt).MM./(r.^6) ...
        - S_MMD8.*D.(Salt).MM./(r.^8);

    Output.XX = - S_XXD6.*C.(Salt).XX./(r.^6)...
        - S_XXD8.*D.(Salt).XX./(r.^8);
elseif strcmpi(Potential_Type,'repulsion')
    Output.MX = B.MX*exp(-alpha.MX.*r);
    
    Output.MM = B.MM*exp(-alpha.MM.*r);

    Output.XX = B.XX*exp(-alpha.XX.*r);
end
Output.r = r;
end
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND
% kJ/mol
% DS = dispersion scaling
% ES = Epsilon scaling (for JC only)
function Data = Pair_Potential_Saver(Startpoint,Endpoint,Spacing,Salt,Watermodel,Dispersion_Scale,Epsilon_Scale)
%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
Ang_per_cm = 1e+8; % Angstroms per cm
nm_per_cm = 1e+7; % nm per cm
Ang_per_m = 1e+10; % Angstroms per m
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1
C_unit = 1e-60; % erg cm^6
D_unit = 1e-76; % erg cm^8
nm_per_Ang = 0.1; % nm per Angstrom

%% Huggins-Mayer Dipole-Dipole Dispersion Parameter C: PM = +-   PP = ++     MM = --
C.LiF.PP  = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.PM  = (0.8*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiF.MM  = (14.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiCl.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.PM = (2.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiCl.MM   = (111*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiBr.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.PM = (2.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiBr.MM = (185*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.LiI.PP = (0.073*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.PM = (3.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.LiI.MM = (378*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaF.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.PM  = (4.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaF.MM  = (16.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaCl.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.PM  = (11.2*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaCl.MM  = (116*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaBr.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.PM  = (14.0*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaBr.MM  = (196*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.NaI.PP  = (1.68*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.PM  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.NaI.MM  = (392*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KF.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.PM  = (19.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KF.MM  = (18.6*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KCl.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.PM  = (48*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KCl.MM  = (124.5*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KBr.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.PM  = (60*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KBr.MM  = (206*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.KI.PP  = (24.3*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.PM  = (82*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.KI.MM  = (403*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbF.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.PM  = (31*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbF.MM  = (18.9*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbCl.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.PM  = (79*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbCl.MM  = (130*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbBr.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.PM  = (99*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbBr.MM  = (215*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.RbI.PP  = (59.4*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.PM  = (135*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.RbI.MM  = (428*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

C.CsF.PP  = (152*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.PM  = (52*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6
C.CsF.MM  = (19.1*C_unit)*kj_per_erg*(nm_per_cm^6)*NA; % kJ/mol nm^6

%% Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: PM = +-   PP = ++     MM = --
D.LiF.PP  = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.PM   = (0.6*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiF.MM    = (17*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiCl.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.PM = (2.4*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiCl.MM   = (223*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiBr.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.PM = (3.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiBr.MM = (423*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.LiI.PP = (0.03*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.PM = (5.3*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.LiI.MM = (1060*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaF.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.PM  = (3.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaF.MM  = (20*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaCl.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.PM  = (13.9*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaCl.MM  = (233*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaBr.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.PM  = (19*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaBr.MM  = (450*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.NaI.PP  = (0.8*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.PM  = (31*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.NaI.MM  = (1100*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KF.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.PM  = (21*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KF.MM  = (22*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KCl.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.PM  = (73*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KCl.MM  = (250*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KBr.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.PM  = (99*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KBr.MM  = (470*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.KI.PP  = (24*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.PM  = (156*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.KI.MM  = (1130*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbF.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.PM  = (40*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbF.MM  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbCl.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.PM  = (134*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbCl.MM  = (260*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbBr.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.PM  = (180*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbBr.MM  = (490*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.RbI.PP  = (82*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.PM  = (280*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.RbI.MM  = (1200*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

D.CsF.PP  = (278*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.PM  = (78*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8
D.CsF.MM  = (23*D_unit)*kj_per_erg*(nm_per_cm^8)*NA; % kJ/mol nm^8

%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
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

%% Generate range (r) in nm
r = Startpoint:Spacing:Endpoint;

%% Split Salt Into Component Metal and Halide
[Metal,Halide] = Separate_Metal_Halide(Salt);

%% Calculate Pauling Coefficients beta: PM = +-   PP = ++     MM = --
beta.PP = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.PM = 1 + q.(Metal)/valence.(Metal) + q.(Halide)/valence.(Halide); % Unitless
beta.MM = 1 + 2*q.(Halide)/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: PM = +-   PP = ++     MM = --
alpha.PP = 1/rho.(Salt); % nm^-1
alpha.PM = 1/rho.(Salt); % nm^-1
alpha.MM = 1/rho.(Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: PM = +-   PP = ++     MM = --
B.PP = beta.PP*b*exp(2*sigma.(Metal)/rho.(Salt));
B.PM = beta.PM*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Salt));
B.MM = beta.MM*b*exp(2*sigma.(Halide)/rho.(Salt));

%% Build PES: Plus-Minus
for i = 1:length(Dispersion_Scale)
    DS = Dispersion_Scale(i);
    if ismembertol(1.0,DS,1e-5)
        Theory = 'TF';
    else
        Theory = strrep(strrep(['TF_D' sprintf('%7.5f',DS)],'.','P'),'-','N');
    end
    
    % Plus - Minus total potential
    U_MetHal.Total = k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r) ...
        + B.PM*exp(-alpha.PM.*r) ...
        - DS.*C.(Salt).PM./(r.^6) ...
        - DS.*D.(Salt).PM./(r.^8);

    % components
    U_MetHal.f = 1./r; % Electrostatics function f(r)
    U_MetHal.g = - DS.*C.(Salt).PM./(r.^6) - DS.*D.(Salt).PM./(r.^8); % Dispersion g(r)
    U_MetHal.h = B.PM*exp(-alpha.PM.*r);% Short range repulsion

    % Plus - Minus total derivative
    U_MetHal.dTotal = -k_0*(e_c^2)*q.(Metal)*q.(Halide)./(r.^2) ...
        - alpha.PM*B.PM*exp(-alpha.PM.*r) ...
        + DS.*6.*C.(Salt).PM./(r.^7) ...
        + DS.*8.*D.(Salt).PM./(r.^9);

    % components
    U_MetHal.df = 1./(r.^2);% Electrostatics function
    U_MetHal.dg = - DS.*6.*C.(Salt).PM./(r.^7) - DS.*8.*D.(Salt).PM./(r.^9) ; % Dispersion
    U_MetHal.dh = alpha.PM*B.PM*exp(-alpha.PM.*r);% Short range repulsion

    % remove infinities
    %U_MetHal = Remove_Infinities(U_MetHal);

    %% Build PES: Plus - Plus

    % Plus - Plus Total Potential
    U_MetMet.Total = k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r) ...
        + B.PP*exp(-alpha.PP.*r) ...
        - DS.*C.(Salt).PP./(r.^6) ...
        - DS.*D.(Salt).PP./(r.^8);

    % components
    U_MetMet.f = 1./r; % Electrostatics function f(r)
    U_MetMet.g = - DS.*C.(Salt).PP./(r.^6) - DS.*D.(Salt).PP./(r.^8); % Dispersion g(r)
    U_MetMet.h = B.PP*exp(-alpha.PP.*r); % Short range repulsion

    % Plus - Plus Total Derivative
    U_MetMet.dTotal = -k_0*(e_c^2)*q.(Metal).*q.(Metal)./(r.^2) ...
        - alpha.PP*B.PP*exp(-alpha.PP.*r)...
        + DS.*6.*C.(Salt).PP./(r.^7) ...
        + DS.*8.*D.(Salt).PP./(r.^9);

    % Components
    U_MetMet.df = 1./(r.^2); % Electrostatics function f(r)
    U_MetMet.dg = - DS.*6.*C.(Salt).PP./(r.^7) - DS.*8.*D.(Salt).PP./(r.^9); % Dispersion g(r)
    U_MetMet.dh = alpha.PP*B.PP*exp(-alpha.PP.*r); % Short range repulsion

    % Remove infinities
    %U_MetMet = Remove_Infinities(U_MetMet);

    %% Build PES: Minus - Minus
    % Minus - Minus Total Potential
    U_HalHal.Total = k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r)...
        + B.MM*exp(-alpha.MM.*r) ...
        - DS.*C.(Salt).MM./(r.^6)...
        - DS.*D.(Salt).MM./(r.^8);

    % components
    U_HalHal.f = 1./r; % Electrostatics function f(r)
    U_HalHal.g = - DS.*C.(Salt).MM./(r.^6) - DS.*D.(Salt).MM./(r.^8); % Dispersion g(r)
    U_HalHal.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion

    % Minus - Minus Total Derivative
    U_HalHal.dTotal = -k_0*(e_c^2)*q.(Halide).*q.(Halide)./(r.^2)...
        - alpha.MM*B.MM*exp(-alpha.MM.*r)...
        + DS.*6.*C.(Salt).MM./(r.^7)...
        + DS.*8.*D.(Salt).MM./(r.^9);

    % components
    U_HalHal.df = 1./(r.^2); % Electrostatics function f(r)
    U_HalHal.dg = - DS.*6.*C.(Salt).MM./(r.^7) - DS.*8.*D.(Salt).MM./(r.^9); % Dispersion g(r)
    U_HalHal.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion

    % Convert r to angstrom
    r_Ang = r./nm_per_Ang;

    %% Place into data structure
    MetMet = [Metal Metal];
    HalHal = [Halide Halide];

    Data.(Salt).Pair.(Theory) = cell(length(r),19);
    Data.(MetMet).Pair.(Theory) = cell(length(r),19);
    Data.(HalHal).Pair.(Theory) = cell(length(r),19);
    for j = 1:length(r)
        % Salt
        Data.(Salt).Pair.(Theory){j,1} = 0; % a
        Data.(Salt).Pair.(Theory){j,2} = 0; % b
        Data.(Salt).Pair.(Theory){j,3} = 0; % c
        Data.(Salt).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(Salt).Pair.(Theory){j,5} = 0; % Volume
        Data.(Salt).Pair.(Theory){j,6} = {Metal 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
        Data.(Salt).Pair.(Theory){j,7} = U_MetHal.Total(j); % Energy
        Data.(Salt).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(Salt).Pair.(Theory){j,9} = 0; % Data type
        Data.(Salt).Pair.(Theory){j,10} = k_0*q.(Metal)*q.(Halide)./(r(j)); % Coulomb real space
        Data.(Salt).Pair.(Theory){j,11} = U_MetHal.g(j) + U_MetHal.h(j); % VDW real space
        Data.(Salt).Pair.(Theory){j,12} = U_MetHal.h(j); % Repulsive Only
        Data.(Salt).Pair.(Theory){j,13} = U_MetHal.g(j); % Dispersion Only
        Data.(Salt).Pair.(Theory){j,14} = k_0*q.(Metal)*q.(Halide)./(r(j)); % Coul-SR: Salt
        Data.(Salt).Pair.(Theory){j,15} = U_MetHal.g(j) + U_MetHal.h(j); % LJ-SR: Salt
        Data.(Salt).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
        Data.(Salt).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
        Data.(Salt).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
        Data.(Salt).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal

        % Metal-Metal
        Data.(MetMet).Pair.(Theory){j,1} = 0; % a
        Data.(MetMet).Pair.(Theory){j,2} = 0; % b
        Data.(MetMet).Pair.(Theory){j,3} = 0; % c
        Data.(MetMet).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(MetMet).Pair.(Theory){j,5} = 0; % Volume
        Data.(MetMet).Pair.(Theory){j,6} = {Metal 0 0 0; Metal 0 0 r_Ang(j)}; % Fractional coords
        Data.(MetMet).Pair.(Theory){j,7} = U_MetMet.Total(j); % Energy
        Data.(MetMet).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(MetMet).Pair.(Theory){j,9} = 0; % Data type
        Data.(MetMet).Pair.(Theory){j,10} = k_0*q.(Metal)*q.(Metal)./(r(j)); % Coulomb real space
        Data.(MetMet).Pair.(Theory){j,11} = U_MetMet.g(j) + U_MetMet.h(j); % VDW real space
        Data.(MetMet).Pair.(Theory){j,12} = U_MetMet.h(j); % Repulsive Only
        Data.(MetMet).Pair.(Theory){j,13} = U_MetMet.g(j); % Dispersion Only
        Data.(MetMet).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
        Data.(MetMet).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
        Data.(MetMet).Pair.(Theory){j,16} = k_0*q.(Metal)*q.(Metal)./(r(j)); % Coul-SR: Met-Met
        Data.(MetMet).Pair.(Theory){j,17} = U_MetMet.g(j) + U_MetMet.h(j); % LJ-SR: Met-Met
        Data.(MetMet).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
        Data.(MetMet).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal

        % Halide-Halide
        Data.(HalHal).Pair.(Theory){j,1} = 0; % a
        Data.(HalHal).Pair.(Theory){j,2} = 0; % b
        Data.(HalHal).Pair.(Theory){j,3} = 0; % c
        Data.(HalHal).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(HalHal).Pair.(Theory){j,5} = 0; % Volume
        Data.(HalHal).Pair.(Theory){j,6} = {Halide 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
        Data.(HalHal).Pair.(Theory){j,7} = U_HalHal.Total(j); % Energy
        Data.(HalHal).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(HalHal).Pair.(Theory){j,9} = 0; % Data type
        Data.(HalHal).Pair.(Theory){j,10} = k_0*q.(Halide)*q.(Halide)./(r(j)); % Coulomb real space
        Data.(HalHal).Pair.(Theory){j,11} = U_HalHal.g(j) + U_HalHal.h(j); % VDW real space
        Data.(HalHal).Pair.(Theory){j,12} = U_HalHal.h(j); % Repulsive Only
        Data.(HalHal).Pair.(Theory){j,13} = U_HalHal.g(j); % Dispersion Only
        Data.(HalHal).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
        Data.(HalHal).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
        Data.(HalHal).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
        Data.(HalHal).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
        Data.(HalHal).Pair.(Theory){j,18} = k_0*q.(Halide)*q.(Halide)./(r(j)); % Coul-SR: Hal-Hal
        Data.(HalHal).Pair.(Theory){j,19} = U_HalHal.g(j) + U_HalHal.h(j); % LJ-SR: Hal-Hal
    end
end


%% JC Potential
%% Conversion factor
kJ_per_kcal = 4.184;

%% JC Ion Parameters in SPC/E water
if strcmp(Watermodel,'SPC/E')
    Param.Li.sigma = (0.791*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.212*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.593*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.737*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (2.021*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.257*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.711*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.751*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.919*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    Param.Li.sigma = (1.025*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.0279896)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.369*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.0874393)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.705*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.1936829)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.813*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.3278219)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.976*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.4065394)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.303*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0033640)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.513*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0355910)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.608*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0586554)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.860*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0536816)*kJ_per_kcal; % kJ/mol

elseif strcmp(Watermodel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    Param.Li.sigma = (0.808*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.1039884)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.226*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.1684375)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.590*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.2794651)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.709*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4331494)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.888*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.3944318)*kJ_per_kcal; % kJ/mol


    Param.F.sigma = (2.538*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0015752)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.760*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0116615)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.768*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0303773)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.952*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0417082)*kJ_per_kcal; % kJ/mol
else
    error(['Unknown Water Model: "' Watermodel ...
        '". Please choose one of "SPC/E", "TIP3P", or "TIP4PEW".'])
end

% Use combining rules
epsilon_MM = Param.(Metal).epsilon; % kJ/mol
epsilon_HH = Param.(Halide).epsilon; % kJ/mol
epsilon_Salt = sqrt(Param.(Metal).epsilon*Param.(Halide).epsilon); % kJ/mol

sigma_MM = Param.(Metal).sigma;
sigma_HH = Param.(Halide).sigma;
sigma_Salt = (1/2)*(Param.(Metal).sigma + Param.(Halide).sigma);

q_M = 1*e_c; % charge of lithium in C
q_H = -1*e_c; % charge of fluoride in C

%% Build JC PES for LiF with dispersion scaling
for i = 1:length(Dispersion_Scale)
    DS = Dispersion_Scale(i);
    if ismembertol(1.0,DS,1e-5)
        Theory = 'JC';
    else
        Theory = strrep(strrep(['JC_D' sprintf('%7.5f',DS)],'.','P'),'-','N');
    end

    U_Salt = k_0.*q_M.*q_H./(r) + 4*epsilon_Salt.*( ( sigma_Salt./r ).^12 - ( DS*sigma_Salt./r ).^6 );

    U_MetMet= k_0.*q_M.*q_M./(r) + 4*epsilon_MM.*( ( sigma_MM./r ).^12 - ( DS*sigma_MM./r ).^6 );

    U_HalHal  = k_0.*q_H.*q_H./(r) + 4*epsilon_HH.*( ( sigma_HH./r ).^12 - ( DS*sigma_HH./r ).^6 );

    %% Place into data structure
    Data.(Salt).Pair.(Theory) = cell(length(r),19);
    Data.(MetMet).Pair.(Theory) = cell(length(r),19);
    Data.(HalHal).Pair.(Theory) = cell(length(r),19);
    for j = 1:length(r)
        % Salt
        Data.(Salt).Pair.(Theory){j,1} = 0; % a
        Data.(Salt).Pair.(Theory){j,2} = 0; % b
        Data.(Salt).Pair.(Theory){j,3} = 0; % c
        Data.(Salt).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(Salt).Pair.(Theory){j,5} = 0; % Volume
        Data.(Salt).Pair.(Theory){j,6} = {Metal 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
        Data.(Salt).Pair.(Theory){j,7} = U_Salt(j); % Total Energy
        Data.(Salt).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(Salt).Pair.(Theory){j,9} = 0; % Data Type
        Data.(Salt).Pair.(Theory){j,10} = k_0.*q_M.*q_H./(r(j)); % Coulomb real space
        Data.(Salt).Pair.(Theory){j,11} = 4*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12 - ( DS*sigma_Salt./r(j) ).^6 ); % VDW real space
        Data.(Salt).Pair.(Theory){j,12} = 4*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12); % Repulsive Only
        Data.(Salt).Pair.(Theory){j,13} = 4*epsilon_Salt.*( - ( DS*sigma_Salt./r(j) ).^6 ); % Dispersion Only
        Data.(Salt).Pair.(Theory){j,14} = k_0.*q_M.*q_H./(r(j)); % Coul-SR: Salt
        Data.(Salt).Pair.(Theory){j,15} = 4*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12 - ( DS*sigma_Salt./r(j) ).^6 ); % LJ-SR: Salt
        Data.(Salt).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
        Data.(Salt).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
        Data.(Salt).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
        Data.(Salt).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal

        % Metal-Metal
        Data.(MetMet).Pair.(Theory){j,1} = 0; % a
        Data.(MetMet).Pair.(Theory){j,2} = 0; % b
        Data.(MetMet).Pair.(Theory){j,3} = 0; % c
        Data.(MetMet).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(MetMet).Pair.(Theory){j,5} = 0; % Volume
        Data.(MetMet).Pair.(Theory){j,6} = {Metal 0 0 0; Metal 0 0 r_Ang(j)}; % Fractional coords
        Data.(MetMet).Pair.(Theory){j,7} = U_MetMet(j); % Energy
        Data.(MetMet).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(MetMet).Pair.(Theory){j,9} = 0; % Data Type
        Data.(MetMet).Pair.(Theory){j,10} = k_0.*q_M.*q_M./(r(j)); % Coulomb real space
        Data.(MetMet).Pair.(Theory){j,11} = 4*epsilon_MM.*( ( sigma_MM./r(j) ).^12 - ( DS*sigma_MM./r(j) ).^6 ); % VDW real space
        Data.(MetMet).Pair.(Theory){j,12} = 4*epsilon_MM.*( ( sigma_MM./r(j) ).^12 ); % Repulsive Only
        Data.(MetMet).Pair.(Theory){j,13} = 4*epsilon_MM.*( - ( DS*sigma_MM./r(j) ).^6 ); % Dispersion Only
        Data.(MetMet).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
        Data.(MetMet).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
        Data.(MetMet).Pair.(Theory){j,16} = k_0.*q_M.*q_M./(r(j)); % Coul-SR: Met-Met
        Data.(MetMet).Pair.(Theory){j,17} = 4*epsilon_MM.*( ( sigma_MM./r(j) ).^12 - ( DS*sigma_MM./r(j) ).^6 ); % LJ-SR: Met-Met
        Data.(MetMet).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
        Data.(MetMet).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal

        % Halide-Halide
        Data.(HalHal).Pair.(Theory){j,1} = 0; % a
        Data.(HalHal).Pair.(Theory){j,2} = 0; % b
        Data.(HalHal).Pair.(Theory){j,3} = 0; % c
        Data.(HalHal).Pair.(Theory){j,4} = r_Ang(j); % BL
        Data.(HalHal).Pair.(Theory){j,5} = 0; % Volume
        Data.(HalHal).Pair.(Theory){j,6} = {Halide 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
        Data.(HalHal).Pair.(Theory){j,7} = U_HalHal(j); % Energy
        Data.(HalHal).Pair.(Theory){j,8} = 0; % Time of calculation
        Data.(HalHal).Pair.(Theory){j,9} = 0; % Data type
        Data.(HalHal).Pair.(Theory){j,10} = k_0.*q_H.*q_H./(r(j)); % Coulomb real space
        Data.(HalHal).Pair.(Theory){j,11} = 4*epsilon_HH.*( ( sigma_HH./r(j) ).^12 - ( DS*sigma_HH./r(j) ).^6 ); % VDW real space
        Data.(HalHal).Pair.(Theory){j,12} = 4*epsilon_HH.*( ( sigma_HH./r(j) ).^12 ); % Repulsive Only
        Data.(HalHal).Pair.(Theory){j,13} = 4*epsilon_HH.*( - ( DS*sigma_HH./r(j) ).^6 ); % Dispersion Only
        Data.(HalHal).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
        Data.(HalHal).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
        Data.(HalHal).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
        Data.(HalHal).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
        Data.(HalHal).Pair.(Theory){j,18} = k_0.*q_H.*q_H./(r(j)); % Coul-SR: Hal-Hal
        Data.(HalHal).Pair.(Theory){j,19} = 4*epsilon_HH.*( ( sigma_HH./r(j) ).^12 - ( DS*sigma_HH./r(j) ).^6 ); % LJ-SR: Hal-Hal
    end
end

% %% Build JC PES for LiF with epsilon scaling
% for i = 1:length(Epsilon_Scale)
%     Gamma = Epsilon_Scale(i); % Epsilon scaling factor
%     if ismembertol(1.0,Gamma,1e-5)
%         Theory = 'JC';
%     else
%         Theory = strrep(strrep(['JC_E' sprintf('%7.5f',Gamma)],'.','P'),'-','N');
%     end
% 
%     U_Salt = k_0.*q_M.*q_H./(r) + 4*Gamma*epsilon_Salt.*( ( sigma_Salt./r ).^12 - ( sigma_Salt./r ).^6 );
% 
%     U_MetMet= k_0.*q_M.*q_M./(r) + 4*Gamma*epsilon_MM.*( ( sigma_MM./r ).^12 - ( sigma_MM./r ).^6 );
% 
%     U_HalHal  = k_0.*q_H.*q_H./(r) + 4*Gamma*epsilon_HH.*( ( sigma_HH./r ).^12 - ( sigma_HH./r ).^6 );
% 
%     %% Place into data structure
%     Data.(Salt).Pair.(Theory) = cell(length(r),19);
%     Data.(MetMet).Pair.(Theory) = cell(length(r),19);
%     Data.(HalHal).Pair.(Theory) = cell(length(r),19);
%     for j = 1:length(r)
%         % Salt
%         Data.(Salt).Pair.(Theory){j,1} = 0; % a
%         Data.(Salt).Pair.(Theory){j,2} = 0; % b
%         Data.(Salt).Pair.(Theory){j,3} = 0; % c
%         Data.(Salt).Pair.(Theory){j,4} = r_Ang(j); % BL
%         Data.(Salt).Pair.(Theory){j,5} = 0; % Volume
%         Data.(Salt).Pair.(Theory){j,6} = {Metal 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
%         Data.(Salt).Pair.(Theory){j,7} = U_Salt(j); % Energy
%         Data.(Salt).Pair.(Theory){j,8} = 0; % Time of calculation
%         Data.(Salt).Pair.(Theory){j,9} = 0; % Data Type
%         Data.(Salt).Pair.(Theory){j,10} = k_0.*q_M.*q_H./(r(j)); % Coulomb real space
%         Data.(Salt).Pair.(Theory){j,11} = 4*Gamma*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12 - ( sigma_Salt./r(j) ).^6 ); % VDW real space
%         Data.(Salt).Pair.(Theory){j,12} = 4*Gamma*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12 ); % Repulsive Only
%         Data.(Salt).Pair.(Theory){j,13} = 4*Gamma*epsilon_Salt.*( - ( sigma_Salt./r(j) ).^6 ); % Dispersion Only
%         Data.(Salt).Pair.(Theory){j,14} = k_0.*q_M.*q_H./(r(j)); % Coul-SR: Salt
%         Data.(Salt).Pair.(Theory){j,15} = 4*Gamma*epsilon_Salt.*( ( sigma_Salt./r(j) ).^12 - ( sigma_Salt./r(j) ).^6 ); % LJ-SR: Salt
%         Data.(Salt).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
%         Data.(Salt).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
%         Data.(Salt).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
%         Data.(Salt).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal
% 
%         % Metal-Metal
%         Data.(MetMet).Pair.(Theory){j,1} = 0; % a
%         Data.(MetMet).Pair.(Theory){j,2} = 0; % b
%         Data.(MetMet).Pair.(Theory){j,3} = 0; % c
%         Data.(MetMet).Pair.(Theory){j,4} = r_Ang(j); % BL
%         Data.(MetMet).Pair.(Theory){j,5} = 0; % Volume
%         Data.(MetMet).Pair.(Theory){j,6} = {Metal 0 0 0; Metal 0 0 r_Ang(j)}; % Fractional coords
%         Data.(MetMet).Pair.(Theory){j,7} = U_MetMet(j); % Energy
%         Data.(MetMet).Pair.(Theory){j,8} = 0; % Time of calculation
%         Data.(MetMet).Pair.(Theory){j,9} = 0; % Data Type
%         Data.(MetMet).Pair.(Theory){j,10} = k_0.*q_M.*q_M./(r(j)); % Coulomb real space
%         Data.(MetMet).Pair.(Theory){j,11} = 4*Gamma*epsilon_MM.*( ( sigma_MM./r(j) ).^12 - ( sigma_MM./r(j) ).^6 ); % VDW real space
%         Data.(MetMet).Pair.(Theory){j,12} = 4*Gamma*epsilon_MM.*( ( sigma_MM./r(j) ).^12 ); % Repulsion Only
%         Data.(MetMet).Pair.(Theory){j,13} = 4*Gamma*epsilon_MM.*( - ( sigma_MM./r(j) ).^6 ); % Dispersion Only
%         Data.(MetMet).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
%         Data.(MetMet).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
%         Data.(MetMet).Pair.(Theory){j,16} = k_0.*q_M.*q_M./(r(j)); % Coul-SR: Met-Met
%         Data.(MetMet).Pair.(Theory){j,17} = 4*Gamma*epsilon_MM.*( ( sigma_MM./r(j) ).^12 - ( sigma_MM./r(j) ).^6 ); % LJ-SR: Met-Met
%         Data.(MetMet).Pair.(Theory){j,18} = 0; % Coul-SR: Hal-Hal
%         Data.(MetMet).Pair.(Theory){j,19} = 0; % LJ-SR: Hal-Hal
% 
%         % Halide-Halide
%         Data.(HalHal).Pair.(Theory){j,1} = 0; % a
%         Data.(HalHal).Pair.(Theory){j,2} = 0; % b
%         Data.(HalHal).Pair.(Theory){j,3} = 0; % c
%         Data.(HalHal).Pair.(Theory){j,4} = r_Ang(j); % BL
%         Data.(HalHal).Pair.(Theory){j,5} = 0; % Volume
%         Data.(HalHal).Pair.(Theory){j,6} = {Halide 0 0 0; Halide 0 0 r_Ang(j)}; % Fractional coords
%         Data.(HalHal).Pair.(Theory){j,7} = U_HalHal(j); % Energy
%         Data.(HalHal).Pair.(Theory){j,8} = 0; % Time of calculation
%         Data.(HalHal).Pair.(Theory){j,9} = 0; % Data type
%         Data.(HalHal).Pair.(Theory){j,10} = k_0.*q_H.*q_H./(r(j)); % Coulomb real space
%         Data.(HalHal).Pair.(Theory){j,11} = 4*Gamma*epsilon_HH.*( ( sigma_HH./r(j) ).^12 - ( sigma_HH./r(j) ).^6 ); % VDW real space
%         Data.(HalHal).Pair.(Theory){j,12} = 4*Gamma*epsilon_HH.*( ( sigma_HH./r(j) ).^12 ); % Repulsion Only
%         Data.(HalHal).Pair.(Theory){j,13} = 4*Gamma*epsilon_HH.*( - ( sigma_HH./r(j) ).^6 ); % Dispersion Only
%         Data.(HalHal).Pair.(Theory){j,14} = 0; % Coul-SR: Salt
%         Data.(HalHal).Pair.(Theory){j,15} = 0; % LJ-SR: Salt
%         Data.(HalHal).Pair.(Theory){j,16} = 0; % Coul-SR: Met-Met
%         Data.(HalHal).Pair.(Theory){j,17} = 0; % LJ-SR: Met-Met
%         Data.(HalHal).Pair.(Theory){j,18} = k_0.*q_H.*q_H./(r(j)); % Coul-SR: Hal-Hal
%         Data.(HalHal).Pair.(Theory){j,19} = 4*Gamma*epsilon_HH.*( ( sigma_HH./r(j) ).^12 - ( sigma_HH./r(j) ).^6 ); % LJ-SR: Hal-Hal
%     end
% end

end
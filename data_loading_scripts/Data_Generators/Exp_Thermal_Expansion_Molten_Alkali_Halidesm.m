%% Molten salt densities
% Data source:
% (1) Physical Properties Relevant to Energy Storage II - https://nvlpubs.nist.gov/nistpubs/Legacy/NSRDS/nbsnsrds61p2.pdf
% (2) Thermodynamic and Transport Properties for Molten salts: Correlation Equations ...
% Used (1) for: LiF, LiCl, LiBr, LiI, NaF, NaCl, KF, KCl
% Densities using equation: rho = a_rho - b_rho*T

% Parameter a, units g / cm^3
Experiment.LiF.Liquid.a_rho		= 2.3581;
Experiment.LiCl.Liquid.a_rho	= 1.8842;
Experiment.LiBr.Liquid.a_rho	= 3.06546;
Experiment.LiI.Liquid.a_rho		= 3.79045;

Experiment.NaF.Liquid.a_rho		= 2.7550;
Experiment.NaCl.Liquid.a_rho	= 2.1393;
Experiment.NaBr.Liquid.a_rho	= 3.1748;
Experiment.NaI.Liquid.a_rho		= 3.6274;

Experiment.KF.Liquid.a_rho		= 2.6464;
Experiment.KCl.Liquid.a_rho		= 2.1359;
Experiment.KBr.Liquid.a_rho		= 2.9583;
Experiment.KI.Liquid.a_rho		= 3.3594;

Experiment.RbF.Liquid.a_rho		= 3.9983;
Experiment.RbCl.Liquid.a_rho	= 3.1210;
Experiment.RbBr.Liquid.a_rho	= 3.7390;
Experiment.RbI.Liquid.a_rho		= 3.9499;

Experiment.CsF.Liquid.a_rho		= 4.8985;
Experiment.CsCl.Liquid.a_rho	= 3.7692;
Experiment.CsBr.Liquid.a_rho	= 4.2449;
Experiment.CsI.Liquid.a_rho		= 4.2410;

% Parameter b
% Units: g / (cm^3 K)
Experiment.LiF.Liquid.b_rho		= 4.902e-4;
Experiment.LiCl.Liquid.b_rho	= 4.328e-4;
Experiment.LiBr.Liquid.b_rho	= 6.5146e-4;
Experiment.LiI.Liquid.b_rho		= 9.1780e-4;

Experiment.NaF.Liquid.b_rho		= 3.63e-4;
Experiment.NaCl.Liquid.b_rho	= 5.430e-4;
Experiment.NaBr.Liquid.b_rho	= 8.169e-4;
Experiment.NaI.Liquid.b_rho		= 9.491e-4;

Experiment.KF.Liquid.b_rho		= 6.515e-4;
Experiment.KCl.Liquid.b_rho		= 5.831e-4;
Experiment.KBr.Liquid.b_rho		= 8.253e-4;
Experiment.KI.Liquid.b_rho		= 9.557e-4;

Experiment.RbF.Liquid.b_rho		= 1.02e-3;
Experiment.RbCl.Liquid.b_rho	= 8.832e-4;
Experiment.RbBr.Liquid.b_rho	= 1.0718e-3;
Experiment.RbI.Liquid.b_rho		= 1.1435e-3;

Experiment.CsF.Liquid.b_rho		= 1.2806e-3;
Experiment.CsCl.Liquid.b_rho	= 1.0650e-3;
Experiment.CsBr.Liquid.b_rho	= 1.2234e-3;
Experiment.CsI.Liquid.b_rho		= 1.1834e-3;

%% Melting Points of Alkali Halides
% Data Source: CRC Handbook (2014)
% Units: K
Experiment.LiF.mp	= 1121.4;
Experiment.LiCl.mp	= 883;
Experiment.LiBr.mp	= 550;
Experiment.LiI.mp	= 743;

Experiment.NaF.mp	= 1269;
Experiment.NaCl.mp	= 1075.168;
Experiment.NaBr.mp	= 1020;
Experiment.NaI.mp	= 1577;

Experiment.KF.mp	= 1131;
Experiment.KCl.mp	= 1044;
Experiment.KBr.mp	= 1007;
Experiment.KI.mp	= 954;

Experiment.RbF.mp	= 1068;
Experiment.RbCl.mp	= 997;
Experiment.RbBr.mp	= 965;
Experiment.RbI.mp	= 929;

Experiment.CsF.mp	= 976;
Experiment.CsCl.mp	= 919;
Experiment.CsBr.mp	= 909;
Experiment.CsI.mp	= 905;

%% Heats of fusion
% Data sources:
% (1) LiF, LiCl, LiBr, LiI, NaF, NaCl, KF, KCl: Physical Properties Relevant to Energy Storage II - https://nvlpubs.nist.gov/nistpubs/Legacy/NSRDS/nbsnsrds61p2.pdf
% (2) RbF, CsF: https://pubs.rsc.org/en/content/articlelanding/1973/f1/f19736902026
% (3) NaBr, NaI, KBr, KI, RbCl, RbBr, RbI, CsCl, CsBr, CsI https://pubs.acs.org/doi/10.1021/j100831a023
% Units kJ/mol
kJ_per_kcal = 4.184;
Experiment.LiF.dH_fus	= 6.43*kJ_per_kcal;
Experiment.LiCl.dH_fus	= 4.76*kJ_per_kcal;
Experiment.LiBr.dH_fus	= 4.22*kJ_per_kcal;
Experiment.LiI.dH_fus	= 3.50*kJ_per_kcal;

Experiment.NaF.dH_fus	= 8.03*kJ_per_kcal;
Experiment.NaCl.dH_fus	= 6.73*kJ_per_kcal;
Experiment.NaBr.dH_fus	= 6.24*kJ_per_kcal;
Experiment.NaI.dH_fus	= 5.64*kJ_per_kcal;

Experiment.KF.dH_fus	= 7.05*kJ_per_kcal;
Experiment.KCl.dH_fus	= 6.34*kJ_per_kcal;
Experiment.KBr.dH_fus	= 6.10*kJ_per_kcal;
Experiment.KI.dH_fus	= 5.74*kJ_per_kcal;

Experiment.RbF.dH_fus	= 6.18*kJ_per_kcal;
Experiment.RbCl.dH_fus	= 5.67*kJ_per_kcal;
Experiment.RbBr.dH_fus	= 5.57*kJ_per_kcal;
Experiment.RbI.dH_fus	= 5.27*kJ_per_kcal;

Experiment.CsF.dH_fus	= 3.35*kJ_per_kcal;
Experiment.CsCl.dH_fus	= 4.84*kJ_per_kcal;
Experiment.CsBr.dH_fus	= 5.64*kJ_per_kcal;
Experiment.CsI.dH_fus	= 5.64*kJ_per_kcal;

%% Element Molar masses
% Unit: g / mol
MM.Li	= 6.9410;
MM.Na	= 22.990;
MM.K	= 39.098;
MM.Rb	= 85.468;
MM.Cs   = 132.91;
MM.F    = 18.998;
MM.Cl   = 35.453;
MM.Br   = 79.904;
MM.I    = 126.90;

%% Molar masses, Molar volume models for MD, and Entropies of fusion
NA =  6.02214076e23; % molecules / mol
cm3_to_Ang3 = 1e-24; % cm^3 / Angstroms^3
DataDir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\ThermExp';

% Unit: g / mol
Salts = fieldnames(Experiment);
for idx = 1:length(Salts)
    Salt = Salts{idx};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    Experiment.(Salt).MM = MM.(Metal) + MM.(Halide);
    
    % Create a model for linear thermal expansion
    b = ( Experiment.(Salt).Liquid.a_rho/Experiment.(Salt).MM )*NA*cm3_to_Ang3; % Y-intercept, in units of [formula unit / Angstrom^3].
    m = (-Experiment.(Salt).Liquid.b_rho/Experiment.(Salt).MM )*NA*cm3_to_Ang3; % Slope, in units of [formula unit / (Angstrom^3 . K)]
    mdl_vol = @(X) (1./(m.*X + b)); % Function that outputs molar volume in units of [Angstrom^3]/[formula unit] for a given temperature
    mdl.predict = mdl_vol;

    % Save the model to file.
    filename = fullfile(DataDir,[Salt '_L_Exp.mat']);
    save(filename,'mdl')
    
    % Calculate entropies of fusion
    % Units: J / mol K
    Experiment.(Salt).dS_fus = Experiment.(Salt).dH_fus*1000/Experiment.(Salt).mp;
end

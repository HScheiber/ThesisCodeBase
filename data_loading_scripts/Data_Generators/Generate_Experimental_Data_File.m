clear;
%% Experimental lattice constants
% Data taken from Alkali Halides handbook
Experiment.LiF.Rocksalt.a = 4.02620;  % Angstrom
Experiment.LiCl.Rocksalt.a = 5.13988; % Angstrom
Experiment.LiBr.Rocksalt.a = 5.501;   % Angstrom
Experiment.LiI.Rocksalt.a = 6.012;    % Angstrom

Experiment.NaF.Rocksalt.a = 4.6329;   % Angstrom
Experiment.NaCl.Rocksalt.a = 5.64009; % Angstrom
Experiment.NaBr.Rocksalt.a = 5.97299; % Angstrom
Experiment.NaI.Rocksalt.a = 6.4728;   % Angstrom

Experiment.KF.Rocksalt.a = 5.344; % Angstrom
Experiment.KCl.Rocksalt.a = 6.29294;  % Angstrom
Experiment.KBr.Rocksalt.a = 6.5982;   % Angstrom
Experiment.KI.Rocksalt.a = 7.06555;   % Angstrom

Experiment.RbF.Rocksalt.a = 5.6516;   % Angstrom
Experiment.RbCl.Rocksalt.a = 6.5898;  % Angstrom
Experiment.RbBr.Rocksalt.a = 6.8908;  % Angstrom
Experiment.RbI.Rocksalt.a = 7.3466;   % Angstrom

Experiment.CsF.Rocksalt.a = 6.020;    % Angstrom
Experiment.CsCl.Rocksalt.a = 7.079;   % Angstrom

Experiment.CsCl.CsCl.a = 4.1200;  % Angstrom
Experiment.CsBr.CsCl.a = 4.2953;  % Angstrom
Experiment.CsI.CsCl.a = 4.5667;   % Angstrom

% Reference temperatures of experimental lattice constant measurements
Experiment.LiF.Rocksalt.Ta = 298.15;  % K
Experiment.LiCl.Rocksalt.Ta = 298.15; % K
Experiment.LiBr.Rocksalt.Ta = 298.15; % K
Experiment.LiI.Rocksalt.Ta = 298.15;  % K

Experiment.NaF.Rocksalt.Ta = 293.15;  % K
Experiment.NaCl.Rocksalt.Ta = 298.15; % K
Experiment.NaBr.Rocksalt.Ta = 298.15; % K
Experiment.NaI.Rocksalt.Ta = 298.15;  % K

Experiment.KF.Rocksalt.Ta = 298.15;   % K
Experiment.KCl.Rocksalt.Ta = 298.15;  % K
Experiment.KBr.Rocksalt.Ta = 298.15;  % K
Experiment.KI.Rocksalt.Ta = 298.15;   % K

Experiment.RbF.Rocksalt.Ta = 293.15;  % K
Experiment.RbCl.Rocksalt.Ta = 293.15; % K
Experiment.RbBr.Rocksalt.Ta = 293.15; % K
Experiment.RbI.Rocksalt.Ta = 300.15;  % K

Experiment.CsF.Rocksalt.Ta = 293.15;  % K
Experiment.CsCl.Rocksalt.Ta = 758.15; % K

Experiment.CsCl.CsCl.Ta = 293.15; % K
Experiment.CsBr.CsCl.Ta = 298.15; % K
Experiment.CsI.CsCl.Ta = 293.15;  % K

%% Experimental thermal expansion data from: Thermal Expansion, nonmetallic solids
% https://apps.dtic.mil/dtic/tr/fulltext/u2/a129116.pdf

% Experimental thermal expansion coefficients [K^{-1}]
Experiment.LiF.Rocksalt.alpha	= [0 0 0.2 0.5 0.7 2.2 5.5 10.0 19.6 26.1 31.0 33.9 34.4 38.6 41.7 44.0 46.5 49.6 53.0 58.2 66.7].*1e-6;
Experiment.LiCl.Rocksalt.alpha	= [0 0 18.6 25.2 33.0 37.8 41.2 43.7 44.0 46.2 47.9 50.1 52.4 54.3 57.4].*1e-6;
Experiment.LiBr.Rocksalt.alpha	= [0 0 24.0 28.6 37.7 42.7 46.4 48.6 48.9 50.8 52.4 54.8 58.5 63.1].*1e-6;
Experiment.LiI.Rocksalt.alpha	= [0 0 0.2 0.5 0.7 2.2 5.5 10.0 19.6 26.1 31.0 33.9 34.4 38.6 41.7 44.0 46.5 49.6 53.0 58.2 66.7].*1e-6.*(60.0/34.0404); % Using scaled version of LiF

Experiment.NaF.Rocksalt.alpha	= [0 0 0.4 0.4 0.5 0.6 3 9.2 14.4 22.5 27.4 30.9 33 33.5 36 37 38 38.5 41.4 44.3 46.4 49.5 53.4 57.5 62.3].*1e-6;
Experiment.NaCl.Rocksalt.alpha	= [0 0 0.007 0.06 0.22 0.56 1.20 9.6 19.0 25.3 31.7 35.4 38.2 40.0 43.5 45.9 49.0 52.7 57.7 63.5 69.5].*1e-6;
Experiment.NaBr.Rocksalt.alpha	= [0 0 0.02 0.2 4 14.2 24.7 29.7 34.5 37.6 40 41.6 41.8 43.4 44.8 46.2 48 50 52.2 54.5 56.7 59.1 61.5 63.8].*1e-6;
Experiment.NaI.Rocksalt.alpha	= [0 0 0.005 0.06 0.7 2.5 4.9 7.3 21 30 34.5 39.1 42 43.7 44.7 45.3 48.5 49 52.2 56.7].*1e-6;

Experiment.KF.Rocksalt.alpha	= [0 0 16.7 19.2 24.3 26.3 27.8 29.3 30.3 31.1 31.2 31.4 33.1 35 37 39.1 41.4 43.7].*1e-6;
Experiment.KCl.Rocksalt.alpha	= [0 0 0.1 0.7 2.5 6.4 10.8 20.7 25.8 30.5 33.3 35.4 36.5 36.6 39.6 42.6 45.7 49.5 53.5 57.6 62.1].*1e-6;
Experiment.KBr.Rocksalt.alpha	= [0 0 0.01 4.3 16.7 25 28.6 33.3 35.7 37.4 38.3 38.5 40.6 43.2 46.7 50 54.3 61.7 72.4].*1e-6;
Experiment.KI.Rocksalt.alpha	= [0 0 0.02 0.04 0.54 2.1 4.7 7.6 21.5 28.5 31.8 35.5 37.5 39 40.2 42 43 48 53.2 59].*1e-6;

Experiment.RbF.Rocksalt.alpha	= [0 0 17.9 23.5 26.8 27.2 27.3744 28.8564 30.5128 32.2564 34.0872 36.0923 38.0974].*1e-6; % Extended beyond 298 K using scaled KF data
Experiment.RbCl.Rocksalt.alpha	= [0 0 -0.0002 -0.0007 0.0001 0.2 1.2 2.2 3.2 10.7 15 21.7 25.8 30.3 33 34.8 36.2 36.6 38.8 42 45 48.3 51.9 55.3].*1e-6;
Experiment.RbBr.Rocksalt.alpha	= [0 0 -0.001 -0.005 0.1 6.2 20.5 26.5 29.1 32.2 34.4 36.2 37.5 37.7 40.1 42.3 44.1 46.4 49.9 55.8].*1e-6;
Experiment.RbI.Rocksalt.alpha	= [0 0 -0.003 -0.006 -0.015 -0.032 -0.042 -0.035 0.02 0.3 2.4 6 10.1 24 30 32.3 34.8 36.5 38 39.1 40.5 41.8 44.5 47.3 50.2 53.2 56.2].*1e-6;

Experiment.CsF.Rocksalt.alpha	= [0 0 0.001 0.004 0.009 0.02 0.042 0.083 0.148 0.25 0.52 0.92 1.4 1.91 2.45 3.03 3.63 4.26 4.9 6.55 14.05 17.35 19.24 21.35 33.85 33.9038 34.1211 35.9684 38.0331 40.2064 42.4884 44.9877 47.4870].*1e-6;
Experiment.CsCl.Rocksalt.alpha	= [0 0 10.5 24.7 31.5 35.9 40.5 42.4 44.4 46.2 46.6 51.1 55.8 61.1 66.3].*1e-6;

Experiment.CsCl.CsCl.alpha		= [0 0 10.5 24.7 31.5 35.9 40.5 42.4 44.4 46.2 46.6 51.1 55.8 61.1 66.3].*1e-6;
Experiment.CsBr.CsCl.alpha		= [0 0 0.01 0.12 1.5 5.52 10 15 31.6 35.8 38.2 41 43.2 45.2 46.8 49.4 51.5 55.9 60.5 65.2 70 73.8].*1e-6;
Experiment.CsI.CsCl.alpha		= [0 0 0.013 0.22 2.5 7.8 14 20 36 39.1 41.5 43.5 44.8 46.5 48.3 50.5 52.7 57.4 62 67 72 74.7].*1e-6;

% Experimental thermal expansion coefficient temperatures [K]
Experiment.LiF.Rocksalt.Talpha	= [-1 0 10  20  25  50  75  100  150  200  250  293  300  400  500  600  700  800  900  1000 1100];
Experiment.LiCl.Rocksalt.Talpha	= [-1 0 75 100 150 200 250 293 300 350 400 450 500 550 600];
Experiment.LiBr.Rocksalt.Talpha	= [-1 0 75 100 150 200 250 293 300 350 400 450 500 550];
Experiment.LiI.Rocksalt.Talpha	= [-1 0 10  20  25  50  75  100  150  200  250  293  300  400  500  600  700  800  900  1000 1100];

Experiment.NaF.Rocksalt.Talpha	= [-1 0 10 15 20 25 50 75 100 150 200 250 293 300 350 400 450 500 600 700 800 900 1000 1100 1200];
Experiment.NaCl.Rocksalt.Talpha	= [-1 0 5 10 15 20 25 50 75 100 150 200 250 293 400 500 600 700 800 900 1000];
Experiment.NaBr.Rocksalt.Talpha	= [-1 0 5 10 25 50 75 100 150 200 250 293 300 350 400 450 500 550 600 650 700 750 800 850];
Experiment.NaI.Rocksalt.Talpha	= [-1 0 2 5 10 15 20 25 50 75 100 150 200 250 293 350 400 450 500 550];

Experiment.KF.Rocksalt.Talpha	= [-1 0 80 100 150 175 200 225 250 275 293 300 350 400 450 500 550 600];
Experiment.KCl.Rocksalt.Talpha	= [-1 0 10 20 30 40 50 75 100 150 200 250 293 300 400 500 600 700 800 900 1000];
Experiment.KBr.Rocksalt.Talpha	= [-1 0 5 25 50 75 100 150 200 250 293 300 400 500 600 700 800 900 1000];
Experiment.KI.Rocksalt.Talpha	= [-1 0 4 5 10 15 20 25 50 75 100 150 200 250 293 350 400 500 600 700];

Experiment.RbF.Rocksalt.Talpha	= [-1 0 131 179 263 293 300 350 400 450 500 550 600]; % Extended beyond 298 K using scaled KF data
Experiment.RbCl.Rocksalt.Talpha	= [-1 0 2 3 5 10 20 25 30 40 50 75 100 150 200 250 293 300 400 500 600 700 800 900];
Experiment.RbBr.Rocksalt.Talpha	= [-1 0 3 5 10 25 50 75 100 150 200 250 293 300 400 500 600 700 800 900];
Experiment.RbI.Rocksalt.Talpha	= [-1 0 2 3 4 5 6 7 8 10 15 20 25 50 75 100 150 200 250 293 350 400 500 600 700 800 900];

Experiment.CsF.Rocksalt.Talpha	= [-1 0 3 4 5 6 7 8 9 10 12 14 16 18 20 22 24 26 28 30 57 65 75 85 283 293 300 350 400 450 500 550 600];
Experiment.CsCl.Rocksalt.Talpha	= [-1 0 25 50 75 100 150 200 250 293 300 400 500 600 700];

Experiment.CsCl.CsCl.Talpha		= [-1 0 25 50 75 100 150 200 250 293 300 400 500 600 700];
Experiment.CsBr.CsCl.Talpha		= [-1 0 2 5 10 15 20 25 50 75 100 150 200 250 293 350 400 500 600 700 800 875];
Experiment.CsI.CsCl.Talpha		= [-1 0 2 5 10 15 20 25 50 75 100 150 200 250 293 350 400 500 600 700 800 850];

%% Wurtzite lattice constants for LiCl, LiBr, and LiI
% From: Bach, A.; Fischer, D.; Jansen, M. Synthesis of a New Modification of 
% Lithium Chloride Confirming Theoretical Predictions. Zeitschrift fur Anorg. und Allg. Chemie 2009, 635, 2406–2409.
Experiment.LiCl.Wurtzite.a = 3.852; % Wurtzite lattice parameter a at T = 223.15 K (Angstrom)
Experiment.LiCl.Wurtzite.da = 0.0001; % Error in Wurtzite lattice parameter a at T = 223.15 K (Angstrom)
Experiment.LiCl.Wurtzite.c = 6.118; % Wurtzite lattice parameter c at T = 223.15 K (Angstrom)
Experiment.LiCl.Wurtzite.dc = 0.0001; % Wurtzite lattice parameter c at T = 223.15 K (Angstrom)
Experiment.LiCl.Wurtzite.c_over_a = 6.118/3.852; % Unitless
Experiment.LiCl.Wurtzite.Ta = 223.15;


% From: Liebold-Ribeiro, Y.; Fischer, D.; Jansen, M. Experimental 
% Substantiation of the “Energy Landscape Concept” for Solids: Synthesis 
% of a New Modification of LiBr. Angew. Chemie - Int. Ed. 2008, 47, 4428–4431.
Experiment.LiBr.Wurtzite.a = 4.1509; % Wurtzite lattice parameter a at T = 223.15 K (Angstrom)
Experiment.LiBr.Wurtzite.da = 0.0005; % Error in Wurtzite lattice parameter a at T = 223.15 K (Angstrom)
Experiment.LiBr.Wurtzite.c = 6.6502; % Wurtzite lattice parameter c at T = 223.15 K (Angstrom)
Experiment.LiBr.Wurtzite.dc = 0.0002; % Error in Wurtzite lattice parameter c at T = 223.15 K (Angstrom)
Experiment.LiBr.Wurtzite.c_over_a = 6.6502/4.1509;
Experiment.LiBr.Wurtzite.Ta = 223.15;

% From: ?an?arevi?, Ž. P.; Schön, J. C.; Fischer, D.; Jansen, M. 
% Theoretical and Experimental Exploration of the Energy Landscape of Lil. 
% In Materials Science Forum; 2005; Vol. 494, pp 61–66.
Experiment.LiI.Wurtzite.a = 4.514; % Wurtzite lattice parameter a at T = 273.15 K (Angstrom)
Experiment.LiI.Wurtzite.da = 0.001; % Assumed Error in Wurtzite lattice parameter a at T = 273.15 K (Angstrom)
Experiment.LiI.Wurtzite.c = 7.311; % Wurtzite lattice parameter c at T = 273.15 K (Angstrom)
Experiment.LiI.Wurtzite.dc = 0.001; % Assumed Error in Wurtzite lattice parameter c at T = 273.15 K (Angstrom)
Experiment.LiI.Wurtzite.c_over_a = 7.311/4.514;
Experiment.LiI.Wurtzite.Ta = 273.15;

%% Errors in lattice constants a
% Data taken from Alkali Halides handbook
% Units Angstrom
Experiment.LiF.Rocksalt.da = 5e-05;
Experiment.LiCl.Rocksalt.da = 4e-05;
Experiment.LiBr.Rocksalt.da = 6e-3;
Experiment.LiI.Rocksalt.da = 7e-3;

Experiment.NaF.Rocksalt.da = 5e-4;
Experiment.NaCl.Rocksalt.da = 3e-5;
Experiment.NaBr.Rocksalt.da = 5e-5;
Experiment.NaI.Rocksalt.da = 5e-4;

Experiment.KF.Rocksalt.da = 3e-3;
Experiment.KCl.Rocksalt.da = 8e-5;
Experiment.KBr.Rocksalt.da = 2e-4;
Experiment.KI.Rocksalt.da = 1.5e-4;

Experiment.RbF.Rocksalt.da = 1e-4;
Experiment.RbCl.Rocksalt.da = 2e-4;
Experiment.RbBr.Rocksalt.da = 2e-4;
Experiment.RbI.Rocksalt.da = 2e-4;

Experiment.CsF.Rocksalt.da = 6e-3;
Experiment.CsCl.Rocksalt.da = 4e-3;

Experiment.CsCl.CsCl.da = 5e-4;
Experiment.CsBr.CsCl.da = 5e-4;
Experiment.CsI.CsCl.da = 5e-4;

%% 0 K lattice constants by Integration of experimental linear thermal expansion coefficients
Salts = fieldnames(Experiment);
for idx = 1:length(Salts)
    Salt = Salts{idx};
    Dat = Experiment.(Salt);
    Structures = fieldnames(Dat);
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        if isfield(Experiment.(Salt).(Structure),'Talpha')
            T = Experiment.(Salt).(Structure).Talpha;
            Alpha = Experiment.(Salt).(Structure).alpha;
            T_ref = Experiment.(Salt).(Structure).Ta;

            % Interpolate
            T_inp = 0:0.001:T_ref;
            Alpha_inp = interp1(T,Alpha,T_inp,'makima');

            % Integrate from 0 to the measured temperature
            Correction = - trapz(T_inp,Alpha_inp);
            Experiment.(Salt).(Structure).a_zero = Experiment.(Salt).(Structure).a + Correction*Experiment.(Salt).(Structure).a;
            if isfield(Experiment.(Salt).(Structure),'c')
                Experiment.(Salt).(Structure).c_zero = Experiment.(Salt).(Structure).c + Correction*Experiment.(Salt).(Structure).c;
                Experiment.(Salt).(Structure).c_over_a_zero = Experiment.(Salt).(Structure).c_zero/Experiment.(Salt).(Structure).a_zero; % Unitless
            end
        elseif isfield(Experiment.(Salt).Rocksalt,'Talpha')
            T = Experiment.(Salt).Rocksalt.Talpha;
            Alpha = Experiment.(Salt).Rocksalt.alpha;
            T_ref = Experiment.(Salt).(Structure).Ta;

            % Interpolate
            T_inp = 0:0.001:T_ref;
            Alpha_inp = interp1(T,Alpha,T_inp,'makima');

            % Integrate from 0 to the measured temperature
            Correction = - trapz(T_inp,Alpha_inp);
            Experiment.(Salt).(Structure).a_zero = Experiment.(Salt).(Structure).a + Correction*Experiment.(Salt).(Structure).a;
            if isfield(Experiment.(Salt).(Structure),'c')
                Experiment.(Salt).(Structure).c_zero = Experiment.(Salt).(Structure).c + Correction*Experiment.(Salt).(Structure).c;
                Experiment.(Salt).(Structure).c_over_a_zero = Experiment.(Salt).(Structure).c_zero/Experiment.(Salt).(Structure).a_zero; % Unitless
            end
        else
            Experiment.(Salt).(Structure).a_zero = nan;
        end
    end
end

%% Melting Points of Alkali Halides
% Data Source: CRC Handbook (2017-2018)
% Units: K
Experiment.LiF.mp	= 1121.4;
Experiment.LiCl.mp	= 883;
Experiment.LiBr.mp	= 823;
Experiment.LiI.mp	= 742;

Experiment.NaF.mp	= 1269;
Experiment.NaCl.mp	= 1075.168;
Experiment.NaBr.mp	= 1020;
Experiment.NaI.mp	= 934;

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

% Add in uncertainties
Salts = fieldnames(Experiment);
for idx = 1:length(Salts)
    Experiment.(Salts{idx}).dmp = 1;
end
Experiment.LiF.dmp = 0.1;
Experiment.NaCl.dmp = 0.001;

%% Molten salt densities function parameters
% Data source:
% Thermodynamic and Transport Properties for Molten salts: Correlation Equations ...
% Densities using equation: rho = a_rho - b_rho*T

% Parameter a, units g / cm^3
Experiment.LiF.Liquid.a_rho		= 2.3581;
Experiment.LiCl.Liquid.a_rho	= 1.8842;
Experiment.LiBr.Liquid.a_rho	= 3.0658;
Experiment.LiI.Liquid.a_rho		= 3.7902;

Experiment.NaF.Liquid.a_rho		= 2.7550;
Experiment.NaCl.Liquid.a_rho	= 2.1389;
Experiment.NaBr.Liquid.a_rho	= 3.1748;
Experiment.NaI.Liquid.a_rho		= 3.6274;

Experiment.KF.Liquid.a_rho		= 2.6464;
Experiment.KCl.Liquid.a_rho		= 2.1359;
Experiment.KBr.Liquid.a_rho		= 2.9583;
Experiment.KI.Liquid.a_rho		= 3.3594;

Experiment.RbF.Liquid.a_rho		= 3.9953;
Experiment.RbCl.Liquid.a_rho	= 3.1210;
Experiment.RbBr.Liquid.a_rho	= 3.7390;
Experiment.RbI.Liquid.a_rho		= 3.9499;

Experiment.CsF.Liquid.a_rho		= 4.8985;
Experiment.CsCl.Liquid.a_rho	= 3.7692;
Experiment.CsBr.Liquid.a_rho	= 4.2449;
Experiment.CsI.Liquid.a_rho		= 4.2550;

% Parameter b
% Units: g / (cm^3 K)
Experiment.LiF.Liquid.b_rho		= 4.902e-4;
Experiment.LiCl.Liquid.b_rho	= 4.328e-4;
Experiment.LiBr.Liquid.b_rho	= 6.52e-4;
Experiment.LiI.Liquid.b_rho		= 9.176e-4;

Experiment.NaF.Liquid.b_rho		= 6.36e-4;
Experiment.NaCl.Liquid.b_rho	= 5.426e-4;
Experiment.NaBr.Liquid.b_rho	= 8.169e-4;
Experiment.NaI.Liquid.b_rho		= 9.491e-4;

Experiment.KF.Liquid.b_rho		= 6.515e-4;
Experiment.KCl.Liquid.b_rho		= 5.831e-4;
Experiment.KBr.Liquid.b_rho		= 8.253e-4;
Experiment.KI.Liquid.b_rho		= 9.557e-4;

Experiment.RbF.Liquid.b_rho		= 1.0211e-3;
Experiment.RbCl.Liquid.b_rho	= 8.832e-4;
Experiment.RbBr.Liquid.b_rho	= 1.0718e-3;
Experiment.RbI.Liquid.b_rho		= 1.1435e-3;

Experiment.CsF.Liquid.b_rho		= 1.2806e-3;
Experiment.CsCl.Liquid.b_rho	= 1.0650e-3;
Experiment.CsBr.Liquid.b_rho	= 1.2234e-3;
Experiment.CsI.Liquid.b_rho		= 1.1833e-3;

%% Molten salt densities at the melting point
% Data source:
% CRC Handbook pp 4-124

% Densities units g / cm^3
Experiment.LiF.Liquid.rho_mp	= 1.81;
Experiment.LiCl.Liquid.rho_mp   = 1.502;
Experiment.LiBr.Liquid.rho_mp   = 2.528;
Experiment.LiI.Liquid.rho_mp    = 3.109;

Experiment.NaF.Liquid.rho_mp    = 1.948;
Experiment.NaCl.Liquid.rho_mp	= 1.556;
Experiment.NaBr.Liquid.rho_mp	= 2.342;
Experiment.NaI.Liquid.rho_mp	= 2.742;

Experiment.KF.Liquid.rho_mp		= 1.910;
Experiment.KCl.Liquid.rho_mp	= 1.527;
Experiment.KBr.Liquid.rho_mp	= 2.127;
Experiment.KI.Liquid.rho_mp		= 2.448;

Experiment.RbF.Liquid.rho_mp	= 2.87;
Experiment.RbCl.Liquid.rho_mp	= 2.248;
Experiment.RbBr.Liquid.rho_mp	= 2.715;
Experiment.RbI.Liquid.rho_mp	= 2.904;

Experiment.CsF.Liquid.rho_mp	= 3.649;
Experiment.CsCl.Liquid.rho_mp	= 2.79;
Experiment.CsBr.Liquid.rho_mp	= 3.133;
Experiment.CsI.Liquid.rho_mp	= 3.197;

%% Fractional Change in Volume during melting (dV / Vs)
% Data sources:
% (1) LiF LiCl LiBr NaF NaCl NaBr KI RbCl RbBr CsCl
% Source: Alkali Halides: A Handbook of Physical Properties

% (2) LiI NaI KF KCl KBr RbI
% Source: https://link.springer.com/article/10.1007/BF00503248

% (3) CsBr CsI
% Source: https://pubs.acs.org/doi/pdf/10.1021/ja01615a015

% unitless
Experiment.LiF.Liquid.dVVs	  = 0.294;
Experiment.LiCl.Liquid.dVVs   = 0.262;
Experiment.LiBr.Liquid.dVVs   = 0.243;
Experiment.LiI.Liquid.dVVs    = 0.202;

Experiment.NaF.Liquid.dVVs    = 0.274;
Experiment.NaCl.Liquid.dVVs   = 0.250;
Experiment.NaBr.Liquid.dVVs	  = 0.224;
Experiment.NaI.Liquid.dVVs	  = 0.186;

Experiment.KF.Liquid.dVVs     = 0.172;
Experiment.KCl.Liquid.dVVs	  = 0.173;
Experiment.KBr.Liquid.dVVs	  = 0.166;
Experiment.KI.Liquid.dVVs	  = 0.159;

Experiment.RbF.Liquid.dVVs	  = 0.156;
Experiment.RbCl.Liquid.dVVs	  = 0.143;
Experiment.RbBr.Liquid.dVVs	  = 0.135;
Experiment.RbI.Liquid.dVVs	  = 0.123;

Experiment.CsF.Liquid.dVVs	  = nan;
Experiment.CsCl.Liquid.dVVs	  = 0.105;
Experiment.CsBr.Liquid.dVVs	  = 0.268;
Experiment.CsI.Liquid.dVVs	  = 0.285;

%% Heats of fusion
% Data sources:
% (1) CRC Handbook for most
% (2) For LiBr: https://pubs.acs.org/doi/10.1021/j100831a023
% Other sources: Physical Properties Relevant to Energy Storage II - https://nvlpubs.nist.gov/nistpubs/Legacy/NSRDS/nbsnsrds61p2.pdf
%                https://pubs.rsc.org/en/content/articlelanding/1973/f1/f19736902026
%                https://pubs.acs.org/doi/10.1021/j100831a023
% Units kJ/mol
Experiment.LiF.dH_fus	= 27.09;
Experiment.LiCl.dH_fus	= 19.8;
Experiment.LiBr.dH_fus	= 17.66;
Experiment.LiI.dH_fus	= 14.6;

Experiment.NaF.dH_fus	= 33.35;
Experiment.NaCl.dH_fus	= 28.16;
Experiment.NaBr.dH_fus	= 26.23;
Experiment.NaI.dH_fus	= 23.7;

Experiment.KF.dH_fus	= 27.2;
Experiment.KCl.dH_fus	= 26.28;
Experiment.KBr.dH_fus	= 25.52;
Experiment.KI.dH_fus	= 24.0;

Experiment.RbF.dH_fus	= 25.8;
Experiment.RbCl.dH_fus	= 24.4;
Experiment.RbBr.dH_fus	= 23.3;
Experiment.RbI.dH_fus	= 22.1;

Experiment.CsF.dH_fus	= 21.7;
Experiment.CsCl.dH_fus	= 20.4;
Experiment.CsBr.dH_fus	= 23.6;
Experiment.CsI.dH_fus	= 25.7;

%% Additional data: Molar masses, Entropies of fusion, Redundant Geometric Data
% Elemental Molar masses
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

N_A =  6.02214076e23; % molecules / mol
cm3_per_Ang3 = 1e-24; % cm^3 / Angstroms^3
DataDir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\ThermExp';

% Unit: g / mol
Salts = fieldnames(Experiment);
for idx = 1:length(Salts)
    Salt = Salts{idx};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    % Molar masses (g/mol)
    Experiment.(Salt).MM = MM.(Metal) + MM.(Halide);
    
    % Entropies of fusion
    % Units: J / mol K
    Experiment.(Salt).dS_fus = Experiment.(Salt).dH_fus*1000/Experiment.(Salt).mp;
    
    % Liquid molar volumes at MP in Ang^3 / Forumla unit
    Experiment.(Salt).Liquid.V_mp = (Experiment.(Salt).MM/Experiment.(Salt).Liquid.rho_mp)*(1/(cm3_per_Ang3*N_A));
    
    % Solid molar volumes at MP in Ang^3 / Forumla unit
    % (Delta V / Vs)
    switch Salt
        case {'CsBr' 'CsI'}
            Experiment.(Salt).CsCl.V_mp = Experiment.(Salt).Liquid.V_mp/(1 + Experiment.(Salt).Liquid.dVVs);
        otherwise
            Experiment.(Salt).Rocksalt.V_mp = Experiment.(Salt).Liquid.V_mp/(1 + Experiment.(Salt).Liquid.dVVs);
    end
    
    % Redundant data
    if isfield(Experiment.(Salt),'Rocksalt')
        Experiment.(Salt).Rocksalt.b = Experiment.(Salt).Rocksalt.a;
        Experiment.(Salt).Rocksalt.b_zero = Experiment.(Salt).Rocksalt.a_zero;
        Experiment.(Salt).Rocksalt.c = Experiment.(Salt).Rocksalt.a;
        Experiment.(Salt).Rocksalt.c_zero = Experiment.(Salt).Rocksalt.a_zero;
        Experiment.(Salt).Rocksalt.c_over_a = 1;
        
        % Volumes and densities
        Experiment.(Salt).Rocksalt.V = (Experiment.(Salt).Rocksalt.a^3)/4; % Ang^3 / Forumla unit
        Experiment.(Salt).Rocksalt.V_zero = (Experiment.(Salt).Rocksalt.a_zero^3)/4; % Ang^3 / Forumla unit
        Experiment.(Salt).Rocksalt.density = Experiment.(Salt).MM/(Experiment.(Salt).Rocksalt.V*cm3_per_Ang3*N_A); % g/cm^3
        Experiment.(Salt).Rocksalt.density_zero = Experiment.(Salt).MM/(Experiment.(Salt).Rocksalt.V_zero*cm3_per_Ang3*N_A); % g/cm^3
    end
    if isfield(Experiment.(Salt),'CsCl')
        Experiment.(Salt).CsCl.b = Experiment.(Salt).CsCl.a;
        Experiment.(Salt).CsCl.b_zero = Experiment.(Salt).CsCl.a_zero;
        Experiment.(Salt).CsCl.c = Experiment.(Salt).CsCl.a;
        Experiment.(Salt).CsCl.c_zero = Experiment.(Salt).CsCl.a_zero;
        Experiment.(Salt).CsCl.c_over_a = 1;
        
        % Volumes and densities
        Experiment.(Salt).CsCl.V = (Experiment.(Salt).CsCl.a^3); % Ang^3 / Forumla unit
        Experiment.(Salt).CsCl.V_zero = (Experiment.(Salt).CsCl.a_zero^3); % Ang^3 / Forumla unit
        
        Experiment.(Salt).CsCl.density = Experiment.(Salt).MM/(Experiment.(Salt).CsCl.V*cm3_per_Ang3*N_A); % g/cm^3
        Experiment.(Salt).CsCl.density_zero = Experiment.(Salt).MM/(Experiment.(Salt).CsCl.V_zero*cm3_per_Ang3*N_A); % g/cm^3
    end
    if isfield(Experiment.(Salt),'Wurtzite')
        Experiment.(Salt).Wurtzite.E = nan;
        Experiment.(Salt).Wurtzite.b = Experiment.(Salt).Wurtzite.a;
        Experiment.(Salt).Wurtzite.b_zero = Experiment.(Salt).Wurtzite.a_zero;
        
        % Volumes and densities
        Experiment.(Salt).Wurtzite.V = (Experiment.(Salt).Wurtzite.a*Experiment.(Salt).Wurtzite.b*...
            cosd(30)*Experiment.(Salt).Wurtzite.c)/2; % Ang^3 / Forumla unit
        Experiment.(Salt).Wurtzite.V_zero = (Experiment.(Salt).Wurtzite.a_zero*Experiment.(Salt).Wurtzite.b_zero*...
            cosd(30)*Experiment.(Salt).Wurtzite.c_zero)/2; % Ang^3 / Forumla unit
        Experiment.(Salt).Wurtzite.density = Experiment.(Salt).MM/(Experiment.(Salt).Wurtzite.V*cm3_per_Ang3*N_A); % g/cm^3
        Experiment.(Salt).Wurtzite.density_zero = Experiment.(Salt).MM/(Experiment.(Salt).Wurtzite.V_zero*cm3_per_Ang3*N_A); % g/cm^3
    end
end

%% Lattice energies with errors for LiX and NaCl
% Data source: Scheiber H. O., Patey G. N. (2021)
Experiment.LiF.Rocksalt.E = -1054.2; % Lattice energy (kJ/mol)
Experiment.LiF.Rocksalt.dE = 1.3000; % Error in Lattice energy (kJ/mol)
Experiment.LiCl.Rocksalt.E = -865.4; % Lattice energy (kJ/mol)
Experiment.LiCl.Rocksalt.dE = 1.5000; % Error in Lattice energy (kJ/mol)
Experiment.LiBr.Rocksalt.E = -821.3; % Lattice energy (kJ/mol)
Experiment.LiBr.Rocksalt.dE = 1.5; % Error in Lattice energy (kJ/mol)
Experiment.LiI.Rocksalt.E = -764.5; % Lattice energy (kJ/mol)
Experiment.LiI.Rocksalt.dE = 1.3; % Error in Lattice energy (kJ/mol)
Experiment.NaCl.Rocksalt.E = -791.3; % Lattice energy (kJ/mol)
Experiment.NaCl.Rocksalt.dE = 0.8; % Error in Lattice energy (kJ/mol)

%% Save data
home = find_home;
Data_Directory = [home filesep 'data'];
save(fullfile(Data_Directory,'Alkali_Halides_Exp.mat'),'Experiment');

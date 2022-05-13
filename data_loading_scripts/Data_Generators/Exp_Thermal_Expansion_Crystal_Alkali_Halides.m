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

%% Calculate Data
DataDir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\ThermExp';
Salts = fieldnames(Experiment);

for idx = 1:length(Salts)
    Salt = Salts{idx};
    Dat = Experiment.(Salt);
    Structures = fieldnames(Dat);
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        T = Experiment.(Salt).(Structure).Talpha;
        Alpha = Experiment.(Salt).(Structure).alpha;
        
        T_inp = -1:0.1:T(end);
        Alpha_inp = interp1(T,Alpha,T_inp,'makima');
        
        % Generate a model
        modelfun = @(b,x)b(1) + b(2)*x + b(3)*x.^2;
        beta0 = [1e-6 1e-7 -1e-10];
        opts = statset('Display','iter','TolFun',1e-20,'TolX',1e-20,'maxiter',1000);
        model = fitnlm(T(end-5:end),Alpha(end-5:end),modelfun,beta0,'Options',opts);

        Tfit = ((T(end)+0.1):0.1:2500)';
        [afit,afitci] = predict(model,Tfit);
        
        h = figure;
        hold on
        plot(T_inp,Alpha_inp)
        scatter(T,Alpha)
        plot(Tfit,afit)
        plot(Tfit,afitci(:,1),':g')
        plot(Tfit,afitci(:,2),':g')
        title([Salt ' ' Structure])
        drawnow;
        waitfor(h)
        
        % Create one continuous linear thermal expansion including extrapolation
        T_fullmdl = [T_inp Tfit'];
        Alpha_fullmdl = [Alpha_inp afit'];
        
        % Calculate the zero a
        a_ref = Experiment.(Salt).(Structure).a;
        T_ref = Experiment.(Salt).(Structure).Ta;
        a_0 = a_ref*(1 - trapz(T_fullmdl(T_fullmdl<=T_ref),Alpha_fullmdl(T_fullmdl<=T_ref)));
        
        % Create a model for linear thermal expansion
        thermal_LiF = @(X) (a_0 + ...
            a_0*trapz(T_fullmdl(T_fullmdl<=X),...
            Alpha_fullmdl(T_fullmdl<=X)));
        
        mdl.predict = thermal_LiF;
        
        % Save the model to file.
        filename = fullfile(DataDir,[Salt '_' StructureLabel(Structure) '_Exp.mat']);
        save(filename,'mdl')
    end
end
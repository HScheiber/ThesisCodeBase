%% Global Calculation settings
clear;
[home,project,computer,slurm] = find_home;

ResultsFile = 'AH_Structure_Classifier_Results.mat';
JsonResultsFile = 'AHStructureDetector_tests.json';
Salt_Struc_Model_T = {...
                     {'LiF' 'Rocksalt' 'Mie' 1600}; ...
                     {'LiI' 'Wurtzite' 'JC' 800}; ...
                     {'LiI' 'Sphalerite' 'JC' 800}; ...
                     {'LiI' 'BetaBeO' 'JC' 800}; ...
                     {'LiI' 'FiveFive' 'JC' 800}; ...
                     {'CsI' 'NiAs' 'BF' 700}; ...
                     {'CsI' 'AntiNiAs' 'BF' 700}; ...
                     {'CsI' 'CsCl' 'BH' 1200}; ...
                     {'CsI' 'Liquid' 'JC' 900} ...
                     };
NN_Times = {[0  0]; ...
            [2  1]; ...
            [4  1]; ...
            [10 1]; ...
            [10 5]; ...
            [20 5]; ...
            [50 5]};
Qlm_Voronoi = {[false false]; ...
               [true false]; ...
               [false true]};

Results = struct;
for idx = 1:length(Salt_Struc_Model_T)
    Salt = Salt_Struc_Model_T{idx}{1};
    Structure = Salt_Struc_Model_T{idx}{2};
    Theory = Salt_Struc_Model_T{idx}{3};
    T0 = Salt_Struc_Model_T{idx}{4};
    Results.(Structure).Salt = Salt;
    Results.(Structure).Theory = Theory;
    Results.(Structure).T0 = T0;
    
    %% Load Settings object and set some shared calculation settings
    Settings = Initialize_MD_Settings;
    Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
    Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
    Settings.Nodes = 1; % Minimum number of cores to request for calculation.
    Settings.Cores = 8; % Minimum number of cores to request for calculation. Set to -1 for entire node
    Settings.MPI_Ranks = 8;
    Settings.OMP_Threads = 1;
    Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
    Settings.SinglePrecision = false; % choose true for single precision mode, false for double
    Settings.Project_Directory_Name = 'AH_Structure_Analyzer_Tests'; % Name of project directory to contain job within the main project folder
    Settings.MinMDP.nsteps_min = 1000;
    Settings.MinMDP.Parallel_Min = false;
    Settings.Use_Conv_cell = false;

    %% NaCl/WBK polarized vs unpolarized Alexandria model at Tm starting in rocksalt structure: NVT
    Settings.Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH, BF
    Settings.Salt = Salt;
    Settings.Structure = Structure;
    Settings.SigmaEpsilon = true;
    %a_RS = 5.857263144; % Angstoms
    %T0 = 909;
    Settings.RunEnergyAnalysis = {}; %{'Potential' 'Total-Energy' 'Temperature' 'Pressure'};

    Settings.MDP.RVDW_Cutoff = 0.5; % nm
    Settings.MDP.RCoulomb_Cutoff = 0.6; % nm
    Settings.MDP.RList_Cutoff = 0.6; % nm
    Settings.MDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings.MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction
    Settings.Output_Energies = 1000; % Number of steps that else between writing energies to energy file.
    Settings.Output_Coords = 1000; % Number of steps between outputting coordinates
    Settings.MDP.dt = 0.001; % Time step in ps for md type calculations
    Settings.PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
    Settings.CheckAmorphousLiquid = false;
    Settings.Annealing = 'no'; % Options: 'no' 'single' 'periodic'
    Settings.MDP.Trajectory_Time = 0.05; % ns
    Settings.MDP.dt = 0.001; % Time step in ps for md type calculations

    Settings.Model = 'Alexandria'; % Name of the current model. Leave blank for the default JC/TF/BH model
    Settings.JobID = 'Test'; % An ID that is tacked onto the folder name of all current jobs
    Settings.N_atoms = 1000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Settings.c_over_a = 1;
    Settings.Geometry = [];
    Settings.Skip_Minimization = false;

    % load the model
    Settings = Alexandria_Potential_Parameters(Settings,'vdW_Type',Theory);
    if strcmp(Theory,'JC')
        SetPot = Initialize_MD_Settings;
        SetPot.Salt = Salt;
        [JCMX,JCMM,JCXX] = JC_Potential_Parameters(SetPot);
        Settings.S.E.MM = Settings.S.E.MM/JCMM.epsilon;
        Settings.S.E.XX = Settings.S.E.XX/JCXX.epsilon;
        Settings.S.E.MX = Settings.S.E.MX/JCMX.epsilon;
        Settings.S.S.MM = Settings.S.S.MM/JCMM.sigma;
        Settings.S.S.XX = Settings.S.S.XX/JCXX.sigma;
        Settings.S.S.MX = Settings.S.S.MX/JCMX.sigma;
    end
    Settings.GaussianCharge = true; % Turn on Gaussian distributed charges when true
    Settings.Polarization = false; % Turn on polarizible Drude model when true
    Settings.niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
    Settings.emtol_polarization = 1e-1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

    % Barostat Options
    Settings.Isotropy = 'isotropic';
    Settings.Target_P = [1 1]; % Bar
    Settings.Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Settings.Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
    Settings.Nstpcouple = Get_nstcouple(Settings.Time_Constant_P,Settings.MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
    Settings.ScaleCompressibility = 1;% 1.128664;

    % Thermostat Options
    Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Settings.Nsttcouple = Get_nstcouple(Settings.Time_Constant_T,Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 
    Settings.Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Settings.MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    
    Settings.RunPostProcessor = false;
    Setup_LiX_Simulation(Settings);
    
    % Postprocessing
    Settings.SavePredictionsImage = true;
    Settings.SaveFeatures = false; % Save structure fraction vs time image when true for each temperature check
    Settings.SavePredictions = false; % Save structure fraction vs time image when true for each temperature check
    Settings.SaveTrajectory = true;
    Settings.TimePerFrame = 1; % ps
    
    [Settings.WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
    for jdx = 1:length(NN_Times)
        Settings.ML_TimeLength = NN_Times{jdx}(1);
        Settings.ML_TimeStep   = NN_Times{jdx}(2);
        for kdx = 1:length(Qlm_Voronoi)
            Settings.Qlm_Average = Qlm_Voronoi{kdx}(1);
            Settings.Voronoi = Qlm_Voronoi{kdx}(2);
            Classifier = ['NN_' num2str(Settings.ML_TimeLength) '_' num2str(Settings.ML_TimeStep)];
            if Settings.Qlm_Average
                Classifier = [Classifier 'a'];
            elseif Settings.Voronoi
                Classifier = [Classifier 'v'];
            end
            Results.(Structure).(Classifier) = Alkali_Halide_Classifier_Tester(Settings);
        end
    end
end
save(fullfile(project,Settings.Project_Directory_Name,ResultsFile),'Results');
result_json = jsonencode(Results,'PrettyPrint',true);
fid = fopen(fullfile(project,Settings.Project_Directory_Name,JsonResultsFile),'wt');
fprintf(fid,result_json);
fclose(fid);

% [Settings.WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
% MD_Postprocessor(Settings)

clear;
%Shared_Settings.Use_Conv_cell = true; % When true, use the conventional unit cell, when false use the primitive unit cell
% 'Rocksalt' 'Wurtzite' 'Sphalerite' 'BetaBeO' 'FiveFive' 'CsCl' 'NiAs' 'AntiNiAs'
Structures = {'AntiNiAs'};

for idx = 1:numel(Structures)
    Structure = Structures{idx};

    % Load Settings object and set some shared calculation settings
    T0 = 0; % Initial temperature

    Shared_Settings = Initialize_MD_Settings;
    Shared_Settings.BatchMode = false; % Sets up batch job when true, or runs immediately when false
    Shared_Settings.MPI_Ranks = 1;
    Shared_Settings.OMP_Threads = 8;
    Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
    Shared_Settings.Find_Min_Params = true;
    Shared_Settings.MinMDP.Parallel_Min = false;
    Shared_Settings.MinMDP.Verbose = true;
    Shared_Settings.Output_Coords = 1000; % Number of steps between outputting coordinates
    Shared_Settings.Calc_Energies = 100; % Number of steps that elapse between calculating the energies.
    Shared_Settings.Output_Energies = 100; % Number of steps that else between writing energies to energy file.
    
    Shared_Settings.Theory = 'JC'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Shared_Settings.Salt = 'LiI'; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
    Shared_Settings.Structure = Structure; % 
    Shared_Settings.JobID = 'Struc_Check'; % An ID that is tacked onto the folder name of all current jobs
    Shared_Settings.N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
    Shared_Settings.c_over_a = 1;
    Shared_Settings.RunPostProcessor = true;

    % Barostat options
    Shared_Settings.Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
    Shared_Settings.Isotropy = 'isotropic';
    Shared_Settings.Thermal_Solid = true;
    Shared_Settings.PreEquilibration = 0;
    Shared_Settings.Target_P = 1.0; % Target pressure in bar
    Shared_Settings.Time_Constant_P = 1.0; % [ps] time constant for coupling P. Should be 20*Nstpcouple*timestep
    Shared_Settings.Nstpcouple = Get_nstcouple(Shared_Settings.Time_Constant_P,Shared_Settings.MDP.dt);

    % Cutoff options
    Shared_Settings.MDP.RVDW_Cutoff = 1.4; % nm
    Shared_Settings.MDP.RCoulomb_Cutoff = 1.4; % nm
    Shared_Settings.MDP.RList_Cutoff = 1.4; % nm
    Shared_Settings.MDP.Disp_Correction = true;

    % Thermostat Options
    Shared_Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
    Shared_Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
    Shared_Settings.Nsttcouple = Get_nstcouple(Shared_Settings.Time_Constant_T,Shared_Settings.MDP.dt); %[ps] The frequency for coupling the temperature. 
    Shared_Settings.Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
    Shared_Settings.MDP.Initial_T = T0; % Initial termpature at which to generate velocities
    Shared_Settings.T0 = T0; % K, Initial temperature


    Shared_Settings.Annealing = 'single'; % Options: 'no' 'single' 'periodic'
    Shared_Settings.Annealing_Times = [0 100]; % [ps] A list with the number of annealing reference/control points used
    Shared_Settings.Annealing_Temps = [0 1000]; % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.
    Shared_Settings.MDP.Trajectory_Time = 0.1; % ns

    Shared_Settings.Expand_a_SC = 1.5;
    Shared_Settings.Expand_b_SC = 1.5;
    Shared_Settings.Expand_c_SC = 1.5;

    % Postprocessing
    Shared_Settings.SavePredictionsImage = true;
    Shared_Settings.ML_TimeLength = 0;
    Shared_Settings.ML_TimeStep = 0;
    Shared_Settings.TimePerFrame = 1; % ps
    Shared_Settings.SaveFeatures = false; % Save structure fraction vs time image when true for each temperature check
    Shared_Settings.SavePredictions = false; % Save structure fraction vs time image when true for each temperature check

    if ~isfield(Shared_Settings,'WorkDir')
        [Shared_Settings.WorkDir,Shared_Settings.JobName,Shared_Settings.Full_Model_Name] = GetMDWorkdir(Shared_Settings);
    end

    Setup_LiX_Simulation(Shared_Settings)
    Shared_Settings.Qlm_Average = false;
    Shared_Settings.Voronoi = true;
    MD_Postprocessor(Shared_Settings)
    Shared_Settings.Qlm_Average = false;
    Shared_Settings.Voronoi = false;
    MD_Postprocessor(Shared_Settings)
end
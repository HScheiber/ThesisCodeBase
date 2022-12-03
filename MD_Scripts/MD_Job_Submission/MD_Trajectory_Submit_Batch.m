%%%%% MD_Trajectory_Submit_Batch %%%%%%

%% Global Calculation settings
clear;
[home,project,computer,slurm] = find_home;

% Set up calculations below
skip_calculations = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs. NOTE: Does not work on sockeye

% Load Settings object and set some shared calculation settings
Shared_Settings = Initialize_MD_Settings;
Shared_Settings.N_Calc = 4; % Number of jobs to link together.
Shared_Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Shared_Settings.BatchMode = true; % Sets up batch job when true, or runs immediately when false
Shared_Settings.Hours = 6; % Max time for each job (hours)
Shared_Settings.Mins = 0; % Max time for job (minutes)
Shared_Settings.Nodes = 1; % Minimum number of cores to request for calculation.
Shared_Settings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Shared_Settings.MPI_Ranks = -1;
Shared_Settings.OMP_Threads = 1;
Shared_Settings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Shared_Settings.SinglePrecision = false; % choose true for single precision mode, false for double
Shared_Settings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Shared_Settings.npme = []; % Number of rank assigned to PME
Shared_Settings.dd = []; % Domain decomposition
Shared_Settings.Project_Directory_Name = 'Molten_Salts_MD'; % Name of project directory to contain job within the main project folder
Shared_Settings.MinMDP.nsteps_min = 1000;
Shared_Settings.TimePerFrame = 1; % post-processing time per frame check in ps
Shared_Settings.Output_Coords = 1000; % output coords every 1 ps
Shared_Settings.MinMDP.Parallel_Min = false;

% Initial calculation index
idx=0;

switch lower(computer)
    case {'unbearabull' 'cedar'}
    case 'graham'
        %% NaCl/WBK polarized Alexandria model at Tm starting in rocksalt vs liquid structure: NPT, checking gromacs version
        Salt = 'NaCl';
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Rocksalt' 'Liquid'};
        NaCl_T0 = 1075.168; % exp melting point
        
        for jdx = 1:numel(Structures)
            Structure = Structures{jdx};
            
            idx = idx+1;
            Settings_array(idx) = Shared_Settings;
            Settings_array(idx).RunEnergyAnalysis = {'Temperature' 'Pressure' 'Volume' 'Enthalpy'};
            Settings_array(idx).Structure = Structure;
            Settings_array(idx).gmx_version = '2019'; % One of: '', '4.6.7' or '2019'.
            Settings_array(idx).N_Calc = 8; % Number of jobs to link together.
            Settings_array(idx).Nodes = 0; % Number of nodes to request for calculation.
            Settings_array(idx).Cores = 1; % Minimum number of cores to request for calculation. Set to -1 for entire node
            Settings_array(idx).Hours = 6; % Max time for each job (hours)
            Settings_array(idx).MPI_Ranks = 1;
            Settings_array(idx).Mempernode = '15000M';

            Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
            Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
            Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
            Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
            Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

            Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
            Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

            Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
            Settings_array(idx).JobID = 'GmxVersionChk'; % An ID that is tacked onto the folder name of all current jobs
            Settings_array(idx).N_atoms = 2000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
            Settings_array(idx).c_over_a = 1;

            % load the model
            Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
            Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
            Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true
            
            % Polarization settings
            Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
            Settings_array(idx).emtol_polarization = 0.1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence
            
            % Misc settings
            Settings_array(idx).Output_Energies = 1; % Number of steps that else between writing energies to energy file.
            Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
            Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
            Settings_array(idx).PreEquilibration = 50; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
            Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
            Settings_array(idx).MDP.Trajectory_Time = 1.0; % ns

            % Barostat Options
            Settings_array(idx).Isotropy = 'isotropic';
            Settings_array(idx).Target_P = 1; % Bar
            Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
            Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
            Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
            Settings_array(idx).UseMoltenCompressibility = false;

            % Thermostat Options
            Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
            Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
            Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
            Settings_array(idx).Target_T = NaCl_T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            Settings_array(idx).MDP.Initial_T = NaCl_T0; % Initial termpature at which to generate velocities
        end
        
        %% LiI and NaCl/WBK polarized Alexandria model at Tm starting in rocksalt vs liquid structure: NVE, 0.5 fs time step
        Salts = {'NaCl' 'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Rocksalt' 'Liquid'};
        Densities_Liq = [16.28021241 13.31010631];
        a_RSs = [5.857263144 6.1675256591]; % Angstoms
        T0s = [1075.168 742]; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_0.5ts'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 0.1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.0005; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 1.0; % ns

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NVE, 0.25 fs time step
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        a_RSs = 6.1675256591; % Angstoms
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_0.25ts'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 0.1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.0001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns
                Settings_array(idx).MDP.dt = 0.00025; % Time step in ps for md type calculations

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NVE, 1e-2 emtol
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        a_RSs = 6.1675256591; % Angstoms
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_1e-2_emtol'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-2; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NVE, 1e-3 emtol
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        a_RSs = 6.1675256591; % Angstoms
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_1e-3_emtol'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-3; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NVE, 1e-3 emtol and 5000 steps
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        a_RSs = 6.1675256591; % Angstoms
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_5k_min'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 5000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-3; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NVE, default: 1e-1 emtol and 1000 steps
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        a_RSs = 6.1675256591; % Angstoms
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            a_RS = a_RSs(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Potential' 'Total-Energy' 'Temperature' 'Pressure'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNVE_def'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 1; % ns
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations

                % Barostat Options
                Settings_array(idx).Barostat = 'no'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)

                % Thermostat Options
                Settings_array(idx).Thermostat = 'no'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NPT: 1e-2 emtol and 1000 steps
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Temperature' 'Pressure' 'Volume' 'Enthalpy'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNPT_1e-2_emtol'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-2; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 10; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 1000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                
                % Barostat Options
                Settings_array(idx).Isotropy = 'isotropic';
                Settings_array(idx).Target_P = 1; % Bar
                Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
                Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
                Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
                Settings_array(idx).UseMoltenCompressibility = false;

                % Thermostat Options
                Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
                Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
            end
        end
        
        %% LiI/WBK polarized Alexandria model at Tm starting in liquid structure: NPT: 1e-1 emtol, 1000 steps, 0.5 fs timestep
        Salts = {'LiI'};
        Theory = 'BF';
        vdW_Type = 'WBK';
        Structures = {'Liquid'};
        Densities_Liq = 13.31010631;
        T0s = 742; % [K] exp MPs
        
        for sdx = 1:numel(Salts)
            Salt = Salts{sdx};
            T0 = T0s(sdx);
            Density_Liq = Densities_Liq(sdx);
            
            for jdx = 1:numel(Structures)
                Structure = Structures{jdx};

                idx = idx+1;
                Settings_array(idx) = Shared_Settings;
                Settings_array(idx).RunEnergyAnalysis = {'Temperature' 'Pressure' 'Volume' 'Enthalpy'};
                Settings_array(idx).Structure = Structure;

                switch Structure
                    case 'Rocksalt'
                        Settings_array(idx).Geometry = Default_Crystal(Settings_array(idx),'Center_Coordinates',true);
                        Settings_array(idx).Geometry.a = a_RS;
                        Settings_array(idx).Geometry.b = a_RS;
                        Settings_array(idx).Geometry.c = a_RS;
                        Settings_array(idx).Skip_Minimization = true;
                    case 'Liquid'
                        Settings_array(idx).Ref_Density = Density_Liq; % ion pairs / nm^3
                        Settings_array(idx).Skip_Minimization = false;
                end

                Settings_array(idx).MDP.RVDW_Cutoff = 1.0; % nm
                Settings_array(idx).MDP.RCoulomb_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.RList_Cutoff = 1.1; % nm
                Settings_array(idx).MDP.Disp_Correction = true; % Adds in long-range dispersion correction
                Settings_array(idx).MinMDP.Disp_Correction = true; % Adds in long-range dispersion correction

                Settings_array(idx).Theory = Theory; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
                Settings_array(idx).Salt = Salt; % Input model(s) to use: JC, JC3P, JC4P, TF, BH

                Settings_array(idx).Model = 'Alexandria_pol'; % Name of the current model. Leave blank for the default JC/TF/BH model
                Settings_array(idx).JobID = 'ChkNPT_0.5ts'; % An ID that is tacked onto the folder name of all current jobs
                Settings_array(idx).N_atoms = 1728; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
                Settings_array(idx).c_over_a = 1;

                % load the model
                Settings_array(idx) = Alexandria_Potential_Parameters(Settings_array(idx),'vdW_Type',vdW_Type);
                Settings_array(idx).GaussianCharge = true; % Turn on Gaussian distributed charges when true
                Settings_array(idx).Polarization = true; % Turn on polarizible Drude model when true

                % Polarization settings
                Settings_array(idx).niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
                Settings_array(idx).emtol_polarization = 1e-1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

                % Misc settings
                Settings_array(idx).Output_Energies = 20; % Number of steps that else between writing energies to energy file.
                Settings_array(idx).Output_Coords = 2000; % Number of steps between outputting coordinates
                Settings_array(idx).MDP.dt = 0.001; % Time step in ps for md type calculations
                Settings_array(idx).PreEquilibration = 10; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
                Settings_array(idx).Annealing = 'no'; % Options: 'no' 'single' 'periodic'
                Settings_array(idx).MDP.Trajectory_Time = 0.5; % ns
                Settings_array(idx).MDP.dt = 0.0005; % Time step in ps for md type calculations
                
                % Barostat Options
                Settings_array(idx).Isotropy = 'isotropic';
                Settings_array(idx).Target_P = 1; % Bar
                Settings_array(idx).Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
                Settings_array(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
                Settings_array(idx).Nstpcouple = Get_nstcouple(Settings_array(idx).Time_Constant_P,Settings_array(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
                Settings_array(idx).UseMoltenCompressibility = false;

                % Thermostat Options
                Settings_array(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
                Settings_array(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
                Settings_array(idx).Nsttcouple = Get_nstcouple(Settings_array(idx).Time_Constant_T,Settings_array(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
                Settings_array(idx).Target_T = T0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
                Settings_array(idx).MDP.Initial_T = T0; % Initial termpature at which to generate velocities
            end
        end
        
    case 'narval'


    otherwise
end

%% Check for already running jobs
if check_running && isunix
    if slurm
        er = 1;
        jdx = 1;
        while er ~= 0
            [er,out] = system('squeue -u $USER -o %200Z | tail -n +2');
            jdx = jdx+1;
            if jdx > 10
                disp('Failed to check current running jobs.')
                out = '';
                break
            end
        end

        % Grab the currently running jobs
        out = strtrim(out);
        c_jobs = strtrim(split(out,newline));

        % Remove the outer directory and find unique jobs
        c_jobs_unique = strrep(unique(c_jobs),fullfile(Shared_Settings.project,Shared_Settings.Project_Directory_Name,filesep),'');

        % Replace file separator with '_'
        c_jobs_unique = strrep(c_jobs_unique,filesep,'_');
    else % PBS without squeue
        er = 1;
        jdx = 1;
        while er ~= 0
            [er,out] = system('qstat -u $USER');
        end
        
        if isempty(out)
            % No jobs running
            c_jobs_unique = {};
        else
            er = 1;
            jdx = 1;
            while er ~= 0
                [er,out] = system('qstatf () { qstat -f $*| tr ''\n'' ''&'' | sed -e ''s/&\t//g'' | tr ''&'' ''\n''; }; qstatfi() { grep -e "$*" | sed -e "s/^\(\S*\).*,$*=\([^,]*\).*/\1 \2/g" ;}; qstatf | qstatfi PBS_O_WORKDIR | grep $USER');
                jdx = jdx+1;
                if jdx > 10
                    disp('Failed to check current running jobs.')
                    out = '';
                    break
                end
            end
            c_jobs = strtrim(splitlines(strtrim(out)));

            % Remove the outer directory and find unique jobs
            c_jobs_unique = strrep(unique(c_jobs),fullfile(Shared_Settings.project,Shared_Settings.Project_Directory_Name,filesep),'');

            % Replace file separator with '_'
            c_jobs_unique = strrep(c_jobs_unique,filesep,'_');
        end
    end
else
    c_jobs_unique = {};
end

%% Prepare batch scripts
% Loop through all models and build their submission file
for idx = 1:length(Settings_array)
    if any(idx == skip_calculations)
        continue
    end
    Settings = Settings_array(idx);
    
    [WorkDir,JobName,~] = GetMDWorkdir(Settings);
    TaskName = [Settings.Salt '_' JobName];
    
    % Check if job is already complete (including postprocessing)
    OutConfFile = fullfile(WorkDir,[JobName '_OutConf.' Settings.CoordType]);
    OutConfFile2 = strrep(OutConfFile,'scratch','project');
    PostProcessFlagFile = fullfile(WorkDir,'POSTPROCESS_COMPLETE');
    PostProcessFlagFile2 = strrep(PostProcessFlagFile,'scratch','project');
    if check_complete && ( isfile(OutConfFile) || isfile(OutConfFile2) )
        Settings.N_Calc = 1;
        if (Settings.RunPostProcessor && (isfile(PostProcessFlagFile) || isfile(PostProcessFlagFile2)) ) || ~Settings.RunPostProcessor
            disp([TaskName ': Job already completed. Skipping Job Submission.'])
            continue
        end
    end
    
    % Check if Model is already running
    if check_running && ismember(TaskName,c_jobs_unique) 
        disp([TaskName ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    % Finish job setup
    Setup_LiX_Simulation(Settings);
end
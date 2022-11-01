function Settings = Initialize_MD_Settings
% Script to generate default molecular dynamics inputs for simulations

% A script for generating MD inputs for Lithium halides salts using GROMACS

% Damp_Types:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% ModelID is a ID to keep track of your input model. Will be added to your
% input model name. Can be a string or a numeric.

% Job name convention: 
%$MAINDIR/Salt/(Date)_(Initial Struc Label)_(Model)_(ModelID)_(NVE or NPT)/(Job Files)

%% Job Resource and Gromacs mdrun Settings
Settings.Submit_Jobs = true; % Set to true to submit MD jobs to batch script or to run locally, otherwise just produce input files.
Settings.BatchMode = true; % Sets up batch job when true, or runs immediately when false
Settings.Continue = false; % Where possible, continue a FINISHED calculation when true
Settings.JobSettings.Hours = 3; % Max time for each job (hours)
Settings.JobSettings.N_Calc = 1; % Number of jobs to link together.
Settings.JobSettings.Mins = 0; % Max time for job (minutes)
Settings.JobSettings.Nodes = 1; % Minimum number of cores to request for calculation.
Settings.JobSettings.Cores = -1; % Minimum number of cores to request for calculation. Set to -1 for entire node
Settings.JobSettings.Mempernode = '0'; % Memory request for server (default = '-1', max per core = '0', eg '3G' for cedar or 3gb for sockeye)
Settings.JobSettings.SinglePrecision = false; % choose true for single precision mode, false for double
Settings.JobSettings.BigNode = false; % For cedar and sockeye, choose the large node types when true.
Settings.JobSettings.MPI_Ranks = -1; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Settings.JobSettings.npme = []; % Number of MPI ranks assigned to PME
Settings.JobSettings.dd = []; % Sets the domain decomposition cells. Requires 3 numbers, X, Y, and Z corresponding to the number of cells in each direction.
% Leave dd empty to not use.
Settings.JobSettings.dds = 0.8; % default = 0.8. Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to provide a margin in which dynamic load balancing can act while preserving the minimum cell size.
% dds may need to be set lower for cluster jobs due to high inhomogeniety in the system
Settings.JobSettings.DLB = true; % Turns on Gromacs dynamic load balancing when true
Settings.JobSettings.TunePME = true; % Optimizes PME load between PP/PME ranks or GPU/CPU
Settings.Verbose = true;

%% Interaction and Structure Settings
% These can be arrays
Settings.Salt = 'LiF'; % 'LiF' 'LiCl' 'LiBr' 'LiI'
Settings.Structure = 'Rocksalt'; % One of: 'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'Liquid' 'Previous'
Settings.Use_Conv_cell = true; % When true, use the conventional unit cell, when false use the primitive unit cell
Settings.Liquid_Interface = false; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
Settings.Liquid_Fraction = 0.5; % Only meaninful when Liquid_Interface = true. Sets the approximate fraction of the total number of atoms that will initialize as Liquid
Settings.Theory = ''; % Input model(s) to use: JC, JC3P, JC4P, TF, BH
Settings.Model = ''; % Name of the current model. Leave blank for the default JC/TF/BH model
Settings.WaterModel = 'SPC/E'; % For use with JC model
Settings.JobID = ''; % An ID that is tacked onto the folder name of all current jobs
Settings.N_atoms = 5000; % Minimum number of atoms to include in box or size of search box for cluster jobs. This will automatically resize as needed
Settings.Cutoff_Buffer = 1.2; % Set this value to some fraction greater than 1 to allow NPT simulation boxes to shrink without making the box size minimal axis smaller than the cutoff. Also allows for dynamic setting of rlist.
Settings.ScaleInitialLiqDensity = 1; % Scale initial guess liquid density by this amount

% Advanced forcefield options
Settings.GaussianCharge = false; % Turn on Gaussian distributed charges when true
Settings.Polarization = false; % Turn on polarizible Drude model when true
Settings.niter_polarization = 1000; % Maximum number of iterations for optimizing the shell positions
Settings.emtol_polarization = 0.1; % [kJ/(mol nm)] A tolerance for self consistent polarization convergence

% Only active for 'previous' structures:
Settings.Prev_geom_loc = '';
Settings.Prev_geom_time = []; % When selecting 'previous' this selects the frame to pick in ns
Settings.Target_Density = []; % [molecules/nm^3] Target density for creating liquid clusters

Settings.TF_Paramset = 0; % choose from 0-3
% TF Parameter sets for C6/C8 coefficients (for when using as a script)
% 0 = default TF Parameters
% 1 = D3 values
% 2 = Best literature values available
% 3 = D4 with C6 and C8 generated on-the-fly

%% Thermostat Options
Settings.Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'v-rescale' (set NO for NVE)
Settings.Target_T = 0; % Target temperature in kelvin. Does not apply when thermostat option 'no' is chosen
Settings.Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be 20*Nsttcouple*timestep
Settings.Nsttcouple = 10; %[ps] The frequency for coupling the temperature. 
Settings.MDP.Initial_T = 0; % Initial termpature at which to generate velocities

%% Barostat Options
Settings.Barostat = 'Berendsen'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
% Note Berendsen equilibrates faster but Parrinello-Rahman produces more realistic density fluctuations
Settings.Isotropy = 'isotropic'; % One of: isotropic, semiisotropic, anisotropic, surface-tension
Settings.Target_P = 1.0; % Target pressure in bar
Settings.Time_Constant_P = 1.0; % 0.2 [ps] time constant for coupling P. Should be 20*Nstpcouple*timestep
Settings.Nstpcouple = 50; % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps
Settings.ScaleCompressibility = 1; % Use this to raise or lower the barostat compressibility from the default experimental isothermal compressibility
Settings.UseMoltenCompressibility = false; % When this is true, uses the experimental molten salt compressibility rather than the solid compressibility

%% Simulated Annealing Options
Settings.Annealing = 'no'; % Options: 'no' 'single' 'periodic'
Settings.Annealing_Times = []; % [ps] A list with the number of annealing reference/control points used
Settings.Annealing_Temps = [];   % [K] A list of temperatures at the annealing reference/control points used. Must be equal in length to previous line.

%% MDP SETTINGS
Settings.MDP.Trajectory_Time = 1; % Trajectory time in nanoseconds. Set to 0 for single point energy calculation.
Settings.MDP.dt = 0.001; % Time step in ps for md type calculations
Settings.MDP.integrator = 'md'; % What type of calculation is run for single point energy calculations (steep = energy min, md = molecular dynamics)
Settings.MDP.LJtol = 1e-5; % When doing PME for VdW-interactions, this is used to control the relative strength of the dispersion potential at rvdw in the same way as ewald-rtol controls the electrostatic potential.
Settings.MDP.CutOffScheme = 'Verlet'; % Either 'group' or 'Verlet' (does NOT apply to tabulated potentials, these are set to group)
Settings.MDP.VerletBT = 0.005; %  (0.005) [kJ mol-1 ps-1]This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets rlist unless set to -1, in which case rlist will be used.
Settings.MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME', 'PME-User', or 'Cut-off'
Settings.MDP.VDWType = 'Cut-off'; % Define the type of van der waals potential used. One of 'PME' or 'Cut-off'
Settings.MDP.RList_Cutoff = 1.5; % 1.5 % nm. This should be larger or equal to RCoulomb/RVDW
Settings.MDP.RCoulomb_Cutoff = 1.5; % nm. if set to less than 0, then Rc = a;
Settings.MDP.RVDW_Cutoff = 1.4; % 1.2 nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1
Settings.MDP.Fourier_Spacing = 0.12; % used 0.1 in minimization. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
Settings.MDP.PME_Order = 4; % Interpolation order for PME. 4 equals cubic interpolation (default).
Settings.MDP.Ewald_rtol = 1e-5; %1e-7 Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
% Also note that Ewald_rtol ≈ erfc(β × rcoulomb) where β (AKA alpha) is the “Ewald width parameter.”
Settings.MDP.continuation = 'no'; % Yes = do not apply constraints to the start configuration and do not reset shells, useful for exact coninuation and reruns
Settings.MDP.Num_Groups = 1;
Settings.MDP.vdw_modifier = 'Potential-shift'; % Potential-shift-Verlet, Potential-shift, None, Force-switch, Potential-switch
Settings.MDP.Disp_Correction = true; % Adds in a long-range dispersion correction when true
Settings.Table_Length = 4.01; % How far should tabulated potentials extend in nm
Settings.CoordType = 'gro'; % Either pdb, gro, or g96 (use g96 for extra precision)

%% Trajectory update/output options: larger number of steps decreases resolution
Settings.Output_Coords = 1000; % Number of steps between outputting coordinates
Settings.Output_Coords_Compressed = 0; % Number of steps between outputting coordinates in compressed format
Settings.Output_Velocity = 0; % Number of steps between outputting velocities
Settings.Output_Forces = 0; % Number of steps between outputting forces
Settings.Calc_Energies = 100; % Number of steps that elapse between calculating the energies.
Settings.Output_Energies = 1000; % Number of steps that else between writing energies to energy file.
Settings.Update_NeighbourList = 10; % Frequency to update the neighbor list. When this is 0, the neighbor list is made only once.

%% Pre-Minimization Settings
Settings.Skip_Minimization = false; % TURNED ON FOR LIQUIDS AUTOMATICALLY. When true, skip the pre-minimization step 
Settings.Find_Min_Params = true; % When true, finds lowest energy parameters for IC based on Data_Types. When false, uses input IC
Settings.Find_Similar_Params = false; % If no exact minimized geometry is found for the input model, find geometry for a similar model
Settings.Data_Types = 1; % Allowed data types for automatic search of initial conditions (0 = normal, 1 = cell optimized, 2 = full optimized, 3 = atom optimized only)
Settings.MinMDP.Parallel_Min = true; % run minimization routine using matlab parallel when true
Settings.MinMDP.Verbose = false; % Sets verbosity of minimization routine
Settings.MinMDP.OptPos = false; % Optimize both position and lattice parameters when true
Settings.MinMDP.Maintain_Symmetry = true; % maintains cell symmetry when true
Settings.MinMDP.emtol = 0.088439; %1e-3 [kJ mol-1 nm-1] The position minimization is converged (per cycle) when the maximum force is smaller than this value
Settings.MinMDP.MaxCycles = 100; % Maximum number of optimization cycles.
Settings.MinMDP.Energy_Tol = 1; % kJ/mol. Convergence criteria 1: must have change in energy between cycles is less than this value
Settings.MinMDP.Gradient_Tol_RMS = 1e-2; %kJ/(mol A) Convergence criteria 2: must have RMS cell parameter gradient less than this value for convergence
Settings.MinMDP.Gradient_Tol_Max = 1e-2; %kJ/(mol A) Convergence criteria 3: must have maximum cell parameter gradient less than this value for convergence
Settings.MinMDP.E_Unphys = -2000; % [kJ/mol] unphysical energy cutoff
Settings.MinMDP.nsteps_point = 0; % Number of steps to perform before ending (should be 0 for single point energy calculations)
Settings.MinMDP.point_integrator = 'md'; % What type of calculation is run for single point energy calculations (steep = energy min, md = molecular dynamics)
Settings.MinMDP.dt = 0.001; % Time step in ps for md type calculations
Settings.MinMDP.min_integrator = 'steep'; % 'steep', 'cg', or 'l-bfgs'
Settings.MinMDP.nsteps_min = 1000; % Number of steps to perform before stopping for energy minimization runs
Settings.MinMDP.CutOffScheme = 'Verlet'; % Either 'group' or 'Verlet' (does NOT apply to tabulated potentials, these are set to group)
Settings.MinMDP.VerletBT = -1; % This sets the maximum allowed error for pair interactions per particle caused by the Verlet buffer, which indirectly sets rlist unless set to -1, in which case rlist will be used.
Settings.MinMDP.VDWType = 'Cut-off'; % Define the type of van der waals potential used. One of 'PME', 'Cut-off', 'Shift', or 'Switch'
Settings.MinMDP.Auto_Cutoff = false; % Set auto or manual cutoff distances in nm (true = automatically half the smallest box length)
Settings.MinMDP.RList_Cutoff = 1.9; %1.9 nm. This should be larger or equal to RCoulomb/RVDW
Settings.MinMDP.RCoulomb_Cutoff = 1.9; %1.9 nm. if set to less than 0, then Rc = a;
Settings.MinMDP.RVDW_Cutoff = 1.9; %1.9 nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1
Settings.MinMDP.rvdw_switch = 1.8; %1.8 nm. where to start switching the LJ potential. Only applies when VDWType = switch
Settings.MinMDP.Fourier_Spacing = 0.10; % nm. Default 0.12 nm. Grid dimensions in PME are controlled with fourierspacing
Settings.MinMDP.PME_Order = 4; % Interpolation order for PME. 4 equals cubic interpolation (default).
Settings.MinMDP.Ewald_rtol = 1e-7; %1e-7 Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
Settings.MinMDP.vdw_modifier = 'Potential-shift'; % Potential-shift-Verlet, Potential-shift, None, Force-switch, Potential-switch
Settings.MinMDP.CoordType = 'g96';
Settings.MinMDP.Disp_Correction = false; % Adds in a long-range dispersion correction when true
Settings.MinMDP.Disp_Correction_Tables = false; % ONLY TURN THIS ON FOR DEBUGGING, can be used to enable disp correction for tables
Settings.MinMDP.Use_Conv_cell = false; % When true, use the conventional unit cell, when false use the primitive unit cell


Settings.PreEquilibration = 0; % [ps] Use this to run an initial rapid equilibration

%% Geometry editing and cluster settings
Settings.pbc = 'on'; % periodic boundary conditions
Settings.GenCluster = false; % Turns on cluster generation when true
Settings.Cluster_Frac = 0.5; % Fraction of the total atoms in the cluster
Settings.Expand_LP = true; % Switch to turn on/off ALL expansions
Settings.Thermal_Solid = false; % Turns on thermal expansion for solids prior to equilibration based on experimental data

% Expand a lattice parameter of the supercell by this factor.
% Creates an empty volume upon expansion. Values less than 1 not defined.
Settings.Expand_a_SC = 1;
Settings.Expand_b_SC = 1;
Settings.Expand_c_SC = 1;

% Expand a lattice parameter and all positions by this factor (FC held fixed)
Settings.Expand_a = 1;
Settings.Expand_b = 1;
Settings.Expand_c = 1;

Settings.c_over_a = 1; % relative length of c vs a (does not apply to cluster calculations)

Settings.Manual_Box = false; % When set to true, rather than setting the number of atoms in a box, user sets the a, b, and c dimensions of the box
Settings.Manual_Box_a = 1; % Box length a in nm
Settings.Manual_Box_b = 1; % Box length b in nm
Settings.Manual_Box_c = 1; % Box length c in nm
Settings.Manual_Radius = 1;% for cluster calculations, this sets the radius of the solid sphere (nm)


%% I/O Job Settings
Settings.Delete_MDPout = false; % Automatically delete MDP out file if true
Settings.Delete_MDlog = false; % Delete MD log file if true
Settings.Delete_ConfOut = false; % Delete the output configuration if true
Settings.Delete_TPR = false; % Delete the tpr file after used if true
Settings.Delete_Backups = true; % Automatically delete any gromacs backup files found if true
Settings.Delete_Minimize = false; % Deletes the pre-minimization folder when true
Settings.Delete_Equil = true; % Deletes the equilibration folder, if it exists, when true
Settings.Delete_outputs = false; % delete stde and stdo after job finishes
Settings.Delete_subms = false; % delete subm scripts after job finishes
Settings.Delete_cpt = true; % delete all but one last checkpoint after job finishes
Settings.MaxWarn = 2; % Maximum allowed warnings for gmx grompp

%% Topology settings
Settings.Top_gen_pairs = 'no'; % Automatically generate pairs
Settings.Top_fudgeLJ = 1.0; % Rescale LJ interaction by this amount for 1-4 bonded atoms
Settings.Top_fudgeQQ = 1.0; % Rescale Coulomb interaction by this amount for 1-4 bonded atoms
Settings.Table_StepSize = 0.0005; % nm, suitable for double precision

[Settings.home,Settings.project,Settings.computer,Settings.slurm,...
    Settings.BO_Models,Settings.qsub,Settings.passlog,Settings.pipe,...
    Settings.wsl,Settings.MLModelDir,Settings.scratch_dir] = find_home;
Settings.Project_Directory_Name = 'Melting_Point_Studies'; % Name of project directory to contain job within the main project folder

%% Settings for melting point calculations
Settings.CheckTime = 25; % ps. Time between checking for melting/freezing
Settings.MaxCheckTime = 1000; % ps. Max time for melting/freezing runs
Settings.T0 = 800; % K, Initial temperature
Settings.lb = 500; % K, lower bound on MP search
Settings.ub = 2000; % K, upper bound on MP search
Settings.MeltFreezeThreshold = 0.25; % Required CHANGE in the fraction of the box either frozen or melted to flag the end of the simulation
% The default SlopeThreshold is equivalent to a 10 % change in fraction of liquid over 1000 ps
Settings.SlopeThreshold = 1e10; % A second requirement: the change in the fraction per unit time must be smaller than the absolute value of this threshold for the system to be considered at the melting point. Units of [Structure Fraction/ps]
Settings.MaxMPIterations = 100; % Maximum iterations for MP calculations
Settings.StepTolerance = 0; % K, tolerence for the melting point step size
Settings.MaxTDiff = 0.01; % K, maximum change in temperature between points before selecting new initial conditions
Settings.FinalDensityProfile = true; % Set to true to run a final density profile along the Z-dimension
Settings.Delete_T_History = false; % deletes the intermediate temperature check files when true
Settings.MinInterfaceWidth = 0.075; % [nm] +- distance from the solid-liquid interface within which to minimize
Settings.RefStructure = 'Liquid'; % Reference structure used for determination of melting or freezing
Settings.QECompressibility = 1e-6; % Compressibility used during the rapid-equilibration stages
Settings.CheckAmorphousLiquid = true; % When true, this enables a check for liquid -> amorphous based on the mean-squared displacement
Settings.CheckAmorphousHalide = false; % (specific to non-MP liquid calc) When false, only check the metal diffusion, when true: check both metal and halide diffusion
Settings.AmorphousDiffThreshold = 1e-6; % [cm^2/s] When CheckAmorphousLiquid is true, this sets the threshold for amorphous vs liquid

% optimizer settings
Settings.Optimizer = 'MPSearcher'; % One of fmincon, patternsearch, MPSearcher, or bayesopt
Settings.ExplorationRatio = 1; % only applies to bayesopt
Settings.Acquisition_Function = 'expected-improvement'; % only applies to bayesopt
Settings.MaxMPIterations = 10; % Maximum iterations for MP calculations, does not apply to MPSearcher optimizer.
Settings.N_Seed_Points = 2; % Number of bayesopt seed points
Settings.KernelFunction = 'ardmatern52';
Settings.InitialMeshSize = 100; % Applies to patternsearch and MPsearcher
Settings.MeshSizeMultiplier = 2; % During the first phase of MPsearcher, multiply the stepsize by this value after each attempt
Settings.FunctionTolerance = 0;
Settings.BracketThreshold = 5; % [K] Sets the target bracket for the melting point
Settings.MinStepSize = 0.25; % [K] Sets the minimum step size between iterations in K
Settings.UseDerivativeWeighting = false; % Set true to use derivative-weighted midpoints
Settings.IgnoreBounds = false; % Ignores recording the MP bounds if set to true. Usually only used on the final iteration.

%% Postprocessor settings
Settings.RunPostProcessor = true; % Turns on or off postprocessing of the trajectory
Settings.SaveTrajectory = true; % Save processed xyz trajectory file for each temperature check when true
Settings.SavePredictionsImage = true; % Save structure fraction vs time image when true for each temperature check
Settings.SaveFeatures = false; % Save structure fraction vs time image when true for each temperature check
Settings.SavePredictions = false; % Save structure fraction vs time image when true for each temperature check
Settings.ML_TimeLength = 21;
Settings.ML_TimeStep = 5;
Settings.TimePerFrame = 5; % ps
Settings.Qlm_Average = true;
Settings.Voronoi = false;


% Default model modification parameters
Settings.S = Init_Scaling_Object;
Settings.CR_Damp = Init_CRDamping_Object; % note: b = steepness of damping, r_d = position in nm.
Settings.C6_Damp = Init_C6Damping_Object;
[Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;

% Finite T settings
Settings.Liquid_Test_Time = 100; % ps. Sets time for primary liquid calculation
Settings.Liquid_Equilibrate_Time = 25; % ps
Settings.Solid_Test_Time = 30; % ps. Sets time for combination of solid equilibration + test time

% MP settings
Settings.MP_Liquid_Test_Time = 100; % ps. Time used for calculation of liquid MSD in melting point calculations.
Settings.MP_Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MP_Equilibrate_Liquid = 20; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
Settings.MinTimeStep = 0.00025; % [ps] minimum time step

end
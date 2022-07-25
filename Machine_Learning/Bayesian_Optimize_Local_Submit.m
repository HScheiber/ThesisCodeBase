%%%%% Bayesian_Optimize_Local_Submit %%%%%%
%% INFO %%
%% Primary loss options
% Model.Loss_Options.regularization = 'L2'; % Set the regularization scheme
% Model.Loss_Options.Rocksalt.LE = 1;
% Model.Loss_Options.Wurtzite.RLE = 1;
% Model.Loss_Options.NiAs.RLE = 1;
% Model.Loss_Options.Sphalerite.RLE = 1;
% Model.Loss_Options.FiveFive.RLE = 1;

% Model.Loss_Options.Rocksalt.a = 1/10;
% Model.Loss_Options.Wurtzite.a = 1/20;
% Model.Loss_Options.Wurtzite.c = 1/20;
% Model.Loss_Options.NiAs.a = 1/20;
% Model.Loss_Options.NiAs.c = 1/20;
% Model.Loss_Options.Sphalerite.a = 1/10;
% Model.Loss_Options.FiveFive.a = 1/20;
% Model.Loss_Options.FiveFive.c = 1/20;

%% Incorporating energy gaps into loss function
% Value: The target gap between reference and current structure
% Gap.Value < 0 -> "reference structure" is favoured
% Gap.Value > 0 -> "current structure" is favoured
% Gap.Value = 0 -> structures are equal in energy
% Model.Loss_Options.(Structure).Gap.Value = 0;

% Weight: The weighting for this gap in the overall loss function: 
% larger values means more weight. Do not use negative weights!
% Default weight is zero (excluded from loss function)
% Model.Loss_Options.(Structure).Gap.Weight = 0;

% Type: Comparison function. Pick one of: lt | gt | eq | ge | le | ne
% Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
% Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
% Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
% Model.Loss_Options.(Structure).Gap.Type = @lt;

% Ref: The reference structure for the gap. The difference between the
% current structure and the reference structure is considered.
% Model.Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

% Defaults:
% Model.Loss_Options.(Structure).Gap.Value = 0;
% Model.Loss_Options.(Structure).Gap.Weight = 0;
% Model.Loss_Options.(Structure).Gap.Type = @lt;
% Model.Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

%% Set up Model
clear;
idx = 0;

%% BH Model (MP Test) XY
Salts = {'LiI'};
Theories = {'BH'};
Replicates = 1;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: FE
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['XY' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end
            
            % Some Job settings
            Models(idx).JobSettings.OMP_Threads = 1;
            Models(idx).JobSettings.MPI_Ranks = 8;
            Models(idx).JobSettings.Cores = 8;
            Models(idx).Cutoff_Buffer = 1.2; % This affects Structure_Minimization as well as other aspects of code
            
            % T=0 Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            
            % Finite T loss
            Models(idx).Loss_Options.Fusion_Enthalpy = 0; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 0; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 0; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 0; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.MP  = 1; % Fitting the experimental MP, using the experimental structure as the solid
            
            % Finite T options (including MP options)
            Models(idx).c_over_a = 2;
            Models(idx).Liquid_Test_Time = 50; % ps. simulation time to sample the liquid (second half averaged for enthalpy / volume)
            Models(idx).Solid_Test_Time = 30; % ps. simulation time to sample the solid (second half averaged for enthalpy / volume)
            Models(idx).Liquid_Interface = true; % When true, creates an system with half STRUCTURE half LIQUID for melting point testing
            Models(idx).Optimizer = 'MPSearcher';
            Models(idx).lb = 0; % K, lower bound on MP search
            Models(idx).ub = 2200; % K, upper bound on MP search
            Models(idx).MeltFreezeThreshold = 0.25; % CHANGE in fraction [0,1] OR Number of atoms (1,inf) of liquid/solid required to establish a phase change
            Models(idx).BracketThreshold = 5; % [K] Sets the target bracket for the melting point
            Models(idx).N_atoms = 2000;
            Models(idx).MaxCheckTime = 5000; % ps. Max time for melting/freezing runs
            Models(idx).Equilibrate_Solid = 15; % number of ps to equilibrate the solid for, use 0 to skip. Only works for flat solid-liquid interface
            Models(idx).Equilibrate_Liquid = 10; % number of ps to equilibrate the liquid for, use 0 to skip. Only works for flat solid-liquid interface
            Models(idx).PreEquilibration = 0.3; % ps. Relax the prepared system for this amount of time at the start with ultrafast relaxation settings.
            Models(idx).InitialMeshSize = 20;
            Models(idx).MeshSizeMultiplier = 5;
            Models(idx).QECompressibility = 1e-7; % sets the compressibility during the system preparation stages
            %Models(idx).MinInterfaceWidth = 0.15; % [nm] +- distance from the solid-liquid interface within which to minimize
            Models(idx).ScaleInitialLiqDensity = 1.25; 
            
            % Barostat Options
            Models(idx).Isotropy = 'semiisotropic';
            Models(idx).Target_P = [1 1]; % Bar
            Models(idx).Barostat = 'Parrinello-Rahman'; % Options: 'no' 'Berendsen' 'Parrinello-Rahman' 'MTTK' (set NO for NVT)
            Models(idx).Time_Constant_P = 1; % 0.2 [ps] time constant for coupling P. Should be at least 20 times larger than (Nstpcouple*timestep)
            Models(idx).Nstpcouple = Get_nstcouple(Models(idx).Time_Constant_P,Models(idx).MDP.dt); % [ps] The frequency for coupling the pressure. The box is scaled every nstpcouple steps. 
            % Thermostat Options
            Models(idx).Thermostat = 'v-rescale'; % Options: 'no' 'berendsen' 'nose-hoover' 'andersen' 'andersen-massive' 'nose-hoover' (set NO for NVE)
            Models(idx).Time_Constant_T = 0.2; %[ps] time constant for coupling T. Should be at least 20*Nsttcouple*timestep
            Models(idx).Nsttcouple = Get_nstcouple(Models(idx).Time_Constant_T,Models(idx).MDP.dt); %[ps] The frequency for coupling the temperature. 
            % Electrostatics
            Models(idx).MDP.CoulombType = 'PME'; % Define the type of coulomb potential used. One of 'PME' or 'Cut-off'
            Models(idx).MDP.Ewald_rtol = 1e-5; % Default (1e-5) The relative strength of the Ewald-shifted direct potential at rcoulomb. Decreasing this will give a more accurate direct sum, but then you need more wave vectors for the reciprocal sum.
            Models(idx).MDP.Fourier_Spacing = 0.12;
            Models(idx).MDP.VerletBT = -1;
            Models(idx).MDP.RList_Cutoff = 1.1; % 1.5 % nm. This should be larger or equal to RCoulomb/RVDW
            Models(idx).MDP.RCoulomb_Cutoff = 1.1; % nm. if set to less than 0, then Rc = a;
            Models(idx).MDP.RVDW_Cutoff = 1.0; % 1.2 nm. note that rlist ? rCoulomb = RVDW when using Verlet and VerletBT = -1
            % End of Finite T options
            
            % Aux options
            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            Models(idx).Parallel_Bayesopt = false;
            Models(idx).Parallel_Struct_Min = true;
        end
    end
end

%% Submit jobs
workdir = pwd;
% Loop through all models and build their submission file
for idx = 1:length(Models)
    Model = Models(idx);
    Calc_Name = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
    Calc_Dir = [workdir filesep Calc_Name];
    Diary_Loc = [Calc_Dir filesep Calc_Name '.log'];
    
    % Move to new directory
    if ~isfolder(Calc_Dir)
        mkdir(Calc_Dir)
    end
    cd(Calc_Dir);
    
    % Turn diary on and submit job
    diary(Diary_Loc);
    Bayesian_Optimize_LiX_Parameters(Model)
    diary('off')
    close all % closes figures
end
cd(workdir);
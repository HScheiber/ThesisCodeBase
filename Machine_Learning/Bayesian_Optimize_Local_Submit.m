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

%% BH Model (Test) XX
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
            Models(idx).Trial_ID = ['XX' Rep];
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
            
            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            
            Models(idx).Loss_Options.Fusion_Enthalpy = 1; % Fitting the experimental enthalpy difference of the liquid and solid at the experimental MP
            Models(idx).Loss_Options.MP_Volume_Change = 1; % Fitting the experimental change in volume due to melting at the experimental MP
            Models(idx).Loss_Options.Liquid_MP_Volume = 1; % Fitting the experimental volume per formula unit at the experimental MP
            Models(idx).Loss_Options.Solid_MP_Volume  = 1; % Fitting the experimental volume of the experimental solid structure at the experimental MP
            Models(idx).Loss_Options.MP  = 1; % Fitting the experimental MP, using the experimental structure as the solid
            Models(idx).Equilibrate_Liquid = 20; % ps
            Models(idx).Output_Coords = 1000;
            
            
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
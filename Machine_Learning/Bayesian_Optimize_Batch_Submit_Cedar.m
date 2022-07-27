%%%%% Bayesian_Optimize_Batch_Submit_Cedar %%%%%%
%% INFO %%
%% Primary loss options
% Models(idx).Loss_Options.regularization = 'L2'; % Set the regularization scheme
% Models(idx).Loss_Options.Rocksalt.LE = 1;
% Models(idx).Loss_Options.Wurtzite.RLE = 1;
% Models(idx).Loss_Options.NiAs.RLE = 1;
% Models(idx).Loss_Options.Sphalerite.RLE = 1;
% Models(idx).Loss_Options.FiveFive.RLE = 1;

% Models(idx).Loss_Options.Rocksalt.a = 1/10;
% Models(idx).Loss_Options.Wurtzite.a = 1/20;
% Models(idx).Loss_Options.Wurtzite.c = 1/20;
% Models(idx).Loss_Options.NiAs.a = 1/20;
% Models(idx).Loss_Options.NiAs.c = 1/20;
% Models(idx).Loss_Options.Sphalerite.a = 1/10;
% Models(idx).Loss_Options.FiveFive.a = 1/20;
% Models(idx).Loss_Options.FiveFive.c = 1/20;

%% Incorporating energy gaps into loss function
% Value: The target gap between reference and current structure
% Gap.Value < 0 -> "reference structure" is favoured
% Gap.Value > 0 -> "current structure" is favoured
% Gap.Value = 0 -> structures are equal in energy
% Models(idx).Loss_Options.(Structure).Gap.Value = 0;

% Weight: The weighting for this gap in the overall loss function: 
% larger values means more weight. Do not use negative weights!
% Default weight is zero (excluded from loss function)
% Models(idx).Loss_Options.(Structure).Gap.Weight = 0;

% Type: Comparison function. Pick one of: lt | gt | eq | ge | le | ne
% Less than:    the Actual_Gap must be less than    the Target_Gap or the loss is non-zero.
% Greater than: the Actual_Gap must be greater than the Target_Gap or the loss is non-zero.
% Equal to:     the Actual_Gap must be equal to     the Target_Gap or the loss is non-zero.
% Models(idx).Loss_Options.(Structure).Gap.Type = @lt;

% Ref: The reference structure for the gap. The difference between the
% current structure and the reference structure is considered.
% Models(idx).Loss_Options.(Structure).Gap.Ref = 'Rocksalt';

% Defaults:
% Models(idx).Loss_Options.(Structure).Gap.Value = 0;
% Models(idx).Loss_Options.(Structure).Gap.Weight = 0;
% Models(idx).Loss_Options.(Structure).Gap.Type = @lt;
% Models(idx).Loss_Options.(Structure).Gap.Ref = 'Rocksalt';


%% Additional_GAdjust options: 'MM' 'MX' 'XX'
% Introduces additional inverted Gaussians with 3 additional parameters each. 
% Leave empty to exclude
% e.g.: Settings.Additional_GAdjust = {'MM'};
%
% Additional_GAdjust_Ranges: Must be a 1D cell array of the same length as Additional_GAdjust
% Each cell Is a 3 (row) x 2 (column) matrix with:
% First row as the range of the Gaussian depth in kJ/mol (generally should be negative)
% Second row as the range of the Gaussian center in nm (must be strictly positive)
% Third row is the Gaussian standard deviation in nm (must be strictly positive and non-zero)
% e.g. Settings.Additional_GAdjust_Ranges = {[-10 0; 0.3 0.6; 0.01 0.1]};

% Models(idx).Fix_Charge = false;
% Models(idx).Additivity = false;
% Models(idx).Additional_MM_Disp = false;
% Models(idx).Additional_GAdjust = {};
% Models(idx).Additional_GAdjust_Ranges = {};

%% Incorporating dispersion "rational" damping
% Models(idx).C6Damp.MM = 1; % Place rational damping on MM dispersion interaction
% Models(idx).C6Damp.MX = 1; % Place rational damping on MX dispersion interaction
% Models(idx).C6Damp.XX = 1; % Place rational damping on XX dispersion interaction

%% Global Calculation settings
clear;
Calculation.N_Calc = 6; % Number of chained calculations
Calculation.hours = 3; % Number of hours for each job in the job chain
Calculation.slurm_memory = '3800M'; % memory per core
Calculation.nodes = 1;
Calculation.cpus_per_node = 12;
Calculation.pbs_memory = [num2str(floor((186/32)*Calculation.cpus_per_node)) 'gb']; % Total memory request: floor(186gb / number of cores)
% Set up models below
skip_models = [];
check_complete = true; % Checks if job is already completed, skips completed jobs
check_running = true; % Checks if a job is already running, skips running jobs
continue_completed = false; % If a job is already complete, but you wish to continue, this will rename the previous *fullopt.mat file and restart. Must be used with check_complete = false

% Initial Model index
idx=0;

%% JC Models EV EW, EX
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Replicates = 1:5;
for sidx = 1:length(Salts)
    Salt = Salts{sidx};
    
    for ridx = 1:length(Replicates)
        Rep = num2str(Replicates(ridx));
        
        %% Model JC: EV
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EV' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EW
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EW' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
        %% Model JC: EW
        idx = idx+1;
        Models(idx) = Initialize_LiX_BO_Settings;
        Models(idx).Salt = Salt;
        Models(idx).Theory = 'JC';
        Models(idx).Trial_ID = ['EW' Rep];
        Models(idx).final_opt_type = 'fminsearchbnd';
        if Replicates(ridx) > 5
            Models(idx).Loss_Convergence = 1e-8;
            Models(idx).Param_Convergence = 1e-5;
        else
            Models(idx).Loss_Convergence = 1e-6;
            Models(idx).Param_Convergence = 1e-3;
        end
        
        % Loss
        Models(idx).Loss_Options.Rocksalt.LE = 1;
        Models(idx).Loss_Options.Rocksalt.a = 1;
        Models(idx).Loss_Options.Wurtzite.RLE = 1;
        Models(idx).Loss_Options.FiveFive.RLE = 1;
        Models(idx).Loss_Options.CsCl.RLE = 1/3;
        Models(idx).Q_value = 0.97; % Value for the charge scale. Only meaningful when Fix_Charge = true

        
        Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
        Models(idx).Fix_Charge = true;
        Models(idx).Additivity = true;
        Models(idx).SDMM_Range = [0 50];
        
    end
end

%% TF & BH Models GA, GB, GC, GD
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:5;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: GA
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GA' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GB
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GB' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GC
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GC' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: GD
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GD' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% TF & BH Models GE, GF, GG, GH, GI
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theories = {'BH' 'TF'};
Replicates = 1:5;

for tidx = 1:length(Theories)
    Theory = Theories{tidx};
    
    for sidx = 1:length(Salts)
        Salt = Salts{sidx};

        for ridx = 1:length(Replicates)
            Rep = num2str(Replicates(ridx));

            %% Model TF & BH: GE
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GE' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GF
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GF' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.Wurtzite.a = 2/3;
            Models(idx).Loss_Options.Wurtzite.c = 1/3;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GG
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GG' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = true;
            
            %% Model TF & BH: GH
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GH' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = true;
            Models(idx).Additivity = false;
            
            %% Model TF & BH: GI
            idx = idx+1;
            Models(idx) = Initialize_LiX_BO_Settings;
            Models(idx).Salt = Salt;
            Models(idx).Theory = Theory;
            Models(idx).Trial_ID = ['GI' Rep];
            Models(idx).final_opt_type = 'fminsearchbnd';
            if Replicates(ridx) > 5
                Models(idx).Loss_Convergence = 1e-8;
                Models(idx).Param_Convergence = 1e-5;
            else
                Models(idx).Loss_Convergence = 1e-6;
                Models(idx).Param_Convergence = 1e-3;
            end

            % Loss
            Models(idx).Loss_Options.Rocksalt.LE = 1;
            Models(idx).Loss_Options.Rocksalt.a = 1;
            Models(idx).Loss_Options.Wurtzite.RLE = 1;
            Models(idx).Loss_Options.FiveFive.RLE = 1;
            Models(idx).Loss_Options.Sphalerite.RLE = 1;
            Models(idx).Loss_Options.BetaBeO.RLE = 1;
            Models(idx).Loss_Options.AntiNiAs.RLE = 1;
            Models(idx).Loss_Options.NiAs.RLE = 1;
            Models(idx).Loss_Options.CsCl.RLE = 1;

            Models(idx).Structures = Auto_Structure_Selection(Models(idx).Loss_Options);
            Models(idx).SigmaEpsilon = true;
            Models(idx).Fix_Charge = false;
            Models(idx).Additivity = true;
            
        end
    end
end

%% Check for already running jobs
if check_running && isunix
    er = 1;
    jdx = 1;
    while er ~= 0
        [er,out] = system('squeue -u $USER -o %100j | tail -n +2');
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
    
    % Remove the job number tags
    c_jobs_unique = unique(regexprep(c_jobs,'-([0-9]){3}',''));
else
    c_jobs_unique = {''};
end

%% Prepare batch scripts
[home,project,computer] = find_home;

switch computer
    case 'cedar'
        pbs = false;
        account = 'rrg-patey-ad';
    case 'graham'
        pbs = false;
        account = 'def-patey';
    case 'sockeye'
        pbs = true;
    case 'patey'
        pbs = false;
    case 'unbearabull'
        pbs = false;
end

% Find the scheduler and get the batch script template file
if pbs
    fn = fullfile(home,'Machine_Learning','templates','PBS_template.subm');
    subm_txt = fileread(fn);
    subm_txt = strrep(subm_txt,'##MEMORY##',Calculation.pbs_memory);
else % slurm
    fn = fullfile(home,'Machine_Learning','templates','SLURM_template.subm');
    subm_txt = fileread(fn);
    subm_txt = strrep(subm_txt,'##MEMORY##',Calculation.slurm_memory);
    subm_txt = strrep(subm_txt,'##ACCOUNT##',account);
end

% Add in some global calculation parameters
subm_txt = strrep(subm_txt,'##HOURS##',num2str(Calculation.hours));
subm_txt = strrep(subm_txt,'##NODES##',num2str(Calculation.nodes));
subm_txt = strrep(subm_txt,'##CPUS##',num2str(Calculation.cpus_per_node));

% Loop through all models and build their submission file
for idx = 1:length(Models)
    if any(idx == skip_models)
        continue
    end
    Model = Models(idx);
    Model_Name = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
    Model_Name_abrv = [Model.Theory '_Model_' Model.Trial_ID];
    subm_txt_fin = strrep(subm_txt,'##MODEL_NAME##',Model_Name);
    
    % Generate and move to the submission directory
    submit_dir = fullfile(project,'Model_Building',Model.Salt,Model_Name_abrv);
    if ~exist(submit_dir, 'dir')
       mkdir(submit_dir)
    end
    cd(submit_dir);
    
    % Check if job is already complete
    if check_complete
        obs = dir(['*Model_' Model.Trial_ID '*fullopt.mat']);
        if ~isempty(obs)
            disp([Model_Name ': Job already completed. Skipping Job Submission.'])
            continue
        end
    end
    
    % Check if Model is already running
    if check_running && ismember(Model_Name,c_jobs_unique) 
        disp([Model_Name ': Job already Running. Skipping Job Submission.'])
        continue
    end
    
    if continue_completed
        obs = dir(['*Model_' Model.Trial_ID '*fullopt.mat']);
        if ~isempty(obs)
            disp([Model_Name ': Job already completed. Continuing completed Job.'])
            
            src = fullfile(obs.folder,obs.name);
            dest = fullfile(obs.folder,strrep(obs.name,'fullopt','oldopt'));
            movefile(src,dest);
        end
    end
    
    % Save input file
    inp_file = fullfile(submit_dir,[Model_Name '.inp']);
    save(inp_file,'Model','-mat')
    
    % Save job submission script
    subm_file = fullfile(submit_dir,[Model_Name '.subm']);
    fid = fopen(subm_file,'wt');
    fwrite(fid,regexprep(subm_txt_fin,{'\r'}',{''}));
    fclose(fid);
    
    % Submit job
    MultiSubmitBayesOpt(Calculation.N_Calc,submit_dir,[Model_Name '.subm'])
    delete(subm_file); % This is used as a template and can be removed
end

% Return
cd(project);
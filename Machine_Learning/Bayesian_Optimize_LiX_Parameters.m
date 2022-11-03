function Bayesian_Optimize_LiX_Parameters(Input_Settings)
    
    %% Load the model parameters
    if isstruct(Input_Settings)
        Settings = Input_Settings;
        batch_subm = false;
    elseif isfile(Input_Settings)
        Settings_dat = load(Input_Settings,'-mat');
        Settings = Settings_dat.Settings;
        batch_subm = true;
    else
        error('Input model does not exist, or is not a compatible data structure.')
    end
    
    % Disable an annoying warning
    warning('off','stats:classreg:learning:impl:GPImpl:GPImpl:SigmaMustBeGreaterThanSigmaLowerBound');
    
    if ~isfield(Settings,'OuterDir')
        Settings.OuterDir = pwd;
    end
    Intermediate_BO_file = fullfile(Settings.OuterDir,'intermediate_bayesian_opt.mat');
    Intermediate_BO_backup = fullfile(Settings.OuterDir,'intermediate_bayesian_opt.mat.PREV');
    Intermediate_Fullopt_file = fullfile(Settings.OuterDir,'intermediate_secondary_opt.mat');
    Intermediate_Fullopt_backup = fullfile(Settings.OuterDir,'intermediate_secondary_opt.mat.PREV');
    if isfield(Settings,'Diary_Loc') && ~isempty(Settings.Diary_Loc)
        diary(Settings.Diary_Loc)
    end
    
    % Initialize some global settings for later
    [Settings.Metal,Settings.Halide] = Separate_Metal_Halide(Settings.Salt);
    Settings.Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]);
    [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts] = MD_Batch_Template(Settings.JobSettings);
    
    Deterministic = ~any([Settings.Loss_Options.Fusion_Enthalpy Settings.Loss_Options.MP_Volume_Change Settings.Loss_Options.Liquid_MP_Volume ...
        Settings.Loss_Options.Solid_MP_Volume Settings.Loss_Options.Liquid_DM_MP Settings.Loss_Options.MP ] > sqrt(eps));
    
    if Settings.UseCoupledConstraint
        NumCC = 1;
        CCDet = true;
    else
        NumCC = 0;
        CCDet = [];
    end
    
    %% Display input summary on first iteration
    if ~isfile(Intermediate_BO_file) && ~isfile(Intermediate_Fullopt_file)
        disp(newline)
        disp(repmat('*',1,30))
        disp('Model Input Parameters')
        disp(repmat('*',1,30))
        disp(Settings)
        disp(repmat('*',1,30))
        disp('Structures:')
        disp(Settings.Structures)
        disp('Loss Function')
        disp(repmat('*',1,30))
        disp(['Regularization: ' Settings.Loss_Options.regularization])
        if Settings.Additional_Function.MM.N >= 0 && ~isempty(Settings.Additional_Function.MM.Range)
            disp(repmat('*',1,30))
            disp('Additional MM Function Included.')
            disp(['Additional MM Function Exponent: ' num2str(Settings.Additional_Function.MM.N)])
            disp(['Additional MM Function Scale Range: ' num2str(Settings.Additional_Function.MM.Range(1)) ...
                ' - ' num2str(Settings.Additional_Function.MM.Range(2))])
            disp(repmat('*',1,30))
        end
        if Settings.Additional_Function.XX.N >= 0 && ~isempty(Settings.Additional_Function.XX.Range)
            disp(repmat('*',1,30))
            disp('Additional XX Function Included.')
            disp(['Additional XX Function Exponent: ' num2str(Settings.Additional_Function.XX.N)])
            disp(['Additional XX Function Scale Range: ' num2str(Settings.Additional_Function.XX.Range(1)) ...
                ' - ' num2str(Settings.Additional_Function.XX.Range(2))])
            disp(repmat('*',1,30))
        end
        if Settings.Additional_Function.MX.N >= 0 && ~isempty(Settings.Additional_Function.MX.Range)
            disp(repmat('*',1,30))
            disp('Additional MX Function Included.')
            disp(['Additional MX Function Exponent: ' num2str(Settings.Additional_Function.MX.N)])
            disp(['Additional MX Function Scale Range: ' num2str(Settings.Additional_Function.MX.Range(1)) ...
                ' - ' num2str(Settings.Additional_Function.MX.Range(2))])
            disp(repmat('*',1,30))
        end
        if Settings.Loss_Options.Experimental_LE
            disp('Reference Rocksalt Energy: Born-Haber Experimental')
        else
            disp('Reference Rocksalt Energy: DFT')
        end
        if Settings.Loss_Options.Experimental_LP
            disp('Reference Rocksalt/Wurtzite Parameters: Experiment')
        else
            disp('Reference Rocksalt/Wurtzite Parameters: DFT')
        end
        tol = sqrt(eps);
        if Settings.Loss_Options.Fusion_Enthalpy > tol
            disp(['Enthalpy of Fusion at Experimental MP - Weight: ' num2str(Settings.Loss_Options.Fusion_Enthalpy)])
        end
        if Settings.Loss_Options.MP_Volume_Change > tol
            disp(['Volume Change of Fusion at Experimental MP - Weight: ' num2str(Settings.Loss_Options.MP_Volume_Change)])
        end
        if Settings.Loss_Options.Liquid_MP_Volume > tol
            disp(['Liquid Volume at Experimental MP - Weight: ' num2str(Settings.Loss_Options.Liquid_MP_Volume)])
        end
        if Settings.Loss_Options.Solid_MP_Volume > tol
            disp(['Solid Volume at Experimental MP - Weight: ' num2str(Settings.Loss_Options.Solid_MP_Volume)])
        end
        if Settings.Loss_Options.Liquid_DM_MP > tol
            disp([Settings.Metal ' Diffusion Constant at Experimental MP - Weight: ' num2str(Settings.Loss_Options.Liquid_DM_MP)])
        end
        
        if Settings.Loss_Options.MP > tol
            disp(['Melting Point - Weight: ' num2str(Settings.Loss_Options.MP)])
        end
        
        for idx = 1:length(Settings.Structures)
            Structure = Settings.Structures{idx};
            if Settings.Loss_Options.(Structure).LE > tol
                disp([Structure ' Absolute Lattice Energy - Weight: ' num2str(Settings.Loss_Options.(Structure).LE)])
            end
            if Settings.Loss_Options.(Structure).RLE > tol
                disp([Structure ' Relative Lattice Energy - Weight: ' num2str(Settings.Loss_Options.(Structure).RLE)])
            end
            if Settings.Loss_Options.(Structure).a > tol
                disp([Structure ' Lattice Parameter a - Weight: ' num2str(Settings.Loss_Options.(Structure).a)])
            end
            if Settings.Loss_Options.(Structure).b > tol
                disp([Structure ' Lattice Parameter b - Weight: ' num2str(Settings.Loss_Options.(Structure).b)])
            end
            if Settings.Loss_Options.(Structure).c > tol
                disp([Structure ' Lattice Parameter c - Weight: ' num2str(Settings.Loss_Options.(Structure).c)])
            end
            if Settings.Loss_Options.(Structure).V > tol
                disp([Structure ' Crystal Volume - Weight: ' num2str(Settings.Loss_Options.(Structure).V)])
            end
            if Settings.Loss_Options.(Structure).Gap.Weight > 0
                
                Ref_Structure = Settings.Loss_Options.(Structure).Gap.Ref;

                disp([Ref_Structure ' - ' Structure ' Energy Gap - Value: ' num2str(Settings.Loss_Options.(Structure).Gap.Value)])
                disp([Ref_Structure ' - ' Structure ' Energy Gap - Weight: ' num2str(Settings.Loss_Options.(Structure).Gap.Weight)])
                
                compare_fun_info = functions(Settings.Loss_Options.(Structure).Gap.Type);
                
                switch compare_fun_info.function % lt | gt | eq | ge | le | ne
                    case 'lt'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must be less than Target Gap.']);
                    case 'gt'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must be greater than Target Gap.']);
                    case 'eq'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must be equal to Target Gap.']);
                    case 'ge'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must be greater than or equal to Target Gap.']);
                    case 'le'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must be less than or equal to Target Gap.']);
                    case 'ne'
                        disp([Ref_Structure ' - ' Structure ' Energy Gap Enforcement: Actual Gap must not be equal to Target Gap.']);
                end
                disp('Otherwise loss function increases.');
            end
        end
        disp(repmat('*',1,30))
        disp(repmat('*',1,30))
    end
    
    %% Filenames
    FileBase = [Settings.Salt '_' Settings.Theory '_Model_' Settings.Trial_ID];
    Results_filename = [FileBase '_bayesopt.mat'];
    Final_point_filename = [FileBase '_final_point.mat'];
    Full_opt_filename = [FileBase '_fullopt.mat'];
    
    % Check for old results files
    obs = dir(['*Model_' Settings.Trial_ID '*bayesopt.mat']);
    if ~isempty(obs)
        Results_filename = obs.name;
    end
    obs = dir(['*Model_' Settings.Trial_ID '*final_point.mat']);
    if ~isempty(obs)
        Final_point_filename = obs.name;
    end
    obs = dir(['*Model_' Settings.Trial_ID '*fullopt.mat']);
    if ~isempty(obs)
        Full_opt_filename = obs.name;
    end
    %% Set the parameter ranges
    params = bayesopt_params(Settings);
    seed_points = length(params)*Settings.Initial_N_Multiplier;

    if Settings.Parallel_Bayesopt
        Settings.Parallel_Struct_Min = false;
        Settings.Parallel_LiX_Minimizer = false;
        Settings.MinMDP.Parallel_Min = false;
    elseif Settings.Parallel_LiX_Minimizer
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_Struct_Min = false;
        Settings.MinMDP.Parallel_Min = false;
    elseif Settings.Parallel_Struct_Min
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_LiX_Minimizer = false;
        Settings.MinMDP.Parallel_Min = true;
    end

    fun = @(x)LiX_Minimizer(Settings,x);
    constraint_fun = @(x)LiX_Constraint_Fcn(Settings,x);
    
    % Set up parallel features
    if Settings.Parallel_Bayesopt
        PrefCores = feature('numcores');
        if ~isempty(gcp('nocreate')) % If there is currently a parallel pool
            Cur_Pool = gcp;
            Cur_Workers = Cur_Pool.NumWorkers;
            
            % Start the parallel pool
            if Cur_Workers ~= PrefCores
                delete(Cur_Pool);
                % Create a "local" cluster object
                local_cluster = parcluster('local');
                
                % Modify the JobStorageLocation to a temporary directory
                [~,~,computer,~] = find_home;
                switch computer
                    case {'cedar' 'narval' 'graham'}
                        tmp = fullfile(getenv('SLURM_TMPDIR'),'local_cluster_jobs');
                    case 'sockeye'
                        tmp = fullfile(getenv('TMPDIR'),'local_cluster_jobs');
                    otherwise
                        tmp = fullfile(tempname,'local_cluster_jobs');
                end
                if ~isfolder(tmp)
                    mkdir(tmp);
                end
                local_cluster.JobStorageLocation = tmp;
                
                parpool(local_cluster,PrefCores);
            end
        else
            % Create a "local" cluster object
            local_cluster = parcluster('local');
            local_cluster.NumWorkers = PrefCores;
            
            % Modify the JobStorageLocation to a temporary directory
            [~,~,computer] = find_home;
            switch computer
                case {'cedar' 'graham' 'narval'}
                    tmp = fullfile(getenv('SLURM_TMPDIR'),'local_cluster_jobs');
                case 'sockeye'
                    tmp = fullfile(getenv('TMPDIR'),'local_cluster_jobs');
                otherwise
                    tmp = fullfile(tempname,'local_cluster_jobs');
            end
            if ~isfolder(tmp)
                mkdir(tmp)
            end
            local_cluster.JobStorageLocation = tmp;
            
            parpool(local_cluster,PrefCores); 
        end
    end
    
    %% Global optimization with bayesian optimization
    % Deal with restart requests
    if Settings.Continue_Calculation && (isfile(Intermediate_Fullopt_file) || isfile(Intermediate_Fullopt_backup))
        run_bayesopt = false;
        continue_bayesopt = false;
        continue_fullopt = true;
    elseif Settings.Continue_Calculation && isfile(Results_filename)
        run_bayesopt = false;
        continue_bayesopt = false;
        continue_fullopt = false;
    elseif Settings.Continue_Calculation && (isfile(Intermediate_BO_file) || isfile(Intermediate_BO_backup))
        run_bayesopt = true;
        continue_bayesopt = true;
        continue_fullopt = false;
    else
        if Settings.Continue_Calculation
            disp('No checkpoint files found. Starting new calculation.');
        end
        run_bayesopt = true;
        continue_bayesopt = false;
        continue_fullopt = false;
    end

    if ~run_bayesopt
        % Skip primary bayesian optimization step
        
        % Find and Load the previous results
        obs = dir(['*Model_' Settings.Trial_ID '*_bayesopt.mat']);         
        res = load(obs.name); % Loads the results variable
        results = res.results;
        
    elseif Settings.Continue_Calculation && continue_bayesopt
        % If you want to resume:
        try
            dat = load(Intermediate_BO_file,'-mat');
        catch
            % Catch corrupt files
            disp('Unable to load intermediate bayesian optimization data.')
            disp('Attempting to load backup step.')
            
            try
                dat = load(Intermediate_BO_backup,'-mat');
            catch
                disp('Unable to load backup step of Bayesian Optimization, restarting calculation.')
                if isfile(Intermediate_BO_file)
                    delete(Intermediate_BO_file)
                end
                if isfile(Intermediate_BO_backup)
                    delete(Intermediate_BO_backup)
                end
                Bayesian_Optimize_LiX_Parameters(Input_Settings)
                return
            end
        end
        
        % Catch corrupt files
        if ~isfield(dat,'BayesoptResults') || ~isprop(dat.BayesoptResults,'TotalElapsedTime')
            % Catch corrupt files
            disp('Unable to load intermediate bayesian optimization data.')
            disp('Attempting to load backup step.')
            try
                dat = load(Intermediate_BO_backup,'-mat');
            catch
                disp('Unable to load backup step of Bayesian Optimization, restarting calculation.')
                if isfile(Intermediate_BO_file)
                    delete(Intermediate_BO_file)
                end
                if isfile(Intermediate_BO_backup)
                    delete(Intermediate_BO_backup)
                end
                Bayesian_Optimize_LiX_Parameters(Input_Settings)
                return
            end
        end
        
        switch Settings.initial_opt_type
            case 'bayesopt'
                if Settings.ShowPlots
                    plotopt = 'all';
                else
                    plotopt = [];
                end
                
                BayesoptResults = dat.BayesoptResults;
                remaining_evals = Settings.Max_Bayesian_Iterations ...
                    + size(BayesoptResults.Options.InitialX,1) ...
                    - BayesoptResults.NumObjectiveEvaluations;
                
                if isfield(BayesoptResults.Options,'KernelFunction')
                    results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Settings.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',1,...
                        'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',remaining_evals,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',2000,...
                        'KernelFunction',Settings.KernelFunction,'XConstraintFcn',constraint_fun,...
                        'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet);
                else
                    results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Settings.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',1,...
                        'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',remaining_evals,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',2000,'XConstraintFcn',constraint_fun,...
                        'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet);
                end
                % Save results
                save(Results_filename,'results')                
        end
    else
        rng('shuffle'); % Reset the random seed to current time
        switch Settings.initial_opt_type
            case 'bayesopt'
                if Settings.ShowPlots
                    plotopt = 'all';
                else
                    plotopt = [];
                end
                
                if ~isempty(Settings.Initialize_From_Model)
                    Prev_Data = LoadPrevModelData(Settings,params);
                    Iterations = size(Prev_Data.InitialObjective,1) + Settings.Max_Bayesian_Iterations;
                    
                    results = bayesopt_priv(fun,params,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Settings.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',1,...
                        'UseParallel',Settings.Parallel_Bayesopt,'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',Iterations,'NumSeedPoints',seed_points,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',2000,...
                        'KernelFunction',Settings.KernelFunction,'XConstraintFcn',constraint_fun,...
                        'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet,...
                        'InitialX',Prev_Data.InitialX,...
                        'InitialObjective',Prev_Data.InitialObjective,...
                        'InitialConstraintViolations',Prev_Data.InitialConstraintViolations,...
                        'InitialErrorValues',Prev_Data.InitialErrorValues,...
                        'InitialUserData',Prev_Data.InitialUserData,...
                        'InitialObjectiveEvaluationTimes',Prev_Data.InitialObjectiveEvaluationTimes,...
                        'InitialIterationTimes',Prev_Data.InitialIterationTimes);
                else
                    results = bayesopt_priv(fun,params,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Settings.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',1,...
                        'UseParallel',Settings.Parallel_Bayesopt,'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',Settings.Max_Bayesian_Iterations,'NumSeedPoints',seed_points,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',2000,...
                        'KernelFunction',Settings.KernelFunction,'XConstraintFcn',constraint_fun,...
                        'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet);
                end
            case 'surrogateopt'
                if Settings.ShowPlots
                    plotopt = 'surrogateoptplot';
                else
                    plotopt = [];
                end
                % surrogateopt options
                options = optimoptions('surrogateopt','Display','iter',...
                    'MaxFunctionEvaluations',Settings.Max_Bayesian_Iterations,...
                    'MinSurrogatePoints',Settings.MinSurrogatePoints,...
                    'MinSampleDistance',Settings.MinSampleDistance,...
                    'PlotFcn',plotopt,...
                    'UseParallel',Settings.Parallel_Bayesopt,...
                    'CheckpointFile','intermediate_surrogate_opt.mat');

                ub = nan(1,length(params));
                lb = nan(1,length(params));
                for idx = 1:length(params)
                    lb(idx) = params(idx).Range(1);
                    ub(idx) = params(idx).Range(2);
                end
                [x0,fval,~,output,trials] = surrogateopt(fun,lb,ub,options);
                results.ObjectiveFcn = fun;
                results.XAtMinObjective = array2table(x0,'VariableNames',{params.Name});
                results.XTrace = array2table(trials.X,'VariableNames',{params.Name});
                results.NumObjectiveEvaluations = output.funccount;
                results.MinEstimatedObjective = fval;
                results.MinObjective = fval;
        end
        % Save results
        save(Results_filename,'results')
    end
    
    %% Secondary optimization with different aqcuisition function focused on locally optimizing result
    if strcmpi(Settings.initial_opt_type,'bayesopt') && strcmpi(Settings.second_opt_type,'bayesopt') && ~isfile('secondary_completed')
        if Settings.ShowPlots
            plotopt = 'all';
        else
            plotopt = [];
        end
        
        % Make file to signal that secondary optimization has already began
        if ~isfile('secondary_started')
            disp('Beginning Intermediate Optimization Step.')
            fclose(fopen('secondary_started', 'w'));
            BayesoptResults = results;
        else
            try
                dat = load(Intermediate_BO_file);
                BayesoptResults = dat.BayesoptResults;
            catch
                % Catch corrupt files
                disp('Unable to load intermediate bayesian optimization data.')
                disp('Attempting to load backup step.')

                try
                    dat = load(Intermediate_BO_backup,'-mat');
                    BayesoptResults = dat.BayesoptResults;
                catch
                    disp('Unable to load backup step of Bayesian Optimization.')
                    disp('Restarting intermediate optimization from finished initial optimization.')
                    if isfile(Intermediate_BO_file)
                        delete(Intermediate_BO_file)
                    end
                    if isfile(Intermediate_BO_backup)
                        delete(Intermediate_BO_backup)
                    end
                    BayesoptResults = results;
                end
            end
        end
        
        % Calculate remaining evaluations for secondary optimization
        remaining_evals = Settings.Max_Secondary_Iterations ...
            + Settings.Max_Bayesian_Iterations ...
            + size(BayesoptResults.Options.InitialX,1) ...
            - BayesoptResults.NumObjectiveEvaluations;
        
        if remaining_evals > 0
            if isfield(BayesoptResults.Options,'KernelFunction')
                results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                    'PlotFcn',plotopt,'ExplorationRatio',Settings.ExplorationRatio_Secondary,...
                    'AcquisitionFunctionName',Settings.Secondary_Acquisition_Function,'Verbose',1,...
                    'ParallelMethod','clipped-model-prediction',...
                    'MaxObjectiveEvaluations',remaining_evals,...
                    'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                    'GPActiveSetSize',2000,...
                    'KernelFunction',Settings.SecondaryKernelFunction,'XConstraintFcn',constraint_fun,...
                    'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet);
            else
                results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                    'PlotFcn',plotopt,'ExplorationRatio',Settings.ExplorationRatio_Secondary,...
                    'AcquisitionFunctionName',Settings.Secondary_Acquisition_Function,'Verbose',1,...
                    'ParallelMethod','clipped-model-prediction','MaxObjectiveEvaluations',remaining_evals,...
                    'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                    'GPActiveSetSize',2000,...
                    'XConstraintFcn',constraint_fun,...
                    'NumCoupledConstraints',NumCC,'AreCoupledConstraintsDeterministic',CCDet);
            end
        else
            results = BayesoptResults;
        end
        
        % overwrite results
        save(Results_filename,'results')
        fclose(fopen('secondary_completed', 'w'));
        delete('secondary_started');
    end
    
    %% Perform final local optimization %%
    if ~isfile(Final_point_filename)
        switch Settings.final_opt_type
            case 'fminsearch'
                % First, reset the parallelization scheme for fminsearch as it has no parallel option
                % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
                if Settings.Parallel_Bayesopt || Settings.Parallel_LiX_Minimizer
                    Settings.Parallel_Bayesopt = false;
                    Settings.Parallel_Struct_Min = false;
                    Settings.Parallel_LiX_Minimizer = true;
                    Settings.MinMDP.Parallel_Min = false;
                elseif Settings.Parallel_Struct_Min
                    Settings.Parallel_Bayesopt = false;
                    Settings.Parallel_Struct_Min = true;
                    Settings.Parallel_LiX_Minimizer = false;
                    Settings.MinMDP.Parallel_Min = true;
                end
                fun = @(x)LiX_Minimizer(Settings,x);

                % Options for the Nelder–Mead search
                optionsNM = optimset('Display','iter','MaxIter',Settings.Max_Local_Iterations,...
                    'TolFun',Settings.Loss_Convergence,'TolX',Settings.Param_Convergence,...
                    'OutputFcn',@outputFcn_secondary_opt);

                if continue_fullopt
                    try
                        dat = load(Intermediate_Fullopt_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;

                        remaining_evals = Settings.Max_Local_Iterations - length(intermediate_data)+2;
                        if remaining_evals <= 0
                            remaining_evals = 1;
                        end

                        optionsNM.MaxIter = remaining_evals;
                    catch
                        disp('Unable to load local optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Fullopt_file)
                            delete(Intermediate_Fullopt_file);
                        end
                        try

                            dat = load(Intermediate_Fullopt_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;

                            remaining_evals = Settings.Max_Local_Iterations - length(intermediate_data)+2;
                            if remaining_evals <= 0
                                remaining_evals = 1;
                            end

                            optionsNM.MaxIter = remaining_evals;
                        catch
                            if isfile(Intermediate_Fullopt_backup)
                                delete(Intermediate_Fullopt_backup);
                            end
                            disp('Unable to load backup secondary optimization checkpoint file, restarting.')
                            [~,midx] = min(results.ObjectiveTrace);
                            x0 = results.XTrace{midx,:};
                        end
                    end
                else
                    [~,midx] = min(results.ObjectiveTrace);
                    x0 = results.XTrace{midx,:};
                end

                [full_opt_point,~,~,full_opt_results] = fminsearch(fun,x0,optionsNM);
            case 'fminsearchbnd'

                % Constraints
                Ranges = [params.Range];
                lb = Ranges(1:2:end);
                ub = Ranges(2:2:end);
                
                % First, reset the parallelization scheme for fminsearchbnd as it has no parallel option
                % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
                if Settings.Parallel_Bayesopt || Settings.Parallel_LiX_Minimizer
                    Settings.Parallel_Bayesopt = false;
                    Settings.Parallel_Struct_Min = false;
                    Settings.Parallel_LiX_Minimizer = true;
                    Settings.MinMDP.Parallel_Min = false;
                elseif Settings.Parallel_Struct_Min
                    Settings.Parallel_Bayesopt = false;
                    Settings.Parallel_Struct_Min = true;
                    Settings.Parallel_LiX_Minimizer = false;
                    Settings.MinMDP.Parallel_Min = true;
                end
                fun = @(x)LiX_Minimizer(Settings,x);

                % Options for the Nelder–Mead search
                optionsNM = optimset('Display','iter','MaxIter',Settings.Max_Local_Iterations,...
                    'MaxFunEvals',Settings.MaxFunEvals,'TolFun',Settings.Loss_Convergence,'TolX',Settings.Param_Convergence,...
                    'OutputFcn',@outputFcn_secondary_opt);

                if continue_fullopt
                    try
                        dat = load(Intermediate_Fullopt_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;

                        remaining_evals = Settings.Max_Local_Iterations - length(intermediate_data)+2;

                        if remaining_evals <= 0
                            remaining_evals = 1;
                        end

                        optionsNM.MaxIter = remaining_evals;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Fullopt_file)
                            delete(Intermediate_Fullopt_file);
                        end
                        try
                            dat = load(Intermediate_Fullopt_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;

                            remaining_evals = Settings.Max_Local_Iterations - length(intermediate_data)+2;

                            if remaining_evals <= 0
                                remaining_evals = 1;
                            end

                            optionsNM.MaxIter = remaining_evals;
                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Fullopt_backup)
                                delete(Intermediate_Fullopt_backup);
                            end
                            remaining_evals = Settings.Max_Local_Iterations;
                            [~,midx] = min(results.ObjectiveTrace);
                            x0 = results.XTrace{midx,:};
                        end
                    end
                else
                    remaining_evals = Settings.Max_Local_Iterations;
                    [~,midx] = min(results.ObjectiveTrace);
                    x0 = results.XTrace{midx,:};
                end

                [full_opt_point,~,exitflag,full_opt_results] = fminsearchbnd(fun,x0,lb,ub,optionsNM);

                % Check to see if the max local iterations is exceeded
                if exitflag == 0 && remaining_evals > 1 && Settings.switch_final_opt % calculation exceeded number of allowed iterations
                    % switch to nealder-mead simplex minimization
                    disp('Swicthing local optimization to patternsearch')
                    Settings.final_opt_type = 'patternsearch';

                    % Overwrite the input file if it exists
                    if batch_subm
                        save(Input_Settings,'Settings','-mat')
                    end

                    Bayesian_Optimize_LiX_Parameters(Settings)
                    return
                end
            case 'patternsearch'

                %% Options for Patternsearch

                % Constraints
                Ranges = [params.Range];
                lb = Ranges(1:2:end);
                ub = Ranges(2:2:end);

                % Load previous data
                if continue_fullopt
                    try
                        dat = load(Intermediate_Fullopt_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;    

                        % Establish initial mesh size
                        if isfield(intermediate_data(end).optimValues,'meshsize')
                            init_meshsize = intermediate_data(end).optimValues.meshsize;
                        else
                            init_meshsize = 0.1;
                        end

                        current_iterations = length(intermediate_data);
                        max_iter = Settings.Max_Local_Iterations - current_iterations;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Fullopt_file)
                            delete(Intermediate_Fullopt_file);
                        end
                        try
                            dat = load(Intermediate_Fullopt_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;    

                            % Establish initial mesh size
                            if isfield(intermediate_data(end).optimValues,'meshsize')
                                init_meshsize = intermediate_data(end).optimValues.meshsize;
                            else
                                init_meshsize = 0.1;
                            end

                            current_iterations = length(intermediate_data);
                            max_iter = Settings.Max_Local_Iterations - current_iterations;

                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Fullopt_backup)
                                delete(Intermediate_Fullopt_backup);
                            end
                            [~,midx] = min(results.ObjectiveTrace);
                            x0 = results.XTrace{midx,:};
                            
                            init_meshsize = 0.1;
                            max_iter = Settings.Max_Local_Iterations;
                        end
                    end
                else
                    [~,midx] = min(results.ObjectiveTrace);
                    x0 = results.XTrace{midx,:};
                    init_meshsize = 0.1;
                    max_iter = Settings.Max_Local_Iterations;
                end

        %         optionsNM = optimset('Display','iter','MaxIter',50,'MaxFunEvals',Model.MaxFunEvals,...
        %             'TolFun',Model.Loss_Convergence,'TolX',Model.Param_Convergence);

                if max_iter <= 0
                    max_iter = 1;
                end

                if Settings.ShowPlots
                    plotopt = 'psplotbestf';
                else
                    plotopt = [];
                end

                options = optimoptions(@patternsearch,'Display','iter','MaxIterations',max_iter,...
                    'UseParallel',Settings.Parallel_Bayesopt,'UseVectorized',false,'PlotFcn',plotopt,...
                    'InitialMeshSize',init_meshsize,'StepTolerance',Settings.Param_Convergence,'FunctionTolerance',Settings.Loss_Convergence,...
                    'PollOrderAlgorithm','Success','Cache','On','UseCompletePoll',false,...
                    'PollMethod','GPSPositiveBasis2N','MaxMeshSize',Inf,'MeshTolerance',1e-8,...
                    'OutputFcn',@outputFcn_patternsearch_opt);%,'MaxFunctionEvaluations',Model.MaxFunEvals);
        %             'SearchFcn',{@searchneldermead,1,optionsNM});

                [full_opt_point,~,~,full_opt_results] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);

                % Check to see if the max local iterations is exceeded
        %         if exitflag == 0  && max_iter > 1 && Model.switch_final_opt % calculation exceeded number of allowed iterations
        %             disp('Swicthing local optimization to fminsearchbnd')
        %             Model.final_opt_type = 'fminsearchbnd';
        %             
        %             % Overwrite the input file if it exists
        %             if batch_subm
        %                 save(Input_Model,'Settings','-mat')
        %             end
        %             
        %             Bayesian_Optimize_LiX_Parameters(Settings)
        %             return
        %         end
            case 'fmincon' % fmincon uses gradients!

                % Constraints
                Ranges = [params.Range];
                lb = Ranges(1:2:end);
                ub = Ranges(2:2:end);

                % Load previous data
                if continue_fullopt
                    try
                        dat = load(Intermediate_Fullopt_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;    

                        current_iterations = length(intermediate_data);
                        max_iter = Settings.Max_Local_Iterations - current_iterations;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Fullopt_file)
                            delete(Intermediate_Fullopt_file);
                        end
                        try
                            dat = load(Intermediate_Fullopt_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;    
                            
                            current_iterations = length(intermediate_data);
                            max_iter = Settings.Max_Local_Iterations - current_iterations;
                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Fullopt_backup)
                                delete(Intermediate_Fullopt_backup);
                            end
                            [~,midx] = min(results.ObjectiveTrace);
                            x0 = results.XTrace{midx,:};
                            max_iter = Settings.Max_Local_Iterations;
                        end
                    end
                else
                    [~,midx] = min(results.ObjectiveTrace);
                    x0 = results.XTrace{midx,:};
                    max_iter = Settings.Max_Local_Iterations;
                end

                if max_iter <= 0
                    max_iter = 1;
                end

                % For unconstrained problem with fmincon, OptimalityTolerance is the max of the gradient
                options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set',...
                    'DiffMaxChange',1e-1,'OptimalityTolerance',1e-1,'UseParallel',Settings.Parallel_Bayesopt,...
                    'MaxIterations',max_iter,'FiniteDifferenceStepSize',sqrt(eps),...
                    'StepTolerance',Settings.Param_Convergence,'FunctionTolerance',Settings.Loss_Convergence,...
                    'FiniteDifferenceType','forward','MaxFunctionEvaluations',Inf,...
                    'OutputFcn',@outputFcn_secondary_opt);

                [full_opt_point,~,~,full_opt_results,~,~,~] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
            case 'none' % if no final optimization is selected, I still want to output the best bayesian optimization result
                full_opt_results = struct;
                full_opt_results.iterations = 0;
                full_opt_results.funccount = 0;
                full_opt_results.algorithm = 'none';
                full_opt_results.message = 'Local optimization skipped. Results are for best global optimization point.';
                [~,midx] = min(results.ObjectiveTrace);
                full_opt_point = results.XTrace{midx,:};
            otherwise
                error('Unknown final optimization scheme. Choose one of "fminsearch", "patternsearch", "fminsearchbnd", or "none"')
        end
        
        if ~Deterministic && ~strcmpi(Settings.final_opt_type,'none')
            try
                dat = load(Intermediate_Fullopt_file).intermediate_data;
            catch
                disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                if isfile(Intermediate_Fullopt_file)
                    delete(Intermediate_Fullopt_file);
                end
                    try
                        dat = load(Intermediate_Fullopt_backup,'-mat');
                    catch
                        disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                        if isfile(Intermediate_Fullopt_backup)
                            delete(Intermediate_Fullopt_backup);
                        end
                        Bayesian_Optimize_LiX_Parameters(Input_Settings);
                        return
                    end
            end
            fvals = zeros(1,length(dat));
            for idx = 1:length(dat)
                fvals(idx) = dat(idx).optimValues.fval;
            end
            [~,minidx] = min(fvals);
            full_opt_point = dat(minidx).x;
        end
        
        % Save a structure containing some calculation properties
        Calculation_properties.Salt = Settings.Salt;
        Calculation_properties.Theory = Settings.Theory;
        Calculation_properties.Fix_Charge = Settings.Fix_Charge;
        Calculation_properties.Additivity = Settings.Additivity;
        Calculation_properties.Additional_MM_Disp = Settings.Additional_MM_Disp;
        Calculation_properties.Additional_GAdjust = Settings.Additional_GAdjust;
        Calculation_properties.SigmaEpsilon = Settings.SigmaEpsilon;
        save(Final_point_filename,'full_opt_results','full_opt_point',...
            'Calculation_properties');
    else
        try
            dat = load(Final_point_filename);
        catch
            disp('Unable to load final point data file, continuing final optimization.')
            if isfile(Final_point_filename)
                delete(Final_point_filename)
            end
            Bayesian_Optimize_LiX_Parameters(Input_Settings)
            return
        end
        Calculation_properties = dat.Calculation_properties;
        full_opt_results = dat.full_opt_results;
        full_opt_point = dat.full_opt_point;
    end
    
    %% Final test of parameters on all structures for output
    % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
    if Settings.Parallel_Bayesopt || Settings.Parallel_LiX_Minimizer
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_Struct_Min = false;
        Settings.Parallel_LiX_Minimizer = true;
        Settings.MinMDP.Parallel_Min = false;
    elseif Settings.Parallel_Struct_Min
        Settings.Parallel_Bayesopt = false;
        Settings.Parallel_LiX_Minimizer = false;
        Settings.MinMDP.Parallel_Min = true;
    end
    Settings.Delete_Equil = false; % save the final MP calculation directories
    Settings.Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
    Settings.Verbose = true;
    [loss,~,UserData] = LiX_Minimizer(Settings,full_opt_point,...
        'Extra_Properties',true,'Therm_Prop_Override',true);
    
    ParNames = {params.Name};
    Pars = table;
    for idx = 1:numel(ParNames)
        Pars.(ParNames{idx}) = full_opt_point(idx);
    end
    
    % Load targets
    DFT = Load_Best_DFT_Data;
    Structures = Settings.Structures;
    N = length(Structures);
    if Settings.Loss_Options.Experimental_LE || Settings.Loss_Options.Experimental_LP
        Exp = Load_Experimental_Data;
        
        if Settings.Loss_Options.Experimental_LE
            E_Correction = Exp.(Settings.Salt).Rocksalt.E - DFT.(Settings.Salt).Rocksalt.Energy;        
            for idx = 1:N
                try
                    DFT.(Settings.Salt).(Structures{idx}).Energy = DFT.(Settings.Salt).(Structures{idx}).Energy + E_Correction;
                catch
                    DFT.(Settings.Salt).(Structures{idx}).Energy = nan;
                    DFT.(Settings.Salt).(Structures{idx}).a = nan;
                    DFT.(Settings.Salt).(Structures{idx}).b = nan;
                    DFT.(Settings.Salt).(Structures{idx}).c = nan;
                end
            end
        end
        if Settings.Loss_Options.Experimental_LP
            DFT.(Settings.Salt).Rocksalt.a = Exp.(Settings.Salt).Rocksalt.a_zero;
            DFT.(Settings.Salt).Rocksalt.b = Exp.(Settings.Salt).Rocksalt.b_zero;
            DFT.(Settings.Salt).Rocksalt.c = Exp.(Settings.Salt).Rocksalt.c_zero;
            DFT.(Settings.Salt).Rocksalt.V = Exp.(Settings.Salt).Rocksalt.V_zero;
            if isfield(Exp.(Settings.Salt),'Wurtzite')
                DFT.(Settings.Salt).Wurtzite.a = Exp.(Settings.Salt).Wurtzite.a_zero;
                DFT.(Settings.Salt).Wurtzite.b = Exp.(Settings.Salt).Wurtzite.b_zero;
                DFT.(Settings.Salt).Wurtzite.c = Exp.(Settings.Salt).Wurtzite.c_zero;
                DFT.(Settings.Salt).Wurtzite.V = Exp.(Settings.Salt).Wurtzite.V_zero;
            end
        end
    end
    
    Minimization_Data = UserData.Minimization_Data;
    Finite_T_Data = UserData.Finite_T_Data;
    
    disp(repmat('*',1,80))
    disp(['Final Results - [Salt: ' Settings.Salt '] - [Potential Form: ' Settings.Theory '] - [Model name: ' Settings.Trial_ID ']'])
    disp(repmat('*',1,80))
    
    format long g
    for idx = 1:length(Settings.Structures)
        disp(repmat('*',1,60))
        disp(Settings.Structures{idx})
        disp('Target (Exp/DFT) Data:')
        disp(repmat('*',1,60))
        disp(DFT.(Settings.Salt).(Settings.Structures{idx}))
        disp(repmat('*',1,60))
        disp(Settings.Structures{idx})
        disp('Model Result:')
        disp(repmat('*',1,60))
        disp(Minimization_Data{idx})
    end
    disp(repmat('*',1,120))
    disp('Finite Temperature Data (Enthalpy in kJ/mol, Entropy in J/(mol K), Volume in A^3/Forumla Unit, Temperature in K)')
    disp(repmat('*',1,120))
    disp(Finite_T_Data)
    disp(repmat('*',1,120))
    
    disp(['Final Optimized Loss: ' num2str(loss,'%.10e')])
    disp('Final Optimized Parameters:')
    disp(Pars)
    
    % Save final results
    save(Full_opt_filename,'full_opt_results','loss','full_opt_point',...
        'Minimization_Data','Finite_T_Data','Pars','Calculation_properties');
    diary off
    
end
function Bayesian_Optimize_LiX_Parameters(Input_Model)
    
    %% Load the model parameters
    if isstruct(Input_Model)
        Model = Input_Model;
        batch_subm = false;
    elseif isfile(Input_Model)
        Model_dat = load(Input_Model,'-mat');
        Model = Model_dat.Model;
        batch_subm = true;
    else
        error('Input model does not exist, or is not a compatible data structure.')
    end
    
    % Disable an annoying warning
    warning('off','stats:classreg:learning:impl:GPImpl:GPImpl:SigmaMustBeGreaterThanSigmaLowerBound');
    
    Model.OuterDir = pwd;
    Intermediate_BO_file = fullfile(Model.OuterDir,'intermediate_bayesian_opt.mat');
    Intermediate_BO_backup = fullfile(Model.OuterDir,'intermediate_bayesian_opt.mat.PREV');
    Intermediate_Secondary_file = fullfile(Model.OuterDir,'intermediate_secondary_opt.mat');
    Intermediate_Seconary_backup = fullfile(Model.OuterDir,'intermediate_secondary_opt.mat.PREV');
    if isfield(Model,'Diary_Loc') && ~isempty(Model.Diary_Loc)
        diary(Model.Diary_Loc)
    end
    
    % Initialize some global settings for later
    [Model.Metal,Model.Halide] = Separate_Metal_Halide(Model.Salt);
    Model.Longest_Cutoff = max([Model.MDP.RList_Cutoff Model.MDP.RCoulomb_Cutoff Model.MDP.RVDW_Cutoff]);
    [~,Model.gmx,Model.gmx_loc,Model.mdrun_opts] = MD_Batch_Template(Model.JobSettings);
    
    Finite_T_Calc = any([Model.Loss_Options.Fusion_Enthalpy Model.Loss_Options.MP_Volume_Change Model.Loss_Options.Liquid_MP_Volume ...
        Model.Loss_Options.Solid_MP_Volume Model.Loss_Options.Liquid_DM_MP Model.Loss_Options.MP ] > sqrt(eps));
    
    if Finite_T_Calc
        Deterministic = false; % Thermal properties are not deterministic
    else
        Deterministic = true; % Lattice energies are deterministic
    end
    
    %% Display input summary on first iteration
    if ~isfile(Intermediate_BO_file) && ~isfile(Intermediate_Secondary_file)
        disp(newline)
        disp(repmat('*',1,30))
        disp('Model Input Parameters')
        disp(repmat('*',1,30))
        disp(Model)
        disp(repmat('*',1,30))
        disp('Structures:')
        disp(Model.Structures)
        disp('Loss Function')
        disp(repmat('*',1,30))
        disp(['Regularization: ' Model.Loss_Options.regularization])
        if Model.Additional_Function.MM.N >= 0 && ~isempty(Model.Additional_Function.MM.Range)
            disp(repmat('*',1,30))
            disp('Additional MM Function Included.')
            disp(['Additional MM Function Exponent: ' num2str(Model.Additional_Function.MM.N)])
            disp(['Additional MM Function Scale Range: ' num2str(Model.Additional_Function.MM.Range(1)) ...
                ' - ' num2str(Model.Additional_Function.MM.Range(2))])
            disp(repmat('*',1,30))
        end
        if Model.Additional_Function.XX.N >= 0 && ~isempty(Model.Additional_Function.XX.Range)
            disp(repmat('*',1,30))
            disp('Additional XX Function Included.')
            disp(['Additional XX Function Exponent: ' num2str(Model.Additional_Function.XX.N)])
            disp(['Additional XX Function Scale Range: ' num2str(Model.Additional_Function.XX.Range(1)) ...
                ' - ' num2str(Model.Additional_Function.XX.Range(2))])
            disp(repmat('*',1,30))
        end
        if Model.Additional_Function.MX.N >= 0 && ~isempty(Model.Additional_Function.MX.Range)
            disp(repmat('*',1,30))
            disp('Additional MX Function Included.')
            disp(['Additional MX Function Exponent: ' num2str(Model.Additional_Function.MX.N)])
            disp(['Additional MX Function Scale Range: ' num2str(Model.Additional_Function.MX.Range(1)) ...
                ' - ' num2str(Model.Additional_Function.MX.Range(2))])
            disp(repmat('*',1,30))
        end
        if Model.Loss_Options.Experimental_LE
            disp('Reference Rocksalt Energy: Born-Haber Experimental')
        else
            disp('Reference Rocksalt Energy: DFT')
        end
        if Model.Loss_Options.Experimental_LP
            disp('Reference Rocksalt/Wurtzite Parameters: Experiment')
        else
            disp('Reference Rocksalt/Wurtzite Parameters: DFT')
        end
        tol = sqrt(eps);
        if Model.Loss_Options.Fusion_Enthalpy > tol
            disp(['Enthalpy of Fusion at Experimental MP - Weight: ' num2str(Model.Loss_Options.Fusion_Enthalpy)])
        end
        if Model.Loss_Options.MP_Volume_Change > tol
            disp(['Volume Change of Fusion at Experimental MP - Weight: ' num2str(Model.Loss_Options.MP_Volume_Change)])
        end
        if Model.Loss_Options.Liquid_MP_Volume > tol
            disp(['Liquid Volume at Experimental MP - Weight: ' num2str(Model.Loss_Options.Liquid_MP_Volume)])
        end
        if Model.Loss_Options.Solid_MP_Volume > tol
            disp(['Solid Volume at Experimental MP - Weight: ' num2str(Model.Loss_Options.Solid_MP_Volume)])
        end
        if Model.Loss_Options.Liquid_DM_MP > tol
            disp([Model.Metal ' Diffusion Constant at Experimental MP - Weight: ' num2str(Model.Loss_Options.Liquid_DM_MP)])
        end
        
        if Model.Loss_Options.MP > tol
            disp(['Melting Point - Weight: ' num2str(Model.Loss_Options.MP)])
        end
        
        for idx = 1:length(Model.Structures)
            Structure = Model.Structures{idx};
            if Model.Loss_Options.(Structure).LE > tol
                disp([Structure ' Absolute Lattice Energy - Weight: ' num2str(Model.Loss_Options.(Structure).LE)])
            end
            if Model.Loss_Options.(Structure).RLE > tol
                disp([Structure ' Relative Lattice Energy - Weight: ' num2str(Model.Loss_Options.(Structure).RLE)])
            end
            if Model.Loss_Options.(Structure).a > tol
                disp([Structure ' Lattice Parameter a - Weight: ' num2str(Model.Loss_Options.(Structure).a)])
            end
            if Model.Loss_Options.(Structure).b > tol
                disp([Structure ' Lattice Parameter b - Weight: ' num2str(Model.Loss_Options.(Structure).b)])
            end
            if Model.Loss_Options.(Structure).c > tol
                disp([Structure ' Lattice Parameter c - Weight: ' num2str(Model.Loss_Options.(Structure).c)])
            end
            if Model.Loss_Options.(Structure).V > tol
                disp([Structure ' Crystal Volume - Weight: ' num2str(Model.Loss_Options.(Structure).V)])
            end
            if Model.Loss_Options.(Structure).Gap.Weight > 0
                
                Ref_Structure = Model.Loss_Options.(Structure).Gap.Ref;

                disp([Ref_Structure ' - ' Structure ' Energy Gap - Value: ' num2str(Model.Loss_Options.(Structure).Gap.Value)])
                disp([Ref_Structure ' - ' Structure ' Energy Gap - Weight: ' num2str(Model.Loss_Options.(Structure).Gap.Weight)])
                
                compare_fun_info = functions(Model.Loss_Options.(Structure).Gap.Type);
                
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
    FileBase = [Model.Salt '_' Model.Theory '_Model_' Model.Trial_ID];
    Results_filename = [FileBase '_bayesopt.mat'];
    Final_point_filename = [FileBase '_final_point.mat'];
    Full_opt_filename = [FileBase '_fullopt.mat'];
    
    % Check for old results files
    obs = dir(['*Model_' Model.Trial_ID '*bayesopt.mat']);
    if ~isempty(obs)
        Results_filename = obs.name;
    end
    obs = dir(['*Model_' Model.Trial_ID '*final_point.mat']);
    if ~isempty(obs)
        Final_point_filename = obs.name;
    end
    obs = dir(['*Model_' Model.Trial_ID '*fullopt.mat']);
    if ~isempty(obs)
        Full_opt_filename = obs.name;
    end
    %% Set the parameter ranges
    params = bayesopt_params(Model);
    seed_points = length(params)*10;

    if Model.Parallel_Bayesopt
        Model.Parallel_Struct_Min = false;
        Model.Parallel_LiX_Minimizer = false;
        Model.MinMDP.Parallel_Min = false;
    elseif Model.Parallel_LiX_Minimizer
        Model.Parallel_Bayesopt = false;
        Model.Parallel_Struct_Min = false;
        Model.MinMDP.Parallel_Min = false;
    elseif Model.Parallel_Struct_Min
        Model.Parallel_Bayesopt = false;
        Model.Parallel_LiX_Minimizer = false;
        Model.MinMDP.Parallel_Min = true;
    end

    fun = @(x)LiX_Minimizer(Model,x);

    % Set up parallel features
    if Model.Parallel_Bayesopt
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
    if Model.Restart_Calculation && (isfile(Intermediate_Secondary_file) || isfile(Intermediate_Seconary_backup))
        run_bayesopt = false;
        continue_bayesopt = false;
        continue_fullopt = true;
    elseif Model.Restart_Calculation && isfile(Results_filename)
        run_bayesopt = false;
        continue_bayesopt = false;
        continue_fullopt = false;
    elseif Model.Restart_Calculation && (isfile(Intermediate_BO_file) || isfile(Intermediate_BO_backup))
        run_bayesopt = true;
        continue_bayesopt = true;
        continue_fullopt = false;
    else
        if Model.Restart_Calculation
            disp('No checkpoint files found. Starting new calculation.');
        end
        run_bayesopt = true;
        continue_bayesopt = false;
        continue_fullopt = false;
    end

    if ~run_bayesopt
        % Skip primary bayesian optimization step
        
        switch Model.initial_opt_type
            case 'bayesopt'
                % Find and Load the previous results
                obs = dir(['*Model_' Model.Trial_ID '*_bayesopt.mat']);         
                res = load(obs.name); % Loads the results variable
                results = res.results;
        end

    elseif Model.Restart_Calculation && continue_bayesopt
        % If you want to resume:
        try
            dat = load(Intermediate_BO_file);
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
                Bayesian_Optimize_LiX_Parameters(Input_Model)
                return
            end
        end
        
        % Catch corrupt files
        if ~isfield(dat,'BayesoptResults') || isnan(dat.BayesoptResults.MinObjective)
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
                Bayesian_Optimize_LiX_Parameters(Input_Model)
                return
            end
        end
        
        switch Model.initial_opt_type
            case 'bayesopt'
                if Model.ShowPlots
                    plotopt = 'all';
                else
                    plotopt = [];
                end
                
                BayesoptResults = dat.BayesoptResults;
                remaining_evals = Model.Max_Bayesian_Iterations - BayesoptResults.NumObjectiveEvaluations;
                
                if isfield(BayesoptResults.Options,'KernelFunction')
                    results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Model.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Model.Acquisition_Function,'Verbose',1,...
                        'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',remaining_evals,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',Model.Max_Bayesian_Iterations,...
                        'KernelFunction',Model.KernelFunction);
                else
                    results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                        'ExplorationRatio',Model.ExplorationRatio,'PlotFcn',plotopt,...
                        'AcquisitionFunctionName',Model.Acquisition_Function,'Verbose',1,...
                        'ParallelMethod','clipped-model-prediction',...
                        'MaxObjectiveEvaluations',remaining_evals,...
                        'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                        'GPActiveSetSize',Model.Max_Bayesian_Iterations);
                end
                % Save results
                save(Results_filename,'results')                
        end
    else
        rng('shuffle'); % Reset the random seed to current time
        switch Model.initial_opt_type
            case 'bayesopt'
                if Model.ShowPlots
                    plotopt = 'all';
                else
                    plotopt = [];
                end
                results = bayesopt_priv(fun,params,'IsObjectiveDeterministic',Deterministic,...
                    'ExplorationRatio',Model.ExplorationRatio,'PlotFcn',plotopt,...
                    'AcquisitionFunctionName',Model.Acquisition_Function,'Verbose',1,...
                    'UseParallel',Model.Parallel_Bayesopt,'ParallelMethod','clipped-model-prediction',...
                    'MaxObjectiveEvaluations',Model.Max_Bayesian_Iterations,'NumSeedPoints',seed_points,...
                    'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                    'GPActiveSetSize',Model.Max_Bayesian_Iterations,...
                    'KernelFunction',Model.KernelFunction);
                
            case 'surrogateopt'
                if Model.ShowPlots
                    plotopt = 'surrogateoptplot';
                else
                    plotopt = [];
                end
                % surrogateopt options
                options = optimoptions('surrogateopt','Display','iter',...
                    'MaxFunctionEvaluations',Model.Max_Bayesian_Iterations,...
                    'MinSurrogatePoints',Model.MinSurrogatePoints,...
                    'MinSampleDistance',Model.MinSampleDistance,...
                    'PlotFcn',plotopt,...
                    'UseParallel',Model.Parallel_Bayesopt,...
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
    if strcmpi(Model.initial_opt_type,'bayesopt') && strcmpi(Model.second_opt_type,'bayesopt') && ~isfile('secondary_completed')
        if Model.ShowPlots
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
        remaining_evals = Model.Max_Secondary_Iterations + Model.Max_Bayesian_Iterations - BayesoptResults.NumObjectiveEvaluations;
        
        if remaining_evals > 0
            if isfield(BayesoptResults.Options,'KernelFunction')
                results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                    'PlotFcn',plotopt,'ExplorationRatio',Model.ExplorationRatio_Secondary,...
                    'AcquisitionFunctionName',Model.Secondary_Acquisition_Function,'Verbose',1,...
                    'ParallelMethod','clipped-model-prediction',...
                    'MaxObjectiveEvaluations',remaining_evals,...
                    'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                    'GPActiveSetSize',Model.Max_Secondary_Iterations + Model.Max_Bayesian_Iterations,...
                    'KernelFunction',Model.SecondaryKernelFunction);
            else
                results = resume(BayesoptResults,'IsObjectiveDeterministic',Deterministic,...
                    'PlotFcn',plotopt,'ExplorationRatio',Model.ExplorationRatio_Secondary,...
                    'AcquisitionFunctionName',Model.Secondary_Acquisition_Function,'Verbose',1,...
                    'ParallelMethod','clipped-model-prediction',...
                    'MaxObjectiveEvaluations',remaining_evals,...
                    'OutputFcn',@saveToFile,'SaveFileName',Intermediate_BO_file,...
                    'GPActiveSetSize',Model.Max_Secondary_Iterations + Model.Max_Bayesian_Iterations);
                
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
        switch Model.final_opt_type
            case 'fminsearch'
                % First, reset the parallelization scheme for fminsearch as it has no parallel option
                % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
                if Model.Parallel_Bayesopt || Model.Parallel_LiX_Minimizer
                    Model.Parallel_Bayesopt = false;
                    Model.Parallel_Struct_Min = false;
                    Model.Parallel_LiX_Minimizer = true;
                    Model.MinMDP.Parallel_Min = false;
                elseif Model.Parallel_Struct_Min
                    Model.Parallel_Bayesopt = false;
                    Model.Parallel_Struct_Min = true;
                    Model.Parallel_LiX_Minimizer = false;
                    Model.MinMDP.Parallel_Min = true;
                end
                fun = @(x)LiX_Minimizer(Model,x);

                % Options for the Nelder–Mead search
                optionsNM = optimset('Display','iter','MaxIter',Model.Max_Local_Iterations,...
                    'TolFun',Model.Loss_Convergence,'TolX',Model.Param_Convergence,...
                    'OutputFcn',@outputFcn_secondary_opt);

                if continue_fullopt
                    try
                        dat = load(Intermediate_Secondary_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;

                        remaining_evals = Model.Max_Local_Iterations - length(intermediate_data)+2;
                        if remaining_evals <= 0
                            remaining_evals = 1;
                        end

                        optionsNM.MaxIter = remaining_evals;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Secondary_file)
                            delete(Intermediate_Secondary_file);
                        end
                        try

                            dat = load(Intermediate_Seconary_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;

                            remaining_evals = Model.Max_Local_Iterations - length(intermediate_data)+2;
                            if remaining_evals <= 0
                                remaining_evals = 1;
                            end

                            optionsNM.MaxIter = remaining_evals;
                        catch
                            if isfile(Intermediate_Seconary_backup)
                                delete(Intermediate_Seconary_backup);
                            end
                            disp('Unable to load backup secondary optimization checkpoint file, restarting.')
                            x0 = results.XAtMinObjective{:,:};
                        end
                    end
                else
                    x0 = results.XAtMinObjective{:,:};
                end

                [full_opt_point,~,~,full_opt_results] = fminsearch(fun,x0,optionsNM);
            case 'fminsearchbnd'

                % Constraints
                Ranges = [params.Range];
                lb = Ranges(1:2:end);
                ub = Ranges(2:2:end);
                
                % First, reset the parallelization scheme for fminsearchbnd as it has no parallel option
                % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
                if Model.Parallel_Bayesopt || Model.Parallel_LiX_Minimizer
                    Model.Parallel_Bayesopt = false;
                    Model.Parallel_Struct_Min = false;
                    Model.Parallel_LiX_Minimizer = true;
                    Model.MinMDP.Parallel_Min = false;
                elseif Model.Parallel_Struct_Min
                    Model.Parallel_Bayesopt = false;
                    Model.Parallel_Struct_Min = true;
                    Model.Parallel_LiX_Minimizer = false;
                    Model.MinMDP.Parallel_Min = true;
                end
                fun = @(x)LiX_Minimizer(Model,x);

                % Options for the Nelder–Mead search
                optionsNM = optimset('Display','iter','MaxIter',Model.Max_Local_Iterations,...
                    'MaxFunEvals',Model.MaxFunEvals,'TolFun',Model.Loss_Convergence,'TolX',Model.Param_Convergence,...
                    'OutputFcn',@outputFcn_secondary_opt);

                if continue_fullopt
                    try
                        dat = load(Intermediate_Secondary_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;

                        remaining_evals = Model.Max_Local_Iterations - length(intermediate_data)+2;

                        if remaining_evals <= 0
                            remaining_evals = 1;
                        end

                        optionsNM.MaxIter = remaining_evals;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Secondary_file)
                            delete(Intermediate_Secondary_file);
                        end
                        try
                            dat = load(Intermediate_Seconary_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;

                            remaining_evals = Model.Max_Local_Iterations - length(intermediate_data)+2;

                            if remaining_evals <= 0
                                remaining_evals = 1;
                            end

                            optionsNM.MaxIter = remaining_evals;
                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Seconary_backup)
                                delete(Intermediate_Seconary_backup);
                            end
                            remaining_evals = Model.Max_Local_Iterations;
                            x0 = results.XAtMinObjective{:,:};
                        end
                    end
                else
                    remaining_evals = Model.Max_Local_Iterations;
                    x0 = results.XAtMinObjective{:,:};
                end

                [full_opt_point,~,exitflag,full_opt_results] = fminsearchbnd(fun,x0,lb,ub,optionsNM);

                % Check to see if the max local iterations is exceeded
                if exitflag == 0 && remaining_evals > 1 && Model.switch_final_opt % calculation exceeded number of allowed iterations
                    % switch to nealder-mead simplex minimization
                    disp('Swicthing local optimization to patternsearch')
                    Model.final_opt_type = 'patternsearch';

                    % Overwrite the input file if it exists
                    if batch_subm
                        save(Input_Model,'Model','-mat')
                    end

                    Bayesian_Optimize_LiX_Parameters(Model)
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
                        dat = load(Intermediate_Secondary_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;    

                        % Establish initial mesh size
                        if isfield(intermediate_data(end).optimValues,'meshsize')
                            init_meshsize = intermediate_data(end).optimValues.meshsize;
                        else
                            init_meshsize = 0.1;
                        end

                        current_iterations = length(intermediate_data);
                        max_iter = Model.Max_Local_Iterations - current_iterations;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Secondary_file)
                            delete(Intermediate_Secondary_file);
                        end
                        try
                            dat = load(Intermediate_Seconary_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;    

                            % Establish initial mesh size
                            if isfield(intermediate_data(end).optimValues,'meshsize')
                                init_meshsize = intermediate_data(end).optimValues.meshsize;
                            else
                                init_meshsize = 0.1;
                            end

                            current_iterations = length(intermediate_data);
                            max_iter = Model.Max_Local_Iterations - current_iterations;

                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Seconary_backup)
                                delete(Intermediate_Seconary_backup);
                            end
                            x0 = results.XAtMinObjective{:,:};
                            init_meshsize = 0.1;
                            max_iter = Model.Max_Local_Iterations;
                        end
                    end
                else
                    x0 = results.XAtMinObjective{:,:};
                    init_meshsize = 0.1;
                    max_iter = Model.Max_Local_Iterations;
                end

        %         optionsNM = optimset('Display','iter','MaxIter',50,'MaxFunEvals',Model.MaxFunEvals,...
        %             'TolFun',Model.Loss_Convergence,'TolX',Model.Param_Convergence);

                if max_iter <= 0
                    max_iter = 1;
                end

                if Model.ShowPlots
                    plotopt = 'psplotbestf';
                else
                    plotopt = [];
                end

                options = optimoptions(@patternsearch,'Display','iter','MaxIterations',max_iter,...
                    'UseParallel',Model.Parallel_Bayesopt,'UseVectorized',false,'PlotFcn',plotopt,...
                    'InitialMeshSize',init_meshsize,'StepTolerance',Model.Param_Convergence,'FunctionTolerance',Model.Loss_Convergence,...
                    'PollOrderAlgorithm','Success','Cache','On','UseCompletePoll',true,...
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
        %                 save(Input_Model,'Model','-mat')
        %             end
        %             
        %             Bayesian_Optimize_LiX_Parameters(Model)
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
                        dat = load(Intermediate_Secondary_file);
                        intermediate_data = dat.intermediate_data;
                        x0 = intermediate_data(end).x;    

                        current_iterations = length(intermediate_data);
                        max_iter = Model.Max_Local_Iterations - current_iterations;
                    catch
                        disp('Unable to load secondary optimization checkpoint file, attempting to load backup.')
                        if isfile(Intermediate_Secondary_file)
                            delete(Intermediate_Secondary_file);
                        end
                        try
                            dat = load(Intermediate_Seconary_backup,'-mat');
                            intermediate_data = dat.intermediate_data;
                            x0 = intermediate_data(end).x;    

                            current_iterations = length(intermediate_data);
                            max_iter = Model.Max_Local_Iterations - current_iterations;
                        catch
                            disp('Unable to load backup secondary optimization checkpoint file, starting new.')
                            if isfile(Intermediate_Seconary_backup)
                                delete(Intermediate_Seconary_backup);
                            end
                            x0 = results.XAtMinObjective{:,:};
                            max_iter = Model.Max_Local_Iterations;
                        end
                    end
                else
                    x0 = results.XAtMinObjective{:,:};
                    max_iter = Model.Max_Local_Iterations;
                end

                if max_iter <= 0
                    max_iter = 1;
                end

                % For unconstrained problem with fmincon, OptimalityTolerance is the max of the gradient
                options = optimoptions(@fmincon,'Display','iter','Algorithm','active-set',...
                    'DiffMaxChange',1e-1,'OptimalityTolerance',1e-1,'UseParallel',Model.Parallel_Bayesopt,...
                    'MaxIterations',max_iter,'FiniteDifferenceStepSize',sqrt(eps),...
                    'StepTolerance',Model.Param_Convergence,'FunctionTolerance',Model.Loss_Convergence,...
                    'FiniteDifferenceType','forward','MaxFunctionEvaluations',Inf,...
                    'OutputFcn',@outputFcn_secondary_opt);

                [full_opt_point,~,~,full_opt_results,~,~,~] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
            case 'none' % if no final optimization is selected, I still want to output the best bayesian optimization result
                full_opt_results = struct;
                full_opt_results.iterations = 0;
                full_opt_results.funccount = 0;
                full_opt_results.algorithm = 'none';
                full_opt_results.message = 'Local optimization skipped. Results are for best global optimization point.';
                full_opt_point = results.XAtMinObjective{:,:};
            otherwise
                error('Unknown final optimization scheme. Choose one of "fminsearch", "patternsearch", "fminsearchbnd", or "none"')
        end
        
        % Save a structure containing some calculation properties
        Calculation_properties.Salt = Model.Salt;
        Calculation_properties.Theory = Model.Theory;
        Calculation_properties.Fix_Charge = Model.Fix_Charge;
        Calculation_properties.Additivity = Model.Additivity;
        Calculation_properties.Additional_MM_Disp = Model.Additional_MM_Disp;
        Calculation_properties.Additional_GAdjust = Model.Additional_GAdjust;
        Calculation_properties.SigmaEpsilon = Model.SigmaEpsilon;
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
            Bayesian_Optimize_LiX_Parameters(Input_Model)
            return
        end
        Calculation_properties = dat.Calculation_properties;
        full_opt_results = dat.full_opt_results;
        full_opt_point = dat.full_opt_point;
    end
    
    %% Final test of parameters on all structures for output
    % If using Parallel_Bayesopt, change it to Parallel_LiX_Minimizer
    if Model.Parallel_Bayesopt || Model.Parallel_LiX_Minimizer
        Model.Parallel_Bayesopt = false;
        Model.Parallel_Struct_Min = false;
        Model.Parallel_LiX_Minimizer = true;
        Model.MinMDP.Parallel_Min = false;
    elseif Model.Parallel_Struct_Min
        Model.Parallel_Bayesopt = false;
        Model.Parallel_LiX_Minimizer = false;
        Model.MinMDP.Parallel_Min = true;
    end
    Model.Delete_Equil = false; % save the final MP calculation directories
    Model.Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'NiAs' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
    Model.Verbose = true;
    [loss,~,UserData] = LiX_Minimizer(Model,full_opt_point,...
        'Extra_Properties',true,'Therm_Prop_Override',true);
    
    if Model.SigmaEpsilon && Model.Additivity && strcmp(Model.Theory,'JC')        
        % Sigma scaling
        sigma_MM = full_opt_point(1);
        sigma_XX = full_opt_point(2);

        % Epsilon scaling
        Epsilon_MM = full_opt_point(3);
        Epsilon_XX = full_opt_point(4);
        
        Sigma_MX   = (sigma_MM + sigma_XX)/2;
        Epsilon_MX = sqrt(Epsilon_MM*Epsilon_XX);
        
        Pars(1:2) = full_opt_point(1:2);
        Pars(3)   = Sigma_MX;
        Pars(4:5) = full_opt_point(3:4);
        Pars(6)   = Epsilon_MX;
        
        if ~Model.Fix_Charge
            Pars(7) = full_opt_point(5);
            if Model.Additional_MM_Disp
                Pars(4) = Pars(4) + full_opt_point(6);
            end
        else
            Pars(7) = Model.Q_value;
            if Model.Additional_MM_Disp
                Pars(4) = Pars(4) + full_opt_point(5);
            end
        end
        
    elseif Model.Additivity && strcmp(Model.Theory,'JC')

        Scale.D.MM = full_opt_point(1);
        Scale.D.XX = full_opt_point(2);

        % Repulsion
        Scale.R.MM = full_opt_point(3);
        Scale.R.XX = full_opt_point(4);

        PotSettings = Initialize_MD_Settings;
        PotSettings.Salt = Model.Salt;
        PotSettings.S = Init_Scaling_Object;
        [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

        % Unscaled
        MX_Epsilon = MXParams.epsilon;
        MX_Sigma   = MXParams.sigma;

        MX_R = 4*MX_Epsilon*MX_Sigma^12;
        MX_D = 4*MX_Epsilon*MX_Sigma^6;

        % Scaled
        MM_Epsilon = MMParams.epsilon*(Scale.D.MM^2)*(1/Scale.R.MM);
        MM_Sigma = MMParams.sigma*(1/(Scale.D.MM^(1/6)))*(Scale.R.MM^(1/6));

        XX_Epsilon = XXParams.epsilon*(Scale.D.XX^2)*(1/Scale.R.XX);
        XX_Sigma = XXParams.sigma*(1/(Scale.D.XX^(1/6)))*(Scale.R.XX^(1/6));

        MX_Epsilon = sqrt(MM_Epsilon*XX_Epsilon);
        MX_Sigma   = (MM_Sigma + XX_Sigma)/2;

        MX_R_scaled = 4*MX_Epsilon*MX_Sigma^12;
        MX_D_scaled = 4*MX_Epsilon*MX_Sigma^6;

        Scale.D.MX = MX_D_scaled/MX_D;
        Scale.R.MX = MX_R_scaled/MX_R;

        Pars(1:2) = full_opt_point(1:2);
        Pars(3)   = Scale.D.MX;
        Pars(4:5) = full_opt_point(3:4);
        Pars(6)   = Scale.R.MX;
        
        if ~Model.Fix_Charge
            Pars(7) = full_opt_point(5);
            if Model.Additional_MM_Disp
                Pars(1) = Pars(1) + full_opt_point(6);
            end
        else
            Pars(7) = Model.Q_value;
            if Model.Additional_MM_Disp
                Pars(1) = Pars(1) + full_opt_point(5);
            end
        end
    elseif strcmp(Model.Theory,'JC')
        Pars(1:6) = full_opt_point(1:6);
        if Model.Fix_Charge
            Pars(7) = Model.Q_value;
        else
            Pars(7) = full_opt_point(7);
        end
    elseif strcmp(Model.Theory,'TF') % TF model
        % Loose form of exp-C6-C8 model
        if Model.SigmaEpsilon
            % Default model parameters: all length-scale units are nm,
            % energy scale units are kJ/mol
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Model.Salt;
            PotSettings.S = Init_Scaling_Object;
            [defTFMX,defTFMM,defTFXX] = TF_Potential_Parameters(PotSettings);
            
            % Input parameters
            r0_MM = full_opt_point(1); % nm
            r0_XX = full_opt_point(2); % nm
           
            if Model.Additivity
                r0_MX = (r0_MM + r0_XX)/2; % nm
                
                epsilon_MM = full_opt_point(3); % kJ/mol
                epsilon_XX = full_opt_point(4); % kJ/mol
                epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                
                gamma_MX = full_opt_point(5); % Unitless
                gamma_MM = gamma_MX; % Unitless
                gamma_XX = gamma_MX; % Unitless
                pidx = 5;
            else
                r0_MX = full_opt_point(3); % nm
                
                epsilon_MM = full_opt_point(4); % kJ/mol
                epsilon_XX = full_opt_point(5); % kJ/mol
                epsilon_MX = full_opt_point(6); % kJ/mol
                
                gamma_MM = full_opt_point(7); % Unitless
                gamma_XX = full_opt_point(8); % Unitless
                gamma_MX = full_opt_point(9); % Unitless
                pidx = 9;
            end
            
            % Convert to Condensed form
            alpha_MM = gamma_MM/r0_MM;
            alpha_XX = gamma_XX/r0_XX;
            alpha_MX = gamma_MX/r0_MX;
            
            B_MM = 48*epsilon_MM*exp(gamma_MM)/(48 - 7*gamma_MM);
            B_XX = 48*epsilon_XX*exp(gamma_XX)/(48 - 7*gamma_XX);
            B_MX = 48*epsilon_MX*exp(gamma_MX)/(48 - 7*gamma_MX);
            
            C_MM = 4*epsilon_MM*gamma_MM*(r0_MM^6)/(48 - 7*gamma_MM);
            C_XX = 4*epsilon_XX*gamma_XX*(r0_XX^6)/(48 - 7*gamma_XX);
            C_MX = 4*epsilon_MX*gamma_MX*(r0_MX^6)/(48 - 7*gamma_MX);
            
            D_MM = 3*epsilon_MM*gamma_MM*(r0_MM^8)/(48 - 7*gamma_MM);
            D_XX = 3*epsilon_XX*gamma_XX*(r0_XX^8)/(48 - 7*gamma_XX);
            D_MX = 3*epsilon_MX*gamma_MX*(r0_MX^8)/(48 - 7*gamma_MX);
            
            % Convert to scaling w.r.t. TF (Pars output)
            % Dispersion scale
            Pars(1) = C_MM/defTFMM.C;
            Pars(2) = C_XX/defTFXX.C;
            Pars(3) = C_MX/defTFMX.C;
            Pars(4) = D_MM/defTFMM.D;
            Pars(5) = D_XX/defTFXX.D;
            Pars(6) = D_MX/defTFMX.D;
            
            % Repulsive wall exponent scale
            Pars(7) = alpha_MM/defTFMM.alpha;
            Pars(8) = alpha_XX/defTFXX.alpha;
            Pars(9) = alpha_MX/defTFMX.alpha;
            
            % Repulsive wall prefactor scale
            Pars(10) = B_MM/defTFMM.B;
            Pars(11) = B_XX/defTFXX.B;
            Pars(12) = B_MX/defTFMX.B;
            
            % Scaling Coulombic Charge
            if Model.Fix_Charge
                Pars(13) = Model.Q_value;
            else
                Pars(13) = full_opt_point(pidx+1);
            end
            
        % Tight form of exp-C6-C8 model
        else
            % C6 Scaling parameters
            Pars(1:3) = full_opt_point(1:3);

            if Model.Fix_C8
                Scale.D6D.MM = Pars(1);
                Scale.D6D.XX = Pars(2);
                Scale.D6D.MX = Pars(3);

                % Calculate value of C8 using recursive relations
                [Metal,Halide] = Separate_Metal_Halide(Model.Salt);

                % Default TF params: C in units of kJ/mol nm^6, D in units of kJ/mol nm^8
                PotSettings = Initialize_MD_Settings;
                PotSettings.Salt = Model.Salt;
                PotSettings.S = Init_Scaling_Object;
                [TF_MX,TF_MM,TF_XX] = TF_Potential_Parameters(PotSettings);

                % Conversion factors
                Bohr_nm = 0.0529177; % a_0 - > Angstrom
                c6conv = 1e-3/2625.4999/((0.052917726)^6); % J/mol nm^6 - > au (from sourcecode)
                J_kJ = 1e-3; % J - > kJ
                Ha_kJmol = 2625.4999; % Ha - > kJ/mol
                c6units = (1/c6conv)*J_kJ; % au - > kJ/mol nm^6
                c8units = (Ha_kJmol)*(Bohr_nm^8); % au - > kJ/mol nm^8

                % Factor used to calculate C8
                sqrt_Q.Li = 5.019869340000000;
                sqrt_Q.Na = 6.585855360000000;
                sqrt_Q.K  = 7.977627530000000;
                sqrt_Q.Rb = 9.554616980000000;
                sqrt_Q.Cs = 11.02204549000000;
                sqrt_Q.F  = 2.388252500000000;
                sqrt_Q.Cl = 3.729323560000000;
                sqrt_Q.Br = 4.590896470000000;
                sqrt_Q.I  = 5.533218150000000;

                % Calculate Scaled C8 using recursion relation from D3 paper
                C8.MM = 3.0*(Scale.D6D.MM*TF_MM.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Metal)*c8units; % in kJ/mol nm^8
                C8.XX = 3.0*(Scale.D6D.XX*TF_XX.C/c6units)*sqrt_Q.(Halide)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8
                C8.MX = 3.0*(Scale.D6D.MX*TF_MX.C/c6units)*sqrt_Q.(Metal)*sqrt_Q.(Halide)*c8units; % in kJ/mol nm^8

                % Update the scaling
                Pars(4) = C8.MM/TF_MM.D;
                Pars(5) = C8.XX/TF_XX.D;
                Pars(6) = C8.MX/TF_MX.D;
                pidx = 3;
            else
                Pars(4:6) = full_opt_point(4:6);
                pidx = 6;
            end

            % Alpha (TF exponential steepness repulsive parameter)
            if Model.Fix_Alpha
                Pars(7:9) = 1;
            else
                Pars(7) = full_opt_point(pidx+1);
                Pars(8) = full_opt_point(pidx+2);
                Pars(9) = full_opt_point(pidx+3);
                pidx = pidx+3;
            end

            % Repulsive wall prefactor
            Pars(10) = full_opt_point(pidx+1);
            Pars(11) = full_opt_point(pidx+2);
            Pars(12) = full_opt_point(pidx+3);
            pidx = pidx+3;

            if Model.Fix_Charge
                Pars(13) = Model.Q_value;
            else
                Pars(13) = full_opt_point(pidx+1);
            end
        end
    elseif strcmp(Model.Theory,'BH') % Buckingham model
        
        % Loose form of exp-C6 model
        if Model.SigmaEpsilon
            
            % Default model parameters: all length-scale units are nm, energy scale units are kJ/mol
            PotSettings = Initialize_MD_Settings;
            PotSettings.Salt = Model.Salt;
            PotSettings.S = Init_Scaling_Object;
            [defBHMX,defBHMM,defBHXX] = BH_Potential_Parameters(PotSettings);
            
            % Input parameters
            r0_MM = full_opt_point(1); % nm
            r0_XX = full_opt_point(2); % nm
            
            if Model.Additivity
                r0_MX = (r0_MM + r0_XX)/2; % nm
                
                epsilon_MM = full_opt_point(3); % kJ/mol
                epsilon_XX = full_opt_point(4); % kJ/mol
                epsilon_MX = sqrt(epsilon_MM*epsilon_XX); % kJ/mol
                
                gamma_MX = full_opt_point(5); % Unitless
                gamma_MM = gamma_MX; % Unitless
                gamma_XX = gamma_MX; % Unitless
                pidx = 5;
            else
                r0_MX = full_opt_point(3); % nm
                
                epsilon_MM = full_opt_point(4); % kJ/mol
                epsilon_XX = full_opt_point(5); % kJ/mol
                epsilon_MX = full_opt_point(6); % kJ/mol
                
                gamma_MM = full_opt_point(7); % Unitless
                gamma_XX = full_opt_point(8); % Unitless
                gamma_MX = full_opt_point(9); % Unitless
                pidx = 9;
            end
            
            % Convert to Condensed form
            alpha_MM = gamma_MM/r0_MM;
            alpha_XX = gamma_XX/r0_XX;
            alpha_MX = gamma_MX/r0_MX;
            
            B_MM = 6*epsilon_MM*exp(gamma_MM)/(gamma_MM - 6);
            B_XX = 6*epsilon_XX*exp(gamma_XX)/(gamma_XX - 6);
            B_MX = 6*epsilon_MX*exp(gamma_MX)/(gamma_MX - 6);
            
            C_MM = epsilon_MM*gamma_MM*(r0_MM^6)/(gamma_MM - 6);
            C_XX = epsilon_XX*gamma_XX*(r0_XX^6)/(gamma_XX - 6);
            C_MX = epsilon_MX*gamma_MX*(r0_MX^6)/(gamma_MX - 6);
            
            % Convert to scaling w.r.t. default BH
            Pars(1) = C_MM/defBHMM.C;
            Pars(2) = C_XX/defBHXX.C;
            Pars(3) = C_MX/defBHMX.C;
            
            Pars(4) = alpha_MM/defBHMM.alpha;
            Pars(5) = alpha_XX/defBHXX.alpha;
            Pars(6) = alpha_MX/defBHMX.alpha;
            
            Pars(7) = B_MM/defBHMM.B;
            Pars(8) = B_XX/defBHXX.B;
            Pars(9) = B_MX/defBHMX.B;
            
            % Scaling Coulombic Charge
            if Model.Fix_Charge
                Pars(10) = Model.Q_value;
            else
                Pars(10) = full_opt_point(pidx+1);
            end
            
        % Tight form of exp-C6 model
        else

            if Model.Additivity
                Scale.D.MM = full_opt_point(1);
                Scale.D.XX = full_opt_point(2);
                Scale.D.MX = sqrt(Scale.D.MM*Scale.D.XX);

                % Repulsion
                Scale.R.MM = full_opt_point(3);
                Scale.R.XX = full_opt_point(4);
                Scale.R.MX = sqrt(Scale.R.MM*Scale.R.XX);

                Pars(1:2) = full_opt_point(1:2);
                Pars(3)   = Scale.D.MX;
                Pars(4:5) = full_opt_point(3:4);
                Pars(6)   = Scale.R.MX;
                pidx = 4;
            else
                Pars(1:6) = full_opt_point(1:6);
                pidx = 6;
            end

            % Alpha (TF exponential steepness repulsive parameter)
            if Model.Fix_Alpha
                Pars(7:9) = 1;
            else
                Pars(7) = full_opt_point(pidx+1);
                Pars(8) = full_opt_point(pidx+2);
                Pars(9) = full_opt_point(pidx+3);
                pidx = pidx+3;
            end

            if ~Model.Fix_Charge
                Pars(10) = full_opt_point(pidx+1);
                pidx = pidx+1;

                if Model.Additional_MM_Disp
                    Pars(1) = Pars(1) + full_opt_point(pidx+1);
                end
            else
                Pars(10) = Model.Q_value;
                if Model.Additional_MM_Disp
                    Pars(1) = Pars(1) + full_opt_point(pidx+1);
                end
            end
        end
        
    end
    
    % Add the Gaussians Params
    if ~isempty(Model.Additional_GAdjust)
        N_Pars = length(Pars);
        N_param = length(full_opt_point);

        N_adj = length(Model.Additional_GAdjust);
        added_params = N_adj*3;

        par_idx = N_Pars+1;
        param_idx = N_param-added_params+1;
        
        Pars(par_idx:par_idx-1+added_params) = full_opt_point(param_idx:N_param);
    end
    
    if Model.Additional_Function.MM.N >= 0
        Ex_fun = 1;
    else
        Ex_fun = 0;
    end

    if Model.Additional_Function.XX.N >= 0
        Ex_fun = Ex_fun + 1;
    end
    if Model.Additional_Function.MX.N >= 0
        Ex_fun = Ex_fun + 1;
    end

    if Ex_fun > 0
        L_pars = length(Pars);
        Pars(L_pars+1:L_pars+Ex_fun) = full_opt_point(end+1-Ex_fun:end);
    end
    
    Minimization_Data = UserData.Minimization_Data;
    Finite_T_Data = UserData.Finite_T_Data;
    format long g
    En = zeros(size(Model.Structures));
    for idx = 1:length(Model.Structures)
        En(idx) = Minimization_Data{idx}.E;
        switch Model.Structures{idx}
            case {'Rocksalt' 'Sphalerite'}
                En(end+1) = Minimization_Data{idx}.a; %#ok<AGROW>
            case {'Wurtzite' 'NiAs' 'FiveFive'}
                En(end+1) = Minimization_Data{idx}.a; %#ok<AGROW>
                En(end+1) = Minimization_Data{idx}.c; %#ok<AGROW>
        end
        disp(Model.Structures{idx})
        disp(Minimization_Data{idx})
        Minimization_Data{idx}.Structure = Model.Structures{idx};
    end
    disp(repmat('*',1,120))
    disp('Finite Temperature Data (Enthalpy in kJ/mol, Entropy in J/(mol K), Volume in A^3/Forumla Unit, Temperature in K)')
    disp(repmat('*',1,120))
    disp(Finite_T_Data)
    disp(repmat('*',1,120))
    
    disp(['Final Optimized Loss: ' num2str(loss,'%.10e')])
    if Model.SigmaEpsilon && strcmp(Model.Theory,'JC')
        disp(['Parameters (Sigma/Epsilon): ' num2str(Pars,'\t%.15e')])
    else
        disp(['Parameters (Dispersion/Repulsion Scale): ' num2str(Pars,'\t%.15e')])
    end
    disp(['0 K Lattice Energies (kJ/mol): ' num2str(En,'\t%.15e')])
    
    % Save final results
    save(Full_opt_filename,'full_opt_results','loss','full_opt_point',...
        'Minimization_Data','Finite_T_Data','Pars','Calculation_properties');
    diary off
    
end
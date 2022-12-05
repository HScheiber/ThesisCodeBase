function [Tm_estimate,WorkDir,Aborted,T_dat] = Find_Melting_Point(Settings)
total_timer = tic;
Aborted = false;

if ~isstruct(Settings)
    Settings = load(Settings,'-mat').Settings;
end

% Grab some additional settings that depend on inputs
if ~isfield(Settings,'WorkDir')
    [Settings.WorkDir,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
end
WorkDir = Settings.WorkDir;
if ~isfolder(WorkDir)
    mkdir(WorkDir)
end
if ~isfield(Settings,'OuterDir')
    Settings.OuterDir = WorkDir;
end
Settings = Update_MD_Settings(Settings);
Settings.Submit_Jobs = false; % should be set to false for MP jobs

Settings.mdrun_opts = regexprep(Settings.mdrun_opts,' -maxh [0-9]+','','once');
Settings.CurrentTFile = fullfile(Settings.WorkDir,[Settings.JobName '_MP.mat']);
Settings.PrevTFile = fullfile(Settings.WorkDir,[Settings.JobName '_MP.mat.PREV']);
Settings.BayesoptCheckpointFile = fullfile(Settings.WorkDir,[Settings.JobName '_MP_bayesopt.mat']);
Settings.RefGeomDir = fullfile(Settings.WorkDir,'InitGeom');
ResultsFile = fullfile(Settings.WorkDir,[Settings.JobName '_MPResults.mat']);

% Make RefGeomDir
if ~exist(Settings.RefGeomDir,'dir')
    mkdir(Settings.RefGeomDir);
end
save(fullfile(Settings.WorkDir,'Calc_Settings.mat'),'Settings')
diary off
diary(fullfile(Settings.WorkDir,'Calculation_diary.log'))

% Check if calculation already completed
if isfile(ResultsFile) && ~Settings.Continue
    try
        T_dat = load(ResultsFile,'T_dat').T_dat;
        Max_T_freeze = max(T_dat.T_Freeze_Trace);
        Min_T_Melt = min(T_dat.T_Melt_Trace);
        T_dat.dT = [Max_T_freeze Min_T_Melt];
        Tm_estimate = mean(T_dat.dT);
        if ~isfield(T_dat,'Alt_Structure')
            T_dat.Alt_Structure = false;
        end
        
        if Settings.Verbose
            disp(repmat('*',1,40));
        end
        if T_dat.Alt_Structure
            if Settings.Verbose
                disp(['Calculation previously aborted at T = ' num2str(T_dat.T,'%.4f') ' K.']);
            end
            Aborted = true;
        else
            if Settings.Verbose
                disp(['Calculation previously completed. Melting Point estimated at Tm = ' num2str(Tm_estimate,'%.4f') ' K. Error bounds: Tm = ' ...
                    num2str(Max_T_freeze,'%.4f') ' - ' num2str(Min_T_Melt,'%.4f') ' K.']);
            end
        end
        if Settings.Verbose
            disp(repmat('*',1,40));
        end
        diary off
        if isfield(Settings,'Diary_Loc') && ~isempty(Settings.Diary_Loc)
            diary(Settings.Diary_Loc)
        end
        return
    catch
    end
elseif isfile(ResultsFile) && Settings.Continue
    delete(ResultsFile)
end

% Check if continuation possible
Loaded_Prev_Data = false;
if isfile(Settings.CurrentTFile)
    try
        T_dat = load(Settings.CurrentTFile,'T_dat').T_dat;
        Loaded_Prev_Data = true;
    catch
        delete(Settings.CurrentTFile)
        if isfile(Settings.PrevTFile)
            try
                T_dat = load(Settings.PrevTFile,'-mat','T_dat').T_dat;
                Loaded_Prev_Data = true;
                copyfile(Settings.PrevTFile,Settings.CurrentTFile)
            catch
                delete(Settings.PrevTFile)
            end
        end
    end
    
    if ~isfield(T_dat,'Alt_Structure')
        T_dat.Alt_Structure = false;
    end
end

if ~Loaded_Prev_Data
    % Create suitable initial conditions (may be re-used for other steps)
    Output = Setup_LiX_Simulation(Settings);
    
    if Output.Aborted
        if Output.SolidMelted
            
            % Solid melted
            time_to_phase_change = 1;
            f = 1/time_to_phase_change;
            df = 5000/(time_to_phase_change);
            
            % Update the current T_mp upper error bound
            T_dat.T = Settings.T0;
            T_dat.dT = [Settings.lb Settings.T0]; % Error bounds on the melting temperature, initialized as the prior bounds
            T_dat.df_bracket = [nan df]; % The gradient for each member of the melting point bounds
            T_dat.T_ref = Settings.T0;
            T_dat.Ref_Density_Trace = Settings.T0;
            T_dat.T_Trace = [Settings.T0]; % Keeps track of temperature at each step of iteration
            T_dat.f_Trace = f; % Keeps track of time to phase change at each step of iteration
            T_dat.df_Trace = df; % Keeps track of gradient at each step of iteration
            T_dat.Freeze_Trace = false; % Keeps track of whether the process froze at each iteration
            T_dat.Melt_Trace = true; % Keeps track of whether the process melted at each iteration
            T_dat.T_Freeze_Trace = []; % This will be a trace of only the temperatures where freezing was detected
            T_dat.T_Melt_Trace = Settings.T0; % This will be a trace of only the temperatures where melting was detected
            T_dat.Alt_Structure = false; % This is a switch to kill the calculation if the wrong structure is freezing out
            
        elseif Output.LiquidFroze || Output.LiquidAmorphous
            
            % Liquid froze or became amorphous
            time_to_phase_change = 1;
            f = 1/time_to_phase_change;
            df = -5000/(time_to_phase_change);
            
            % Update the current T_mp lower error bound
            T_dat.T = Settings.T0;
            T_dat.dT = [Settings.T0 Settings.ub]; % Error bounds on the melting temperature, initialized as the prior bounds
            T_dat.df_bracket = [df nan]; % The gradient for each member of the melting point bounds
            T_dat.T_ref = Settings.T0;
            T_dat.Ref_Density_Trace = Settings.T0;
            T_dat.T_Trace = Settings.T0; % Keeps track of temperature at each step of iteration
            T_dat.f_Trace = f; % Keeps track of time to phase change at each step of iteration
            T_dat.df_Trace = df; % Keeps track of gradient at each step of iteration
            T_dat.Freeze_Trace = true; % Keeps track of whether the process froze at each iteration
            T_dat.Melt_Trace = false; % Keeps track of whether the process melted at each iteration
            T_dat.T_Freeze_Trace = Settings.T0; % This will be a trace of only the temperatures where freezing was detected
            T_dat.T_Melt_Trace = []; % This will be a trace of only the temperatures where melting was detected
            T_dat.Alt_Structure = false; % This is a switch to kill the calculation if the wrong structure is freezing out
            
        else
            % Solid or liquid crystallized to an unwanted structure
            % homogeneously, or calculation failed due to bad potential
            % In this case, abort the calculation...
            f = -1;
            df = 0;
            T_dat.Alt_Structure = true;
            
            T_dat.T = Settings.T0;
            T_dat.dT = [Settings.lb Settings.ub]; % Error bounds on the melting temperature, initialized as the prior bounds
            T_dat.df_bracket = [nan nan]; % The gradient for each member of the melting point bounds
            T_dat.T_ref = Settings.T0;
            T_dat.Ref_Density_Trace = Settings.T0;
            T_dat.T_Trace = Settings.T0; % Keeps track of temperature at each step of iteration
            T_dat.f_Trace = f; % Keeps track of time to phase change at each step of iteration
            T_dat.df_Trace = df; % Keeps track of gradient at each step of iteration
            T_dat.Freeze_Trace = []; % Keeps track of whether the process froze at each iteration
            T_dat.Melt_Trace = []; % Keeps track of whether the process melted at each iteration
            T_dat.T_Freeze_Trace = []; % This will be a trace of only the temperatures where freezing was detected
            T_dat.T_Melt_Trace = []; % This will be a trace of only the temperatures where melting was detected
            T_dat.Alt_Structure = true; % This is a switch to kill the calculation if the wrong structure is freezing out
        end
        save(Settings.CurrentTFile,'T_dat')
    else
        % Make a calculation history file with with current calculation state and
        % history, accessible to all previous functions
        T_dat.T = Settings.T0;
        T_dat.dT = [Settings.lb Settings.ub]; % Error bounds on the melting temperature, initialized as the prior bounds
        T_dat.df_bracket = [nan nan]; % The gradient for each member of the melting point bounds
        T_dat.T_ref = Settings.T0;
        T_dat.Ref_Density_Trace = Settings.T0;
        T_dat.T_Trace = []; % Keeps track of temperature at each step of iteration
        T_dat.f_Trace = []; % Keeps track of time to phase change at each step of iteration
        T_dat.df_Trace = []; % Keeps track of gradient at each step of iteration
        T_dat.Freeze_Trace = []; % Keeps track of whether the process froze at each iteration
        T_dat.Melt_Trace = []; % Keeps track of whether the process melted at each iteration
        T_dat.T_Freeze_Trace = []; % This will be a trace of only the temperatures where freezing was detected
        T_dat.T_Melt_Trace = []; % This will be a trace of only the temperatures where melting was detected
        T_dat.Alt_Structure = false; % This is a switch to kill the calculation if the wrong structure is freezing out
        save(Settings.CurrentTFile,'T_dat')

        % Rename the gro file at the given reference temperature
        Strucure_Ref_File_old = fullfile(Settings.WorkDir,[Settings.JobName '.' Settings.CoordType]);
        Strucure_Ref_File = fullfile(Settings.RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '.' Settings.CoordType]);
        movefile(Strucure_Ref_File_old, Strucure_Ref_File)

        % Move the solid info file
        Sol_Ref_File_old = fullfile(Settings.WorkDir,[Settings.JobName '_SolInfo.mat']);
        Sol_Ref_File = fullfile(Settings.RefGeomDir,[Settings.JobName '_' num2str(T_dat.T_ref,'%.4f') '_SolInfo.mat']);
        movefile(Sol_Ref_File_old, Sol_Ref_File)
    end
else
    if Settings.Verbose
        disp('Calculation previously started. Restoring optimizer to previous state...')
    end
end

% Define the function that checks the time to reach melting/freezing point
fun = @(T)Melting_Point_Check(T,Settings);

switch lower(Settings.Optimizer)
case 'fmincon'

    % Options for fmincon
    if Settings.Verbose
        display_opt = 'iter';
    else
        display_opt = 'none';
    end
    options = optimoptions(@fmincon,'Display',display_opt,'Algorithm','trust-region-reflective',...
        'OptimalityTolerance',1/Settings.MaxCheckTime,'ObjectiveLimit',0,...
        'UseParallel',false,'MaxIterations',Settings.MaxMPIterations,'StepTolerance',Settings.StepTolerance,...
        'FunctionTolerance',Settings.FunctionTolerance,'MaxFunctionEvaluations',Settings.MaxMPIterations,...
        'SpecifyObjectiveGradient',true);

    % Run optimization algorithm
    [Tm,~,~,~,~,~,~] = fmincon(fun,Settings.T0,[],[],[],[],Settings.lb,Settings.ub,[],options);

case 'bayesopt'
    if Settings.Verbose
        display_opt = 1;
    else
        display_opt = 0;
    end
    if Loaded_Prev_Data
        Loaded_checkpoint = true;
        try
            dat = load(Settings.BayesoptCheckpointFile,'-mat');
            results = dat.BayesoptResults;
        catch
            % Catch corrupt files
            if Settings.Verbose
                disp('Unable to load bayesian optimization checkpoint file.')
                disp('Attempting to load backup checkpoint.')
            end

            try
                dat = load([Settings.BayesoptCheckpointFile '.PREV'],'-mat');
                results = dat.BayesoptResults;
            catch
                disp('Warning: Unable to load backup checkpoint of Bayesian Optimization - Restarting Calculation.')
                if isfile(Settings.BayesoptCheckpointFile)
                    delete(Settings.BayesoptCheckpointFile)
                end
                if isfile([Settings.BayesoptCheckpointFile '.PREV'])
                    delete([Settings.BayesoptCheckpointFile '.PREV'])
                end
                Loaded_checkpoint = false;
            end
        end
        
    else
        Loaded_checkpoint = false;
    end
    
    if Loaded_checkpoint
        remaining_evals = Settings.MaxMPIterations - results.NumObjectiveEvaluations;
        results = resume(results,'IsObjectiveDeterministic',false,...
            'PlotFcn','all','ExplorationRatio',Settings.ExplorationRatio,...
            'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',display_opt,...
            'MaxObjectiveEvaluations',remaining_evals,...
            'OutputFcn',@saveToFile,'SaveFileName',Settings.BayesoptCheckpointFile,...
            'GPActiveSetSize',Settings.MaxMPIterations,...
            'KernelFunction',Settings.KernelFunction);
    else
        Tparam = optimizableVariable('T',[Settings.lb Settings.ub],'Type','real'); % Units: K
        results = bayesopt_priv(fun,Tparam,'IsObjectiveDeterministic',false,...
            'ExplorationRatio',Settings.ExplorationRatio,'PlotFcn','all',...
            'AcquisitionFunctionName',Settings.Acquisition_Function,'Verbose',display_opt,...
            'UseParallel',false,'MaxObjectiveEvaluations',Settings.MaxMPIterations,...
            'NumSeedPoints',Settings.N_Seed_Points,'OutputFcn',@saveToFile,...
            'SaveFileName',Settings.BayesoptCheckpointFile,...
            'GPActiveSetSize',Settings.MaxMPIterations,...
            'KernelFunction',Settings.KernelFunction);
    end
    Tm = results.XAtMinObjective{:,:};
% case 'patternsearch'
% 
%     options = optimoptions(@patternsearch,'Display','iter',...
%         'MaxIterations',Settings.MaxMPIterations,'UseParallel',false,...
%         'PlotFcn','psplotbestf','InitialMeshSize',Settings.InitialMeshSize,...
%         'StepTolerance',Settings.StepTolerance,'FunctionTolerance',Settings.FunctionTolerance,...
%         'PollOrderAlgorithm','Success','Cache','On','UseCompletePoll',false,...
%         'PollMethod','GPSPositiveBasis2N','MaxMeshSize',Inf,'MeshTolerance',Settings.StepTolerance,...
%         'OutputFcn',@(x,y,z)outputFcn_patternsearch_opt(x,y,z,Settings),'MaxFunctionEvaluations',Settings.MaxMPIterations);
% 
%     [Tm,~,~,results] = patternsearch(fun,Settings.T0,...
%         [],[],[],[],Settings.lb,Settings.ub,[],options);
case 'mpsearcher'
    Tm = MPSearcher(fun,Settings.T0,Settings.lb,Settings.ub,Settings);
end

% Copy the final melting temperature and history into the results file
copyfile(Settings.CurrentTFile, ResultsFile)
T_dat = load(ResultsFile,'T_dat').T_dat;

% Copy over contents from the T = Tm simulation folder
Tmtxt = num2str(Tm,'%.4f');
T_folder = fullfile(Settings.WorkDir,['T_' Tmtxt]);
Tm_folder = fullfile(Settings.WorkDir,['Tm_' Tmtxt]);
if ~exist(Tm_folder,'dir')
    mkdir(Tm_folder)
end

T_contents = dir(T_folder);
T_contents = T_contents(cellfun(@isempty,regexp({T_contents.name},...
    ['(^\.)|(' Tmtxt '.gro)|(Grompplog.log)'],'once')));
for idx = 1:length(T_contents)
    copyfile(fullfile(T_folder,T_contents(idx).name), fullfile(Tm_folder,T_contents(idx).name))
    if contains(T_contents(idx).name,'.tpr')
        tpr_file = fullfile(Tm_folder,T_contents(idx).name);
    end
end

% Optional: Calculate the density profile along the Z dimension between 100
if Settings.FinalDensityProfile && ~T_dat.Alt_Structure && isfile(tpr_file)
    trr_file = fullfile(Tm_folder,[Settings.JobName '.trr']);
    density_file = fullfile(Settings.WorkDir,[Settings.JobName '_Density_Profile.xvg']);
    gmx_density_cmd = [Settings.wsl 'echo "0 0" ' Settings.pipe ...
        ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_density ...
        ' -f ' windows2unix(trr_file) ' -s ' ...
        windows2unix(tpr_file) ' -o ' windows2unix(density_file) ...
        ' -center -b ' num2str(Settings.MaxCheckTime/2) ...
        ' -e ' num2str(Settings.MaxCheckTime) ' -sl 200'];
    [errcode,outp] = system(gmx_density_cmd);
    if errcode ~= 0
        disp(outp);
        warning(['Error running gmx density. Problem command: ' newline gmx_density_cmd]);
    end
end

% Optional: delete the temporary Temperature check folders
if Settings.Delete_T_History
    T_dirs = dir(Settings.WorkDir);
    T_dirs = T_dirs([T_dirs.isdir]);
    T_dirs = T_dirs(~cellfun(@isempty,regexp({T_dirs.name},'T_([0-9]|\.)+','once')));
    for idx = 1:length(T_dirs)
        rmdir(fullfile(Settings.WorkDir,T_dirs(idx).name),'s')
    end  
end

Tm_estimate = mean(T_dat.dT);
if T_dat.Alt_Structure
    Tm_estimate = nan;
    Aborted = true;
end
if Settings.Verbose
    disp(['Calculation complete. Epalsed Time: ' datestr(seconds(toc(total_timer)),'HH:MM:SS')])
end
diary off
if isfield(Settings,'Diary_Loc') && ~isempty(Settings.Diary_Loc)
    diary(Settings.Diary_Loc)
end

end
function Prev_Data = LoadPrevModelData(Settings,params,subsample)
Prev_Data = [];
% To do:
% 0. Gather all models of interest from saved model folder
% 1. load all X Points
% 2. load all user data
% 3. Load all iteration times
% 4. Load all error values
% 5. if using coupled constraints, load all constraint violation points
% 6a. Check if any X points fall outside of current ranges. If so, exclude al data from those points
% 6b. Check for any X that are duplicate, remove dupes.
% 7. Re-calculate new objective function values based on the user data.
% Check to make sure user data is sufficient to actually calculate the
% current objective function, if not return an error.
% 8. Pass all data to output structure Prev_Data

Saved_Models_dir = fullfile(Settings.project,Settings.Project_Directory_Name,'Completed',Settings.Salt);

for idx = 1:length(Settings.Initialize_From_Model)
    Prev_model = Settings.Initialize_From_Model{idx};
    files = {dir(fullfile(Saved_Models_dir,[Settings.Salt '_' Settings.Theory '_Model_' Prev_model '*_data.mat'])).name};
    Models{idx} = unique(regexp(files,[Prev_model '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models{idx});
    Models{idx} = Models{idx}(nonempty_idx);
end

Models = unique(horzcat(Models{:}));
N_Models = numel(Models);

% Initialize data arrays
InitialX                        = cell(N_Models,1);
InitialConstraintViolations     = cell(N_Models,1);
InitialErrorValues              = cell(N_Models,1);
InitialUserData                 = cell(N_Models,1);
InitialObjectiveEvaluationTimes = cell(N_Models,1);
InitialIterationTimes           = cell(N_Models,1);

for idx = N_Models:-1:1
    Model = Models{idx};
    
    % Find the fully optimized model
    dat_file = fullfile(Saved_Models_dir,[Settings.Salt '_' Settings.Theory '_Model_' Model '_data.mat']);
    
    if ~isfile(dat_file)
        disp(['No data available for for: ' Settings.Salt ', ' Settings.Theory ', Model ' Model '.']);
        InitialX(idx) = [];
        InitialConstraintViolations(idx) = [];
        InitialErrorValues(idx) = [];
        InitialUserData(idx) = [];
        InitialObjectiveEvaluationTimes(idx) = [];
        InitialIterationTimes(idx) = [];
        continue
    end
    
    % Load data
    try
        data = load(dat_file).full_data;
        InitialX{idx} = data.bayesopt_results.XTrace;
        
        Within_range = true(size(InitialX{idx},1),1);
        for jdx = 1:length(params)
            Within_range = Within_range & (InitialX{idx}{:,jdx} >= params(jdx).Range(1)) & (InitialX{idx}{:,jdx} <= params(jdx).Range(2));
        end
        InitialX{idx} = InitialX{idx}(Within_range,:);
        InitialConstraintViolations{idx} = data.bayesopt_results.ConstraintsTrace;
        if isempty(InitialConstraintViolations{idx})
            InitialConstraintViolations{idx} = double.empty(size(Within_range,1),0);
        else
            InitialConstraintViolations{idx} = InitialConstraintViolations{idx}(Within_range,:);
        end
        InitialErrorValues{idx} = data.bayesopt_results.ErrorTrace(Within_range);
        InitialUserData{idx} = data.bayesopt_results.UserDataTrace(Within_range);
        InitialObjectiveEvaluationTimes{idx} = data.bayesopt_results.ObjectiveEvaluationTimeTrace(Within_range);
        InitialIterationTimes{idx} = data.bayesopt_results.IterationTimeTrace(Within_range);
        
    catch
        disp(['Could not load data for: ' Settings.Salt ', ' Settings.Theory ', Model ' Model '.']);
        InitialX(idx) = [];
        InitialConstraintViolations(idx) = [];
        InitialErrorValues(idx) = [];
        InitialUserData(idx) = [];
        InitialObjectiveEvaluationTimes(idx) = [];
        InitialIterationTimes(idx) = [];
        continue
    end
end

% Combine arrays
Prev_Data.InitialX = vertcat(InitialX{:});
Prev_Data.InitialConstraintViolations = vertcat(InitialConstraintViolations{:});
Prev_Data.InitialErrorValues = vertcat(InitialErrorValues{:});
Prev_Data.InitialUserData = vertcat(InitialUserData{:});
Prev_Data.InitialObjectiveEvaluationTimes = vertcat(InitialObjectiveEvaluationTimes{:});
Prev_Data.InitialIterationTimes = vertcat(InitialIterationTimes{:});

% Check for and remove duplicates in X
datafields = fields(Prev_Data);
[~,idxUnique,~] = unique(Prev_Data.InitialX,'rows');
for idx = 1:numel(datafields)
    field = datafields{idx};
    Prev_Data.(field) = Prev_Data.(field)(idxUnique,:);
end

% Re-calculate objective function based on user data
Prev_Data.InitialObjective = Loss_Recalculate(Settings,Prev_Data.InitialX,...
    Prev_Data.InitialUserData,data.Settings.Structures);

% Tale a Random subsample
Ndata = numel(Prev_Data.InitialObjective);
if subsample < Ndata
    SampleX = randsample(Ndata,1000);
    datafields = fields(Prev_Data);
    for idx = 1:numel(datafields)
        field = datafields{idx};
        Prev_Data.(field) = Prev_Data.(field)(SampleX,:);
    end
end

% Check if coupled constrains has changed...
if Settings.



end
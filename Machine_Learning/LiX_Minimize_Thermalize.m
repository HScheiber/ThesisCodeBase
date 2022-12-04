%% Inputs
function UserData = LiX_Minimize_Thermalize(Input_Settings,varargin)
% Grab input structure
if isstruct(Input_Settings)
    Settings = Input_Settings;
elseif isfile(Input_Settings)
    Settings_dat = load(Input_Settings,'-mat');
    Settings = Settings_dat.Settings;
else
    error('Input model does not exist, or is not a compatible data structure.')
end

% Optional inputs
p = inputParser;
p.FunctionName = 'LiX_Minimizer';
addOptional(p,'Extra_Properties',Settings.Extra_Properties,@(x)validateattributes(x,{'logical'},{'nonempty'}))

parse(p,varargin{:});
Extra_Properties = p.Results.Extra_Properties;

Settings.MinMDP.Verbose = Settings.Verbose;
if Settings.Parallel_LiX_Minimizer && Settings.Parallel_Struct_Min
    Settings.Parallel_Struct_Min = false;
end

Settings.Minimization_Data = Initialize_Minimization_Data(Settings);
Settings.Finite_T_Data = Initialize_Finite_T_Data(Settings);

%% Parallel Setup
N = length(Settings.Structures);
if Settings.Parallel_LiX_Minimizer
    % Set up matlab parallel features
    Parcores = feature('numcores');
    PrefCores = min(Parcores,N);
    if ~isempty(gcp('nocreate'))
        Cur_Pool = gcp;
        Cur_Workers = Cur_Pool.NumWorkers;

        % Start the parallel pool
        if Cur_Workers < PrefCores
            delete(Cur_Pool);
            % Create a "local" cluster object
            local_cluster = parcluster('local');

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
            ppool = parpool(local_cluster,min(Parcores,N));
        else
            ppool = Cur_Pool;
        end
    else
        ppool = parpool(min(Parcores,N)); 
    end

    f(1:N) = parallel.FevalFuture;

    % Run in parallel with parfavel
    for idx = 1:N
        Settings.Structure = Settings.Structures{idx};
        f(idx) = parfeval(ppool,@Structure_Minimization,1,Settings,...
            'Extra_Properties',Extra_Properties);
    end

    wait(f); % Wait for parallel jobs to finish

    % Collect outputs into cell array
    for idx = 1:N
        Settings.Minimization_Data{idx} = f(idx).fetchOutputs;
    end
%% Serial mode    
else
    for idx = 1:N
        Settings.Structure = Settings.Structures{idx};
        Settings.Minimization_Data{idx} = Structure_Minimization(Settings,...
            'Extra_Properties',Extra_Properties);
    end
end

% Initialize Finite T Data structure and update Settings
Settings.Structure = Settings.Finite_T_Data.Structure;
Settings.Geometry = Default_Crystal(Settings,'Center_Coordinates',true);

strmatch = strcmp(Settings.Finite_T_Data.Structure,Settings.Structures);
Settings.Ref_Density = 1/(Settings.Minimization_Data{strmatch}.V*(0.1^3)); % molecules / nm^3
Settings.Geometry.a = Settings.Minimization_Data{strmatch}.a;
Settings.Geometry.b = Settings.Minimization_Data{strmatch}.b;
Settings.Geometry.c = Settings.Minimization_Data{strmatch}.c;

%% Melting point
if ~isempty(gcp('nocreate')) % close the current ppool, it will likely close on its own anyway, causing issues
    delete(gcp);
end

[~,Settings.JobName,Settings.Full_Model_Name] = GetMDWorkdir(Settings);
Settings.WorkDir = fullfile(Settings.OuterDir,'BestPoint_Thermal','Melting_Point');

Settings.BatchMode = false;
Settings.Submit_Jobs = false;
Settings.Skip_Minimization = true; % Skip the automatic geometry minimization
Settings.RefStructure = Settings.Finite_T_Data.Structure;
Verbose = Settings.Verbose;
Settings.Verbose = true;
[Tm_estimate,~,Aborted,T_dat] = Find_Melting_Point(Settings);
Settings.Verbose = Verbose;

Settings.Finite_T_Data.T_dat = T_dat;
if Aborted
    Settings.Finite_T_Data.MP = nan;
else
    Settings.Finite_T_Data.MP = Tm_estimate;
end

%% High T liquid properties
Settings.WorkDir = fullfile(Settings.OuterDir,'BestPoint_Thermal','Liq_Properties_at_MP');

dd = Settings.dd;
npme = Settings.npme;
Settings.dd = [];
Settings.npme = [];
[~,Settings] = MD_Batch_Template(Settings);
Verbose = Settings.Verbose;
Settings.Verbose = true;
Liq_Output = Calc_Liquid_Properties_at_MP(Settings); % Output is nan if liquid converts to >0.85 solid
Settings.Verbose = Verbose;
Settings.dd = dd;
Settings.npme = npme;
[~,Settings] = MD_Batch_Template(Settings);

Settings.Finite_T_Data.Liquid_V_MP = Liq_Output.Liquid_V_MP;
Settings.Finite_T_Data.Liquid_H_MP = Liq_Output.Liquid_H_MP;
Settings.Finite_T_Data.Liquid_DM_MP = Liq_Output.Liquid_DM_MP; % cm^2 / s
Settings.Finite_T_Data.Liquid_DX_MP = Liq_Output.Liquid_DX_MP;

%% High T solid properties
Settings.WorkDir = fullfile(Settings.OuterDir,'BestPoint_Thermal','Sol_Properties_at_MP');

dd = Settings.dd;
npme = Settings.npme;
Settings.dd = [];
Settings.npme = [];
[~,Settings] = MD_Batch_Template(Settings);
Verbose = Settings.Verbose;
Settings.Verbose = true;
Sol_Output = Calc_Solid_Properties_at_MP(Settings);
Settings.Verbose = Verbose;
Settings.dd = dd;
Settings.npme = npme;
[~,Settings] = MD_Batch_Template(Settings);

Settings.Finite_T_Data.Solid_V_MP = Sol_Output.Solid_V_MP;
Settings.Finite_T_Data.Solid_H_MP = Sol_Output.Solid_H_MP;

Settings.Finite_T_Data.Fusion_dH = Settings.Finite_T_Data.Liquid_H_MP - ...
    Settings.Finite_T_Data.Solid_H_MP;

Settings.Finite_T_Data.Fusion_dV = Settings.Finite_T_Data.Liquid_V_MP - ...
    Settings.Finite_T_Data.Solid_V_MP;


% Delete previous calculations that did not complete
files = dir(Settings.scratch_dir);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name};
prev_calcs = subFolderNames(cellfun(@(x) ~isempty(x),regexp(subFolderNames,'.+?_[S|M|L|O]P','once')));
for idx = 1:length(prev_calcs)
    try
        rmdir(fullfile(Settings.scratch_dir,prev_calcs{idx}),'s')
    catch
        if Settings.Verbose
            disp(['Unable to remove failed calculation directory: ' fullfile(Settings.scratch_dir,prev_calcs{idx})])
        end
    end
end

Full_opt_filename = fullfile(Settings.OuterDir,[Settings.JobName '_fullopt.mat']);
UserData.Finite_T_Data = Settings.Finite_T_Data;
UserData.Minimization_Data = Settings.Minimization_Data;
Minimization_Data = UserData.Minimization_Data;
Finite_T_Data = UserData.Finite_T_Data;
save(Full_opt_filename,'Minimization_Data','Finite_T_Data');
end
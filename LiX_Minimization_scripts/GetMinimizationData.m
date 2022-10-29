Settings = Initialize_MD_Settings;

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Theories = {'BF'};
Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};

Settings.JobSettings.Cores = 8;
Settings.JobSettings.MPI_Ranks = 8; % Sets the number of MPI ranks (distributed memory parallel processors). -1 for auto
Settings.JobSettings.OMP_Threads = 1; % Set the number of OMP threads per MPI rank
Settings.MinMDP.Energy_Tol = 1e-3; % kJ/mol
Settings.MinMDP.Gradient_Tol_RMS = 1e-3; % kJ/(mol A)
Settings.MinMDP.Gradient_Tol_Max = 1e-3; % kJ/(mol A)
Settings.Cutoff_Buffer = 1.02;
Settings.MinMDP.RList_Cutoff = 2;
Settings.MinMDP.RCoulomb_Cutoff = 2;
Settings.MinMDP.RVDW_Cutoff = 2;
Settings.GaussianCharge = true;
Settings.initial_opt_type = true;

[~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,~] = MD_Batch_Template(Settings.JobSettings);
Settings.MinMDP.Parallel_Min = false;
setenv('OMP_NUM_THREADS','1');
setenv('GMX_PME_NUM_THREADS','1');
setenv('GMX_PME_NTHREADS','1');
setenv('GMX_OPENMP_MAX_THREADS','1');
setenv('KMP_AFFINITY','disabled');
Settings.mdrun_opts = ' -pin on -ntmpi 1 -ntomp 1';
Settings.gmx = Settings.gmx_loc;

% Set up matlab parallel features
N = length(Theories)*length(Salts)*length(Structures);
Parcores = feature('numcores');
PrefCores = min(Parcores,N);
if ~isempty(gcp('nocreate'))
    Cur_Pool = gcp;
    Cur_Workers = Cur_Pool.NumWorkers;

    % Start the parallel pool
    if Cur_Workers ~= PrefCores
        delete(Cur_Pool);
        % Create a "local" cluster object
        local_cluster = parcluster('local');

        % Modify the JobStorageLocation to a temporary directory
        tmp = tempname;
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
indexes = combvec(1:length(Salts),1:length(Theories),1:length(Structures));
for idx = 1:N
    Settings.Salt = Salts{indexes(1,idx)};
    Settings.Theory = Theories{indexes(2,idx)};
    Settings.Structure = Structures{indexes(3,idx)};
    Settings = Alexandria_Potential_Parameters(Settings,'vdW_Type','WBK');
    
    if strcmp(Settings.Structure,'FiveFive')
        Find_Min_Params = false;
    else
        Find_Min_Params = true;
    end
    
    f(idx) = parfeval(ppool,@Structure_Minimization,1,Settings,'Extra_Properties',true,...
        'Find_Min_Params',Find_Min_Params);
end

wait(f); % Wait for parallel jobs to finish

% Collect outputs into cell array
Data = struct();
for idx = 1:N

    Salt = Salts{indexes(1,idx)};
    Theory = Theories{indexes(2,idx)};
    Structure = Structures{indexes(3,idx)};
    if strcmp(f(idx).State,'finished')
        Data.(Salt).(Theory).(Structure) = f(idx).fetchOutputs;
    else
        Data.(Salt).(Theory).(Structure).E = nan;
    end
end

save(fullfile(Settings.home,'data','MX_Alexandria_Min_Data.mat'),'Data')

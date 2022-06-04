% Tm_vs_final_amt
Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Set60_TvsXliq.mat';
DataKeyword = 'Set60';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salt = 'NaCl';

% Want to gather data on:
% Fraction of solid vs liquid in the final T_LB and T_UB
% and also collect T_UB, T_LB

% Salt directory
SaltDataDir = fullfile(ProjectDir,Salt);

% List all inner directories
d_jobs = dir(SaltDataDir);
d_jobs = d_jobs([d_jobs(:).isdir]);
d_jobs = d_jobs(~ismember({d_jobs(:).name},{'.','..'}));
d_jobs = d_jobs(strncmp(DataKeyword,{d_jobs.name},length(DataKeyword)));

Data.T_UBs = [];
Data.T_LBs = [];
Data.T_Xs = [];

Data.Y_UBs = [];
Data.Y_LBs = [];
Data.Y_Xs = [];

for job_idx = 1:length(d_jobs)
    JobName = d_jobs(job_idx).name;
    disp(repmat('*',1,30))
    disp(['Current Job:' JobName])
    disp(repmat('*',1,30))
    
    JobDir = fullfile(d_jobs(job_idx).folder,JobName);
    res_file = fullfile(JobDir,[JobName '_MPResults.mat']);
    
    if ~isfile(res_file)
        continue
    end
    
    T_dat = load(res_file).T_dat;
    
    T_LB = max(T_dat.T_Trace(logical(T_dat.Freeze_Trace)));
    
    T_UB = min(T_dat.T_Trace(logical(T_dat.Melt_Trace)));
    T_xs = T_dat.T_Trace((~logical(T_dat.Melt_Trace) & ~logical(T_dat.Freeze_Trace)));
    
    
    %% T_LB first
    T_LB_dir = fullfile(JobDir,['T_' num2str(T_LB,'%0.4f')]);
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(T_LB_dir, Settings.Salt, ...
        pyargs('SystemName',JobName,...
        'RefStructure','Liquid',...
        'CheckFullTrajectory',false,...
        'SaveTrajectory',false,...
        'SavePredictionsImage',false,...
        'InitialRefFrac',0.5,...
        'RefChangeThreshold',0.25,...
        'FileType','gro',...
        'T_Ref',T_LB,...
        'T',T_LB,...
        'Verbose',true));
    
    Data.T_LBs(end+1) = T_LB;
    Data.Y_LBs(end+1) = PyOut{4};
    
    %% T_UB second
    T_UB_dir = fullfile(JobDir,['T_' num2str(T_UB,'%0.4f')]);
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(T_UB_dir, Settings.Salt, ...
        pyargs('SystemName',JobName,...
        'RefStructure','Liquid',...
        'CheckFullTrajectory',false,...
        'SaveTrajectory',false,...
        'SavePredictionsImage',false,...
        'InitialRefFrac',0.5,...
        'RefChangeThreshold',0.25,...
        'FileType','gro',...
        'T_Ref',T_UB,...
        'T',T_UB,...
        'Verbose',true));
    
    Data.T_UBs(end+1) = T_UB;
    Data.Y_UBs(end+1) = PyOut{4};
    
    %% Now any T_x's last
    for Txidx = 1:length(T_xs)
        T_x = T_xs(Txidx);
        
        T_x_dir = fullfile(JobDir,['T_' num2str(T_x,'%0.4f')]);
        PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(T_x_dir, Settings.Salt, ...
            pyargs('SystemName',JobName,...
            'RefStructure','Liquid',...
            'CheckFullTrajectory',false,...
            'SaveTrajectory',false,...
            'SavePredictionsImage',false,...
            'InitialRefFrac',0.5,...
            'RefChangeThreshold',0.25,...
            'FileType','gro',...
            'T_Ref',T_x,...
            'T',T_x,...
            'Verbose',true));
        
        Data.T_Xs(end+1) = T_x;
        Data.Y_Xs(end+1) = PyOut{4};
    end
    
    disp(['Job Complete:' JobName])
    disp([num2str(100*job_idx/length(d_jobs),'%0.2f') ' % complete.'])
    disp(repmat('*',1,30))
end

save(SaveDataDir,'Data')
Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Melting_Point_Data.mat';
DataKeyword = 'Set61';
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

for job_idx = 1:length(d_jobs)
    JobName = d_jobs(job_idx).name;
    JobDir = fullfile(d_jobs(job_idx).folder,JobName);
    res_file = fullfile(JobDir,[JobName '_MPResults.mat']);
    
    if ~isfile(res_file)
        continue
    end
    
    T_dat = load(res_file).T_dat;
    
    T_LB = max(T_dat.T_Trace(logical(T_dat.Freeze_Trace)));
    T_UB = min(T_dat.T_Trace(logical(T_dat.Melt_Trace)));
    T_x = T_dat.T_Trace((~logical(T_dat.Melt_Trace) & ~logical(T_dat.Freeze_Trace)));
    
    T_LB_dir = fullfile(JobDir,['T_' num2str(T_LB,'%0.4f')]);
    
    
    
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(T_LB_dir, Settings.Salt, ...
        pyargs('SystemName',JobName,...
        'RefStructure','liquid',...
        'CheckFullTrajectory',true,...
        'SaveTrajectory',false,...
        'SavePredictionsImage',false,...
        'RefChangeThreshold',Settings.MeltFreezeThreshold,...
        'SlopeThreshold',Settings.SlopeThreshold,...
        'FileType',Settings.CoordType,...
        'T_Ref',T_dat.T_ref,...
        'T',T,...
        'TimePerFrame',10));
    
    
    
    
end
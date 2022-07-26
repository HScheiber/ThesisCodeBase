function [WorkDir,JobName,Full_Model_Name] = GetMDWorkdir(Settings)

% If this is part of a bayesopt calculation, return a temporary folder as
% the working directory
if isfield(Settings,'initial_opt_type')
    Hashobj.S = Settings.S;
    Hashobj.CR_Damp = Settings.CR_Damp;
    Hashobj.Model = Settings.Model;
    Hashobj.Theory = Settings.Theory;
    Hashobj.C6_Damp = Settings.C6_Damp;
    Hashobj.GAdjust_MM = Settings.GAdjust_MM;
    Hashobj.GAdjust_MX = Settings.GAdjust_MX;
    Hashobj.GAdjust_XX = Settings.GAdjust_XX;
    
    WorkDir = fullfile(pwd,DataHash(Hashobj));
    JobName = 'Test_MP';
    Full_Model_Name = Settings.Model;
    return
end

% Otherwise, choose a systematic name
if isempty(Settings.Model)
    Settings.S = Init_Scaling_Object;
    Settings.CR_Damp = Init_CRDamping_Object;
    Settings.C6_Damp = Init_C6Damping_Object;
    [Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;
    
    Full_Model_Name = ModelName(Settings);
else
    Full_Model_Name = [Settings.Theory '_Model_' Settings.Model];
end

% Establish ensemble name
if strcmpi(Settings.Thermostat,'no')
    if strcmpi(Settings.Barostat,'no')
        Ensemble = 'NVE';
    else
        Ensemble = 'NPH';
    end
else
    if strcmpi(Settings.Barostat,'no')
        Ensemble = 'NVT';
    else
        Ensemble = 'NPT';
    end
end

% Assign calculation directory
JobName = [Settings.JobID '_' StructureLabel(Settings.Structure) '_' Full_Model_Name '_' Ensemble];
WorkDir = fullfile(Settings.project,Settings.Project_Directory_Name,Settings.Salt,JobName);


end
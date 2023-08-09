function [Settings,ModelFound] = Load_Model_Params(Settings)

ModelFound = true;
if strcmp(Settings.Theory,'JC') && strcmp(Settings.Salt,'LiI') && strcmp(Settings.Model,'BF1b')
    damp_MM = true;
    Settings.Model = 'BF1';
else
    damp_MM = false;
end

[~,~,~,~,Model_Data_Loc] = find_home;

Model_data_Filename = fullfile(Model_Data_Loc,Settings.Salt,[Settings.Salt '_' Settings.Theory '_Model_' Settings.Model '_data.mat']);

if ~isfield(Settings,'S')
    Settings.S = Init_Scaling_Object;
    Settings.C6_Damp = Init_C6Damping_Object;
    Settings.CR_Damp = Init_CRDamping_Object;
    [Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;
end
    
% Special case
if isempty(Settings.Model)
    Settings.S = Init_Scaling_Object;
    Settings.C6_Damp = Init_C6Damping_Object;
    Settings.CR_Damp = Init_CRDamping_Object;
    [Settings.GAdjust_MX,Settings.GAdjust_MM,Settings.GAdjust_XX] = Init_GAdjust_Object;
    
    if strcmp(Settings.Theory,'TF')
        [~,~,~,Settings.S] = TF_Potential_Parameters(Settings);
    end
    
    return
elseif strcmp(Settings.Theory,'JC') && strcmp(Settings.Salt,'LiI') && strcmp(Settings.Model,'A')
    disp(['Model found. Loading ' Settings.Theory ' Model ' Settings.Model])
    Settings.CR_Damp.MM.r_d = 0.26; % liCl = 0.21, LiI = 0.26
    Settings.CR_Damp.MM.b  = 75;
    Settings.S.D.MM = 910;
    Settings.S.D.XX = 0.23;
    return
end

if isfile(Model_data_Filename)
    dat = load(Model_data_Filename);
    data = dat.full_data;

    if isfield(data,'secondary_result')
        params = data.full_opt_point;
    else
        [~,midx] = min(data.bayesopt_results.ObjectiveTrace);
        params = data.bayesopt_results.XTrace{midx,:};
    end
    
    Settings = data.Settings;
    Settings.S = Init_Scaling_Object;
    if ~isfield(Settings,'Comb_rule')
        Settings.Comb_rule = 'Lorentz-Berthelot';
    end
    Ptypes = bayesopt_params(Settings);
    Param_names = {Ptypes.Name};
    Param = table();
    for jdx = 1:numel(Param_names)
        Param.(Param_names{jdx}) = params(jdx);
    end
    Settings = Get_Scaling_Params(Settings,Param);
else
    warning('Model not found.')
    ModelFound = false;
    return
end

if damp_MM
    Settings.C6_Damp.MM = true;
end

% Perturb the potential with Gaussians
if isempty(Settings.Additional_GAdjust)
    Settings.GAdjust_MX = [0 0 1];  %[0 0 1];
    Settings.GAdjust_MM = [0 0 1];
    Settings.GAdjust_XX = [0 0 1];
else
    Settings.GAdjust_MX = [0 0 1];
    Settings.GAdjust_MM = [0 0 1];
    Settings.GAdjust_XX = [0 0 1];
    mx = 1;
    mm = 1;
    xx = 1;
    for idx = 1:length(Settings.Additional_GAdjust)
        int = [Settings.Additional_GAdjust{idx} '_' num2str(idx)];
        switch int
            case ['MM' '_' num2str(idx)]
                Settings.GAdjust_MM(mm,1) = Param(pidx + 1);
                Settings.GAdjust_MM(mm,2) = Param(pidx + 2);
                Settings.GAdjust_MM(mm,3) = Param(pidx + 3);
                mm = mm+1;
            case ['XX' '_' num2str(idx)]
                Settings.GAdjust_XX(xx,1) = Param(pidx + 1);
                Settings.GAdjust_XX(xx,2) = Param(pidx + 2);
                Settings.GAdjust_XX(xx,3) = Param(pidx + 3);
                xx = xx+1;
            case ['MX' '_' num2str(idx)]
                Settings.GAdjust_MX(mx,1) = Param(pidx + 1);
                Settings.GAdjust_MX(mx,2) = Param(pidx + 2);
                Settings.GAdjust_MX(mx,3) = Param(pidx + 3);
                mx = mx+1;
        end
        pidx = pidx + 3;
    end
end

% Additional functions
ints = {'MM' 'XX' 'MX'};
for idx = 1:length(ints)
    int = ints{idx};
    
    if Settings.Additional_Function.(int).N >= 0
        Settings.S.N.(int).Value = Settings.Additional_Function.(int).N;
        Settings.S.N.(int).Scale = Param(pidx + 1);
        pidx = pidx + 1;
    end
end

end
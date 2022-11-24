function Settings = Fit_WBK_to_BH(Settings,Model_BH)

[~,~,~,~,ML_results_dir] = find_home;

Salt = Settings.Salt;
Theory = 'BH';
Model = Model_BH;

% Find the fully optimized model
dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);

if ~isfile(dat_file)
    error(['No model results file found for: ' Salt ', ' Theory ', Model ' Model '.']);
end

% Load data
data = load(dat_file);
data = data.full_data;

if isfield(data,'secondary_result')
    optimvals = nan(1,length(data.secondary_result));
    for jdx = 1:length(data.secondary_result)
        optimvals(jdx) = [data.secondary_result(jdx).optimValues.fval];
    end
    params = data.full_opt_point;
else
    [~,midx] = min(data.bayesopt_results.ObjectiveTrace);
    params = data.bayesopt_results.XTrace{midx,:};
end

pSettings = data.Settings;
Ptypes = bayesopt_params(pSettings);
Param_names = {Ptypes.Name};

Param = table();
for jdx = 1:numel(Param_names)
    Param.(Param_names{jdx}) = params(jdx);
end
pSettings = Get_Scaling_Params(pSettings,Param);
Settings_WBK = pSettings;
Settings_WBK.S = Init_Scaling_Object;
Settings_WBK.SigmaEpsilon = true;

if pSettings.Additivity
    Param.r0_MX = (Param.r0_MM + Param.r0_XX)/2;
    Param.epsilon_MX = sqrt(Param.epsilon_MM * Param.epsilon_XX);
    Param.gamma_MM = Param.gamma_MX;
    Param.gamma_XX = Param.gamma_MX;
    Param_names{end+1} = 'r0_MX';
    Param_names{end+1} = 'epsilon_MX';
    Param_names{end+1} = 'gamma_MM';
    Param_names{end+1} = 'gamma_XX';
end

ints = {'MM' 'XX' 'MX'};
for kdx = 1:numel(ints)
    int = ints{kdx};
    r0_BH = Param.(['r0_' int]);
    gamma_BH = Param.(['gamma_' int]);
    if gamma_BH < 6
        epsilon_BH = -Param.(['epsilon_' int]);
    else
        epsilon_BH = Param.(['epsilon_' int]);
    end
    k = 1/(gamma_BH - 6);

    f_BH = @(r) 6*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) - epsilon_BH*gamma_BH*k*((r0_BH./r).^6);
    df_BH = @(r) -6*(gamma_BH/r0_BH)*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) + 6*epsilon_BH*gamma_BH*k*(r0_BH^6)./(r.^7);
    df_BH_neg = @(r) 6*(gamma_BH/r0_BH)*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) - 6*epsilon_BH*gamma_BH*k*(r0_BH^6)./(r.^7);

    % Find minima and inflection points
    if gamma_BH < 7 % find the well minimum
        sigma_WBK = fzero(df_BH,r0_BH+1);
        eps_WBK = abs(f_BH(sigma_WBK));
    else % In this domain, r0 and epsilon describe the well minima
        sigma_WBK = r0_BH;
        eps_WBK = abs(epsilon_BH);
    end

    % find the inflection points
    r_infl = fminsearchbnd(df_BH,0.99*sigma_WBK,0,sigma_WBK);
    r_infl_outer = fminsearchbnd(df_BH_neg,1.01*sigma_WBK,sigma_WBK,10);

    % Define the WBK function with unknown gamma
    f_WBK_gam = @(r,gamma_WBK) (2*eps_WBK./(1 - (3./(gamma_WBK+3)))).*((sigma_WBK^6)./((sigma_WBK^6) + (r.^6))).*((3./(gamma_WBK+3)).*exp(gamma_WBK.*(1 - (r./sigma_WBK))) - 1);

    % Domain to compare BK and WBK functions
    r_check = r_infl:0.01:r_infl_outer;

    % Define a new function that computes the squared difference between BK and WBK functions for a given input gamma
    diff_BH_WBK = @(gamma_WBK) sum((f_BH(r_check) - f_WBK_gam(r_check,gamma_WBK)).^2);

    % Optimize the value of WBK gamma
    gamma_WBK = fminsearchbnd(diff_BH_WBK,gamma_BH,0,100);
    Settings_WBK.S.S.(int) = sigma_WBK; % nm
    Settings_WBK.S.E.(int) = eps_WBK;   % kJ/mol
    Settings_WBK.S.G.(int) = gamma_WBK; % unitless
end

Settings.S = Settings_WBK.S;
Settings.SigmaEpsilon = true;

end
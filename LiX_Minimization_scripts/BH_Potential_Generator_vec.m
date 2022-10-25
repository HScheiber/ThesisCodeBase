function U = BH_Potential_Generator_vec(Settings)

%% Load BH Parameters
dat = load(fullfile(Settings.home,'data','BH_Default_Param.mat'));
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.M =  Settings.S.Q; % atomic
q.X = -Settings.S.Q; % atomic

if Settings.SigmaEpsilon
    %% Parameters of interest for vdW potential
    sigma.MX = Settings.S.S.MX; % nm
    sigma.MM = Settings.S.S.MM;
    sigma.XX = Settings.S.S.XX;

    epsilon.MX = Settings.S.E.MX; % kJ/mol
    epsilon.MM = Settings.S.E.MM;
    epsilon.XX = Settings.S.E.XX;

    gamma.MX = Settings.S.G.MX; % unitless
    gamma.MM = Settings.S.G.MM;
    gamma.XX = Settings.S.G.XX;

    %% Convert to B*exp - C6/r^6 form
    B.MX = (6*epsilon.MX./(gamma.MX - 6)).*exp(gamma.MX);
    B.MM = (6*epsilon.MM./(gamma.MM - 6)).*exp(gamma.MM);
    B.XX = (6*epsilon.XX./(gamma.XX - 6)).*exp(gamma.XX);

    C.MX = (epsilon.MX.*gamma.MX.*(sigma.MX.^6))./(gamma.MX - 6);
    C.MM = (epsilon.MM.*gamma.MM.*(sigma.MM.^6))./(gamma.MM - 6);
    C.XX = (epsilon.XX.*gamma.XX.*(sigma.XX.^6))./(gamma.XX - 6);

    alpha.MX = gamma.MX./sigma.MX; % nm^(-1)
    alpha.MM = gamma.MM./sigma.MM; % nm^(-1)
    alpha.XX = gamma.XX./sigma.XX; % nm^(-1)
else
    %% B*exp - C6/r^6 form
    B.MX = Settings.S.R.All.*Settings.S.R.MX;
    B.MM = Settings.S.R.All.*Settings.S.R.MM;
    B.XX = Settings.S.R.All.*Settings.S.R.XX;

    C.MX = Settings.S.D.All.*Settings.S.D.MX;
    C.MM = Settings.S.D.All.*Settings.S.D.MM;
    C.XX = Settings.S.D.All.*Settings.S.D.XX;

    alpha.MX = Settings.S.A.All.*Settings.S.A.MX; % nm^(-1)
    alpha.MM = Settings.S.A.All.*Settings.S.A.MM; % nm^(-1)
    alpha.XX = Settings.S.A.All.*Settings.S.A.XX; % nm^(-1)
end

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    % Build PES
    U.(int) = QQ_prefactor.*q.(int(1)).*q.(int(2))./U.r + B.(int).*exp(-alpha.(int).*U.r) - C.(int)./(U.r.^6);
end

end
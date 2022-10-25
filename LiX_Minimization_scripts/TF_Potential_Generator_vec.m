function U = TF_Potential_Generator_vec(Settings)

%% Load TF Parameters
dat = load(fullfile(Settings.home,'data','TF_Default_Param.mat'));
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.M =  Settings.S.Q; % atomic
q.X  = -Settings.S.Q; % atomic

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

if Settings.SigmaEpsilon
    %% Parameters of interest for vdW potential
    sigma.MX = Settings.S.S.MX; % (also called R0) nm
    sigma.MM = Settings.S.S.MM;
    sigma.XX = Settings.S.S.XX;

    gamma.MX = Settings.S.G.MX; % unitless
    gamma.MM = Settings.S.G.MM;
    gamma.XX = Settings.S.G.XX;

    epsilon.MX = Settings.S.E.MX; % kJ/mol
    epsilon.MM = Settings.S.E.MM;
    epsilon.XX = Settings.S.E.XX;

    %% Convert to B*exp(-alpha*r) - C6/r^6 - C8/r^8 form
    kappa.MX = 1./(2 - (gamma.MX./6) - (gamma.MX./8));
    kappa.MM = 1./(2 - (gamma.MM./6) - (gamma.MM./8));
    kappa.XX = 1./(2 - (gamma.XX./6) - (gamma.XX./8));

    alpha.MX = gamma.MX./sigma.MX; % nm^(-1)
    alpha.MM = gamma.MM./sigma.MM; % nm^(-1)
    alpha.XX = gamma.XX./sigma.XX; % nm^(-1)

    B.MX = 2.*epsilon.MX.*kappa.MX.*exp(gamma.MX);
    B.MM = 2.*epsilon.MM.*kappa.MM.*exp(gamma.MM);
    B.XX = 2.*epsilon.XX.*kappa.XX.*exp(gamma.XX);

    C.MX = epsilon.MX.*gamma.MX.*kappa.MX.*(sigma.MX.^6)./6;
    C.MM = epsilon.MM.*gamma.MM.*kappa.MM.*(sigma.MM.^6)./6;
    C.XX = epsilon.XX.*gamma.XX.*kappa.XX.*(sigma.XX.^6)./6;

    D.MX = epsilon.MX.*gamma.MX.*kappa.MX.*(sigma.MX.^8)./8;
    D.MM = epsilon.MM.*gamma.MM.*kappa.MM.*(sigma.MM.^8)./8;
    D.XX = epsilon.XX.*gamma.XX.*kappa.XX.*(sigma.XX.^8)./8;
else
    
    alpha.MX = Settings.S.A.All.*Settings.S.A.MX; % nm^(-1)
    alpha.MM = Settings.S.A.All.*Settings.S.A.MM; % nm^(-1)
    alpha.XX = Settings.S.A.All.*Settings.S.A.XX; % nm^(-1)

    B.MX = Settings.S.R.All.*Settings.S.R.MX;
    B.MM = Settings.S.R.All.*Settings.S.R.MM;
    B.XX = Settings.S.R.All.*Settings.S.R.XX;

    C.MX = Settings.S.D6D.All.*Settings.S.D6D.MX.*Settings.S.D.All.*Settings.S.D.MX;
    C.MM = Settings.S.D6D.All.*Settings.S.D6D.MM.*Settings.S.D.All.*Settings.S.D.MM;
    C.XX = Settings.S.D6D.All.*Settings.S.D6D.XX.*Settings.S.D.All.*Settings.S.D.XX;

    D.MX = Settings.S.D8D.All.*Settings.S.D8D.MX.*Settings.S.D.All.*Settings.S.D.MX;
    D.MM = Settings.S.D8D.All.*Settings.S.D8D.MM.*Settings.S.D.All.*Settings.S.D.MM;
    D.XX = Settings.S.D8D.All.*Settings.S.D8D.XX.*Settings.S.D.All.*Settings.S.D.XX;
end

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    %% Build PES
    
    U.(int) = QQ_prefactor.*q.(int(1)).*q.(int(2))./U.r ...
        + B.(int).*exp(-alpha.(int).*U.r) ...
        - C.(Settings.Salt).(int)./(U.r.^6) ...
        - D.(Settings.Salt).(int)./(U.r.^8);
end

end
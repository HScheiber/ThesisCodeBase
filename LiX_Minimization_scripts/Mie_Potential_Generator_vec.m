function U = Mie_Potential_Generator_vec(Settings)

%% Load JC Parameters
dat = load(fullfile(Settings.home,'data','JC_Default_Param.mat'));
Param = dat.Param;
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.(Settings.Metal) =  Settings.S.Q; % atomic
q.(Settings.Halide)= -Settings.S.Q; % atomic

%% n repulsive exponent parameter
n.MM = Settings.S.n.MM;
n.XX = Settings.S.n.XX;
n.MX = Settings.S.n.MX;

%% Calculate prefactor
P_MM = (n.MM./(n.MM - 6)).*((n.MM./6).^(6./(n.MM - 6)));
P_XX = (n.XX./(n.XX - 6)).*((n.XX./6).^(6./(n.XX - 6)));
P_MX = (n.MX./(n.MX - 6)).*((n.MX./6).^(6./(n.MX - 6)));

%% Calculate parameters of interest for LJ potential
sigma_MM = Settings.S.S.All.*Settings.S.S.MM.*Param.(Settings.Metal).sigma;
sigma_XX = Settings.S.S.All.*Settings.S.S.XX.*Param.(Settings.Halide).sigma;
sigma_MX = Settings.S.S.All.*Settings.S.S.MX.*( Param.(Settings.Metal).sigma + Param.(Settings.Halide).sigma )./2;

epsilon_MM = Settings.S.E.All.*Settings.S.E.MM.*Param.(Settings.Metal).epsilon;
epsilon_XX = Settings.S.E.All.*Settings.S.E.XX.*Param.(Settings.Halide).epsilon;
epsilon_MX = Settings.S.E.All.*Settings.S.E.MX.*sqrt(Param.(Settings.Metal).epsilon.*Param.(Settings.Halide).epsilon);

% Change parameteters into A/r^n - C/r^6 format
A.MM = Settings.S.R.All.*Settings.S.R.MM.*P_MM.*epsilon_MM.*(sigma_MM.^n.MM);
C.MM = Settings.S.D.All.*Settings.S.D.MM.*P_MM.*epsilon_MM.*(sigma_MM.^6);

A.XX = Settings.S.R.All.*Settings.S.R.XX.*P_XX.*epsilon_XX.*(sigma_XX.^n.XX);
C.XX = Settings.S.D.All.*Settings.S.D.XX.*P_XX.*epsilon_XX.*(sigma_XX.^6);

A.MX = Settings.S.R.All.*Settings.S.R.MX.*P_MX.*epsilon_MX.*(sigma_MX.^n.MX);
C.MX = Settings.S.D.All.*Settings.S.D.MX.*P_MX.*epsilon_MX.*(sigma_MX.^6);

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    switch int
        case 'MX'
            Y1 = Settings.Metal;
            Y2 = Settings.Halide;
        case 'MM'
            Y1 = Settings.Metal;
            Y2 = Settings.Metal;
        case 'XX'
            Y1 = Settings.Halide;
            Y2 = Settings.Halide;
    end
    
    % Build PES
    U.(int) = QQ_prefactor.*q.(Y1).*q.(Y2)./U.r + A.(int)./(U.r.^n.(int)) - C.(int)./(U.r.^6);
end

end
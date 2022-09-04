% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETERS AND kJ/mol


function U = TF_Potential_Generator_vec(Settings)


%% Load TF Parameters
dat = load(fullfile(Settings.home,'data','TF_Default_Param.mat'));
C = dat.C; % Huggins-Mayer Dipole-Dipole Dispersion Parameter C: MX = +-   MM = ++     XX = --
D = dat.D; % Huggins-Mayer Dipole-Quadrupole Dispersion Parameter D: MX = +-   MM = ++     XX = --
sigma = dat.sigma; % TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
valence = dat.valence; % TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
rho = dat.rho; % TF Hardness Parameter Rho
b = dat.b; % Huggins-Mayer potential parameter b (same for all salts)
QQ_prefactor = dat.QQ_prefactor;


%% TF Parameter: q (charge)
q.Li =  Settings.S.Q; % atomic
q.Na =  Settings.S.Q; % atomic
q.K  =  Settings.S.Q; % atomic
q.Rb =  Settings.S.Q; % atomic
q.Cs =  Settings.S.Q; % atomic

q.F  = -Settings.S.Q; % atomic
q.Cl = -Settings.S.Q; % atomic
q.Br = -Settings.S.Q; % atomic
q.I  = -Settings.S.Q; % atomic

%% Generate range (r) in nm
U.r = Settings.Table_StepSize:Settings.Table_StepSize:Settings.Table_Length;

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2/valence.(Settings.Metal); % Unitless
beta.MX = 1 + 1/valence.(Settings.Metal) - 1/valence.(Settings.Halide); % Unitless
beta.XX = 1 - 2/valence.(Settings.Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All.*Settings.S.A.MM./rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All.*Settings.S.A.MX./rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All.*Settings.S.A.XX./rho.(Settings.Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All.*Settings.S.R.MM.*beta.MM.*b.*exp(2*sigma.(Settings.Metal)./rho.(Settings.Salt));
B.MX = Settings.S.R.All.*Settings.S.R.MX.*beta.MX.*b.*exp((sigma.(Settings.Metal) + sigma.(Settings.Halide))./rho.(Settings.Salt));
B.XX = Settings.S.R.All.*Settings.S.R.XX.*beta.XX.*b.*exp(2*sigma.(Settings.Halide)./rho.(Settings.Salt));

%% Scale Dispersion
C.(Settings.Salt).MM = Settings.S.D6D.All.*Settings.S.D6D.MM.*Settings.S.D.All.*Settings.S.D.MM.*C.(Settings.Salt).MM;
D.(Settings.Salt).MM = Settings.S.D8D.All.*Settings.S.D8D.MM.*Settings.S.D.All.*Settings.S.D.MM.*D.(Settings.Salt).MM;

C.(Settings.Salt).XX = Settings.S.D6D.All.*Settings.S.D6D.XX.*Settings.S.D.All.*Settings.S.D.XX.*C.(Settings.Salt).XX;
D.(Settings.Salt).XX = Settings.S.D8D.All.*Settings.S.D8D.XX.*Settings.S.D.All.*Settings.S.D.XX.*D.(Settings.Salt).XX;

C.(Settings.Salt).MX = Settings.S.D6D.All.*Settings.S.D6D.MX.*Settings.S.D.All.*Settings.S.D.MX.*C.(Settings.Salt).MX;
D.(Settings.Salt).MX = Settings.S.D8D.All.*Settings.S.D8D.MX.*Settings.S.D.All.*Settings.S.D.MX.*D.(Settings.Salt).MX;

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
    
    %% Build PES
    
    U.(int) = QQ_prefactor.*q.(Y1).*q.(Y2)./U.r ...
        + B.(int).*exp(-alpha.(int).*U.r) ...
        - C.(Settings.Salt).(int)./(U.r.^6) ...
        - D.(Settings.Salt).(int)./(U.r.^8);
end

end
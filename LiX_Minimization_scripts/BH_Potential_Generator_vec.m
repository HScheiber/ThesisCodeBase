% Coulomb-Buckingham model parameters with
% Lorenz-Berthelot combining rules
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETER AND
% kJ/mol

% C6Damping:
% 0 = no (default) damping. This is default of JC model.
% 1 = BJ/rational damping (same as in D3(BJ))
% 2 = Tang Damping (Essentially removes dispersion in JC)
% 3 = MMDRE Damping function (Very weak damping, damps mainly at mid range)
% 4 = PAMoC Damping function (fairly weak damping, damps mainly at mid range)
% 5 = EHFSK Damping function (strong damping)
% 6 = WY damping function (strongest damping)

% GAdjust are N x 3 arrays of gaussian parameters
% (i , 1) is the Gaussian height of the ith adjustment in kj/mol (may be negative or positive)
% (i , 2) is the center point of the ith Gaussian in nm (should be positive)
% (i , 3) is the standard deviation or width in nm (negative and positive values are the same)
% are the same)

% CRDamping = Close Range Damping
function U = BH_Potential_Generator_vec(Settings)

%% Load BH Parameters
dat = load(fullfile(Settings.home,'data','BH_Default_Param.mat'));
Param = dat.Param;
b = dat.b;
rho = dat.rho;
sigma = dat.sigma;
valence = dat.valence;
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.(Settings.Metal) =  Settings.S.Q; % atomic
q.(Settings.Halide)= -Settings.S.Q; % atomic

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2/valence.(Settings.Metal); % Unitless
beta.XX = 1 - 2/valence.(Settings.Halide); % Unitless
beta.MX = sqrt(beta.MM*beta.XX); % Unitless

%% Calculate Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All.*Settings.S.A.MM./rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All.*Settings.S.A.MX./rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All.*Settings.S.A.XX./rho.(Settings.Salt); % nm^-1

%% Calculate Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All.*Settings.S.R.MM.*beta.MM.*b.*exp(2.*sigma.(Settings.Metal)./rho.(Settings.Salt));
B.XX = Settings.S.R.All.*Settings.S.R.XX.*beta.XX.*b.*exp(2.*sigma.(Settings.Halide)./rho.(Settings.Salt));
B.MX = Settings.S.R.All.*Settings.S.R.MX.*beta.MX.*b.*exp((sigma.(Settings.Metal) + sigma.(Settings.Halide))./rho.(Settings.Salt));

%% Calculate parameters of interest for LJ potential: change parameteters into C6/r6 format and apply mixing rules
CMM_pre = 4*Param.(Settings.Metal).epsilon*(Param.(Settings.Metal).sigma^6);
CXX_pre = 4*Param.(Settings.Halide).epsilon*(Param.(Settings.Halide).sigma^6);
C.MM = Settings.S.D.All.*Settings.S.D.MM.*CMM_pre;
C.XX = Settings.S.D.All.*Settings.S.D.XX.*CXX_pre;
C.MX = Settings.S.D.All.*Settings.S.D.MX.*sqrt(CMM_pre.*CXX_pre);

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
    U.(int) = QQ_prefactor.*q.(Y1).*q.(Y2)./U.r + B.(int).*exp(-alpha.(int).*U.r) - C.(int)./(U.r.^6);
end

end
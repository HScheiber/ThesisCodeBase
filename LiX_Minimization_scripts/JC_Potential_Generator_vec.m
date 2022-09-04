% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', or 'TIP4PEW'
% Generates pair potential energy surfaces from 'Startpoint' up to a given
% length input 'Endpoint'. Plotswitch is a logical true or false to
% determine whether to plot the resulting potentials
% Recommended spacing is 0.005 angstroms or 0.0005 nm
% INPUT UNITS MUST BE ALL IN NANOMETERS. OUTPUTS ARE IN NANOMETER AND
% kJ/mol

% CRDamping = Close Range Damping
function U = JC_Potential_Generator_vec(Settings)


%% Load JC Parameters
dat = load(fullfile(Settings.home,'data','JC_Default_Param.mat'));
Param = dat.Param;
QQ_prefactor = dat.QQ_prefactor;

%% Parameter: q (charge)
q.(Settings.Metal) =  Settings.S.Q; % atomic
q.(Settings.Halide)= -Settings.S.Q; % atomic

%% Calculate parameters of interest for LJ potential
sigma_MM = Settings.S.S.All.*Settings.S.S.MM.*Param.(Settings.Metal).sigma;
sigma_XX = Settings.S.S.All.*Settings.S.S.XX.*Param.(Settings.Halide).sigma;
sigma_MX = Settings.S.S.All.*Settings.S.S.MX.*( Param.(Settings.Metal).sigma + Param.(Settings.Halide).sigma )./2;

epsilon_MM = Settings.S.E.All.*Settings.S.E.MM.*Param.(Settings.Metal).epsilon;
epsilon_XX = Settings.S.E.All.*Settings.S.E.XX.*Param.(Settings.Halide).epsilon;
epsilon_MX = Settings.S.E.All.*Settings.S.E.MX.*sqrt(Param.(Settings.Metal).epsilon.*Param.(Settings.Halide).epsilon);

% Change parameteters into A/r12 - C/r6 format
A.MM = Settings.S.R.All.*Settings.S.R.MM.*4.*epsilon_MM.*(sigma_MM.^12);
C.MM = Settings.S.D.All.*Settings.S.D.MM.*4.*epsilon_MM.*(sigma_MM.^6);

A.XX = Settings.S.R.All.*Settings.S.R.XX.*4.*epsilon_XX.*(sigma_XX.^12);
C.XX = Settings.S.D.All.*Settings.S.D.XX.*4.*epsilon_XX.*(sigma_XX.^6);

A.MX = Settings.S.R.All.*Settings.S.R.MX.*4.*epsilon_MX.*(sigma_MX.^12);
C.MX = Settings.S.D.All.*Settings.S.D.MX.*4.*epsilon_MX.*(sigma_MX.^6);


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
    U.(int) = QQ_prefactor.*q.(Y1).*q.(Y2)./U.r + A.(int)./(U.r.^12) - C.(int)./(U.r.^6);
end

end
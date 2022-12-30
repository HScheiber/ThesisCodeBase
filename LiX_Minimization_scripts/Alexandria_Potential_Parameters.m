% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', 'TIP4PEW', or 'SD' (for Smith-Dang NaCl only)
% Output units: Sigma is nanometers, epsilon is kJ/mol
% DS is the dispersion scaling factor (should only affect the r6 term)
% ES is the epsilon scaling factor (increases well depth)
function [Settings,AltParam] = Alexandria_Potential_Parameters(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'Alexandria_Potential_Parameters';
addOptional(p,'vdW_Type','WBK');
addOptional(p,'Coulomb_Only',false);
parse(p,varargin{:});
vdW_Type = p.Results.vdW_Type;
Coulomb_Only = p.Results.Coulomb_Only;
% Allowed vdW types: 'WBK', 'BK', 'LJ_12-6', 'LJ_8-6'

% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
AltParam = struct();

% Switch vdW_Type to Alexandria labels if Vancouver labels used
switch lower(vdW_Type)
    case 'bh'
        vdW_Type = 'BK';
    case 'bf'
        vdW_Type = 'WBK';
    case 'jc'
        vdW_Type = 'LJ_12-6';
    case 'mie'
        vdW_Type = 'LJ_8-6';
end

if ~Coulomb_Only
    switch lower(vdW_Type)
        case 'wbk'
            Settings.SigmaEpsilon = true;

            % Sigma in nm
            sigma.Li = 0.3043;
            sigma.Na = 0.3426;
            sigma.K  = 0.4170;
            sigma.Rb = 0.4478;
            sigma.Cs = 0.4752;
            sigma.F  = 0.3904;
            sigma.Cl = 0.4519;
            sigma.Br = 0.4730;
            sigma.I  = 0.4962;

            % epsilon in kJ/mol
            epsilon.Li = 0.03785;
            epsilon.Na = 0.10013;
            epsilon.K  = 0.31051;
            epsilon.Rb = 0.38827;
            epsilon.Cs = 0.34274;
            epsilon.F  = 0.03312;
            epsilon.Cl = 0.28844;
            epsilon.Br = 0.45204;
            epsilon.I  = 0.89077;

            % gamma (unitless)
            gamma.Li = 16.696;
            gamma.Na = 18.348;
            gamma.K  = 16.896;
            gamma.Rb = 16.833;
            gamma.Cs = 18.629;
            gamma.F  = 14.942;
            gamma.Cl = 15.073;
            gamma.Br = 15.101;
            gamma.I  = 15.104;

            Settings.S.G.MM = gamma.(Metal);
            Settings.S.G.XX = gamma.(Halide);
            Settings.S.G.MX = (gamma.(Metal) + gamma.(Halide))/2;

            Settings.S.E.MM = epsilon.(Metal);
            Settings.S.E.XX = epsilon.(Halide);
            Settings.S.E.MX = 2*epsilon.(Metal)*epsilon.(Halide)/(epsilon.(Metal) + epsilon.(Halide));

            Settings.S.S.MM = sigma.(Metal);
            Settings.S.S.XX = sigma.(Halide);
            Settings.S.S.MX = ( sqrt( ( epsilon.(Metal)*epsilon.(Halide)*gamma.(Metal)*gamma.(Halide)*(sigma.(Metal)*sigma.(Halide))^6 )...
                /((gamma.(Metal) - 6)*(gamma.(Halide) - 6)) )*(Settings.S.G.MX - 6)/(Settings.S.E.MX*Settings.S.G.MX) )^(1/6);

        case 'bk'
            Settings.SigmaEpsilon = false;

            % Exponential prefactor in kJ/mol
            A.Li = 1.2508e11;
            A.Na = 6.8785e13;
            A.K  = 1.8778e10;
            A.Rb = 1.8478e9;
            A.Cs = 1.4744e8;
            A.F  = 1.2372e4;
            A.Cl = 2.7188e4;
            A.Br = 4.4101e4;
            A.I  = 1.1228e5;

            % Exponent (alpha) in nm^-1;
            B.Li = 141.74;
            B.Na = 133.25;
            B.K  = 71.24;
            B.Rb = 57.79;
            B.Cs = 44.38;
            B.F  = 25.61;
            B.Cl = 21.50;
            B.Br = 21.35;
            B.I  = 21.93;

            % C6 dispersion (kj/mol nm^6)
            C.Li = 6.4257e-4;
            C.Na = 1.6807e-3;
            C.K  = 9.7224e-3;
            C.Rb = 1.5399e-2;
            C.Cs = 3.9186e-2;
            C.F  = 1.4416e-3;
            C.Cl = 1.0589e-2;
            C.Br = 1.7476e-2;
            C.I  = 3.0169e-2;

            % Alpha (nm^-1)
            Settings.S.A.MM = B.(Metal);
            Settings.S.A.XX = B.(Halide);
            Settings.S.A.MX = 2*B.(Metal)*B.(Halide)/(B.(Metal) + B.(Halide));

            % Repulsive prefactor (kj/mol)
            Settings.S.R.MM = A.(Metal);
            Settings.S.R.XX = A.(Halide);
            Settings.S.R.MX = (1/2)*( A.(Metal)*(A.(Metal)*B.(Metal)/(A.(Halide)*B.(Halide)))^(-B.(Metal)/(B.(Metal) + B.(Halide))) ...
                +  A.(Halide)*(A.(Halide)*B.(Halide)/(A.(Metal)*B.(Metal)))^(-B.(Halide)/(B.(Metal) + B.(Halide))) );

            % C6 parameter
            Settings.S.D6D.MM = C.(Metal);
            Settings.S.D6D.XX = C.(Halide);
            Settings.S.D6D.MX = sqrt(C.(Metal)*C.(Halide));

            % Gamma
            AltParam.G.MM = -7*lambertw((-1/7)*(6*Settings.S.D6D.MM*(Settings.S.A.MM^6)/Settings.S.R.MM)^(1/7));
            AltParam.G.XX = -7*lambertw((-1/7)*(6*Settings.S.D6D.XX*(Settings.S.A.XX^6)/Settings.S.R.XX)^(1/7));
            AltParam.G.MX = -7*lambertw((-1/7)*(6*Settings.S.D6D.MX*(Settings.S.A.MX^6)/Settings.S.R.MX)^(1/7));

            % Sigma (nm)
            AltParam.S.MM = AltParam.G.MM/Settings.S.A.MM;
            AltParam.S.XX = AltParam.G.XX/Settings.S.A.XX;
            AltParam.S.MX = AltParam.G.MX/Settings.S.A.MX;

            % Epsilon (kJ/mol)
            AltParam.E.MM = Settings.S.D6D.MM*(AltParam.G.MM - 6)/(AltParam.G.MM*(AltParam.S.MM^6));
            AltParam.E.XX = Settings.S.D6D.XX*(AltParam.G.XX - 6)/(AltParam.G.XX*(AltParam.S.XX^6));
            AltParam.E.MX = Settings.S.D6D.MX*(AltParam.G.MX - 6)/(AltParam.G.MX*(AltParam.S.MX^6));

        case 'lj_12-6'
            Settings.SigmaEpsilon = true;

            % Sigma in nm
            sigma.Li = 0.2157;
            sigma.Na = 0.2873;
            sigma.K  = 0.4075;
            sigma.Rb = 0.4608;
            sigma.Cs = 0.5161;
            sigma.F  = 0.3365;
            sigma.Cl = 0.4655;
            sigma.Br = 0.4868;
            sigma.I  = 0.5052;

            % epsilon in kJ/mol
            epsilon.Li = 0.01640;
            epsilon.Na = 0.05477;
            epsilon.K  = 0.02863;
            epsilon.Rb = 0.01975;
            epsilon.Cs = 0.01472;
            epsilon.F  = 0.07646;
            epsilon.Cl = 0.05707;
            epsilon.Br = 0.08802;
            epsilon.I  = 0.20108;

            % Exponent
            Settings.S.n.MM = 12;
            Settings.S.n.XX = 12;
            Settings.S.n.MX = 12;

            % Epsilon in kJ/mol
            Settings.S.E.MM = epsilon.(Metal);
            Settings.S.E.XX = epsilon.(Halide);
            Settings.S.E.MX = sqrt(epsilon.(Metal)*epsilon.(Halide));

            % Sigma in nm
            Settings.S.S.MM = sigma.(Metal);
            Settings.S.S.XX = sigma.(Halide);
            Settings.S.S.MX = (sigma.(Metal) + sigma.(Halide))/2;
        case 'lj_8-6'
            Settings.SigmaEpsilon = true;

            % Sigma in nm
            sigma.Li = 0.2142;
            sigma.Na = 0.2812;
            sigma.K  = 0.3829;
            sigma.Rb = 0.4169;
            sigma.Cs = 0.4830;
            sigma.F  = 0.2739;
            sigma.Cl = 0.3789;
            sigma.Br = 0.3988;
            sigma.I  = 0.4328;

            % epsilon in kJ/mol
            epsilon.Li = 0.22306;
            epsilon.Na = 0.46614;
            epsilon.K  = 0.28882;
            epsilon.Rb = 0.29529;
            epsilon.Cs = 0.12771;
            epsilon.F  = 1.30440;
            epsilon.Cl = 1.37999;
            epsilon.Br = 0.08802;
            epsilon.I  = 2.95471;

            % Exponent
            Settings.S.n.MM = 8;
            Settings.S.n.XX = 8;
            Settings.S.n.MX = 8;

            % Epsilon in kJ/mol
            Settings.S.E.MM = epsilon.(Metal);
            Settings.S.E.XX = epsilon.(Halide);
            Settings.S.E.MX = sqrt(epsilon.(Metal)*epsilon.(Halide));

            % Sigma in nm
            Settings.S.S.MM = sigma.(Metal);
            Settings.S.S.XX = sigma.(Halide);
            Settings.S.S.MX = (sigma.(Metal) + sigma.(Halide))/2;
        otherwise
            error(['Unknown vdW type for Alexandria potential: ' vdW_Type])
    end
end 
%% Coulomb parameters
% Q_core in electron charge units
Q_core.Li = 3.391;
Q_core.Na = 4.007;
Q_core.K  = 5.464;
Q_core.Rb = 6.402;
Q_core.Cs = 8.542;
Q_core.F  = 5.590;
Q_core.Cl = 5.713;
Q_core.Br = 5.778;
Q_core.I  = 7.613;

% Charge distribution parameter (nm^-1)
beta.Li = 13.772;
beta.Na = 12.609;
beta.K  = 12.412;
beta.Rb = 12.360;
beta.Cs = 11.880;
beta.F  = 9.139;
beta.Cl = 6.967;
beta.Br = 6.465;
beta.I  = 5.745;

% Polarizability parameter (nm^3)
alpha.Li = 0.029e-3;
alpha.Na = 0.18e-3;
alpha.K  = 0.81e-3;
alpha.Rb = 1.32e-3;
alpha.Cs = 2.02e-3;
alpha.F  = 1.3e-3;
alpha.Cl = 3.5e-3;
alpha.Br = 4.6e-3;
alpha.I  = 7.5e-3;

% Charge distribution parameters
Settings.S.beta.MM = beta.(Metal)^2./sqrt(beta.(Metal)^2 + beta.(Metal)^2);
Settings.S.beta.XX = beta.(Halide)^2./sqrt(beta.(Halide)^2 + beta.(Halide)^2);
Settings.S.beta.MX = beta.(Metal)*beta.(Halide)./sqrt(beta.(Metal)^2 + beta.(Halide)^2);

% Polarization parameters
Settings.S.PM = alpha.(Metal); % polarizibility parameter alpha for the metal ions in nm^3
Settings.S.PX = alpha.(Halide); % polarizibility parameter alpha for the halide ions in nm^3
Settings.S.QcoreM = Q_core.(Metal); % Sets metal core charge. Total charge is given by S.Q such that Q_core + Q_shell = S.Q for metals
Settings.S.QcoreX = Q_core.(Halide); % Sets halide core charge. Total charge is given by S.Q such that Q_core + Q_shell = -S.Q for halides

end
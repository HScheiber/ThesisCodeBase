% Output units: energies are kJ/mol and distances are nm
function [OutputMX,OutputMM,OutputXX] = BH_Potential_Parameters(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'BH_Potential_Parameters';
addOptional(p,'PlotType','full')
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

%% Conversion factors and fundamental constants
kj_per_erg = 1e-10; % kJ per erg
NA = 6.0221409e23; % Molecules per mole
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

%% JC Ion sigma/epsilon Parameters in SPC/E water
Param.Li.sigma = (0.791*nm_per_Ang)*(2^(5/6)); % nm
Param.Li.epsilon = (0.3367344)*kJ_per_kcal; % kJ/mol

Param.Na.sigma = (1.212*nm_per_Ang)*(2^(5/6)); % nm
Param.Na.epsilon = (0.3526418)*kJ_per_kcal; % kJ/mol

Param.K.sigma = (1.593*nm_per_Ang)*(2^(5/6)); % nm
Param.K.epsilon = (0.4297054)*kJ_per_kcal; % kJ/mol

Param.Rb.sigma = (1.737*nm_per_Ang)*(2^(5/6)); % nm
Param.Rb.epsilon = (0.4451036)*kJ_per_kcal; % kJ/mol

Param.Cs.sigma = (2.021*nm_per_Ang)*(2^(5/6)); % nm
Param.Cs.epsilon = (0.0898565)*kJ_per_kcal; % kJ/mol

Param.F.sigma = (2.257*nm_per_Ang)*(2^(5/6)); % nm
Param.F.epsilon = (0.0074005)*kJ_per_kcal; % kJ/mol

Param.Cl.sigma = (2.711*nm_per_Ang)*(2^(5/6)); % nm
Param.Cl.epsilon = (0.0127850)*kJ_per_kcal; % kJ/mol

Param.Br.sigma = (2.751*nm_per_Ang)*(2^(5/6)); % nm
Param.Br.epsilon = (0.0269586)*kJ_per_kcal; % kJ/mol

Param.I.sigma = (2.919*nm_per_Ang)*(2^(5/6)); % nm
Param.I.epsilon = (0.0427845)*kJ_per_kcal; % kJ/mol


%% Parameter: q (charge)
q.Li =  Settings.S.Q; % atomic
q.Na =  Settings.S.Q; % atomic
q.K  =  Settings.S.Q; % atomic
q.Rb =  Settings.S.Q; % atomic
q.Cs =  Settings.S.Q; % atomic

q.F  = -Settings.S.Q; % atomic
q.Cl = -Settings.S.Q; % atomic
q.Br = -Settings.S.Q; % atomic
q.I  = -Settings.S.Q; % atomic


%% TF Repulsive Size Parameter sigma (AKA r+/-): P = +   M = -
% Metals
sigma.Li = 0.816*nm_per_Ang; % nm
sigma.Na = 1.170*nm_per_Ang; % nm
sigma.K  = 1.463*nm_per_Ang; % nm
sigma.Rb = 1.587*nm_per_Ang; % nm
sigma.Cs = 1.720*nm_per_Ang; % nm

% Halides
sigma.F  = 1.179*nm_per_Ang; % nm
sigma.Cl = 1.585*nm_per_Ang; % nm
sigma.Br = 1.716*nm_per_Ang; % nm
sigma.I  = 1.907*nm_per_Ang; % nm

%% TF Parameter: Number of Valence electrons (for Pauling Coefficient Calculation)
% Metals
valence.Li = 2;
valence.Na = 8;
valence.K = 8;
valence.Rb = 8;
valence.Cs = 8;

% Halides
valence.F = 8;
valence.Cl = 8;
valence.Br = 8;
valence.I = 8;

%% TF Hardness Parameter Rho
rho.LiF = 0.299*nm_per_Ang; % nm
rho.LiCl = 0.342*nm_per_Ang; % nm
rho.LiBr = 0.353*nm_per_Ang; % nm
rho.LiI = 0.430*nm_per_Ang; % nm

rho.NaF = 0.330*nm_per_Ang; % nm
rho.NaCl = 0.317*nm_per_Ang; % nm
rho.NaBr = 0.340*nm_per_Ang; % nm
rho.NaI = 0.386*nm_per_Ang; % nm

rho.KF = 0.338*nm_per_Ang; % nm
rho.KCl = 0.337*nm_per_Ang; % nm
rho.KBr = 0.335*nm_per_Ang; % nm
rho.KI = 0.355*nm_per_Ang; % nm

rho.RbF = 0.328*nm_per_Ang; % nm
rho.RbCl = 0.318*nm_per_Ang; % nm
rho.RbBr = 0.335*nm_per_Ang; % nm
rho.RbI = 0.337*nm_per_Ang; % nm

rho.CsF = 0.282*nm_per_Ang; % nm
rho.CsCl = 0.272*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsBr = 0.289*nm_per_Ang; % nm THIS IS ESTIMATED
rho.CsI = 0.291*nm_per_Ang; % nm THIS IS ESTIMATED

%% Huggins-Mayer potential parameter b (same for all salts)
b = (0.338e-12)*kj_per_erg*NA; % kJ/mol

%% Calculate Pauling Coefficients beta: MX = +-   MM = ++     XX = --
beta.MM = 1 + 2*q.(Metal)/valence.(Metal); % Unitless
beta.XX = 1 + 2*q.(Halide)/valence.(Halide); % Unitless
beta.MX = 1 + 1/valence.(Metal) - 1/valence.(Halide); % Unitless

%% Calculate TF Repulsive Exponential Parameter alpha: MX = +-   MM = ++     XX = --
alpha.MM = Settings.S.A.All*Settings.S.A.MM/rho.(Settings.Salt); % nm^-1
alpha.MX = Settings.S.A.All*Settings.S.A.MX/rho.(Settings.Salt); % nm^-1
alpha.XX = Settings.S.A.All*Settings.S.A.XX/rho.(Settings.Salt); % nm^-1

%% Calculate TF Repulsive Scaling Parameter B: MX = +-   MM = ++     XX = -- (Including scaling)
B.MM = Settings.S.R.All*Settings.S.R.MM*beta.MM*b*exp(2*sigma.(Metal)/rho.(Settings.Salt));
B.XX = Settings.S.R.All*Settings.S.R.XX*beta.XX*b*exp(2*sigma.(Halide)/rho.(Settings.Salt));
B.MX = Settings.S.R.All*Settings.S.R.MX*beta.MX*b*exp((sigma.(Metal) + sigma.(Halide))/rho.(Settings.Salt));

%% Calculate parameters of interest for LJ potential
sigma_MM = Param.(Metal).sigma;
sigma_XX = Param.(Halide).sigma;

epsilon_MM = Param.(Metal).epsilon;
epsilon_XX = Param.(Halide).epsilon;

% Change parameteters into C6/r6 format and apply mixing rules
CMM_pre = 4*epsilon_MM*(sigma_MM^6);
CXX_pre = 4*epsilon_XX*(sigma_XX^6);
C.MM = Settings.S.D.All*Settings.S.D.MM*CMM_pre;
C.XX = Settings.S.D.All*Settings.S.D.XX*CXX_pre;
C.MX = Settings.S.D.All*Settings.S.D.MX*sqrt(CMM_pre*CXX_pre);

OutputMM.B = B.MM; % Repulsion parameter
OutputMM.C = C.MM; % C6 dispersion constant
OutputMM.alpha = alpha.MM; % Exponential steepness parameter

OutputXX.B = B.XX; % Repulsion parameter
OutputXX.C = C.XX; % C6 dispersion constant
OutputXX.alpha = alpha.XX; % Exponential steepness parameter

OutputMX.B = B.MX; % Repulsion parameter
OutputMX.C = C.MX; % C6 dispersion constant
OutputMX.alpha = alpha.MX; % Exponential steepness parameter

if Plotswitch
    figure;
    % Options
    lw=2;
    fs=25;
    
    %% Conversion Factors / Constants
    m_per_nm = 1e+9; % nm per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    J_per_kJ = 1000;
    epsilon_0 = (8.854187817620e-12)*J_per_kJ/(m_per_nm*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in C^-2 nm kJ mol^-1
    
    q.(Metal) = Settings.S.Q; % charge of M in C
    q.(Halide) = -Settings.S.Q; % charge of X in C
    
    r = Startpoint:Settings.Table_StepSize:Settings.Table_Length; % r vector
    
    U.MX.f = 1./r; % Electrostatics function f(r)
    U.MX.g = - (C.MX)./(r.^6); % Dispersion g(r)
    U.MX.h = B.MX*exp(-alpha.MX.*r); % Short range repulsion h(r)
    
    U.MX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MX.dg = - 6*(C.MX)./(r.^7); % Dispersion -g(r)
    U.MX.dh = alpha.MX*B.MX*exp(-alpha.MX.*r); % Short range repulsion -h(r)
    
    U.XX.f = 1./r; % Electrostatics function f(r)
    U.XX.g = - (C.XX)./(r.^6); % Dispersion g(r)
    U.XX.h = B.XX*exp(-alpha.XX.*r); % Short range repulsion h(r)
    
    U.XX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.XX.dg = - 6*(C.XX)./(r.^7); % Dispersion -g(r)
    U.XX.dh = alpha.XX*B.XX*exp(-alpha.XX.*r); % Short range repulsion -h(r)
    
    U.MM.f = 1./r; % Electrostatics function f(r)
    U.MM.g = - (C.MM)./(r.^6); % Dispersion g(r)
    U.MM.h = B.MM*exp(-alpha.MM.*r); % Short range repulsion h(r)
    
    U.MM.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MM.dg = - 6*(C.MM)./(r.^7); % Dispersion -g(r)
    U.MM.dh = alpha.MM*B.MM*exp(-alpha.MM.*r); % Short range repulsion -h(r)
    
    h = cell(1,3);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Halide).*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.(Metal)*q.(Metal).*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.(Halide)*q.(Halide).*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'lj'
            h{1} = plot(r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g + U.MM.h,'Color','g','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g + U.XX.h,'Color','b','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Repulsive Potential';
    end
    
    title(['Plot of ' ttxt ' for ' Settings.Salt ' BH Model'],...
       'Interpreter','latex','fontsize',fs)
    
    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    xlabel('Separation [\AA]','fontsize',fs,'Interpreter','latex');
    ylabel('Potential Energy [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');
    
    ylim(yl);
    xlim([Startpoint Settings.Table_Length]);
    
    % Blank line
    hline = refline([0 0]);
    hline.Color = 'k';
    hline.LineWidth = lw-1;
    hline.LineStyle = '--';
    leg1 = legend([h{:}],{[Metal '$^{+}$' ' - ' Halide '$^{-}$'] [Metal '$^{+}$' ' - ' Metal '$^{+}$'] [Halide '$^{-}$' ' - ' Halide '$^{-}$']});
    leg1.Interpreter = 'latex';
end

end
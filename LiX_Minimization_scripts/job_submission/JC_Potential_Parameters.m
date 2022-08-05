% Jung-Cheatham model parameters adapted for three water models with
% Lorenz-Berthelot combining rules
% Use Watermodel = 'SPC/E', 'TIP3P', 'TIP4PEW', or 'SD' (for Smith-Dang NaCl only)
% Output units: Sigma is nanometers, epsilon is kJ/mol
% DS is the dispersion scaling factor (should only affect the r6 term)
% ES is the epsilon scaling factor (increases well depth)
function [OutputMX,OutputMM,OutputXX] = JC_Potential_Parameters(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'JC_Potential_Parameters';
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

%% Conversion factors
nm_per_Ang = 0.1; % nm per Angstrom
kJ_per_kcal = 4.184; % kj per kcal

%% JC Ion Parameters in SPC/E water
if strcmp(Settings.WaterModel,'SPC/E')
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

elseif strcmp(Settings.WaterModel,'TIP3P')
    %% JC Ion Parameters in TIP3P water
    Param.Li.sigma = (1.025*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.0279896)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.369*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.0874393)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.705*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.1936829)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.813*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.3278219)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.976*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.4065394)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.303*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0033640)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.513*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0355910)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.608*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0586554)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.860*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0536816)*kJ_per_kcal; % kJ/mol

elseif strcmp(Settings.WaterModel,'TIP4PEW')
    %% JC Ion Parameters in TIP4P water
    Param.Li.sigma = (0.808*nm_per_Ang)*(2^(5/6)); % nm
    Param.Li.epsilon = (0.1039884)*kJ_per_kcal; % kJ/mol

    Param.Na.sigma = (1.226*nm_per_Ang)*(2^(5/6)); % nm
    Param.Na.epsilon = (0.1684375)*kJ_per_kcal; % kJ/mol

    Param.K.sigma = (1.590*nm_per_Ang)*(2^(5/6)); % nm
    Param.K.epsilon = (0.2794651)*kJ_per_kcal; % kJ/mol

    Param.Rb.sigma = (1.709*nm_per_Ang)*(2^(5/6)); % nm
    Param.Rb.epsilon = (0.4331494)*kJ_per_kcal; % kJ/mol

    Param.Cs.sigma = (1.888*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cs.epsilon = (0.3944318)*kJ_per_kcal; % kJ/mol

    Param.F.sigma = (2.538*nm_per_Ang)*(2^(5/6)); % nm
    Param.F.epsilon = (0.0015752)*kJ_per_kcal; % kJ/mol

    Param.Cl.sigma = (2.760*nm_per_Ang)*(2^(5/6)); % nm
    Param.Cl.epsilon = (0.0116615)*kJ_per_kcal; % kJ/mol

    Param.Br.sigma = (2.768*nm_per_Ang)*(2^(5/6)); % nm
    Param.Br.epsilon = (0.0303773)*kJ_per_kcal; % kJ/mol

    Param.I.sigma = (2.952*nm_per_Ang)*(2^(5/6)); % nm
    Param.I.epsilon = (0.0417082)*kJ_per_kcal; % kJ/mol
    
elseif strcmp(Settings.WaterModel,'SD') && strcmpi(Settings.Salt,'NaCl')
    Param.Na.sigma = (2.350*nm_per_Ang); % nm
    Param.Na.epsilon = (0.1300)*kJ_per_kcal; % kJ/mol
    
    Param.Cl.sigma = (4.400*nm_per_Ang); % nm
    Param.Cl.epsilon = (0.1000)*kJ_per_kcal; % kJ/mol
else
    error(['Unknown Water Model: "' Settings.WaterModel ...
        '". Please choose one of "SPC/E", "TIP3P", "TIP4PEW", or "SD" (with NaCl).'])
end

% Calculate cross terms using mixing rules on unscaled terms

Param.(Settings.Salt).sigma = (1/2)*(Param.(Metal).sigma + Param.(Halide).sigma);
Param.(Settings.Salt).epsilon = sqrt(Param.(Metal).epsilon*Param.(Halide).epsilon);

%% Apply scaling

% Generate the scaling factors
Gamma_All = Settings.S.E.All*(Settings.S.D.All^2)*(1/Settings.S.R.All); % Epsilon scaling factor
Beta_All = Settings.S.S.All*(1/(Settings.S.D.All^(1/6)))*(Settings.S.R.All^(1/6)); % Sigma scaling factor

Gamma_MM = Gamma_All*Settings.S.E.MM*(Settings.S.D.MM^2)*(1/Settings.S.R.MM);
Beta_MM = Beta_All*Settings.S.S.MM*(1/(Settings.S.D.MM^(1/6)))*(Settings.S.R.MM^(1/6));

Gamma_XX = Gamma_All*Settings.S.E.XX*(Settings.S.D.XX^2)*(1/Settings.S.R.XX);
Beta_XX = Beta_All*Settings.S.S.XX*(1/(Settings.S.D.XX^(1/6)))*(Settings.S.R.XX^(1/6));

Gamma_MX = Gamma_All*Settings.S.E.MX*(Settings.S.D.MX^2)*(1/Settings.S.R.MX);
Beta_MX = Beta_All*Settings.S.S.MX*(1/(Settings.S.D.MX^(1/6)))*(Settings.S.R.MX^(1/6));

Param.(Metal).sigma = Beta_MM*Param.(Metal).sigma;
Param.(Metal).epsilon = Gamma_MM*Param.(Metal).epsilon;

Param.(Halide).sigma = Beta_XX*Param.(Halide).sigma;
Param.(Halide).epsilon = Gamma_XX*Param.(Halide).epsilon;

Param.(Settings.Salt).sigma = Beta_MX*Param.(Settings.Salt).sigma;
Param.(Settings.Salt).epsilon = Gamma_MX*Param.(Settings.Salt).epsilon;

OutputMM = Param.(Metal);
OutputXX = Param.(Halide);
OutputMX = Param.(Settings.Salt);

if Plotswitch
    figure;
    % Options
    lw=3;
    fs=35;
    
    %% Conversion Factors / Constants
    m_per_nm = 1e+9; % nm per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    J_per_kJ = 1000;
    epsilon_0 = (8.854187817620e-12)*J_per_kJ/(m_per_nm*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in C^-2 nm kJ mol^-1

    % Use combining rules
    epsilon_MM = Param.(Metal).epsilon; % kJ/mol
    epsilon_XX = Param.(Halide).epsilon; % kJ/mol
    epsilon_MX = Param.(Settings.Salt).epsilon; % kJ/mol
    
    sigma_MM = Param.(Metal).sigma;
    sigma_XX = Param.(Halide).sigma;
    sigma_MX = Param.(Settings.Salt).sigma;

    q.(Metal) = Settings.S.Q; % charge of M in C
    q.(Halide) = -Settings.S.Q; % charge of X in C
    
    r = Startpoint:Settings.Table_StepSize:Settings.Table_Length; % r vector
    
    U.MX.f = 1./r; % Electrostatics function f(r)
    U.MX.g = - 4*epsilon_MX.*(( sigma_MX./r ).^6 ); % Dispersion g(r)
    U.MX.h = 4*epsilon_MX.*(( sigma_MX./r ).^12 ); % Short range repulsion h(r)
    
    U.MX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MX.dg = - 6*4*epsilon_MX.*(( sigma_MX./r ).^7 ); % Dispersion -g(r)
    U.MX.dh = 12*4*epsilon_MX.*(( sigma_MX./r ).^13 ); % Short range repulsion -h(r)
    
    U.XX.f = 1./r; % Electrostatics function f(r)
    U.XX.g = - 4*epsilon_XX.*(( sigma_XX./r ).^6 ); % Dispersion g(r)
    U.XX.h = 4*epsilon_XX.*(( sigma_XX./r ).^12 ); % Short range repulsion h(r)
    
    U.XX.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.XX.dg = - 6*4*epsilon_XX.*(( sigma_XX./r ).^7 ); % Dispersion -g(r)
    U.XX.dh = 12*4*epsilon_XX.*(( sigma_XX./r ).^13 ); % Short range repulsion -h(r)
    
    U.MM.f = 1./r; % Electrostatics function f(r)
    U.MM.g = - 4*epsilon_MM.*(( sigma_MM./r ).^6 ); % Dispersion g(r)
    U.MM.h = 4*epsilon_MM.*(( sigma_MM./r ).^12 ); % Short range repulsion h(r)
    
    U.MM.df = 1./(r.^2); % Derivative of Electrostatics function -f(r)
    U.MM.dg = - 6*4*epsilon_MM.*(( sigma_MM./r ).^7 ); % Dispersion -g(r)
    U.MM.dh = 12*4*epsilon_MM.*(( sigma_MM./r ).^13 ); % Short range repulsion -h(r)
    
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
    
    title(['Plot of ' ttxt ' for ' Settings.Salt ' JC Model'],...
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
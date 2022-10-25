% Output units: energies are kJ/mol and distances are nm
function [OutputMX,OutputMM,OutputXX] = BH_Potential_Parameters(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'BH_Potential_Parameters';
addOptional(p,'PlotType','full')
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'TightForm',false);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
TightForm = p.Results.TightForm;
% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

if isfield(Settings,'SigmaEpsilon')
    TightForm = ~Settings.SigmaEpsilon;
end

%% Parameter: q (charge)
q.M =  Settings.S.Q; % atomic
q.X  = -Settings.S.Q; % atomic

if TightForm % Input as tight form: B, C, D, alpha parameters
    alpha.MX = Settings.S.A.All*Settings.S.A.MX; % nm^(-1)
    alpha.MM = Settings.S.A.All*Settings.S.A.MM; % nm^(-1)
    alpha.XX = Settings.S.A.All*Settings.S.A.XX; % nm^(-1)
    
    B.MX = Settings.S.R.All*Settings.S.R.MX; % kj/mol
    B.MM = Settings.S.R.All*Settings.S.R.MM;
    B.XX = Settings.S.R.All*Settings.S.R.XX;
    
    C.MX = Settings.S.D6D.All*Settings.S.D6D.MX*Settings.S.D.All*Settings.S.D.MX; % kJ/mol nm^(-6)
    C.MM = Settings.S.D6D.All*Settings.S.D6D.MM*Settings.S.D.All*Settings.S.D.MM;
    C.XX = Settings.S.D6D.All*Settings.S.D6D.XX*Settings.S.D.All*Settings.S.D.XX;
else
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

    %% Convert to B*exp - C6/r^6 form
    B.MX = (6*epsilon.MX/(gamma.MX - 6))*exp(gamma.MX);
    B.MM = (6*epsilon.MM/(gamma.MM - 6))*exp(gamma.MM);
    B.XX = (6*epsilon.XX/(gamma.XX - 6))*exp(gamma.XX);

    C.MX = (epsilon.MX*gamma.MX*(sigma.MX^6))/(gamma.MX - 6);
    C.MM = (epsilon.MM*gamma.MM*(sigma.MM^6))/(gamma.MM - 6);
    C.XX = (epsilon.XX*gamma.XX*(sigma.XX^6))/(gamma.XX - 6);

    alpha.MX = gamma.MX/sigma.MX; % nm^(-1)
    alpha.MM = gamma.MM/sigma.MM; % nm^(-1)
    alpha.XX = gamma.XX/sigma.XX; % nm^(-1)
end

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
            h{1} = plot(r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
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
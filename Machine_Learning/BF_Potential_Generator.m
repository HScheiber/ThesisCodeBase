function [U,C6_out] = BF_Potential_Generator(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'BF_Potential_Generator';
addOptional(p,'PlotType','full');
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'MDP_Minimize',false);
addOptional(p,'Include_Dispersion_Scale',true);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
Startpoint = p.Results.Startpoint;
Plotswitch = p.Results.Plotswitch;
Incl_Disp = p.Results.Include_Dispersion_Scale;
% Allowed plot types: 'full', 'lj', 'full-derivative', 'lj-derivative',
% 'dispersion', 'dispersion-derivative', 'repulsive',
% 'repulsive-derivative'

if p.Results.MDP_Minimize
    MDP = 'MinMDP';
else
    MDP = 'MDP';
end

%% Conversion factors and fundamental constants
nm_per_m = 1e+9; % nm per m
NA = 6.0221409e23; % Molecules per mole
e_c = 1.60217662e-19; % Elementary charge in Coulombs
epsilon_0 = (8.854187817620e-12)*1000/(nm_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 nm^-1
k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ nm C^-2 mol^-1

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

%% Parameter: q (charge)
q.M =  Settings.S.Q; % atomic
q.X  = -Settings.S.Q; % atomic

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

%% Generate range (r) in nm
U.r = Startpoint:Settings.Table_StepSize:Settings.Table_Length;

%% Generate PES
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    Kappa = 2*epsilon.(int)/(1 - (3/(gamma.(int) + 3)));
    Prefactor = Kappa*(3/(gamma.(int) + 3));
    D = (sigma.(int)^6)./((sigma.(int)^6) + (U.r.^6));
    dD = (6*(sigma.(int)^6).*(U.r.^5))./(((sigma.(int)^6) + (U.r.^6)).^2);
    
    % C6 values should be in terms on nm and kJ/mol
    if Incl_Disp
        C6 = Kappa*(sigma.(int)^6);
        C6_out.(int) = 1;
    else
        C6 = 1;
        C6_out.(int) = Kappa*(sigma.(int)^6);
    end
    
    %% Build PES
    
    % Components of potential
    [U.(int).f0,U.(int).df0] = Coulomb_Potential(Settings,U.r,int);
    U.(int).g = - C6./(sigma.(int)^6 + U.r.^6); % Dispersion g(r)
    U.(int).h = Prefactor.*D.*exp(gamma.(int).*(1 - (U.r./sigma.(int)))); % Short range repulsion (with possible close-range coulomb damping)
    
    % Negative components of derivative
    U.(int).dg = (-C6.*6.*(U.r.^5))./((sigma.(int)^6 + U.r.^6).^2); % Dispersion -dg(r)/dr
    U.(int).dh = Prefactor.*dD.*exp(gamma.(int).*(1 - (U.r./sigma.(int)))) ...
               + Prefactor.*D.*(gamma.(int)./sigma.(int)).*exp(gamma.(int).*(1 - (U.r./sigma.(int)))); % Short range repulsion
    
    % Shift the potential to zero at the cutoff
    if contains(Settings.(MDP).vdw_modifier,'potential-shift','IgnoreCase',true)
        D_VDW = (sigma.(int)^6)./((sigma.(int)^6) + (Settings.(MDP).RVDW_Cutoff.^6));
        C6_VDW = Kappa*(sigma.(int)^6);
        EVDW_Cutoff = Prefactor.*D_VDW.*exp(gamma.(int).*(1 - (Settings.(MDP).RVDW_Cutoff./sigma.(int)))) ...
                    - C6_VDW./(sigma.(int)^6 + Settings.(MDP).RVDW_Cutoff.^6);
        
        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U.(int).g = U.(int).g - EVDW_Cutoff./C6_out.(int);
    end
    
    % remove infinities
    U.(int) = Remove_Infinities(U.(int));
    
    % Print
    U.(int).f  = k_0*(e_c^2)*q.(int(1))*q.(int(2)).*U.(int).f0;
    U.(int).df = k_0*(e_c^2)*q.(int(1))*q.(int(2)).*U.(int).df0;
    U.(int).Total = U.(int).f + C6_out.(int).*U.(int).g + U.(int).h;
    U.(int).dTotal = -(U.(int).df + C6_out.(int).*U.(int).dg + U.(int).dh);
end

%% PLOT if plotswitch chosen
if Plotswitch
    figure;
    % Options
    lw=2;
    fs=25;

    h = cell(1,8);
    hold on
    switch lower(PlotType)
        case 'full'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.f + U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.f + U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.f + U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Full Potential';
        case 'full-derivative'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.df + U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.df + U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.df + U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Full Potential';
        case 'coulomb'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.f,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.f,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.f,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Coulomb Potential';
        case 'coulomb-derivative'
            h{1} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.X.*U.MX.df,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,k_0*(e_c^2).*q.M*q.M.*U.MM.df,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,k_0*(e_c^2).*q.X*q.X.*U.XX.df,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-600 1000];
            ttxt = 'Derivative of Coulomb Potential';
        case 'lj'
            h{1} = plot(U.r.*10,U.MX.g + U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Lennard-Jones Potential';
        case 'lj-derivative'
            h{1} = plot(U.r.*10,U.MX.dg + U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dg + U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dg + U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Lennard-Jones Potential';
        case 'dispersion'
            h{1} = plot(U.r.*10,U.MX.g,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.g,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.g,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Dispersion Potential';
        case 'dispersion-derivative'
            h{1} = plot(U.r.*10,U.MX.dg,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dg,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dg,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Derivative of Dispersion Potential';
        case 'repulsive'
            h{1} = plot(U.r.*10,U.MX.h,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
            yl = [-50 10];
            ttxt = 'Repulsive Potential';
        case 'repulsive-derivative'
            h{1} = plot(U.r.*10,U.MX.dh,'Color','r','LineWidth',lw,'LineStyle','-');
            h{2} = plot(U.r.*10,U.MM.dh,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(U.r.*10,U.XX.dh,'Color','g','LineWidth',lw,'Linestyle','-');
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
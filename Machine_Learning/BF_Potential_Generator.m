function [U_MX_out, U_MM_out, U_XX_out,C6_out] = BF_Potential_Generator(Settings,varargin)

% Optional inputs
p = inputParser;
p.FunctionName = 'BF_Potential_Generator';
addOptional(p,'PlotType','full');
addOptional(p,'ReturnAsStructure',false);
addOptional(p,'Startpoint',0);
addOptional(p,'Plotswitch',false);
addOptional(p,'MDP_Minimize',false);
addOptional(p,'Include_Dispersion_Scale',true);
parse(p,varargin{:});
PlotType = p.Results.PlotType;
ReturnAsStructure = p.Results.ReturnAsStructure;
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
r = Startpoint:Settings.Table_StepSize:Settings.Table_Length;

%% If Damping at close range, affects all attractive interactions
for interaction = {'MX' 'XX' 'MM'}
    int = interaction{1};
    
    Kappa = 2*epsilon.(int)/(1 - (3/(gamma.(int) + 3)));
    Beta = Kappa*(3/(gamma.(int) + 3));
    D = (sigma.(int)^6)./((sigma.(int)^6) + (r.^6));
    dD = (6*(sigma.(int)^6).*(r.^5))./(((sigma.(int)^6) + (r.^6)).^2);
    
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
    U.(int).f = 1./r; % Electrostatics function f(r)
    U.(int).g = - C6./(sigma.(int)^6 + r.^6); % Dispersion g(r)
    U.(int).h = Beta.*D.*exp(gamma.(int).*(1 - (r./sigma.(int)))); % Short range repulsion (with possible close-range coulomb damping)
    
    % Negative components of derivative
    U.(int).df = 1./(r.^2); % Electrostatics function (not including Coulomb constant or charges)
    U.(int).dg = (-C6.*6.*(r.^5))./((sigma.(int)^6 + r.^6).^2); % Dispersion -dg(r)/dr
    U.(int).dh = Beta.*dD.*exp(gamma.(int).*(1 - (r./sigma.(int)))) ...
               + Beta.*D.*(gamma.(int)./sigma.(int)).*exp(gamma.(int).*(1 - (r./sigma.(int)))); % Short range repulsion
    
    % Shift the potential to zero at the cutoff
    if contains(Settings.(MDP).vdw_modifier,'potential-shift','IgnoreCase',true)
        D_VDW = (sigma.(int)^6)./((sigma.(int)^6) + (Settings.(MDP).RVDW_Cutoff.^6));
        C6_VDW = Kappa*(sigma.(int)^6);
        EVDW_Cutoff = Beta.*D_VDW.*exp(gamma.(int).*(1 - (Settings.(MDP).RVDW_Cutoff./sigma.(int)))) ...
                    - C6_VDW./(sigma.(int)^6 + Settings.(MDP).RVDW_Cutoff.^6);

        % Shift by the dispersion energy at vdw cutoff radius. only affects one
        % energy component, not derivatives (i.e. forces)
        U.(int).g = U.(int).g - EVDW_Cutoff./C6_out.(int);
    end
    
    % remove infinities
    U.(int) = Remove_Infinities(U.(int));
    
    % Print
    U_out = [r ; U.(int).f ; U.(int).df ; U.(int).g ; U.(int).dg ; U.(int).h ; U.(int).dh];
    U.(int).out = deblank( sprintf(['%16.10e   %16.10e %16.10e   %16.10e %16.10e   %16.10e %16.10e' newline],U_out(:)) );
end

if ReturnAsStructure
	U.MX.f0 = U.MX.f;
    U.MM.f0 = U.MM.f;
    U.XX.f0 = U.XX.f;
    
	U.MX.df0 = U.MX.df;
    U.MM.df0 = U.MM.df;
    U.XX.df0 = U.XX.df;
    
	U.MX.f = k_0*(e_c^2).*q.M*q.X.*U.MX.f;
    U.MM.f = k_0*(e_c^2).*q.M*q.M.*U.MM.f;
    U.XX.f = k_0*(e_c^2).*q.X*q.X.*U.XX.f;
    
	U.MX.df = k_0*(e_c^2).*q.M*q.X.*U.MX.df;
    U.MM.df = k_0*(e_c^2).*q.M*q.M.*U.MM.df;
    U.XX.df = k_0*(e_c^2).*q.X*q.X.*U.XX.df;
    
    U.MX.Total = U.MX.f + C6_out.MX.*U.MX.g + U.MX.h;
    U.MM.Total = U.MM.f + C6_out.MM.*U.MM.g + U.MM.h;
    U.XX.Total = U.XX.f + C6_out.XX.*U.XX.g + U.XX.h;
    
    U.MX.dTotal = -(U.MX.df + C6_out.MX.*U.MX.dg + U.MX.dh);
    U.MM.dTotal = -(U.MM.df + C6_out.MM.*U.MM.dg + U.MM.dh);
    U.XX.dTotal = -(U.XX.df + C6_out.XX.*U.XX.dg + U.XX.dh);
    
    U_MX_out = U.MX;
    U_MM_out = U.MM;
    U_XX_out = U.XX;

    U_MX_out.r = r;
    U_MM_out.r = r;
    U_XX_out.r = r;
else
    U_MX_out = U.MX.out;
    U_MM_out = U.MM.out;
    U_XX_out = U.XX.out;
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
            h{2} = plot(r.*10,U.MM.g + U.MM.h,'Color','b','LineWidth',lw,'Linestyle','-');
            h{3} = plot(r.*10,U.XX.g + U.XX.h,'Color','g','LineWidth',lw,'Linestyle','-');
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
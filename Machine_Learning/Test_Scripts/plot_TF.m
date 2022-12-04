PlotType = 'LJ';
epsilon = 1;
sigma = 0.5;
rmin = 0;
rmax = 10;

Settings = Initialize_MD_Settings;
Settings.Salt = 'LiF';
Settings.MDP.vdw_modifier = 'None';
[~,~,TF_Params] = TF_Potential_Parameters(Settings);

% Calculate Tight form Parameters
A_ij = TF_Params.B;
rho_ij = 1/TF_Params.alpha;
C_ij = TF_Params.C;
D_ij = TF_Params.D;

switch PlotType
case 'LJ'
    ffun = @(r) A_ij.*exp(-r./rho_ij) - C_ij./(r.^6) - D_ij./(r.^8);

    r = rmin:0.001:rmax;
    xmax = rmax;
    
    % Calculate potential
    f = ffun(r);

    % Derivative of potential
    ylims = [-0.5 0.5];
case 'Coulomb'

    q1 = 1;
    q2 = -1;
    Ang_per_m = 1e+10; % Ang per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Ang^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ Ang C^-2 mol^-1

    ffun =  @(r)  k_0*(e_c^2)*q1*q2./(r);
    
    q2 = 1;
    ffun_rep =  @(r)  k_0*(e_c^2)*q1*q2./(r);
    
    % Check if rt is larger than the x limits
    r = rmin:0.001:rmax;
    xmax = rmax;

    % Calculate potential
    f = ffun(r);
    f_rep = ffun_rep(r);

    % Derivative of potential
    ylims = [-8000 8000];
    
case {'Full Attractive' 'Full Repulsive'}

    if strcmp(PlotType,'Full Attractive')
        q1 = 1;
        q2 = -1;
    else
        q1 = 1;
        q2 = 1;
    end
    Ang_per_m = 1e+10; % Ang per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Ang^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ Ang C^-2 mol^-1

    ffun = @(r) 4*epsilon.*( ((sigma./r).^12) - ((sigma./r).^6) ) + ...
        k_0*(e_c^2)*q1*q2./(r);
    fpfun = @(r) -48*epsilon.*( ((sigma./r).^13) - (1/2).*((sigma./r).^7) ) - ...
        k_0*(e_c^2)*q1*q2./(r.^2);

    rt = fzero(fpfun,sigma);

    fmin = ffun(rt);
    fpmin = fpfun(rt);

    % Check if rt is larger than the x limits
    if rt > 18.19
        r = rmin:0.001:rt+rt*1.1;
        xmax = rt+rt*1.1;
    else
        r = rmin:0.001:rmax;
        xmax = rmax;
    end

    % Calculate potential
    f = ffun(r);

    % Derivative of potential
    fprime = fpfun(r);
    ylims = [-1000 1000];
end



% Update plot
figh = figure('WindowState','Maximized');
ax = axes;
cla(ax)
hold(ax,'on')
%yline(ax,0,':k',"LineWidth",1)

col = cbrewer('seq','Blues',3);
col_rep = cbrewer('seq','Reds',3);


plot(ax,[0 xmax],[0 0],'--k',"LineWidth",3)

if strcmp(PlotType,'LJ')
    plot(ax,r.*10,f,'-',"LineWidth",5,'Color',col(end,:))
elseif strcmp(PlotType,'Coulomb')
    plot(ax,r,f,'-',"LineWidth",5,'Color',col(end,:))
    plot(ax,r,f_rep,'-',"LineWidth",5,'Color',col_rep(end,:))
else
    plot(ax,r,f,'-',"LineWidth",5,'Color',col(end,:))
end

ylim(ax,ylims)
xlim(ax,[rmin,xmax])

fs=32;
set(gca,'box','on','TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);

%ylabel(gca,'$u^{LJ}$','fontsize',fs,'Interpreter','latex');

if strcmp(PlotType,'LJ')
    xlabel(gca,'$r$ [\AA]','fontsize',fs,'Interpreter','latex');
%     xticks([0 sigma])
%     yticks([-epsilon 0])
%     xticklabels(gca,{'0' '$\sigma_{ij}$'})
%     yticklabels(gca,{'$-\varepsilon_{ij}$' '0'})
%     
else
    xlabel(gca,'$r$ [\AA]','fontsize',fs,'Interpreter','latex');
    ylabel(gca,'$u_{ij}(r)$ [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');
end

exportgraphics(gca,'TF_Pot.png')
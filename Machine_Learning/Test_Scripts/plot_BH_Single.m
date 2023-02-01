PlotType = 'LJ';
epsilon = 3;
r0 = 1;
gamma = 8.2;

titles = {'$0 \leq \gamma < 6$ with $\varepsilon \leq 0$' '$6 < \gamma  < 7$ with $\varepsilon \geq 0$' ...
    '$\gamma = 7$ with $\varepsilon \geq 0$' '$7 < \gamma <\infty$ with $\varepsilon \geq 0$'};
labels = {'\textbf{a})' '\textbf{b})' '\textbf{c})' '\textbf{d})'};
rmin = 0;
rmax = 3;
ylims = [-3.5 3.5];
fs=32;

% Update plot
figh = figure('WindowState','Maximized');
T = tiledlayout(figh,2,2);

%yline(ax,0,':k',"LineWidth",1)


if gamma < 6 
    epsil = -epsilon;
else
    epsil = epsilon;
end

k = 1/(gamma - 6);

switch PlotType
case 'LJ'
    ffun = @(r) 6*epsil*k*exp(gamma*(1 - r./r0)) - epsil*gamma*k*((r0./r).^6);
    fpfun = @(r) -6*(gamma/r0)*epsil*k*exp(gamma*(1 - r./r0)) + 6*epsil*gamma*k*(r0^6)./(r.^7);
    if gamma < 7
        rt = fzero(fpfun,r0+10);
    else
        rt = fzero(fpfun,r0*0.6);
    end

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
case 'Full Attractive'
    q1 = 1;
    q2 = -1;
    Ang_per_m = 1e+10; % Ang per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Ang^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ Ang C^-2 mol^-1

    ffun = @(r) 6*epsil*k*exp(gamma*(1 - r./r0)) - epsil*gamma*k*((r0./r).^6) + ...
        k_0*(e_c^2)*q1*q2./(r);
    fpfun = @(r) -6*(gamma/r0)*epsil*k*exp(gamma*(1 - r./r0)) + 6*epsil*gamma*k*(r0^6)./(r.^7) - ...
        k_0*(e_c^2)*q1*q2./(r.^2);
    if gamma < 7
        rt = fzero(fpfun,r0+10);
    else
        rt = fzero(fpfun,r0*0.6);
    end

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
case 'Full Repulsive'
    q1 = 1;
    q2 = 1;
    Ang_per_m = 1e+10; % Ang per m
    NA = 6.0221409e23; % Molecules per mole
    e_c = 1.60217662e-19; % Elementary charge in Coulombs
    epsilon_0 = (8.854187817620e-12)*1000/(Ang_per_m*NA); % Vacuum Permittivity C^2 mol kJ^-1 Ang^-1
    k_0 = 1/(4*pi*epsilon_0); % Coulomb constant in kJ Ang C^-2 mol^-1

    ffun = @(r) 6*epsil*k*exp(gamma*(1 - r./r0)) - epsil*gamma*k*((r0./r).^6) + ...
        k_0*(e_c^2)*q1*q2./(r);
    fpfun = @(r) -6*(gamma/r0)*epsil*k*exp(gamma*(1 - r./r0)) + 6*epsil*gamma*k*(r0^6)./(r.^7) - ...
        k_0*(e_c^2)*q1*q2./(r.^2);
    if gamma < 7
        rt = fzero(fpfun,r0+10);
    else
        rt = fzero(fpfun,r0*0.6);
    end

    fmin = ffun(rt);
    fpmin = fpfun(rt);

    % Check if rt is larger than the x limits
    if rt > 18.19
        r = 0:0.001:rt+rt*1.1;
        xmax = rt+rt*1.1;
    else
        r = 0:0.001:rmax;
        xmax = rmax;
    end

    % Calculate potential
    f = ffun(r);

    % Derivative of potential
    fprime = fpfun(r);
end

ax = gca;

plot(ax,r,f,'-b',"LineWidth",3,'color','k')

hold(gca,'on')
plot(ax,[0; xmax],[0; 0],'--k',"LineWidth",1)

%plot(ax,[r0 r0],[ffun(r0) ylims(1)],':r',"LineWidth",3)
%plot(ax,[rt rt],[ffun(rt) ylims(1)],':r',"LineWidth",3)

%plot(ax,[0 r0],[ffun(r0) ffun(r0)],':r',"LineWidth",3)
%plot(ax,[0 rt],[fmin fmin],':r',"LineWidth",3)

%scatter(ax,rt,fmin,80,'r','filled','o',"LineWidth",1,'MarkerEdgeColor','k')
%scatter(ax,r0,ffun(r0),80,'r','filled','o',"LineWidth",1,'MarkerEdgeColor','k')
ylim(ax,ylims)
xlim(ax,[rmin,xmax])


set(ax,'box','on','TickLabelInterpreter','latex');
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(ax,'$r$','fontsize',fs,'Interpreter','latex');
ylabel(ax,'$u_{ij}(r)$','fontsize',fs,'Interpreter','latex');


xticks(ax,0)
yticks(ax,0)

%exportgraphics(figh,'BH_pot.eps')
exportgraphics(figh,'BH_pot.png','resolution',300)
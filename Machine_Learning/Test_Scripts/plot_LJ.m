PlotType = 'LJ';
epsilon = 1;
sigma = 0.5;
rmin = 0;
rmax = 2;

% Calculate Tight form Parameters
B = 4*epsilon*(sigma^12); % Repulsive prefactor
C = 4*epsilon*(sigma^6); % 1/r6 dispersion 

switch PlotType
case 'LJ'
    ffun = @(r) 4*epsilon.*( ((sigma./r).^12) - ((sigma./r).^6) );
    fpfun = @(r) -48*epsilon.*( ((sigma./r).^13) - (1/2).*((sigma./r).^7) );
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
    ylims = [-1.5 1];
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




plot(ax,[0 xmax],[0 0],':k',"LineWidth",3)
plot(ax,[sigma sigma],[0 ylims(1)],':k',"LineWidth",3)
plot(ax,[0 rt],[fmin fmin],':k',"LineWidth",3)
plot(ax,r,f,'-b',"LineWidth",4)
scatter(ax,rt,fmin,50,'k','filled','o',"LineWidth",3,'MarkerEdgeColor','k')
scatter(ax,sigma,ffun(sigma),50,'k','filled','o',"LineWidth",3,'MarkerEdgeColor','k')
ylim(ax,ylims)
xlim(ax,[rmin,xmax])

fs=32;
set(gca,'box','on','TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(gca,'$r$','fontsize',fs,'Interpreter','latex');
ylabel(gca,'$u^{LJ}$','fontsize',fs,'Interpreter','latex');

xticks([0 sigma])
yticks([-epsilon 0])
xticklabels(gca,{'0' '$\sigma$'})
yticklabels(gca,{'$-\varepsilon$' '0'})

exportgraphics(gca,'LJ_pot.eps')
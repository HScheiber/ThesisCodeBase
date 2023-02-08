epsilon = 1;
sigma = 0.5;
rmin = 0;
rmax = 2;


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


% Update plot
figh = figure('WindowState','Maximized');
ax = axes;
cla(ax)
hold(ax,'on')
%yline(ax,0,':k',"LineWidth",1)

col = cbrewer('seq','Blues',3);
col_rep = cbrewer('seq','Reds',3);


plot(ax,[0 xmax],[0 0],'--k',"LineWidth",3)

plot(ax,r,f,'-',"LineWidth",5,'Color',col(end,:))
plot(ax,r,f_rep,'-',"LineWidth",5,'Color',col_rep(end,:))


ylim(ax,ylims)
xlim(ax,[rmin,xmax])

fs=45;
set(gca,'box','on','TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);

%ylabel(gca,'$u^{LJ}$','fontsize',fs,'Interpreter','latex');


xlabel(gca,'$r$ [\AA]','fontsize',fs,'Interpreter','latex');
ylabel(gca,'$u_{ij}(r)$ [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');

exportgraphics(gca,'Coulomb_Pot.png')
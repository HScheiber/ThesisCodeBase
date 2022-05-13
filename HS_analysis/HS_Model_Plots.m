%% Load Data
clear;
Structures = {'Rocksalt' 'Wurtzite'};
Scaling = Init_Scaling_Object;

% Repulsive Wall Height (HS only)
% Units: kJ/mol
Scaling.Z.All = 1e10;
Scaling.Z.MX = 1;
Scaling.Z.MM = 1;
Scaling.Z.XX = 1;

% Repulsive Steepness Parameter b (HS only)
% Units: inverse nanometers
Scaling.b.All = Inf;
Scaling.b.MX = 1;
Scaling.b.MM = 1;
Scaling.b.XX = 1;

% Hard Sphere Distance Parameter r_d (HS only)
% Units are nanometers
Scaling.r_d.All = 1;
Scaling.r_d.MM = 0.09;
rM_over_rX = 0.1:0.005:1.0;
XX_Size = (Scaling.r_d.MM./rM_over_rX);

load('Rocksalt_vs_NiAs_HS_Model_M0.09.mat','-mat')

%% Crystal ionic radii and radius ratios
Li = 76/1000; % nm
Na = 102/1000; % nm
F = 133/1000; % nm
Cl = 181/1000; % nm
Br = 196/1000; % nm
I = 220/1000; % nm

Sa.LiF = Li/F;
Sa.LiCl = Li/Cl;
Sa.LiBr = Li/Br;
Sa.LiI = Li/I;
Sa.NaCl = Na/Cl;


% Sa.LiF = 0.62214;
% Sa.LiCl = 0.52374;
% Sa.LiBr = 0.48214;
% Sa.LiI = 0.44001;
% Sa.NaCl = 0.69437;

salt_colors = cbrewer('qual','Set1',5);
sa_txt_shift = 0.008;
cn_txt_shift = 0;


%% Relative Energy
N = 4;
lw = 3;
fs = 30;
leg_text = {'Coulomb' 'Coul:SR M-M' ...
    'Coul:SR M-X' 'Coul:SR X-X'};
colours = cbrewer('qual','Dark2',N);

x = Scaling.r_d.MM./XX_Size;
y = cell(1,N);
y{1} = (E.(Structures{1}).CoulSR    + E.(Structures{1}).CoulLR) ./ (E.(Structures{2}).CoulSR    + E.(Structures{2}).CoulLR);
y{2} = (E.(Structures{1}).CoulSR_MM) ./ (E.(Structures{2}).CoulSR_MM);
y{3} = (E.(Structures{1}).CoulSR_MX) ./ (E.(Structures{2}).CoulSR_MX);
y{4} = (E.(Structures{1}).CoulSR_XX) ./ (E.(Structures{2}).CoulSR_XX);

fg1 = figure('WindowState','maximized');
axh1 = axes(fg1);
hold(axh1,'on')
RRs = [1.0 0.732 0.414 0.225];
RRt = {'CN = 8' 'CN = 6' 'CN = 4'};

for idx = 1:length(RRs)
    xline(axh1,RRs(idx),'Color','k','LineWidth',3,'LineStyle','-.')
end

% RRcs = [ 0.717 0.610 0.326];
% for idx = 1:length(RRcs)
%     xline(axh1,RRcs(idx),'Color','k','LineWidth',1,'LineStyle','-')
% end

yrrs = 1.03;
for idx = 1:length(RRt)
    xrrs = mean([RRs(idx) RRs(idx+1)]);
    text(axh1,xrrs+cn_txt_shift,yrrs,RRt{idx},'Interpreter','latex','fontsize',fs,...
        'HorizontalAlignment','center')
end

Salt_label = 0.905;
xline(axh1,Sa.LiF,'Color',salt_colors(1,:),'LineWidth',2,'LineStyle',':')
text(axh1,Sa.LiF+sa_txt_shift,Salt_label,'LiF','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(1,:))
xline(axh1,Sa.LiCl,'Color',salt_colors(2,:),'LineWidth',2,'LineStyle',':')
text(axh1,Sa.LiCl+sa_txt_shift,Salt_label,'LiCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(2,:))
xline(axh1,Sa.LiBr,'Color',salt_colors(3,:),'LineWidth',2,'LineStyle',':')
text(axh1,Sa.LiBr+sa_txt_shift,Salt_label,'LiBr','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(3,:))
xline(axh1,Sa.LiI,'Color',salt_colors(4,:),'LineWidth',2,'LineStyle',':')
text(axh1,Sa.LiI+sa_txt_shift,Salt_label,'LiI','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(4,:))
xline(axh1,Sa.NaCl,'Color',salt_colors(5,:),'LineWidth',2,'LineStyle',':')
text(axh1,Sa.NaCl+sa_txt_shift,Salt_label,'NaCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(5,:))


hold on
p = gobjects(1,N);
for idx = 1:N
    p(idx) = plot(axh1,x,y{idx},'Color',colours(idx,:),'LineWidth',lw,'LineStyle','-');
end

title(axh1,['Charged Hard Spheres: ' Structures{1} '/' Structures{2} ...
    ' Energies'],'Interpreter','latex','fontsize',fs)

set(axh1,'box','on','TickLabelInterpreter','latex');
set(axh1,'XMinorTick','on','YMinorTick','on','FontSize',fs);
set(axh1, 'Xdir', 'reverse')
xlabel(axh1,'r(M) / r(X)','fontsize',fs,'Interpreter','latex');
ylabel(axh1,['$E_{' Structures{1} '}$/$E_{' Structures{2} '}$'],'fontsize',fs,'Interpreter','latex');
yline(axh1,1,'Color','k','LineWidth',2,'LineStyle',':')

legend(axh1,p,leg_text,'Interpreter','latex','location','northwest');

ylim(axh1,[0.9 1.1]);
xlim(axh1,[0.1 1]);

%% Structure 2 c/a
fg2 = figure('WindowState','maximized');
axh2 = axes(fg2);
hold(axh2,'on')

yy = E.(Structures{2}).c./E.(Structures{2}).a;
q = plot(axh2,x,yy,'Color','k','LineWidth',lw,'LineStyle','-');

for idx = 1:length(RRs)
    xline(axh2,RRs(idx),'Color','k','LineWidth',3,'LineStyle','-.')
end

yrrs = 1.85;
for idx = 1:length(RRt)
    xrrs = mean([RRs(idx) RRs(idx+1)]);
    text(axh2,xrrs+cn_txt_shift,yrrs,RRt{idx},'Interpreter','latex','fontsize',fs,...
        'HorizontalAlignment','center')
end

Salt_label = 1.51;
xline(axh2,Sa.LiF,'Color',salt_colors(1,:),'LineWidth',2,'LineStyle',':')
text(axh2,Sa.LiF+sa_txt_shift,Salt_label,'LiF','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(1,:))
xline(axh2,Sa.LiCl,'Color',salt_colors(2,:),'LineWidth',2,'LineStyle',':')
text(axh2,Sa.LiCl+sa_txt_shift,Salt_label,'LiCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(2,:))
xline(axh2,Sa.LiBr,'Color',salt_colors(3,:),'LineWidth',2,'LineStyle',':')
text(axh2,Sa.LiBr+sa_txt_shift,Salt_label,'LiBr','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(3,:))
xline(axh2,Sa.LiI,'Color',salt_colors(4,:),'LineWidth',2,'LineStyle',':')
text(axh2,Sa.LiI+sa_txt_shift,Salt_label,'LiI','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(4,:))
xline(axh2,Sa.NaCl,'Color',salt_colors(5,:),'LineWidth',2,'LineStyle',':')
text(axh2,Sa.NaCl+sa_txt_shift,Salt_label,'NaCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(5,:))

title(axh2,['Charged Hard Spheres: ' Structures{2} ' c/a'],'Interpreter','latex','fontsize',fs)

set(axh2,'box','on','TickLabelInterpreter','latex');
set(axh2,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(axh2,'r(M) / r(X)','fontsize',fs,'Interpreter','latex');
ylabel(axh2,'c/a','fontsize',fs,'Interpreter','latex');
yline(axh2,0,'Color','k','LineWidth',2,'LineStyle',':')
ylim(axh2,[1.5 1.9]);
set(axh2, 'Xdir', 'reverse')

%% Lattice Params
N = 3;
fg3 = figure('WindowState','maximized');
axh3 = axes(fg3);
leg_text3 = {[Structures{1} ' a'] [Structures{2} ' a'] ...
    [Structures{2} ' c']};
hold(axh3,'on')

colours = cbrewer('qual','Dark2',N);

yyy = cell(1,N);
yyy{1} = E.(Structures{1}).a./(Scaling.r_d.MM*10);
yyy{2} = E.(Structures{2}).a./(Scaling.r_d.MM*10);
yyy{3} = E.(Structures{2}).c./(Scaling.r_d.MM*10);
p3 = gobjects(1,N);
for idx = 1:N
    p3(idx) = plot(axh3,x,yyy{idx},'Color',colours(idx,:),'LineWidth',lw,'LineStyle','-');
end

for idx = 1:length(RRs)
    xline(axh3,RRs(idx),'Color','k','LineWidth',3,'LineStyle','-.')
end

yrrs = 6;
for idx = 1:length(RRt)
    xrrs = mean([RRs(idx) RRs(idx+1)]);
    text(axh3,xrrs+cn_txt_shift,yrrs,RRt{idx},'Interpreter','latex','fontsize',fs,...
        'HorizontalAlignment','center')
end

Salt_label = 0.2;
xline(axh3,Sa.LiF,'Color',salt_colors(1,:),'LineWidth',2,'LineStyle',':')
text(axh3,Sa.LiF+sa_txt_shift,Salt_label,'LiF','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(1,:))
xline(axh3,Sa.LiCl,'Color',salt_colors(2,:),'LineWidth',2,'LineStyle',':')
text(axh3,Sa.LiCl+sa_txt_shift,Salt_label,'LiCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(2,:))
xline(axh3,Sa.LiBr,'Color',salt_colors(3,:),'LineWidth',2,'LineStyle',':')
text(axh3,Sa.LiBr+sa_txt_shift,Salt_label,'LiBr','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(3,:))
xline(axh3,Sa.LiI,'Color',salt_colors(4,:),'LineWidth',2,'LineStyle',':')
text(axh3,Sa.LiI+sa_txt_shift,Salt_label,'LiI','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(4,:))
xline(axh3,Sa.NaCl,'Color',salt_colors(5,:),'LineWidth',2,'LineStyle',':')
text(axh3,Sa.NaCl+sa_txt_shift,Salt_label,'NaCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(5,:))

title(axh3,'Charged Hard Spheres: Relative Lattice Parameters','Interpreter','latex','fontsize',fs)

set(axh3,'box','on','TickLabelInterpreter','latex');
set(axh3,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(axh3,'r(M) / r(X)','fontsize',fs,'Interpreter','latex');
ylabel(axh3,'Lattice Param / r(M)','fontsize',fs,'Interpreter','latex');
legend(axh3,p3,leg_text3,'Interpreter','latex','location','northwest');

ylim(axh3,[0 8.5]);
set(axh3, 'Xdir', 'reverse')


%% Densities
fg4 = figure('WindowState','maximized');
axh4 = axes(fg4);
hold(axh4,'on')

V_S1 = (E.(Structures{1}).a.^3)./4;
V_S2 = (E.(Structures{2}).a.*E.(Structures{2}).b.*sind(60).*E.(Structures{2}).c)./2;
y4 = V_S1./V_S2;

q4 = plot(axh4,x,y4,'Color','k','LineWidth',lw,'LineStyle','-');

for idx = 1:length(RRs)
    xline(axh4,RRs(idx),'Color','k','LineWidth',3,'LineStyle','-.')
end

yrrs = 1.005;
for idx = 1:length(RRt)
    xrrs = mean([RRs(idx) RRs(idx+1)]);
    text(axh4,xrrs+cn_txt_shift,yrrs,RRt{idx},'Interpreter','latex','fontsize',fs,...
        'HorizontalAlignment','center')
end

Salt_label = 0.9905;
xline(axh4,Sa.LiF,'Color',salt_colors(1,:),'LineWidth',2,'LineStyle',':')
text(axh4,Sa.LiF+sa_txt_shift,Salt_label,'LiF','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(1,:))
xline(axh4,Sa.LiCl,'Color',salt_colors(2,:),'LineWidth',2,'LineStyle',':')
text(axh4,Sa.LiCl+sa_txt_shift,Salt_label,'LiCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(2,:))
xline(axh4,Sa.LiBr,'Color',salt_colors(3,:),'LineWidth',2,'LineStyle',':')
text(axh4,Sa.LiBr+sa_txt_shift,Salt_label,'LiBr','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(3,:))
xline(axh4,Sa.LiI,'Color',salt_colors(4,:),'LineWidth',2,'LineStyle',':')
text(axh4,Sa.LiI+sa_txt_shift,Salt_label,'LiI','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(4,:))
xline(axh4,Sa.NaCl,'Color',salt_colors(5,:),'LineWidth',2,'LineStyle',':')
text(axh4,Sa.NaCl+sa_txt_shift,Salt_label,'NaCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(5,:))


title(axh4,['Charged Hard Spheres: Relative Density of ' Structures{1} ' and ' Structures{2}],'Interpreter','latex','fontsize',fs)

set(axh4,'box','on','TickLabelInterpreter','latex');
set(axh4,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(axh4,'r(M) / r(X)','fontsize',fs,'Interpreter','latex');
ylabel(axh4,['$\frac{d_{' Structures{1} '}}{d_{' Structures{2} '}}$'],'fontsize',fs,'Interpreter','latex');
yline(axh4,0,'Color','k','LineWidth',2,'LineStyle',':')
ylim(axh4,[0.99 1.01]);
set(axh4, 'Xdir', 'reverse')

%% Energy Differences
N = 4;
leg_text = {'Coulomb' 'Coul:SR M-M' ...
    'Coul:SR M-X' 'Coul:SR X-X'};
colours = cbrewer('qual','Dark2',N);

y = cell(1,N);
y{1} = (E.(Structures{1}).CoulSR + E.(Structures{1}).CoulLR) - ( E.(Structures{2}).CoulSR + E.(Structures{2}).CoulLR);
y{2} = (E.(Structures{1}).CoulSR_MM)      - (E.(Structures{2}).CoulSR_MM);
y{3} = (E.(Structures{1}).CoulSR_MX)      - (E.(Structures{2}).CoulSR_MX);
y{4} = (E.(Structures{1}).CoulSR_XX)      - (E.(Structures{2}).CoulSR_XX);
%y{5} = E.(Structures{1}).Pot            - E.(Structures{2}).Pot;
%y{6} = E.(Structures{1}).CoulLR         - E.(Structures{2}).CoulLR;

fg5 = figure('WindowState','maximized');
axh5 = axes(fg5);
hold(axh5,'on')
RRs = [1.0 0.732 0.414 0.225];
RRt = {'CN = 8' 'CN = 6' 'CN = 4'};

for idx = 1:length(RRs)
    xline(axh5,RRs(idx),'Color','k','LineWidth',3,'LineStyle','-.')
end

yrrs = 5;
for idx = 1:length(RRt)
    xrrs = mean([RRs(idx) RRs(idx+1)]);
    text(axh5,xrrs+cn_txt_shift,yrrs,RRt{idx},'Interpreter','latex','fontsize',fs,...
        'HorizontalAlignment','center')
end

Salt_label = -65;
xline(axh5,Sa.LiF,'Color',salt_colors(1,:),'LineWidth',2,'LineStyle',':')
text(axh5,Sa.LiF+sa_txt_shift,Salt_label,'LiF','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(1,:))
xline(axh5,Sa.LiCl,'Color',salt_colors(2,:),'LineWidth',2,'LineStyle',':')
text(axh5,Sa.LiCl+sa_txt_shift,Salt_label,'LiCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(2,:))
xline(axh5,Sa.LiBr,'Color',salt_colors(3,:),'LineWidth',2,'LineStyle',':')
text(axh5,Sa.LiBr+sa_txt_shift,Salt_label,'LiBr','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(3,:))
xline(axh5,Sa.LiI,'Color',salt_colors(4,:),'LineWidth',2,'LineStyle',':')
text(axh5,Sa.LiI+sa_txt_shift,Salt_label,'LiI','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(4,:))
xline(axh5,Sa.NaCl,'Color',salt_colors(5,:),'LineWidth',2,'LineStyle',':')
text(axh5,Sa.NaCl+sa_txt_shift,Salt_label,'NaCl','Interpreter','latex','fontsize',fs,...
    'HorizontalAlignment','left','Rotation',90,'Color',salt_colors(5,:))

hold on
p = gobjects(1,N);
for idx = 1:N
    p(idx) = plot(axh5,x,y{idx},'Color',colours(idx,:),'LineWidth',lw,'LineStyle','-');
end

title(axh5,['Charged Hard Spheres [r(M) = 0.9 \AA]: ' Structures{1} ...
    ' - ' Structures{2} ' Energies'],'Interpreter','latex','fontsize',fs)

set(axh5,'box','on','TickLabelInterpreter','latex');
set(axh5,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(axh5,'r(M) / r(X)','fontsize',fs,'Interpreter','latex');
ylabel(axh5,['$E_{' Structures{1} '}$ - $E_{' Structures{2} '}$ (kJ/mol)'],'fontsize',fs,'Interpreter','latex');
yline(axh5,0,'Color','k','LineWidth',2,'LineStyle',':')

legend(axh5,p,leg_text,'Interpreter','latex','location','northeast');

ylim(axh5,[-70 100]);
set(axh5, 'Xdir', 'reverse')
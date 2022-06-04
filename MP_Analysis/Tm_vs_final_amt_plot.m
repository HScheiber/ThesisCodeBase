Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Set61_TvsXliq.mat';
DataKeyword = 'Set61';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salt = 'NaCl';

Data = load(SaveDataDir,'Data').Data;
sz = 80;
fs = 28;

cols = cbrewer('qual', 'Set2', 3);


figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'FontSize',fs,'TickLabelInterpreter','latex',...
    'XMinorGrid','off','YLimitMethod','padded',...
    'YMinorGrid','off','Position',[0.08 0.13 0.9 0.84]);
hold(axh,'on')

p(1) = scatter(Data.T_Xs,100*Data.Y_Xs,sz,'MarkerEdgeColor','k',...
              'MarkerFaceColor',cols(1,:),...
              'LineWidth',1);
p(2) = scatter(Data.T_UBs,100*Data.Y_UBs,sz,'MarkerEdgeColor','k',...
              'MarkerFaceColor',cols(2,:),...
              'LineWidth',1);
p(3) = scatter(Data.T_LBs,100*Data.Y_LBs,sz,'MarkerEdgeColor','k',...
              'MarkerFaceColor',cols(3,:),...
              'LineWidth',1);
legtxt = {'$T^{x}$' '$T_{UB}$' '$T_{LB}$'};

legend(p,legtxt,'interpreter','latex','location','Southeast')

set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
xlim(axh,'padded')
ylim(axh,'padded')
xlabel(axh,'$T$ [K]','Interpreter','latex');
ylabel(axh,'$x_{\mathrm{liq}}$ [\%]','Interpreter','latex')
grid(axh,'minor')
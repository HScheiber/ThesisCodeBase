Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Set60_TvsXliq.mat';
DataKeyword = 'Set60';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salt = 'NaCl';

Data = load(SaveDataDir,'Data').Data;
sz = 80;
fs = 32;


w = 1;
dist_start = 1280;
dist_end = 1300;
bin_edges = dist_start:w:dist_end;
bin_centers = bin_edges(1:end-1)+w/2;
xlimits = [bin_edges(1) bin_edges(end)];

N_tot = length(Data.T_UBs);
Y_Xs_Dat = histcounts(Data.T_Xs,bin_edges);
Y_UB_Dat = histcounts(Data.T_UBs,bin_edges);
Y_LB_Dat = histcounts(Data.T_LBs,bin_edges);

Colours = cbrewer('qual', 'Set1', 3);
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'FontSize',fs,'TickLabelInterpreter','latex',...
    'XMinorGrid','off','YLimitMethod','padded',...
    'YMinorGrid','off','Position',[0.08 0.13 0.9 0.84]);
hold(axh,'on')



p(1) = area(axh,bin_centers,Y_Xs_Dat./N_tot,'FaceColor',Colours(1,:),...
    'EdgeColor',Colours(1,:),'Linewidth',2,'FaceAlpha',0.25);
p(2) = area(axh,bin_centers,Y_UB_Dat./N_tot,'FaceColor',Colours(2,:),...
    'EdgeColor',Colours(2,:),'Linewidth',2,'FaceAlpha',0.25);
p(3) = area(axh,bin_centers,Y_LB_Dat./N_tot,'FaceColor',Colours(3,:),...
    'EdgeColor',Colours(3,:),'Linewidth',2,'FaceAlpha',0.25);

legtxt = {'$T^{x}$' '$T_{UB}$' '$T_{LB}$'};

legend(axh,p,legtxt,'interpreter','latex','location','Southeast','FontSize',fs)

set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
xlim(axh,'padded')
ylim(axh,'padded')
xlabel(axh,'$T$ [K]','Interpreter','latex');
ylabel(axh,'$\rho(T)$','Interpreter','latex')
grid(axh,'minor')

if strcmp(DataKeyword,'Set61')
    exportgraphics(axh ,'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Manuscript_4\SI_Figures\TUB_TLB_Hist_N11664.pdf',...
        'ContentType','vector','BackgroundColor','none')
elseif strcmp(DataKeyword,'Set60')
    exportgraphics(axh ,'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Manuscript_4\SI_Figures\TUB_TLB_Hist_N2000.pdf',...
        'ContentType','vector','BackgroundColor','none')
end
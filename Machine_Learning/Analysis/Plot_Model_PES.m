
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
Settings = Initialize_MD_Settings;
Settings.Salt = 'LiI';
Settings.Theory = 'BH';
ModelID = 'NB2';
PlotType = 'full';
fs = 28; % font size
lw = 3; % line width
savefile = false; % switch to save the final plots to file
Startpoint = 0.08; % nm

% Find reps of models
if length(ModelID) > 2
    Models = [{''} ModelID];
else
    files = {dir(fullfile(ML_results_dir,Settings.Salt,['*' ModelID '*'])).name};
    Models = unique(regexp(files,[ModelID '([0-9]|[a-z])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models);
    Models = Models(nonempty_idx);
    if ~isempty(Models)
        [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
        Models = Models(idx);
    end
    Models = [{''} Models];
end

N_Models = length(Models);
% Col_MX = [0.2 0.2 0.2; cbrewer('seq','Blues',(N_Models-1)/2); cbrewer('seq','Reds',(N_Models-1)/2) ];
%Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
Col_MX = [0 0 0; cbrewer('qual','Set1',N_Models,'spline')];
Col_MX = max(min(Col_MX,1),0);


%% Other parameters
Settings.WaterModel = 'SPC/E';
Settings.Table_Length = 3; % nm
Settings.Table_StepSize = 0.001;
Settings.MDP.vdw_modifier = 'potential-shift';
Settings.MDP.RVDW_Cutoff = 1; % nm

% Pre allocate graphics objects arrays
h_MX = gobjects(N_Models,1);
h_MM = gobjects(N_Models,1);
h_XX = gobjects(N_Models,1);
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);

% Loop through models
figh = figure('WindowState','Maximized');
ax = axes(figh,'Position',[0.100 0.1100 0.850 0.8150]); % [left bottom width height]
hold(ax,'On')
for idx = N_Models:-1:1
    ML_Model_Name = Models{idx};
    Settings.Model = ML_Model_Name;
    [Settings,ModelFound] = Load_Model_Params(Settings);
    if ~ModelFound
        h_MM(idx) = [];
        h_XX(idx) = [];
        h_MX(idx) = [];
        continue
    end
    
    if strcmp(Settings.Theory,'JC')
        [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings,...
            'Startpoint',Startpoint,'ReturnAsStructure',true);
        disp([Settings.Salt ' ' ML_Model_Name ': Q = ' num2str(Settings.S.Q)])
    elseif strcmp(Settings.Theory,'TF')
        
        [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings,...
            'Startpoint',Startpoint,'ReturnAsStructure',true);
    elseif strcmp(Settings.Theory,'BH')
        
        [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings,...
            'Startpoint',Startpoint,'ReturnAsStructure',true);
    elseif strcmp(Settings.Theory,'BD')
        [U_MX, U_MM, U_XX] = BD_Potential_Generator(Settings,...
            'Startpoint',Startpoint,'ReturnAsStructure',true);
    elseif strcmp(Settings.Theory,'Mie')
        [U_MX, U_MM, U_XX] = Mie_Potential_Generator(Settings,...
            'Startpoint',Startpoint,'ReturnAsStructure',true);
    end
    
    switch lower(PlotType)
        case 'full'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.f + U_MX.g + U_MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.f + U_MM.g + U_MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.f + U_XX.g + U_XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'full-derivative'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.df + U_MX.dg + U_MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.df + U_MM.dg + U_MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.df + U_XX.dg + U_XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'lj'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.g + U_MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.g + U_MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.g + U_XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'lj-derivative'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.dg + U_MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.dg + U_MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.dg + U_XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'dispersion'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.g,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.g,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.g,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'dispersion-derivative'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'repulsive'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
        case 'repulsive-derivative'
            h_MX(idx) = plot(ax,U_MX.r.*10,U_MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
            h_MM(idx) = plot(ax,U_MM.r.*10,U_MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
            h_XX(idx) = plot(ax,U_XX.r.*10,U_XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
    end
end

switch lower(PlotType)
    case 'full'
        yl = [-800 2000];
        ttxt = 'Full Potential';
    case 'full-derivative'
        yl = [-600 1000];
        ttxt = 'Derivative of Full Potential';
    case 'lj'
        yl = [-10 10];
        ttxt = 'Lennard-Jones Potential';
    case 'lj-derivative'
        yl = [-50 30];
        ttxt = 'Derivative of Lennard-Jones Potential';
    case 'dispersion'
        yl = [-50 10];
        ttxt = 'Dispersion Potential';
    case 'dispersion-derivative'
        yl = [-50 10];
        ttxt = 'Derivative of Dispersion Potential';
    case 'repulsive'
        yl = [-20 50];
        ttxt = 'Repulsive Potential';
    case 'repulsive-derivative'
        yl = [-20 50];
        ttxt = 'Derivative of Repulsive Potential';
end

if ~savefile
    title(ax,['Plot of ' ttxt ' for ' Settings.Salt ' ' Settings.Theory ' Models'],...
       'Interpreter','latex','fontsize',fs)
end

set(ax,'box','on','TickLabelInterpreter','latex');
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel(ax,'Separation [\AA]','fontsize',fs,'Interpreter','latex');
ylabel(ax,'Pair Potential [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');

ylim(yl);
xlim([0 6.5]);

% Blank line
hline = refline([0 0]);
hline.Color = 'k';
hline.LineWidth = lw-1;
hline.LineStyle = '-.';

if savefile
    Model_Labels = [{[Metal '$^{+}$' Halide '$^{-}$']}  ...
        {[Metal '$^{+}$' Metal '$^{+}$']} ...
        {[Halide '$^{-}$' Halide '$^{-}$']}];
    legend(ax,[h_MX(2:end); h_MM(2:end); h_XX(2:end)],Models{2:end},'Interpreter','latex','location','northeast',...
        'NumColumns',1);
    
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
else
    Model_Labels = [{[Metal '$^{+}$' Halide '$^{-}$']}  ...
        {[Metal '$^{+}$' Metal '$^{+}$']} ...
        {[Halide '$^{-}$' Halide '$^{-}$']}];
    legh = legend(ax,h_MX(2:end),Models{2:end},'Interpreter','latex',...
        'NumColumns',1,'Location','NorthEast');
    legend('boxoff')
    
    
    %% Copy the axes and plot the second legned
    ah = copyobj( ax, figh);
    ah.Visible = 'off';
    lh = legend(ah, [h_MX(1); h_MM(1); h_XX(1)],Model_Labels,'Interpreter','latex',...
        'NumColumns',1,'Location','SouthEast');
    legend('boxoff')
    
end
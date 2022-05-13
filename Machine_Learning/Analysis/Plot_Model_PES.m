% Script for plotting an error comparison for multiple models, to show how
% each performs.
%#ok<*UNRCH>

% %% Analysis Parameters: Comparison LiBr JC C3 vs C6 set
% Salt = 'LiI';
% Theory = 'JC';
% Models = {'' 'CBa' 'BC6' 'CDa'};
% Model_Labels = {'JC' 'Target $E$[RS]' 'Target $E$[RS] - $E$[WZ]' 'Target $a_{RS}$'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_PES.png'];
% 
% PlotType = 'lj';
% 
% N_Models = length(Models);
% Col_MX = [0.2 0.2 0.2; cbrewer('qual','Set2',N_Models)];
% %Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
% %Col_MX = [0 0 0; cbrewer('qual','Dark2',N_Models)];


%% Analysis Parameters: Comparison LiBr JC additivity vs no additivity
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'' 'C3a' 'C3b' 'C3c' 'C3d' 'C3e' 'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5'};
% Model_Labels = {'JC' 'No Additivity' 'No Additivity' 'No Additivity' 'No Additivity' 'No Additivity' ...
%     'Additivity' 'Additivity' 'Additivity' 'Additivity' 'Additivity'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_PES.png'];
% 
% PlotType = 'lj';
% 
% N_Models = length(Models);
% Col_MX = [0.2 0.2 0.2; cbrewer('seq','Blues',(N_Models-1)/2); cbrewer('seq','Reds',(N_Models-1)/2) ];
% %Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
% %Col_MX = [0 0 0; cbrewer('qual','Dark2',N_Models)];
% 
% Col_MX = min(Col_MX,1);

%% Analysis Parameters: Comparison LiBr JC charge scaling vs no charge scaling
% Salt = 'LiBr';
% Theory = 'JC';
% Models = {'' 'C5a' 'C5b' 'C5c' 'C5d' 'C5e' 'C6a1' 'C6a2' 'C6a3' 'C6a4' 'C6a5'};
% Model_Labels = {'JC' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' 'Parameter Charge' ...
%     'Fixed Charge' 'Fixed Charge' 'Fixed Charge' 'Fixed Charge' 'Fixed Charge'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = RS/WZ properties}'];
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_PES.png'];
% 
% PlotType = 'lj';
% 
% N_Models = length(Models);
% Col_MX = [0.2 0.2 0.2; cbrewer('seq','Blues',(N_Models-1)/2); cbrewer('seq','Reds',(N_Models-1)/2) ];
% %Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
% %Col_MX = [0 0 0; cbrewer('qual','Dark2',N_Models)];
% 
% Col_MX = min(Col_MX,1);

%% Analysis Parameters: Comparison LiBr JC charge scaling vs no charge scaling
% Salt = 'LiBr';
% Theory = 'TF';
% Models = {'' 'B1'};
% Model_Labels = {'JC' 'Model B1'};
% Title_txt = ['\bf{Comparison of ' Salt ' ' Theory ' Models (Fixed Charge). Target = Everything}'];
% filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_PES.png'];
% 
% PlotType = 'full';
% 
% N_Models = length(Models);
% % Col_MX = [0.2 0.2 0.2; cbrewer('seq','Blues',(N_Models-1)/2); cbrewer('seq','Reds',(N_Models-1)/2) ];
% %Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
% Col_MX = [0 0 0; cbrewer('qual','Dark2',N_Models)];
% 
% Col_MX = min(Col_MX,1);

ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
Settings = Initialize_MD_Settings;
Settings.Salt = 'LiI';
Settings.Theory = 'BH';
Basenum = 'F';
Midnum = 'A';

% Find reps of models
files = {dir(fullfile(ML_results_dir,['*' Basenum Midnum '*'])).name};
Models = unique(regexp(files,[Basenum Midnum '([0-9]|[a-z])+'],'match','once'));
nonempty_idx = ~cellfun(@isempty,Models);
Models = Models(nonempty_idx);
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end
Models = [{''} Models];

Title_txt = ['\bf{Comparison of ' Settings.Salt ' ' Settings.Theory ' Models (Fixed Charge). Target = RS Properties}'];
filename = [Settings.Salt '_' Settings.Theory '_' Models{1} '-' Models{end} '_PES.png'];

PlotType = 'full';

N_Models = length(Models);
% Col_MX = [0.2 0.2 0.2; cbrewer('seq','Blues',(N_Models-1)/2); cbrewer('seq','Reds',(N_Models-1)/2) ];
%Col_MX = [0.5 0.5 0.5; cbrewer('qual','Pastel2',N_Models)];
Col_MX = [0 0 0; cbrewer('qual','Set1',N_Models,'spline')];
Col_MX = max(min(Col_MX,1),0);


%% Other parameters
savefile = true; % switch to save the final plots to file
fs = 28; % font size
lw = 3; % line width
Startpoint = 0.001;
Settings.WaterModel = 'SPC/E';
Settings.Table_Length = 3;
Settings.Table_StepSize = 0.001;
Settings.MDP.vdw_modifier = 'potential-shift';
Settings.MDP.RVDW_Cutoff = 1.9;

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
    if ~isempty(ML_Model_Name)
        [Settings,ModelFound] = Load_Model_Params(Settings);
        if ~ModelFound
            h_MM(idx) = [];
            h_XX(idx) = [];
            h_MX(idx) = [];
            continue
        end
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
    legend(ax,[h_MX(1); h_MM(1); h_XX(1)],Model_Labels,'Interpreter','latex','location','northeast',...
        'NumColumns',1);
    
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
else
    Model_Labels = [{[Metal '$^{+}$' Halide '$^{-}$']}  ...
        {[Metal '$^{+}$' Metal '$^{+}$']} ...
        {[Halide '$^{-}$' Halide '$^{-}$']} ...
        Model_Labels{2:end}];
    legend(ax,[h_MX(1); h_MM(1); h_XX(1); h_MX(2:end)],Model_Labels,'Interpreter','latex','location','northeast',...
        'NumColumns',2);
end
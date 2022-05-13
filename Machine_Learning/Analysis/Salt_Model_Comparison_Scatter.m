% Script for plotting an error comparison for multiple models, to show how
% each performs.
%#ok<*UNRCH>

%% Analysis Parameters: Comparison Model A1

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'JC';
Models = {'D1' 'D1' 'D1' 'BD1'}; % Must be same length as number of salts

%Title_txt = 'Fixed Charge';
Title_txt = '';
Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
filename = ['Salt_Compare_JC_Model_' Models{1} '.png'];

Plot_Energies = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
Energies_ylim = [-Inf Inf];

Plot_Rel_Energies = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
Rel_Energies_ylim = [-Inf Inf];

Plot_a = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive'};
a_ylim = [-Inf Inf];

Plot_c = {'Wurtzite' 'NiAs' 'FiveFive'};
c_ylim = [-Inf Inf];

Plot_ca = {};
ca_ylim = [-Inf Inf];

Plot_loss = {};

%% Other parameters
set_ylim = true;
savefile = true; % switch to save the final plots to file
show_as_percent_error = false; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)

% Calculates error in lattice energy with respect to experimental lattice energies when available
Target_Experimental_Energies = true; 
% Calculates error in lattice parameter with respect to experimental lattice parameters when available
Target_Experimental_latpar = true; 

fs = 22; % font size

%% Pre-allocate data arrays
N_Salts = length(Salts); % number of models plus the DFT results
N_Structures = length(Structures);

Energies = nan(N_Salts,N_Structures);
Rel_Energies = nan(N_Salts,N_Structures);
Rel_Energies_Target = nan(N_Salts,N_Structures);
LatPars_a = nan(N_Salts,N_Structures);
LatPars_c = nan(N_Salts,N_Structures);
LatPars_ca = nan(N_Salts,N_Structures);
Total_loss = nan(N_Salts,1);

% Plot color scheme
Colours = cbrewer('qual','Set3',N_Salts);
%Colours = [0 0 0; cbrewer('qual','Paired',N_Models)];


% Load target data
DFT = Load_Best_DFT_Data;

if Target_Experimental_Energies || Target_Experimental_latpar
    EXP = Load_Experimental_Data;
    
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        
        if Target_Experimental_Energies
            Correction = EXP.(Salt).Rocksalt.E - DFT.(Salt).Rocksalt.Energy;
            for jdx = 1:length(Structures)
                DFT.(Salt).(Structures{jdx}).Energy = DFT.(Salt).(Structures{jdx}).Energy + Correction;
            end
        end

        if Target_Experimental_latpar
            DFT.(Salt).Rocksalt.a = EXP.(Salt).Rocksalt.a;
            DFT.(Salt).Rocksalt.b = EXP.(Salt).Rocksalt.b;
            DFT.(Salt).Rocksalt.c = EXP.(Salt).Rocksalt.c;
            if ~isnan(EXP.(Salt).Wurtzite.a)
                DFT.(Salt).Wurtzite.a = EXP.(Salt).Wurtzite.a;
                DFT.(Salt).Wurtzite.b = EXP.(Salt).Wurtzite.b;
                DFT.(Salt).Wurtzite.c = EXP.(Salt).Wurtzite.c;
            end
        end
    end
    clear('EXP')
end

% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

% Loop over chosen models
for idx = 1:N_Salts
    Salt = Salts{idx};
    Model = Models{idx};
    
    % Find the fully optimized model
    sres = dir([ML_results_dir filesep '*' Salt '*' Theory '*' 'Model_' Model '_*fullopt.mat']);
    
    if length(sres) > 1
        warning(['Multiple results found for: ' Salt ', ' Theory ', Model ' Model '. Using first result.']);
        sres = sres(1);
    elseif isempty(sres)
        error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
    end
    
    model_filename = fullfile(sres.folder,sres.name);
    
    % Load data
    data = load(model_filename);
    try
        Minimization_Data = data.Minimization_Data;
        Loss = data.loss;
        clear('data')
    catch
        warning(['No saved crystal data for ' Salt ' ' Theory ' Model ' Model])
        continue
    end
    
    % Grab available structures from data
    Min_Dat_Structs = cell(1,length(Minimization_Data));
    for jdx = 1:length(Minimization_Data)
        Min_Dat_Structs{jdx} = Minimization_Data{jdx}.Structure;
    end
    
    % Loop through structures
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        
        % Absolute energies
        if contained_in_cell(Structure,Plot_Energies)
            if show_as_percent_error
                Energies(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                    DFT.(Salt).(Structure).Energy)*100/ ...
                    DFT.(Salt).(Structure).Energy;
            else
                Energies(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                    DFT.(Salt).(Structure).Energy; 
            end
        end
        
        % Relative Energies
        if contained_in_cell(Structure,Plot_Rel_Energies)
            
            Target_gap = DFT.(Salt).(Structure).Energy - DFT.(Salt).Rocksalt.Energy;
            
            Model_gap  = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.E - ...
                         Minimization_Data{strcmpi(Min_Dat_Structs,'Rocksalt')}.E;
            
            Rel_Energies(idx,jdx) = Model_gap;
            Rel_Energies_Target(idx,jdx) = Target_gap;
        end
        
        % lattice param a
        if contained_in_cell(Structure,Plot_a)
            if show_as_percent_error
                LatPars_a(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a - DFT.(Salt).(Structure).a)*100/ ...
                    DFT.(Salt).(Structure).a;
            else
                LatPars_a(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a - DFT.(Salt).(Structure).a;
            end
        end
        
        % lattice param c
        if contained_in_cell(Structure,Plot_c)
            if show_as_percent_error
                LatPars_c(idx,jdx) = (Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c - DFT.(Salt).(Structure).c)*100/ ...
                    DFT.(Salt).(Structure).c;
            else
                LatPars_c(idx,jdx) = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c - DFT.(Salt).(Structure).c;
            end
        end
        
        % lattice param c/a
        if contained_in_cell(Structure,Plot_ca)
            Target_ca = DFT.(Salt).(Structure).c / DFT.(Salt).(Structure).a;
            
            Model_ca  = Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.c / ...
                         Minimization_Data{strcmpi(Min_Dat_Structs,Structure)}.a;
            
            if show_as_percent_error
                LatPars_ca(idx,jdx) = (Model_ca - Target_ca)*100 / Target_ca;
            else
                LatPars_ca(idx,jdx) = Model_ca - Target_ca;
            end
        end
    end
    
    % total loss
    if contained_in_cell(Structure,Plot_loss)
        Total_loss(idx) = Loss;
    end
end

N_compares = 0;
if ~isempty(Plot_Energies)
    N_compares = N_compares+1;
end
if ~isempty(Plot_Rel_Energies)
    N_compares = N_compares+1;
end
if ~isempty(Plot_a)
    N_compares = N_compares+1;
end
if ~isempty(Plot_c)
    N_compares = N_compares+1;
end
if ~isempty(Plot_ca)
    N_compares = N_compares+1;
end
if ~isempty(Plot_loss)
    N_compares = N_compares+1;
end

all_y = {Total_loss, Energies, Rel_Energies, LatPars_a, LatPars_c, LatPars_ca};
all_ylim = {[],Energies_ylim, Rel_Energies_ylim, a_ylim, c_ylim, ca_ylim};
plot_types = {'loss' 'LE' 'RLE' 'a' 'c' 'ca'};

bar_x = 1:N_Structures;
Titles = {'Final Optimized Loss' ...
    'Error in $E_L$' ...
    '$E_{L}$[Structure] - $E_{L}$[Rocksalt]' ...
    'Error in $a$' ...
    'Error in $c$' ...
    'Error in $c/a$'};
if show_as_percent_error
    ylabs = {'Loss' ...
        '[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
else
    ylabs = {'Loss' ...
        '[kJ mol$^{-1}$]' ...
        '[kJ mol$^{-1}$]' ...
        '[\AA]' ...
        '[\AA]' ...
        '$(c/a)$'};
end

pobj = gobjects(1,N_Salts);

% Create figure and axis
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
t = tiledlayout(figh,N_compares,1,'TileSpacing','tight');

for idx = 1:length(all_y)
    bar_y = all_y{idx};
    if all(isnan(bar_y),'all')
        continue
    end
    plot_type = plot_types{idx};
    cTitle = Titles{idx};
    ylab = ylabs{idx};
    cYLIM = all_ylim{idx};
    
    axobj = nexttile;
    if strcmp(plot_type,'loss')
        p = bar(axobj,bar_y,'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        p.CData = Colours;
        
        xticks(axobj,1:N_Salts)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Salts+0.5])
        xticklabels(axobj,[])
        set(axobj, 'YScale', 'log')
    elseif strcmp(plot_type,'RLE')
        
        for jdx = 1:size(bar_y,2)
            c_bar_y = bar_y(:,jdx);
            c_bar_x = linspace(jdx-0.3,jdx+0.3,N_Salts);
            
            hold on
            
            plot(axobj,c_bar_x,Rel_Energies_Target(:,jdx),'-','linewidth',2,'color','r')
            for kdx = 1:size(bar_y,1)
                scatter(axobj,c_bar_x(kdx),Rel_Energies_Target(kdx,jdx),100,'MarkerEdgeColor','r',...
                    'MarkerFaceColor',Colours(kdx,:),'linewidth',2)
            end
            
            plot(axobj,c_bar_x,c_bar_y,'--k','linewidth',2);
            for kdx = 1:size(bar_y,1)
                scatter(axobj,c_bar_x(kdx),c_bar_y(kdx),100,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(kdx,:),'linewidth',2)
            end
        end
        
        yline(axobj,0,':k','linewidth',1)
           
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        
        xticks(axobj,1:N_Structures)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        xticklabels(axobj,[])
        if set_ylim
            ylim(axobj,cYLIM);
        else
            ylim(axobj,'padded')
        end
    else

        for jdx = 1:size(bar_y,2)
            c_bar_y = bar_y(:,jdx);
            c_bar_x = linspace(jdx-0.3,jdx+0.3,N_Salts);
            
            hold on
            
            plot(axobj,c_bar_x,c_bar_y,'--k','linewidth',2);
            for kdx = 1:size(bar_y,1)
                pobj(kdx) = scatter(axobj,c_bar_x(kdx),c_bar_y(kdx),100,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(kdx,:),'linewidth',2);
            end
        end
        
        yline(axobj,0,':k','linewidth',1)
        xticks(axobj,1:N_Structures)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        xticklabels(axobj,[])
        if set_ylim
            ylim(axobj,cYLIM);
        else
            ylim(axobj,'padded')
        end
    end
    
    title(axobj, cTitle,'Fontsize',fs,'Interpreter','latex')
    grid(axobj,'minor');
    axobj.YGrid = 'on';
    axobj.XGrid = 'off';
    
    ylabel(axobj, ylab,'Fontsize',fs-1,'Interpreter','latex')
    
end

title(t, Title_txt,'Fontsize',fs+5,'Interpreter','latex')
xticklabels(axobj,strrep(Structures,'FiveFive','5-5'))

h(1) = plot(nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiF',...
    'MarkerFaceColor',Colours(1,:),'linewidth',2,'Color','k');
h(2) = plot(nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiCl',...
    'MarkerFaceColor',Colours(2,:),'linewidth',2,'Color','k');
h(3) = plot(nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiBr',...
    'MarkerFaceColor',Colours(3,:),'linewidth',2,'Color','k');
h(4) = plot(nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiI',...
    'MarkerFaceColor',Colours(4,:),'linewidth',2,'Color','k');
h(5) = plot(nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'DFT',...
    'MarkerFaceColor','w','linewidth',2,'Color','r');


legh = legend(h,'Location','EastOutside','Orientation','Horizontal',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', 1);

legh.Layout.Tile = 'East';
title(legh,'\bf{\underline{Salts}}','Fontsize',fs,'Interpreter','latex')

if savefile
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
end
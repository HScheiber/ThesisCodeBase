% INPUTS: Label should be a string. One of: 'LiF', 'LiLi', 'FF'
% Basis_set should be something like 'Def2QZVPP' or 'AUG-cc-pVDZ'
function Plot_PES_compare_single_basis_full(Label,Basis_set)
% Label='LiF_W';
% Basis_set='pob-TZVP';
%% Plot Options
fs = 18; % Font size
lw = 3; % Line width
ep = 8; % y axis end point in angstrom
home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
%Basis_set = 'Def2QZVPP';
%Basis_set = 'AUG-cc-pVDZ';
%linetypes = {'-' '--' ':' '-.'};
linetypes = {'-' '-' '-' '-'};
plotTime = true; % 1 = include bar plot of CPU times, 0 = exclude it
x_axis = "bl"; % One of "a" or "bl" (bond length)
%% Reference LiF_Data

data_dir = [home filesep 'data'];
if strcmp(Label,'LiF')
    load([data_dir filesep 'LiF_Pair_Potential.mat'],'Data')
elseif strcmp(Label,'LiLi')
    load([data_dir filesep 'LiLi_Ion_Pair_Potential.mat'],'Data')
elseif strcmp(Label,'FF')
    load([data_dir filesep 'FF_Ion_Pair_Potential.mat'],'Data')
elseif strcmp(Label,'LiF_W')
    load([data_dir filesep 'LiF_Wurtzite_Lattice_Energies.mat'],'Data')
elseif strcmp(Label,'LiF_R')
    load([data_dir filesep 'LiF_Rocksalt_Lattice_Energies.mat'],'Data')
else
    error('Error: Unknown Label Selected.')
end

% Unit Conversions
Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole

fields = fieldnames(Data);

%% First Loop: determine number of objects to plot: x
x = 0;
datList = {''};
% Loop through levels of theory
for i = 1:numel(fields)
    M = size(Data.(fields{i}){1,2},1) - 1; % Number of different sub-methods
    N = size(Data.(fields{i}),1); % Number of different basis sets
    % Loop through sub-methods
    for j=1:M
        % Loop through basis sets
        for k=1:N
            % Current Basis set
            Basis = Data.(fields{i}){k,1};
            datType = Data.(fields{i}){k,2}{j,1};
            if strcmp(Basis,Basis_set) %&& ~ismember(datType,datList) 
                x=x+1;
                datList{x} = datType;
            end
        end
    end
end

% Pre-assign colors and sub_type arrays
sub_type = cell(x,1);
time_list = zeros(x,1);
colors = cbrewer('qual','Set1',x,'PCHIP');
figure('WindowState','maximized');
if plotTime
    subplot('Position',[0.056470588235294 0.574718668268892 0.783340336134454 0.382728140241747])
end

%% Second Loop: Plot
q=0;
datList = {''};
p = cell(1,1);
% Loop through levels of theory
for i = numel(fields):-1:1
    
    M = size(Data.(fields{i}){1,2},1) - 1; % Number of different sub-methods
    N = size(Data.(fields{i}),1); % Number of different basis sets
    
    % Loop through sub-types
    for j=1:M
        
        % Loop through basis sets
        for k=1:N
            
            % Current Basis set
            Basis = Data.(fields{i}){k,1};
            
            % If current basis set matches desired, plot it
            datType = Data.(fields{i}){1,2}{j,1};
            if strcmp(Basis,Basis_set) %&& ~ismember(datType,datList) 
                
                q=q+1;
                datList{q} = datType;
                time_list(q) = Data.(fields{i}){k,3};
                % Current sub-type
                %temptext = strrep(fields{i},[Label '_'],'');
                %temptext2 = regexprep(temptext,'_(.)','\($1\)');
                %sub_type{q} = [temptext2 ' - ' Data.(fields{i}){1,2}{j,1}];
                sub_type{q} = [regexprep(regexprep(fields{i},[Label '_'],''),'_','-') ...
                    ' - ' Data.(fields{i}){k,2}{j,1}];
                
                if ismember(Data.(fields{i}){k,4},[1 2 3 4]) 
                    sub_type{q} = regexprep(sub_type{q},'HF','DFT');
                else
                end
                
                hold on
                p{q} = plot(Data.(fields{i}){k,2}{j,2},Data.(fields{i}){k,2}{j,3}.*Ha_kJ_NA,...
                    'Color',colors(q,:),'LineWidth',lw,...
                    'LineStyle',linetypes{randi([1 4],1,1)},...
                    'Visible','on');
            elseif strcmp(Basis,Basis_set) && ismember(datType,datList)
                [~,idx] = ismember(datType,datList);
                if time_list(idx) > Data.(fields{i}){k,3}
                    time_list(idx) = Data.(fields{i}){k,3};
                end
            end
        end   
    end
end

%% Plot TF potential
[TF_r, TF_U_LiF, TF_U_LiLi, TF_U_FF] = LiF_TF_Potential_Generator(false,10);
hold on
if strcmp(Label,'LiF')
    p{end+1} = plot(TF_r,TF_U_LiF,'Color','k','LineWidth',lw,'LineStyle','-',...
    'Visible','on');
    sub_type{end+1} = 'TF (Empirical)';
elseif strcmp(Label,'LiLi')
    p{end+1} = plot(TF_r,TF_U_LiLi,'Color','k','LineWidth',lw,'LineStyle','-',...
    'Visible','on');
    sub_type{end+1} = 'TF (Empirical)';
elseif strcmp(Label,'FF')
    p{end+1} = plot(TF_r,TF_U_FF,'Color','k','LineWidth',lw,'LineStyle','-',...
    'Visible','on');
    sub_type{end+1} = 'TF (Empirical)';
end

%% Plot JC potential
[JC_r, JC_U_LiF, JC_U_LiLi, JC_U_FF] = LiF_JC_Potential_Generator(false,10);
hold on
if strcmp(Label,'LiF')
    p{end+1} = plot(JC_r,JC_U_LiF,'Color','b','LineWidth',lw,'LineStyle','-',...
        'Visible','on');
    sub_type{end+1} = 'JC - SPC/E Water (Empirical)';
elseif strcmp(Label,'LiLi')
    p{end+1} = plot(JC_r,JC_U_LiLi,'Color','b','LineWidth',lw,'LineStyle','-',...
        'Visible','on');
    sub_type{end+1} = 'JC - SPC/E Water (Empirical)';
elseif strcmp(Label,'FF')
    hold on
    p{end+1} = plot(JC_r,JC_U_FF,'Color','b','LineWidth',lw,'LineStyle','-',...
        'Visible','on');
    sub_type{end+1} = 'JC - SPC/E Water (Empirical)';
end

%% Plot options
title(['Comparison of ' strrep(Label,'_','-') ' Potential Energy Surfaces with ' Basis_set ' Basis Set' ],...
    'Interpreter','latex','fontsize',fs)

set(gca,'box','on','TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel('LiF Bond Length (\AA)','fontsize',fs,'Interpreter','latex');
ylabel('Total Energy (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

legend1 = clickableLegend(sub_type,'Interpreter','latex');

if strcmp(Label,'LiF')
    ylim([-900 0]);
elseif strcmp(Label,'LiF_W')
    ylim([-1200 0]);
elseif strcmp(Label,'LiLi')
    ylim([0 1000]);
elseif strcmp(Label,'FF')
    ylim([0 1000]);
end

xlim([0 ep]);
    axis manual;
if plotTime
    set(legend1,'Position',[0.839285714285714 0.574468085106383 0.152310924369746 0.382978723404256]);
else
    set(gca,'Position',[0.0625 0.11318795430945 0.730042016806722 0.81619937694704]); %#ok<*UNRCH>
    set(legend1,'Position',[0.7984375 0.00311526479750779 0.199479166666667 0.995846313603323]);
end

% for i=1:length(p)-2
%     set(p{i},'Visible','off')
% end

%% Bar plot
if plotTime
    subplot('Position',[0.0567226890756303 0.0972644376899706 0.934873949579832 0.360688956433637]);
    hold on
    bar(time_list./(3600),'FaceColor','flat','CData',colors)
    set(gca, 'XTick',1:1:length(datList))
    set(gca,'xticklabel',datList)

    title(['Total CPU Walltime for ' strrep(Label,'_','-') ...
        ' Calculations with ' Basis_set ' Basis Set'],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','off','YMinorTick','off','FontSize',fs);
    xlabel('Theory Type','fontsize',fs,'Interpreter','latex');
    ylabel('CPU Time (hours)','fontsize',fs,'Interpreter','latex');
end
end
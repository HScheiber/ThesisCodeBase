%% Plot Options
fs = 20; % Font size
lw = 3; % Line width
ep = 10; % y axis end point in angstrom
%Basis_set = 'AUG-cc-pVTZ';

Basis_set = 'AUG-cc-pVTZ';
Label = 'LiF';
%linetypes = {'-' '--' ':' '-.'};
linetypes = {'-' '-' '-' '-'};
plotTime = 0; % 1 = include bar plot of CPU times, 0 = exclude it
%% Reference LiF_Data

% Unit Conversion
Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole


%% IMPORT DATA
%load('LiF_PES_Lattice.mat') % Creates variable called DATA

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
colors = cbrewer('qual','Set1',3,'PCHIP');
figure('units','normalized','outerposition',[0 0 1 1]); % New figure
if plotTime == 1
    subplot(2,1,1)
end

%% Second Loop: Plot
q=0;
datList = {''};
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
                sub_type{q} = [regexprep(regexprep(fields{i},'_','-'),'LiF-','') ...
                    ' - ' Data.(fields{i}){k,2}{j,1}];
                
                hold on
                plot(Data.(fields{i}){k,2}{j,2},Data.(fields{i}){k,2}{j,3}.*Ha_kJ*NA,...
                    'Color',colors(Data.(fields{i}){k,4},:),'LineWidth',lw,...
                    'LineStyle',linetypes{randi([1 4],1,1)})
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
[TF_r, TF_U_LiF,~,~] = LiF_TF_Potential_Generator(false,10);
plot(TF_r,TF_U_LiF,'Color','k','LineWidth',lw,'LineStyle','-')
sub_type{end+1} = 'TF (Empirical)';

%% Plot JC potential
[JC_r, JC_U_LiF,~,~] = LiF_JC_Potential_Generator(false,10);
plot(JC_r,JC_U_LiF,'Color','b','LineWidth',lw,'LineStyle','-')
sub_type{end+1} = 'JC - SPC/E Water (Empirical)';

%% Plot options
title(['Comparison of ' Label ' Calculations with ' Basis_set ' Basis Set' ],...
    'Interpreter','latex','fontsize',fs)

set(gca,'box','on','TickLabelInterpreter','latex');
set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
xlabel('LiF Bond Length (\AA)','fontsize',fs,'Interpreter','latex');
ylabel('Total Energy (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');

legend(sub_type)
ylim([-900 0]);
xlim([0 10]); 


%% Bar plot
if plotTime == 1
    subplot(2,1,2);
    bar(time_list./3600,'FaceColor','flat','CData',colors)
    set(gca,'xticklabel',datList)

    title(['Total CPU Time Accross 32 Cores for ' Label ' Calculations with ' Basis_set ' Basis Set'],...
        'Interpreter','latex','fontsize',fs)

    set(gca,'box','on','TickLabelInterpreter','latex');
    set(gca,'XMinorTick','off','YMinorTick','off','FontSize',fs);
    xlabel('Theory Type','fontsize',fs,'Interpreter','latex');
    ylabel('CPU Time (hours)','fontsize',fs,'Interpreter','latex');
end
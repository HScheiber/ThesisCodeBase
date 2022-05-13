%% Plot Options
fs = 18; % Font size
lw = 2; % Line width
ep = 8; % y axis end point in angstrom
plotTime=false;
Label = 'LiF_R';
home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';

%% Reference Data

% Unit Conversion
Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole


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

fields = fieldnames(Data);

% Loop through levels of theory
for i = 1:numel(fields)
    
    figure('WindowState','maximized'); % New figure
    if plotTime
        subplot('Position',[0.056470588235294 0.574718668268892 0.783340336134454 0.382728140241747])
    end
    
    M = size(Data.(fields{i}){1,2},1) - 1; % Number of different sub-methods
    N = size(Data.(fields{i}),1); % Number of different basis sets
    % Loop through sub-types
    for j=1:M
        
        % Current sub-type
        sub_type = Data.(fields{i}){1,2}{j,1};
        
        % Create subplot within current figure
%         if M > 1
%             subplot(2,round(M/2),j)
%         end
        
        % Loop through basis sets
        Basis = cell(N,1);
        colors = cbrewer('qual','Set1',N,'PCHIP');
        time_list = zeros(N,1);
        for k=1:N
            
            % Current Basis set
            Basis{k} = Data.(fields{i}){k,1};
            time_list(k) = Data.(fields{i}){k,3};
            
            hold on
            plot(Data.(fields{i}){k,2}{j,2},Data.(fields{i}){k,2}{j,3}.*Ha_kJ_NA,...
                'Color',colors(k,:),'LineWidth',lw)
        end
    

        %% Plot options
        titletext = strrep(fields{i},'_',' ');
        title([titletext ' Calc with ' sub_type ' Theory'],...
            'Interpreter','latex','fontsize',fs)
        
        set(gca,'box','on','TickLabelInterpreter','latex');
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
        xlabel('Unit Cell Parameter $a$ (\AA)','fontsize',fs,'Interpreter','latex');
        ylabel('Total Energy (kJ mol$^{-1}$)','fontsize',fs,'Interpreter','latex');
        
        legend1 = clickableLegend(Basis,'Location','SouthEast');
        %ylim([-1200 0]);
        xlim([0 ep]);    
        
        if plotTime
            set(legend1,'Position',[0.839285714285714 0.574468085106383 0.152310924369746 0.382978723404256]);
        else
            set(gca,'Position',[0.0625 0.11318795430945 0.730042016806722 0.81619937694704]); %#ok<*UNRCH>
            set(legend1,'Position',[0.7984375 0.00311526479750779 0.199479166666667 0.995846313603323]);
        end
        
    end
    
    %% Bar plot
    if plotTime
        subplot('Position',[0.0567226890756303 0.0972644376899706 0.934873949579832 0.360688956433637]);
        hold on
        bar(time_list./(3600),'FaceColor','flat','CData',colors)
        set(gca, 'XTick',1:1:length(Basis))
        set(gca,'xticklabel',Basis)

        title(['Total CPU Walltime for ' strrep(Label,'_','-') ...
            ' Calculations using ' strrep(strrep(fields{i},[Label '_'],''),'_','-') ' DFT'],...
            'Interpreter','latex','fontsize',fs)

        set(gca,'box','on','TickLabelInterpreter','latex');
        set(gca,'XMinorTick','off','YMinorTick','off','FontSize',fs);
        xlabel('Basis Set','fontsize',fs,'Interpreter','latex');
        ylabel('CPU Time (hours)','fontsize',fs,'Interpreter','latex');
    end
end


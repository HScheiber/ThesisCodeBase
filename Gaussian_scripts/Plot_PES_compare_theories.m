%% Plot Options
fs = 10; % Font size
lw = 2; % Line width
ep = 10; % y axis end point in angstrom


%% Reference Data

% Unit Conversion
Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole

%% IMPORT DATA
load('LiF_PES_Lattice_Energies.mat') % Creates variable called DATA

fields = fieldnames(LiF_Data);

%f = cell(numel(fields),1);

% Loop through levels of theory
for i = 1:numel(fields)
    
    figure('units','normalized','outerposition',[0 0 1 1]); % New figure
    
    % Loop through basis sets
    for j = 1:size(Data.(fields{i}),1)
        
        subplot(2, round(size(Data.(fields{i}),1)/2),j)
        
        Basis = Data.(fields{i}){j,1};
        
        % Loop through output types
        N = size(Data.(fields{i}){j,2},1);
        type = cell(N,1);
        colors = cbrewer('qual','Set1',N-1,'PCHIP');
        for k = 1:N
            type{k} = Data.(fields{i}){j,2}{k,1};
            
            if ~strcmp(type{k},'RMSD')
                hold on
                plot(Data.(fields{i}){j,2}{k,2},...
                    Data.(fields{i}){j,2}{k,3},...
                	'Color',colors(k,:),'LineWidth',lw)
            else
                type(k) = [];
            end
        end
        
        %% Plot options
        titletext = strrep(fields{i},'_',' ');
        title([titletext ' Calc with ' Basis ' Basis Set'],...
            'Interpreter','latex','fontsize',fs)
        
        set(gca,'box','on','TickLabelInterpreter','latex');
        set(gca,'XMinorTick','on','YMinorTick','on','FontSize',fs);
        xlabel('LiF Bond Length (\AA)','fontsize',fs,'Interpreter','latex');
        ylabel('Total Energy (E$_{H}$)','fontsize',fs,'Interpreter','latex');
        
        legend(type)
        ylim([-107.4 -106.6]);
        xlim([0 10]);
    end

end
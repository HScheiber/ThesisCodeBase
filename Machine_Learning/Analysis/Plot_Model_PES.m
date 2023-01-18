
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';
Settings = Initialize_MD_Settings;
Settings = Update_MD_Settings(Settings);
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; % 'LiF' 'LiCl' 'LiBr' 'LiI'

Settings.Theory = 'BH';
ModelID = 'MF';
PlotTypes = 'lj';
fs = 34; % font size
lw = 3; % line width
savefile = false; % switch to save the final plots to file
Startpoint = 0.08; % nm

figh = figure('WindowState','Maximized');
if strcmp(PlotTypes,'both')
    np = 4;
    plts = {'lj' 'full'};
else
    np = 2;
    plts = {PlotTypes};
end
T = tiledlayout(figh,2,np,'TileSpacing','compact','padding','tight');

for pts = 1:numel(plts)
    PlotType = plts{pts};
    for sidx = 1:numel(Salts)
        Settings.Salt = Salts{sidx};

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
        %ax = axes(figh,'Position',[0.100 0.1100 0.850 0.8150]); % [left bottom width height]
        ax = nexttile(T);
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

            if isempty(Models{idx}) && strcmp(Settings.Theory,'BH')
                U = TF_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'JC')
                U = JC_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
                disp([Settings.Salt ' ' ML_Model_Name ': Q = ' num2str(Settings.S.Q)])
            elseif strcmp(Settings.Theory,'LJ')
                U = LJ_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'TF')
                U = TF_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'BH')
                U = BH_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'BD')
                U = BD_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'BE')
                U = BE_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            elseif strcmp(Settings.Theory,'Mie')
                U = Mie_Potential_Generator(Settings,...
                    'Startpoint',Startpoint);
            end

            switch lower(PlotType)
                case 'full'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.f + U.MX.g + U.MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.f + U.MM.g + U.MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.f + U.XX.g + U.XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'full-derivative'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.df + U.MX.dg + U.MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.df + U.MM.dg + U.MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.df + U.XX.dg + U.XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'lj'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.g + U.MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.g + U.MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.g + U.XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'lj-derivative'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.dg + U.MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.dg + U.MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.dg + U.XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'dispersion'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.g,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.g,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.g,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'dispersion-derivative'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.dg,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'repulsive'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.h,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
                case 'repulsive-derivative'
                    h_MX(idx) = plot(ax,U.r.*10,U.MX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'LineStyle','-');
                    h_MM(idx) = plot(ax,U.r.*10,U.MM.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle','--');
                    h_XX(idx) = plot(ax,U.r.*10,U.XX.dh,'Color',Col_MX(idx,:),'LineWidth',lw,'Linestyle',':');
            end
        end

        switch lower(PlotType)
            case 'full'
                yl = [-800 1500];
                yt = -800:400:1400;
                ttxt = 'Full Potential';
            case 'full-derivative'
                yl = [-600 1000];
                ttxt = 'Derivative of Full Potential';
            case 'lj'
                yl = [-2.2 2];
                yt = -2:1:2;
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
        
        set(ax,'box','on','TickLabelInterpreter','latex');
        set(ax,'XMinorTick','on','YMinorTick','on','FontSize',fs);
        grid(ax,'on')
        grid(ax,'minor')
        xticks(ax,0:1:6);
        if pts == 1
            title(ax,Settings.Salt,'Interpreter','latex','fontsize',fs)
            xticklabels(ax,[])
            %xlabel(ax,'Separation [\AA]','fontsize',fs,'Interpreter','latex');
        end
        ylim(ax,yl);
        
        if pts == 1 && sidx == 1
            yticks(ax,yt)
            ylabel(ax,'$u_{ij}^{\textrm{LJ}}$ [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');
        elseif pts == 2 && sidx == 1
            yticks(ax,yt)
            ylabel(ax,'$u_{ij}^{\textrm{CLJ}}$ [kJ mol$^{-1}$]','fontsize',fs,'Interpreter','latex');
        else
            yticks(ax,yt)
            yticklabels(ax,[])
        end
        xlim([1 6.5]);

        % Blank line
        hline = refline([0 0]);
        hline.Color = 'k';
        hline.LineWidth = lw-1;
        hline.LineStyle = '-.';
    end
end


%ylabel(T,'$u_{ij}(r)$ [kJ mol$^{-1}$]','Interpreter','latex','fontsize',fs)
xlabel(T,'Separation $r$ [\AA]','Interpreter','latex','fontsize',fs)

if savefile
    Model_Labels = [{[Metal '$^{+}$' 'X$^{-}$']}  ...
        {[Metal '$^{+}$' Metal '$^{+}$']} ...
        {['X$^{-}$' 'X$^{-}$']}];
    legh = legend(ax,[h_MX(1); h_MM(1); h_XX(1)],Model_Labels,'Interpreter','latex','location','NorthOutside',...
        'NumColumns',3);
    legend('boxoff')
    
    legh.Layout.Tile = 'North';
    
    %set(figh,'renderer','opengl')
    filename = fullfile('C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Thesis\Thesis_Draft\BO_Figures',...
        [Settings.Theory '_Model_' ModelID '_PES.png']);
    exportgraphics(figh,filename,'resolution',600)
    
else
    title(T,['Plot of ' ttxt ' for ' Settings.Theory ' Models'],...
       'Interpreter','latex','fontsize',fs)
    
    Model_Labels = [{[Metal '$^{+}$' Halide '$^{-}$']}  ...
        {[Metal '$^{+}$' Metal '$^{+}$']} ...
        {[Halide '$^{-}$' Halide '$^{-}$']}];
    legh = legend(ax,h_MX(2:end),Models{2:end},'Interpreter','latex',...
        'NumColumns',1,'Location','EastOutside');
    legend('boxoff')
    
    legh.Layout.Tile = 'East';
%     
%     %% Copy the axes and plot the second legned
%     ah = copyobj( ax, figh);
%     ah.Visible = 'off';
%     lh = legend(ah, [h_MX(1); h_MM(1); h_XX(1)],Model_Labels,'Interpreter','latex',...
%         'NumColumns',1,'Location','SouthEast');
%     legend('boxoff')
    
end
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'BF';
ModelID = 'MK';
Reps = 1:5;
Show = []; % Sort by loss function and only show the lowest-loss results.
Show_init = [];
fs = 34; % font size
markers = {'o' 's' '^' 'v' 'd' '>' '<' 'p' 'h' 'x'};
show_as_C6 = false; % When true, plot parameters as B/C6/C8 than the sigma/epsilon parameter value itself
include_loss_panel = false;
loss_panel_height = 0.3;
plot_ref_model = true;
savefile = true; % switch to save the final plots to file

%% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';
% Find reps of models
N_Salts = numel(Salts);
Models = cell(1,N_Salts);
for idx = 1:N_Salts
    files = {dir(fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' ModelID '*_data.mat'])).name};
    Models{idx} = unique(regexp(files,[ModelID '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models{idx});
    Models{idx} = Models{idx}(nonempty_idx);
end
Models = unique(horzcat(Models{:}));
N_Reps = numel(Reps);
Models_of_Interest = cell(1,N_Reps);
for idx = 1:N_Reps
    Models_of_Interest{idx} = [ModelID num2str(Reps(idx))];
end
Models = intersect(Models_of_Interest,Models);
N_Models = numel(Models);

% Make sure the models are sorted
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end

%% Pre-allocate data arrays
switch Theory
    case 'BH'
        if show_as_C6
            ParNames = {'C_MX' 'C_MM' 'C_XX' 'B_MX' 'B_MM' 'B_XX' 'alpha_MX' 'alpha_MM' 'alpha_XX'};
            NPars = 9;
            NY = 2;
            NX = 5;
        else
            ParNames = {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'r0_MX' 'r0_MM' 'r0_XX' 'gamma_MX'};
            NPars = 7;
            NY = 2;
            NX = 4;
        end
    case 'BF'
        if show_as_C6
            ParNames = {'C_MX' 'C_MM' 'C_XX' 'B_MX' 'B_MM' 'B_XX' 'alpha_MX' 'alpha_MM' 'alpha_XX'};
            NPars = 9;
            NY = 2;
            NX = 5;
        else
            ParNames = {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'sigma_MX' 'sigma_MM' 'sigma_XX' 'gamma_MX'};
            NPars = 7;
            NY = 2;
            NX = 4;
        end
    case 'JC'
        if show_as_C6
            ParNames = {'C_MX' 'C_MM' 'C_XX' 'B_MX' 'B_MM' 'B_XX'};
            NPars = 6;
            NY = 2;
            NX = 3;
        else
            ParNames = {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'sigma_MX' 'sigma_MM' 'sigma_XX'};
            NPars = 6;
            NY = 2;
            NX = 3;
        end
        
    case 'Mie'
        if show_as_C6
            ParNames = {'C_MX' 'C_MM' 'C_XX' 'B_MX' 'B_MM' 'B_XX' 'n_MX'};
            NPars = 7;
            NY = 2;
            NX = 3;
        else
            ParNames = {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'sigma_MX' 'sigma_MM' 'sigma_XX' 'n_MX'};
            NPars = 7;
            NY = 2;
            NX = 3;
        end
end

Param_PlotData      = nan(N_Salts,N_Models,NPars);
RefModeData         = nan(N_Salts,NPars);
Total_loss          = nan(N_Salts,N_Models);

% Plot color scheme
Colours = cbrewer('qual','Set3',max(N_Salts,3));

% Load Data
for idx = 1:N_Salts
    Salt = Salts{idx};
    
    RefSettings = Initialize_MD_Settings;
    RefSettings.Salt = Salt;
    RefSettings.Additivity = false;
    Param = struct;
    if strcmp(Theory,'JC')
        [OutputMX,OutputMM,OutputXX] = JC_Potential_Parameters(RefSettings);
        
        % {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'sigma_MX' 'sigma_MM' 'sigma_XX'};
        % {'C_MX' 'C_MM' 'C_XX' 'A_MX' 'A_MM' 'A_XX'};
        
        Param.sigma_MX = OutputMX.sigma;
        Param.sigma_MM = OutputMM.sigma;
        Param.sigma_XX = OutputXX.sigma;
        Param.epsilon_MX = OutputMX.epsilon;
        Param.epsilon_MM = OutputMM.epsilon;
        Param.epsilon_XX = OutputXX.epsilon;
        
        Param.B_MX = 4*Param.epsilon_MX*(Param.sigma_MX^12);
        Param.B_MM = 4*Param.epsilon_MM*(Param.sigma_MM^12);
        Param.B_XX = 4*Param.epsilon_XX*(Param.sigma_XX^12);
        
        Param.C_MX = 4*Param.epsilon_MX*(Param.sigma_MX^6);
        Param.C_MM = 4*Param.epsilon_MM*(Param.sigma_MM^6);
        Param.C_XX = 4*Param.epsilon_XX*(Param.sigma_XX^6);
        
    elseif strcmp(Theory,'Mie')
        [OutputMX,OutputMM,OutputXX] = JC_Potential_Parameters(RefSettings);
        
        % {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'sigma_MX' 'sigma_MM' 'sigma_XX'};
        % {'C_MX' 'C_MM' 'C_XX' 'A_MX' 'A_MM' 'A_XX'};
        
        Param.sigma_MX = OutputMX.sigma;
        Param.sigma_MM = OutputMM.sigma;
        Param.sigma_XX = OutputXX.sigma;
        Param.epsilon_MX = OutputMX.epsilon;
        Param.epsilon_MM = OutputMM.epsilon;
        Param.epsilon_XX = OutputXX.epsilon;
        
        Param.B_MX = 4*Param.epsilon_MX*(Param.sigma_MX^12);
        Param.B_MM = 4*Param.epsilon_MM*(Param.sigma_MM^12);
        Param.B_XX = 4*Param.epsilon_XX*(Param.sigma_XX^12);
        
        Param.C_MX = 4*Param.epsilon_MX*(Param.sigma_MX^6);
        Param.C_MM = 4*Param.epsilon_MM*(Param.sigma_MM^6);
        Param.C_XX = 4*Param.epsilon_XX*(Param.sigma_XX^6);
        
        Param.n_MX = 12;
    else
        [OutputMX,OutputMM,OutputXX] = TF_Potential_Parameters(RefSettings);
        if show_as_C6
            ParNamesTF = ParNames;
        else
            ParNamesTF = {'epsilon_MX' 'epsilon_MM' 'epsilon_XX' 'r0_MX' 'r0_MM' 'r0_XX' 'gamma_MX'};
        end
        Param.C_MX = OutputMX.C;
        Param.C_MM = OutputMM.C;
        Param.C_XX = OutputXX.C;
        Param.B_MX = OutputMX.B;
        Param.B_MM = OutputMM.B;
        Param.B_XX = OutputXX.B;
        Param.alpha_MX = OutputMX.alpha;
        Param.alpha_MM = OutputMM.alpha;
        Param.alpha_XX = OutputXX.alpha;
        
        Param.gamma_MX = -7*lambertw((-1/7)*(6*Param.C_MX*(Param.alpha_MX^6)/Param.B_MX)^(1/7));
        Param.gamma_MM = -7*lambertw((-1/7)*(6*Param.C_MM*(Param.alpha_MM^6)/Param.B_MM)^(1/7));
        Param.gamma_XX = -7*lambertw((-1/7)*(6*Param.C_XX*(Param.alpha_XX^6)/Param.B_XX)^(1/7));
        
        Param.r0_MX = Param.gamma_MX/Param.alpha_MX;
        Param.r0_MM = Param.gamma_MM/Param.alpha_MM;
        Param.r0_XX = Param.gamma_XX/Param.alpha_XX;
        
        Param.epsilon_MX = Param.C_MX*(Param.gamma_MX - 6)/(Param.gamma_MX*(Param.r0_MX^6));
        Param.epsilon_MM = Param.C_MM*(Param.gamma_MM - 6)/(Param.gamma_MM*(Param.r0_MM^6));
        Param.epsilon_XX = Param.C_XX*(Param.gamma_XX - 6)/(Param.gamma_XX*(Param.r0_XX^6));
        
        if Param.gamma_MX < 6
            Param.epsilon_MX = -Param.epsilon_MX;
        end
        
        if Param.gamma_MM < 6
            Param.epsilon_MM = -Param.epsilon_MM;
        end
        
        if Param.gamma_XX < 6
            Param.epsilon_XX = -Param.epsilon_XX;
        end
    end
    for kdx = 1:NPars
        RefModeData(idx,kdx) = Param.(ParNamesTF{kdx});
    end
    
    for jdx = 1:N_Models
        Model = Models{jdx};
        
        % Find the fully optimized model
        dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);
        
        if ~isfile(dat_file)
            disp(['Could not load results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end
        % Load data
        try
            data = load(dat_file).full_data;
            
            if isfield(data,'secondary_result')
                
                optimvals = nan(1,length(data.secondary_result));
                for kdx = 1:length(data.secondary_result)
                    optimvals(kdx) = [data.secondary_result(kdx).optimValues.fval];
                end
                [Total_loss(idx,jdx),midx] = min(optimvals);
            else
                optimvals = data.bayesopt_results.ObjectiveTrace;
                [Total_loss(idx,jdx),midx] = min(optimvals);
                full_opt_point = data.bayesopt_results.XTrace{midx,:};
            end
            
            Param = data.full_opt_point;
            ParTable = Gen_Param_Table(data.Settings,Param);
            
            for kdx = 1:NPars
                Param_PlotData(idx,jdx,kdx) = ParTable.(ParNames{kdx});
            end
            
        catch
            disp(['Could not obtain data for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end
    end
end

% Plot the data
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
switch Theory
    case 'JC'
        PubTheoryName = 'CLJ';
    case 'Mie'
        PubTheoryName = 'Coulomb-Mie';
    case 'BH'
        PubTheoryName = 'CBH';
    case 'BF'
        PubTheoryName = 'Coulomb Wang-Buckingham';
    case 'TF'
        PubTheoryName = 'CBHM';
end

% Plot the loss data at the top
if include_loss_panel
    loss_panel = uipanel(figh,'FontSize',fs,...
        'BorderType','none',...
        'Position',[0 1-loss_panel_height 1 loss_panel_height],...
        'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
    axobj_loss = gobjects(N_Salts,1);

    ww0 = 1 - 0.04;
    bb0 = 0.05;
    ww1 = ww0/N_Salts - ww0/(N_Salts+0.5);
    ww2 = ww1/(max(N_Salts-1,1));
    ww = ww0/(N_Salts+0.5);
    hh = 1 - bb0 - 0.3;
    bb = bb0;
    for idx = 1:N_Salts
        ll = (1-ww0) + (ww0/N_Salts+ww2)*(idx-1);
        axobj_loss(idx) = subplot('Position',[ll bb ww hh],...
            'parent',loss_panel); % [left bottom width height]
        hold(axobj_loss(idx),'on');


        % Gather data
        plot_data = squeeze(Total_loss(idx,:));
        expon = floor(log10(max(plot_data)));

        p = bar(axobj_loss(idx),1:N_Models,plot_data./(10.^expon),'FaceColor',Colours(idx,:),'Visible','on','BarWidth',1,...
            'LineWidth',2);
        for mdx = 1:N_Models
            xpos = p.XEndPoints(mdx);
            scatter(axobj_loss(idx),xpos,plot_data(mdx)./(10.^expon),100,'MarkerEdgeColor','k',...
                'MarkerFaceColor',Colours(idx,:),'linewidth',2,'marker',markers{mdx});
        end

        % Set plot properties
        set(axobj_loss(idx),'box','on','TickLabelInterpreter','latex');
        set(axobj_loss(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
        xticks(axobj_loss(idx),1:N_Models)
        xticklabels(axobj_loss(idx),[])
        axobj_loss(idx).XAxis.TickLength = [0,0];
        ylim(axobj_loss(idx),'padded')
        if ~savefile
            sgtitle(loss_panel,['Minimized Objective Function: ' PubTheoryName ' Model ' ModelID],...
                'Fontsize',fs-2,'Interpreter','latex')
        end
        xlim(axobj_loss(idx),[0 N_Models+1])
        set(axobj_loss(idx), 'YTickMode', 'auto');
        set(axobj_loss(idx),'xminorgrid','off','yminorgrid','on');
        axobj_loss(idx).YAxis.Exponent = 0;
        ylabel(axobj_loss(idx), ['$f\left(\mathbf{x}^{*}\right)/10^{' num2str(expon) '}$'],'Fontsize',fs-4,'Interpreter','latex')
    end
else
    loss_panel_height=0;
end

% Plot parameters
bspace = 0.05;
par_panel = uipanel(figh,'FontSize',fs,...
    'BorderType','none',...
    'Position',[0 0 1 1-loss_panel_height],...
    'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]

T = tiledlayout(par_panel,'flow', 'Padding', 'loose', 'TileSpacing', 'loose'); 
axh = gobjects(1,NPars);
xx = 1:N_Salts;
for idx = 1:NPars
    axh(idx) = nexttile(T);
    hold(axh(idx),'on');
    
    % Gather data
    plot_data = squeeze(Param_PlotData(:,:,idx));

    for kdx = 1:N_Models
        yy = plot_data(:,kdx);
        scatter(axh(idx),xx,yy,100,Colours(xx,:),'filled','MarkerEdgeColor','k',...
            'linewidth',1,'marker',markers{kdx});
    end
    
    if plot_ref_model
        yy = RefModeData(:,idx);
        scatter(axh(idx),xx,yy,200,'r','filled','MarkerEdgeColor','k',...
            'linewidth',3,'marker','x');
    end
    
    xlim(axh(idx),[0.5 N_Salts + 0.5])
    ylim(axh(idx),'padded')
    xticks(axh(idx),[])
    
    [name,units,logscale] = param_name_map(ParNames{idx});
    ylabel(axh(idx),units,'interpreter','latex','fontsize',fs);
    xlabel(axh(idx),name,'interpreter','latex','fontsize',fs);
    set(axh(idx),'box','on','TickLabelInterpreter','latex',...
        'XMinorGrid','off','YMinorGrid','on','YGrid','on',...
        'FontSize',fs)
    if logscale
        set(axh(idx),'YScale','log')
    end
end


h = bar(gca,nan,nan(1,N_Salts),'FaceColor','flat','Visible','on','BarWidth',1,...
    'LineWidth',2);
for kdx = 1:N_Salts
    h(kdx).CData = Colours(kdx,:);
end
legtxt = strcat({'\ \ '},Salts,{'\quad'});
legh = legend(h,legtxt,'Orientation','Horizontal',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts);
legh.Layout.Tile = 'south';

if savefile
    filename = fullfile('C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Thesis\Thesis_Draft\BO_Figures',...
        [Theory '_Model_' ModelID '_Params.png']);
    print(figh,filename,'-dpng','-r0','-noui')
end

function [name,units,logscale] = param_name_map(p_name)
    logscale = false;
    switch p_name
        case 'C_MX'
            name = 'C$_{6,\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$ nm$^{6}$]';
        case 'C_MM'
            name = 'C$_{6,\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[kJ mol$^{-1}$ nm$^{6}$]';
        case 'C_XX'
            name = 'C$_{6,\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$ nm$^{6}$]';
        case 'A_MX'
            name = 'B$_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$ nm$^{12}$]';
        case 'A_MM'
            name = 'B$_{\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[kJ mol$^{-1}$ nm$^{12}$]';
        case 'A_XX'
            name = 'B$_{\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$ nm$^{12}$]';
        case 'B_MX'
            name = 'B$_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$]';
        case 'B_MM'
            name = 'B$_{\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[kJ mol$^{-1}$]';
        case 'B_XX'
            name = 'B$_{\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$]';
        case 'n_MX'
            name = 'Born Exponent $n$';
            units = '';
        case 'epsilon_MX'
            name = '$\epsilon_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$]';
            logscale = true;
        case 'epsilon_MM'
            name = '$\epsilon_{\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[kJ mol$^{-1}$]';
            logscale = true;
        case 'epsilon_XX'
            name = '$\epsilon_{\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[kJ mol$^{-1}$]';
            logscale = true;
        case 'sigma_MX'
            name = '$\sigma_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[nm]';
        case 'sigma_MM'
            name = '$\sigma_{\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[nm]';
        case 'sigma_XX'
            name = '$\sigma_{\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[nm]';
        case 'alpha_MX'
            name = '$\alpha_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[nm$^{-1}$]';
        case 'alpha_MM'
            name = '$\alpha_{\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[nm$^{-1}$]';
        case 'alpha_XX'
            name = '$\alpha_{\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[nm$^{-1}$]';
        case 'r0_MX'
            name = '$r_{0,\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '[nm]';
        case 'r0_MM'
            name = '$r_{0,\textrm{Li}^{+}\textrm{Li}^{+}}$';
            units = '[nm]';
        case 'r0_XX'
            name = '$r_{0,\textrm{X}^{-}\textrm{X}^{-}}$';
            units = '[nm]';
        case 'gamma_MX'
            name = '$\gamma_{\textrm{Li}^{+}\textrm{X}^{-}}$';
            units = '';
        otherwise
            name = p_name;
    end
end
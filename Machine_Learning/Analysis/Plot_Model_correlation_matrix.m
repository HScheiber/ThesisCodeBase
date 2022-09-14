% Script for plotting bayesian optimization results: model parameters correlated
% with various target properties
Salt = 'LiI';
Theory = 'BH';
Model = 'MC1';
Show = []; % Sort by loss function and only show the lowest-loss results.
Show_init = [];
fs = 10; % fontsize
scat_size = 50;
show_as_percent_error = true; % Plot as percent error. Does not apply to 'Loss'. Plot results as error w.r.t. DFT results
show_as_scale = true; % When true, plot parameters as the scaling w.r.t. default model rather than the parameter value itself
show_additivity = false; % When true, shows the distribution of cross interactions
show_as_sigma_epsilon = false; % When true, transform JC parameters into sigma-epsilon form

% Properties to correlate with model parameters
Structure = 'Rocksalt';
% Possible properties to choose from:
% 'CoulSR'   'CoulLR'   'Pot'   'CoulSR_MM'   'CoulSR_MX'   'CoulSR_XX'
% 'a'   'b'   'c'   'E'   'c_over_a' 'Loss' (Loss may depend on multiple structures)
Property = 'Loss';

%% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
model_filename = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);

if ~isfile(model_filename)
    error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
end

if strcmp(Theory,'TF') % Cannot show additivity for TF model
    show_additivity = false;
end

% Load data
data = load(model_filename).full_data;
Settings = data.Settings;
results = data.bayesopt_results;
params = bayesopt_params(data.Settings);
clear('data')

if isempty(Show_init)
    % Get search points
    search_pnts = results.XTrace;
    losses = results.ObjectiveTrace;

    % Get user data: crystal properties at each point
    Crystal_data_full = results.UserDataTrace;
else
    % Get search points
    search_pnts = results.XTrace(1:Show_init,:);
    losses = results.ObjectiveTrace(1:Show_init);

    % Get user data: crystal properties at each point
    Crystal_data_full = results.UserDataTrace(1:Show_init);
end
lossnan = isnan(losses);
losses(lossnan) = [];
search_pnts_nan =  search_pnts(lossnan,:);
search_pnts(lossnan,:) = [];

% Get the Structures (same order at each point)
N = length(Crystal_data_full{1});
Structures = cell(1,N);
for idx = 1:N
    Structures{idx} = Crystal_data_full{1}.Minimization_Data{idx}.Structure;
end

% Get the property of interest
M = length(Crystal_data_full);
structure_idx = strcmpi(Structure,Structures);
Crystal_data = nan(1,M);

switch lower(Property)
    case 'c_over_a'
        for idx = 1:M
            Crystal_data(idx) = Crystal_data_full{idx}.Minimization_Data{structure_idx}.c / Crystal_data_full{idx}.Minimization_Data{structure_idx}.a;
        end
        Heatmap_Label = [Structure ' c/a'];
        clog = false;
    case 'loss'
        Crystal_data = losses;
        Heatmap_Label = 'Loss';
        clog = false;
    otherwise
        for idx = 1:M
            Crystal_data(idx) = Crystal_data_full{idx}.Minimization_Data{structure_idx}.(Property);
        end
        Heatmap_Label = [Structure ' ' Property];
        
        switch Property
            case {'CoulSR' 'CoulLR' 'Pot' 'CoulSR_MM' 'CoulSR_MX' 'CoulSR_XX' 'E'}
                Heatmap_Label = [Heatmap_Label newline '[kJ mol$^{-1}$]'];
            case {'a' 'b' 'c'}
                Heatmap_Label = [Heatmap_Label newline '[\AA]'];
        end
        clog = false;
end

% Sort everything by the loss function
[losses_srt,idx] = sort(losses);
search_pnts_srt = search_pnts(idx,:);
Crystal_data_srt = Crystal_data(idx);

% Extract only best M points
if isempty(Show)
    losses_srt_cut = losses_srt;
    search_pnts_cut = search_pnts_srt;
    Crystal_data_cut = Crystal_data_srt;
else
    losses_srt_cut = losses_srt(1:Show);
    search_pnts_cut = search_pnts_srt(1:Show,:);
    Crystal_data_cut = Crystal_data_srt(1:Show);
end

% Generate the MX parameters by additivity if they aren't present
if ~contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames) && show_additivity
    
    PotSettings = Initialize_MD_Settings;
    PotSettings.Salt = Salt;
    PotSettings.S = Init_Scaling_Object;
    [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

    % Unscaled
    MX_Epsilon = MXParams.epsilon;
    MX_Sigma   = MXParams.sigma;

    MX_R = 4*MX_Epsilon*MX_Sigma^12;
    MX_D = 4*MX_Epsilon*MX_Sigma^6;

    % Scaled
    MM_Epsilon = MMParams.epsilon.*(search_pnts_cut.SDMM.^2).*(1./search_pnts_cut.SRMM);
    MM_Sigma = MMParams.sigma.*(1./(search_pnts_cut.SDMM.^(1/6))).*(search_pnts_cut.SRMM.^(1/6));

    XX_Epsilon = XXParams.epsilon.*(search_pnts_cut.SDXX.^2).*(1./search_pnts_cut.SRXX);
    XX_Sigma = XXParams.sigma.*(1./(search_pnts_cut.SDXX.^(1/6))).*(search_pnts_cut.SRXX.^(1/6));

    MX_Epsilon = sqrt(MM_Epsilon.*XX_Epsilon);
    MX_Sigma   = (MM_Sigma + XX_Sigma)./2;

    MX_R_scaled = 4.*MX_Epsilon.*MX_Sigma.^12;
    MX_D_scaled = 4.*MX_Epsilon.*MX_Sigma.^6;
    
    search_pnts_cut.SDMX = MX_D_scaled./MX_D;
    search_pnts_cut.SRMX = MX_R_scaled./MX_R;
end

% Add the non-addivitive M-M dispersion force back in.
if contained_in_cell('SDMM2',search_pnts_cut.Properties.VariableNames)
    search_pnts_cut.SDMM = search_pnts_cut.SDMM + search_pnts_cut.SDMM2;
    search_pnts_cut.SDMM2 = [];
end

if ~show_as_scale && strcmp(Theory,'JC')
    
    PotSettings = Initialize_MD_Settings;
    PotSettings.Salt = Salt;
    PotSettings.S = Init_Scaling_Object;
    [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

    % Unscaled
    MX_Epsilon = MXParams.epsilon;
    MX_Sigma   = MXParams.sigma;

    MM_R = 4*(MMParams.epsilon)*MMParams.sigma^12;
    MM_D = 4*(MMParams.epsilon)*MMParams.sigma^6;
    
    XX_R = 4*(XXParams.epsilon)*XXParams.sigma^12;
    XX_D = 4*(XXParams.epsilon)*XXParams.sigma^6;
    
    MX_R = 4*MX_Epsilon*MX_Sigma^12;
    MX_D = 4*MX_Epsilon*MX_Sigma^6;
    
    search_pnts_cut.DMM = MM_D.*search_pnts_cut.SDMM;
    search_pnts_cut.DXX = XX_D.*search_pnts_cut.SDXX;
    search_pnts_cut.RMM = MM_R.*search_pnts_cut.SRMM;
    search_pnts_cut.RXX = XX_R.*search_pnts_cut.SRXX;
    
    if contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.DMX = MX_D.*search_pnts_cut.SDMX;
        search_pnts_cut.RMX = MX_R.*search_pnts_cut.SRMX;
    end
    
elseif strcmp(Theory,'BH')
    
    if show_as_scale
        Settings = Get_Scaling_Params(Settings,search_pnts_cut);
        search_pnts_cut.SDMM = Settings.S.D.MM;
        search_pnts_cut.SDXX = Settings.S.D.XX;
        search_pnts_cut.SRMM = Settings.S.R.MM;
        search_pnts_cut.SRXX = Settings.S.R.XX;
        search_pnts_cut.SAMX = Settings.S.A.MX;
        search_pnts_cut.r0_MM = [];
        search_pnts_cut.r0_XX = [];
        search_pnts_cut.epsilon_MM = [];
        search_pnts_cut.epsilon_XX = [];
        search_pnts_cut.gamma_MX = [];
        
        Settings = Get_Scaling_Params(Settings,search_pnts_nan);
        search_pnts_nan.SDMM = Settings.S.D.MM;
        search_pnts_nan.SDXX = Settings.S.D.XX;
        search_pnts_nan.SRMM = Settings.S.R.MM;
        search_pnts_nan.SRXX = Settings.S.R.XX;
        search_pnts_nan.SAMX = Settings.S.A.MX;
        search_pnts_nan.r0_MM = [];
        search_pnts_nan.r0_XX = [];
        search_pnts_nan.epsilon_MM = [];
        search_pnts_nan.epsilon_XX = [];
        search_pnts_nan.gamma_MX = [];
        
    end
end

if show_as_sigma_epsilon && strcmp(Theory,'JC')
    
    PotSettings = Initialize_MD_Settings;
    PotSettings.Salt = Salt;
    PotSettings.S = Init_Scaling_Object;
    [MXParams,MMParams,XXParams] = JC_Potential_Parameters(PotSettings);

    % Unscaled
    MX_Epsilon = MXParams.epsilon;
    MX_Sigma   = MXParams.sigma;
    
    % Scaled
    MM_Epsilon_scaled = MMParams.epsilon.*(search_pnts_cut.SDMM.^2).*(1./search_pnts_cut.SRMM);
    MM_Sigma_scaled = MMParams.sigma.*(1./(search_pnts_cut.SDMM.^(1/6))).*(search_pnts_cut.SRMM.^(1/6));

    XX_Epsilon_scaled = XXParams.epsilon.*(search_pnts_cut.SDXX.^2).*(1./search_pnts_cut.SRXX);
    XX_Sigma_scaled = XXParams.sigma.*(1./(search_pnts_cut.SDXX.^(1/6))).*(search_pnts_cut.SRXX.^(1/6));
    
    if ~contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames)
        MX_Epsilon_scaled = sqrt(MM_Epsilon.*XX_Epsilon);
        MX_Sigma_scaled   = (MM_Sigma + XX_Sigma)./2;
    else
        MX_Epsilon_scaled = MX_Epsilon.*(search_pnts_cut.SDMX.^2).*(1./search_pnts_cut.SRMX);
        MX_Sigma_scaled = MX_Sigma.*(1./(search_pnts_cut.SDMX.^(1/6))).*(search_pnts_cut.SRMX.^(1/6));
    end
    
    if show_as_scale
        search_pnts_cut.SSMM = MM_Sigma_scaled / MMParams.sigma;
        search_pnts_cut.SEMM = MM_Epsilon_scaled / MMParams.epsilon;
        
        search_pnts_cut.SSXX = XX_Sigma_scaled / XXParams.sigma;
        search_pnts_cut.SEXX = XX_Epsilon_scaled / XXParams.epsilon;
        
        if (~contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames) && show_additivity) || ...
                contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames)
            search_pnts_cut.SSMX = MX_Sigma_scaled ./ MX_Sigma;
            search_pnts_cut.SEMX = MX_Epsilon_scaled ./ MX_Epsilon;
        end
    else
        search_pnts_cut.sigma_MM = MM_Sigma_scaled;
        search_pnts_cut.epsilon_MM = MM_Epsilon_scaled;
        
        search_pnts_cut.sigma_XX = XX_Sigma_scaled;
        search_pnts_cut.epsilon_XX = XX_Epsilon_scaled;
        
        if (~contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames) && show_additivity ) || ...
                contained_in_cell('SDMX',search_pnts_cut.Properties.VariableNames)
            search_pnts_cut.sigma_MX = MX_Sigma_scaled;
            search_pnts_cut.epsilon_MX = MX_Epsilon_scaled;
        end
    end
elseif show_as_sigma_epsilon && strcmp(Theory,'TF')
    error('Cannot show TF model in Sigma-Epsilon form.')
end

% Remove entires that are not used
if show_as_sigma_epsilon
    
    % Remove old table entries
    search_pnts_cut.SDMM = [];
    search_pnts_cut.SDXX = [];
    search_pnts_cut.SRMM = [];
    search_pnts_cut.SRXX = [];
    if ismember('SDMX',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.SDMX = [];
        search_pnts_cut.SRMX = [];
    end
    
    if ismember('DMM',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.DMM = [];
        search_pnts_cut.DXX = [];
        search_pnts_cut.RMM = [];
        search_pnts_cut.RXX = [];
    end
    
    if ismember('DMX',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.DMX = [];
        search_pnts_cut.RMX = [];
    end
    
elseif ~show_as_scale && strcmp(Theory,'JC')
    search_pnts_cut.SDMM = [];
    search_pnts_cut.SDXX = [];
    search_pnts_cut.SRMM = [];
    search_pnts_cut.SRXX = [];
    if ismember('SDMX',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.SDMX = [];
        search_pnts_cut.SRMX = [];
    end
    
elseif ~show_as_scale && strcmp(Theory,'TF')
    search_pnts_cut.SD6MM = [];
    search_pnts_cut.SD6XX = [];
    search_pnts_cut.SD6MX = [];
    search_pnts_cut.SD8MM = [];
    search_pnts_cut.SD8XX = [];
    search_pnts_cut.SD8MX = [];
    search_pnts_cut.SAMM = [];
    search_pnts_cut.SAXX = [];
    search_pnts_cut.SAMX = [];
    search_pnts_cut.SRMM = [];
    search_pnts_cut.SRXX = [];
    search_pnts_cut.SRMX = [];
end


% Gather variable names
[~,N] = size(search_pnts_cut);
param_names = search_pnts_cut.Properties.VariableNames;
param_transform = {params.Transform};
pnames = cell(size(param_names));
for idx = 1:N
    switch param_names{idx}
        case 'SDMM'
            pnames{idx} = 'M-M Disp.';
        case 'SDMM2'
            pnames{idx} = 'M-M Added Disp.';
        case 'SDXX'
            pnames{idx} = 'X-X Disp.';
        case 'SDMX'
            pnames{idx} = 'M-X Disp.';
        case 'SRMM'
            pnames{idx} = 'M-M Rep.';
        case 'SRXX'
            pnames{idx} = 'X-X Rep.';
        case 'SRMX'
            pnames{idx} = 'M-X Rep.';
        case 'SQ'
            pnames{idx} = 'Charge';
        case 'SD6MM'
            pnames{idx} = 'M-M $r^{6}$ Disp.';
        case 'SD6XX'
            pnames{idx} = 'X-X $r^{6}$ Disp.';
        case 'SD6MX'
            pnames{idx} = 'M-X $r^{6}$ Disp.';
        case 'SD8MM'
            pnames{idx} = 'M-M $r^{8}$ Disp.';
        case 'SD8XX'
            pnames{idx} = 'X-X $r^{8}$ Disp.';
        case 'SD8MX'
            pnames{idx} = 'M-X $r^{8}$ Disp.';
        case 'SAMM'
            pnames{idx} = 'M-M Rep. Exp.';
        case 'SAXX'
            pnames{idx} = 'X-X Rep. Exp.';
        case 'SAMX'
            pnames{idx} = 'M-X Rep. Exp.';
        case 'sigma_MM'
            pnames{idx} = '$\sigma_{MM}$ [nm]';
        case 'epsilon_MM'
            pnames{idx} = '$\epsilon_{MM}$ [kJ mol$^{-1}$]';
        case 'sigma_XX'
            pnames{idx} = '$\sigma_{XX}$ [nm]';
        case 'epsilon_XX'
            pnames{idx} = '$\epsilon_{XX}$ [kJ mol$^{-1}$]';
        case 'sigma_MX'
            pnames{idx} = '$\sigma_{MX}$ [nm]';
        case 'epsilon_MX'
            pnames{idx} = '$\epsilon_{MX}$ [kJ mol$^{-1}$]';
        case 'SSMM'
            pnames{idx} = '$\sigma_{MM}$';
        case 'SEMM'
            pnames{idx} = '$\epsilon_{MM}$';
        case 'SSXX'
            pnames{idx} = '$\sigma_{XX}$';
        case 'SEXX'
            pnames{idx} = '$\epsilon_{XX}$';
        case 'SSMX'
            pnames{idx} = '$\sigma_{MX}$';
        case 'SEMX'
            pnames{idx} = '$\epsilon_{MX}$';
        case 'D6MM'
            pnames{idx} = 'M-M $r^{6}$ Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'D6XX'
            pnames{idx} = 'X-X $r^{6}$ Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'D6MX'
            pnames{idx} = 'M-X $r^{6}$ Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'D8MM'
            pnames{idx} = 'M-M $r^{8}$ Disp. [kJ mol$^{-1}$ nm$^{8}$]';
        case 'D8XX'
            pnames{idx} = 'X-X $r^{8}$ Disp. [kJ mol$^{-1}$ nm$^{8}$]';
        case 'D8MX'
            pnames{idx} = 'M-X $r^{8}$ Disp. [kJ mol$^{-1}$ nm$^{8}$]';
        case 'AMM'
            pnames{idx} = 'M-M Rep. Exp. [nm$^{-1}$]';
        case 'AXX'
            pnames{idx} = 'X-X Rep. Exp. [nm$^{-1}$]';
        case 'AMX'
            pnames{idx} = 'M-X Rep. Exp. [nm$^{-1}$]';
        case 'RMM'
            pnames{idx} = 'M-M Rep. [kJ mol$^{-1}$]';
        case 'RXX'
            pnames{idx} = 'X-X Rep. [kJ mol$^{-1}$]';
        case 'RMX'
            pnames{idx} = 'M-X Rep. [kJ mol$^{-1}$]';
        case 'DMM'
            pnames{idx} = 'M-M Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'DXX'
            pnames{idx} = 'X-X Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'DMX'
            pnames{idx} = 'M-X Disp. [kJ mol$^{-1}$ nm$^{6}$]';
        case 'r0_MM'
            pnames{idx} = 'M-M $r_{0}$ [nm]';
        case 'r0_XX'
            pnames{idx} = 'X-X $r_{0}$ [nm]';
        case 'r0_MX'
            pnames{idx} = 'M-X $r_{0}$ [nm]';
        case 'gamma_MM'
            pnames{idx} = 'M-M $\gamma$';
        case 'gamma_XX'
            pnames{idx} = 'X-X $\gamma$';
        case 'gamma_MX'
            pnames{idx} = 'M-X $\gamma$';
            
            
        otherwise
            tname = strrep(param_names{idx},'_','-');
            tname = strrep(tname,'GA-','G. Dep. ');
            tname = strrep(tname,'GB-','G. Pos. ');
            tname = strrep(tname,'GC-','G. Wid. ');
            tname = strrep(tname,'-1','');
            
            if ~isempty(regexp(param_names{idx},'GA','once'))
                tname = [tname ' [kJ mol$^{-1}$]'];
            elseif ~isempty(regexp(param_names{idx},'G(B|C)','once'))
                tname = [tname ' [nm]'];
            end
            
            pnames{idx} = tname;
    end
end

% Plot scatter matrix
figh = figure('WindowState','Maximized');

% Make colormap
if show_as_percent_error && ~strcmpi(Property,'loss')
    % Make sure the source map has an odd number of values
    map = flipud(max(min(cbrewer('div','RdBu',11,'spline'),1),0));
    
    % Center and expand color map
    map = center_colormap(map,Crystal_data_cut,0);
    
    colormap(figh,map)
else
    colormap(figh,flipud(hot))
end

x = 0;
P = gobjects(N,N);
T = tiledlayout(figh,N,N, 'Padding', 'loose', 'TileSpacing', 'tight'); 
for idx = 1:N
    for jdx = 1:N
        x = x+1;
        nexttile;
        if idx == jdx
            
            % Do histogram  
            if strcmp(param_transform{idx},'log') || show_as_scale
                [~,edges] = histcounts(log(search_pnts_cut{:,idx}),50);
                P(idx,jdx) = histogram(search_pnts_cut{:,idx},exp(edges));
                
            else
                P(idx,jdx) = histogram(search_pnts_cut{:,idx},50);
            end
            
            %P(idx,jdx) = scatter(search_pnts_cut{:,idx},losses_srt_cut);
            cur_ax = P(idx,jdx).Parent;
            
            if idx ~= N
                cur_ax.XTickLabel = [];
            else
                xlabel(cur_ax,pnames{jdx},'Interpreter','latex','FontSize',fs);
            end
            if jdx ~= 1
                cur_ax.YTickLabel = [];
            else
                ylabel(cur_ax,pnames{idx},'Interpreter','latex','FontSize',fs);
            end
            
            % some axis properties
            set(cur_ax,'XGrid','On','YGrid','On','GridLineStyle','-','Layer','Top',...
                'TickLength',[0 0],'FontSize',fs,'TickLabelInterpreter','latex')
            
            if strcmp(param_transform{idx},'log') || show_as_scale
                set(cur_ax, 'XScale', 'log')
            end
            
        else
            % do scatter plot
            X = search_pnts_cut{:,jdx};
            Y = search_pnts_cut{:,idx};
            
            Xnan = search_pnts_nan{:,jdx};
            Ynan = search_pnts_nan{:,idx};
            
            P(idx,jdx) = scatter(X,Y,10,Crystal_data_cut,'filled','SizeData',scat_size,...
                'MarkerEdgeColor','k');
            cur_ax = P(idx,jdx).Parent;
            hold(cur_ax,'on')
            if isempty(Show)
                scatter(Xnan,Ynan,10,'g','filled','SizeData',scat_size/5,...
                    'MarkerEdgeColor','k');
            end
            
            if idx ~= N
                cur_ax.XTickLabel = [];
            else
                xlabel(cur_ax,pnames{jdx},'Interpreter','latex','FontSize',fs);
            end
            if strcmp(param_transform{jdx},'log') || show_as_scale
                set(cur_ax, 'XScale', 'log')
            end
            
            if jdx ~= 1
                cur_ax.YTickLabel = [];
            else
                ylabel(cur_ax,pnames{idx},'Interpreter','latex','FontSize',fs);
            end
            if strcmp(param_transform{idx},'log') || show_as_scale
                set(cur_ax, 'YScale', 'log')
            end
            
            % some axis properties
            set(cur_ax,'XGrid','On','YGrid','On','GridLineStyle','-','Layer','Top',...
                'FontSize',fs,'TickLabelInterpreter','latex')
        end
        if clog
            set(cur_ax,'ColorScale','log')
        end
    end
end
drawnow

% Set the proper X-limits on the histograms
for idx = 1:N
    ax11 = P(idx,idx).Parent;
    
    if idx == 1
        ax21 = P(2,idx).Parent;
    else
        ax21 = P(idx-1,idx).Parent;
    end
    
    axlim_21 = ax21.XLim;
    ax11.XLim = axlim_21;
end

drawnow

% Set pseudo Y-ticks on the histograms
for idx = 1:N
    ax11 = P(idx,idx).Parent;
    ax11.YGrid='off';
%     if ~strcmp(param_transform{idx},'log')
%     
%         ax11 = P(idx,idx).Parent;
% 
%         if idx == 1
%             ax12 = P(idx,2).Parent;
%         else
%             ax12 = P(idx,idx-1).Parent;
%         end
% 
%         axlabs_12 = ax12.YTick;
%         axlim_11 = ax11.YLim;
% 
%         ax11.YTick = linspace(axlim_11(1),axlim_11(end),length(axlabs_12));
% 
%         bot_scl = ax12.YLim(1)./ax12.YTick(1);
%         if isnan(bot_scl)
%             bot_scl = 1;
%         end
%         top_scl = ax12.YLim(2)./ax12.YTick(end);
%         if isnan(top_scl)
%             top_scl = 1;
%         end
% 
%         ax11.YLim = ax11.YLim.*[bot_scl top_scl];
%     else
%         ax11 = P(idx,idx).Parent;
%         ax11.YGrid='off';
%     end
end

% Correct the Y-label at the 1,1 position
ax11 = P(1,1).Parent;
ax12 = P(1,2).Parent;

axlabs_12 = ax12.YTick;
axlabs_12_txt = arrayfun(@num2str,axlabs_12,'UniformOutput',false);
axlabs_11 = ax11.YTick;

ax11.YTick = linspace(axlabs_11(1),axlabs_11(end),length(axlabs_12));
ax11.YTickLabel = axlabs_12_txt;


% Add colorbar
sig = 2;
mincb = round(min(Crystal_data_cut),sig,'Significant');
maxcb = round(max(Crystal_data_cut),sig,'Significant');
while mincb == maxcb
    sig = sig+1;
    mincb = round(min(Crystal_data_cut),sig,'Significant');
    maxcb = round(max(Crystal_data_cut),sig,'Significant');
end

if show_as_percent_error && ~strcmpi(Property,'loss')
    abs_max = max([abs(mincb) abs(maxcb)]);
    mincb = -abs_max;
    maxcb = abs_max;
end

cb = colorbar(P(1,2).Parent,'Position',[0.935 0.11 0.02 0.815],...
    'TickLabelInterpreter','latex','FontSize',fs,'Limits',[mincb maxcb]); % [left, bottom, width, height]

cb.Title.String = Heatmap_Label;
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = fs+5;

% Add title
title(T,[Salt ' ' Theory ' Model ' Model ': Parameter Correlation Matrix'],...
    'Interpreter','Latex','FontSize',fs+16);
xlabel(T,'Parameter 1','Interpreter','Latex','FontSize',fs+5)
ylabel(T,'Parameter 2','Interpreter','Latex','FontSize',fs+5)
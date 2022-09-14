% Script for plotting bayesian optimization results as a correlation matrix
% Converts JC parameters to sigma/epsilon form
Salt = 'LiI';
Model = 'BH';
Trial = 'MC4';
Show = []; % leave empty to show all
Show_init = [];
fs = 14;
show_as_percent_error = true; % Plot as percent error. Does not apply to 'Loss'. Plot results as error w.r.t. DFT results
show_as_scale = true; % When true, plot parameters as the scaling w.r.t. default model rather than the parameter value itself
show_additivity = false; % When true, shows the distribution of cross interactions
show_as_sigma_epsilon = false; % When true, transform JC parameters into sigma-epsilon form

% Possible properties to choose from:
% 'CoulSR'   'CoulLR'   'Pot'   'CoulSR_MM'   'CoulSR_MX'   'CoulSR_XX'
% 'a'   'b'   'c'   'E'   'c_over_a' 'Loss' (Loss may depend on multiple structures)
Y_Structure = 'Rocksalt';
Y_Property = 'E'; % property to plot along Y-axis
Col_Structure = 'Rocksalt';
Col_Property = 'Loss'; % property to plot on color bar

%% Find model in results
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
sres = dir([ML_results_dir filesep '*' Salt '*' Model '*' 'Model_' Trial '*bayesopt.mat']);

if length(sres) > 1
    warning(['Multiple results found for: ' Salt ', ' Model ', Model ' Trial '. Using first result.']);
    sres = sres(1);
elseif isempty(sres)
    error(['No results found for: ' Salt ', ' Model ', Model ' Trial '.']);
end

model_filename = fullfile(sres.folder,sres.name);

if strcmp(Model,'TF') % Cannot show additivity for TF model
    show_additivity = false;
end

% Load data
data = load(model_filename);
results = data.results;
clear('data')

% Get search points
if ~isempty(Show_init)
    search_pnts = results.XTrace(1:Show_init);
    losses = results.ObjectiveTrace(1:Show_init);
    
    % Get user data: crystal properties at each point
    Crystal_data_full = results.UserDataTrace{1:Show_init};
else
    search_pnts = results.XTrace;
    losses = results.ObjectiveTrace;
    
    % Get user data: crystal properties at each point
    Crystal_data_full = results.UserDataTrace;
end

% Get the Structures (same order at each point)
N = length(Crystal_data_full{1});
Structures = cell(1,N);
for idx = 1:N
    Structures{idx} = Crystal_data_full{1}{idx}.Structure;
end

%% Get the properties of interest: Y-axis
M = length(Crystal_data_full);
structure_idx = strcmpi(Y_Structure,Structures);
Y_data = nan(1,M);

switch lower(Y_Property)
    case 'c_over_a'
        for idx = 1:M
            Y_data(idx) = Crystal_data_full{idx}{structure_idx}.c / Crystal_data_full{idx}{structure_idx}.a;
        end
        Y_Label = [Structure ' c/a'];
        Y_log = false;
    case 'loss'
        Y_data = losses;
        Y_Label = 'Loss';
        Y_log = false;
    otherwise
        for idx = 1:M
            Y_data(idx) = Crystal_data_full{idx}{structure_idx}.(Y_Property);
        end
        Y_Label = [Y_Structure ' ' Y_Property];
        
        switch Y_Property
            case {'CoulSR' 'CoulLR' 'Pot' 'CoulSR_MM' 'CoulSR_MX' 'CoulSR_XX' 'E'}
                Y_Label = [Y_Label ' [kJ mol$^{-1}$]'];
            case {'a' 'b' 'c'}
                Y_Label = [Y_Label ' [\AA]'];
        end
        Y_log = false;
end

%% Get the properties of interest: Color-axis
structure_idx = strcmpi(Col_Structure,Structures);
Col_data = nan(1,M);

switch lower(Col_Property)
    case 'c_over_a'
        for idx = 1:M
            Col_data(idx) = Crystal_data_full{idx}{structure_idx}.c / Crystal_data_full{idx}{structure_idx}.a;
        end
        Col_Label = [Structure ' c/a'];
        Col_log = false;
    case 'loss'
        Col_data = losses;
        Col_Label = 'Loss';
        Col_log = false;
    otherwise
        for idx = 1:M
            Col_data(idx) = Crystal_data_full{idx}{structure_idx}.(Col_Property);
        end
        Col_Label = [Col_Structure ' ' Col_Property];
        
        switch Col_Property
            case {'CoulSR' 'CoulLR' 'Pot' 'CoulSR_MM' 'CoulSR_MX' 'CoulSR_XX' 'E'}
                Col_Label = [Col_Label newline '[kJ mol$^{-1}$]'];
            case {'a' 'b' 'c'}
                Col_Label = [Col_Label newline '[\AA]'];
        end
        Col_log = false;
end

%% Convert to % error
if show_as_percent_error && contained_in_cell(Y_Property,{'Pot' 'E' 'a' 'b' 'c'})
    DFT = Load_Best_DFT_Data;
    DOI = DFT.(Salt).(Y_Structure);
    clear('DFT')
    
    switch Y_Property
        case {'Pot' 'E'}
            DFT_res = DOI.Energy;
        case {'a' 'b' 'c'}
            DFT_res = DOI.(Y_Property);
    end
    
    Y_data = (Y_data - DFT_res).*100./DFT_res;
    Y_Label = strrep(Y_Label,'[\AA]','[\% Error]');
    Y_Label = strrep(Y_Label,'[kJ mol$^{-1}$]','[\% Error]');
    
elseif show_as_percent_error
    warning(['Cannot plot ' Y_Property ' as error w.r.t. DFT results: No DFT comparison available.'])
end

if show_as_percent_error && contained_in_cell(Col_Property,{'Pot' 'E' 'a' 'b' 'c'})
    DFT = Load_Best_DFT_Data;
    DOI = DFT.(Salt).(Col_Structure);
    clear('DFT')
    
    switch Col_Property
        case {'Pot' 'E'}
            DFT_res = DOI.Energy;
        case {'a' 'b' 'c'}
            DFT_res = DOI.(Col_Property);
    end
    
    Col_data = (Col_data - DFT_res).*100./DFT_res;
    Col_Label = strrep(Col_Label,'[\AA]','[\% Error]');
    Col_Label = strrep(Col_Label,'[kJ mol$^{-1}$]','[\% Error]');
    
elseif show_as_percent_error
    warning(['Cannot plot ' Col_Property ' as error w.r.t. DFT results: No DFT comparison available.'])
end

%% Sort everything by the loss function
[losses_srt,idx] = sort(losses);
search_pnts_srt = search_pnts(idx,:);
Y_data_srt = Y_data(idx);
Col_data_srt = Col_data(idx);

% Extract only best M points
if isempty(Show)
    losses_srt_cut = losses_srt;
    search_pnts_cut = search_pnts_srt;
    Y_data_cut = Y_data_srt;
    Col_data_cut = Col_data_srt;
else
    losses_srt_cut = losses_srt(1:Show);
    search_pnts_cut = search_pnts_srt(1:Show,:);
    Y_data_cut = Y_data_srt(1:Show);
    Col_data_cut = Col_data_srt(1:Show);
end

%% Generate the MX parameters by additivity if they aren't present
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

%% Add the non-addivitive M-M dispersion force back in.
if contained_in_cell('SDMM2',search_pnts_cut.Properties.VariableNames)
    search_pnts_cut.SDMM = search_pnts_cut.SDMM + search_pnts_cut.SDMM2;
    search_pnts_cut.SDMM2 = [];
end

%% Convert to absolute values if requested
if ~show_as_scale && strcmp(Model,'JC')
    
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
    
elseif ~show_as_scale && strcmp(Model,'TF')
    PotSettings = Initialize_MD_Settings;
    PotSettings.Salt = Salt;
    PotSettings.S = Init_Scaling_Object;
    [OutputMX,OutputMM,OutputXX] = TF_Potential_Parameters(PotSettings);
    
    search_pnts_cut.D6MM = search_pnts_cut.SD6MM*OutputMM.C;
    search_pnts_cut.D6XX = search_pnts_cut.SD6XX*OutputXX.C;
    search_pnts_cut.D6MX = search_pnts_cut.SD6MX*OutputMX.C;
    
    search_pnts_cut.D8MM = search_pnts_cut.SD8MM*OutputMM.D;
    search_pnts_cut.D8XX = search_pnts_cut.SD8XX*OutputXX.D;
    search_pnts_cut.D8MX = search_pnts_cut.SD8MX*OutputMX.D;
    
    search_pnts_cut.AMM = search_pnts_cut.SAMM*OutputMM.alpha;
    search_pnts_cut.AXX = search_pnts_cut.SAXX*OutputXX.alpha;
    search_pnts_cut.AMX = search_pnts_cut.SAMX*OutputMX.alpha;
    
    search_pnts_cut.RMM = search_pnts_cut.SRMM*OutputMM.B;
    search_pnts_cut.RXX = search_pnts_cut.SRXX*OutputXX.B;
    search_pnts_cut.RMX = search_pnts_cut.SRMX*OutputMX.B;
end

%% Convert to sigma-epsilon form if requested
if show_as_sigma_epsilon && strcmp(Model,'JC')
    
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
elseif show_as_sigma_epsilon && strcmp(Model,'TF')
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
    
elseif ~show_as_scale && strcmp(Model,'JC')
    search_pnts_cut.SDMM = [];
    search_pnts_cut.SDXX = [];
    search_pnts_cut.SRMM = [];
    search_pnts_cut.SRXX = [];
    if ismember('SDMX',search_pnts_cut.Properties.VariableNames)
        search_pnts_cut.SDMX = [];
        search_pnts_cut.SRMX = [];
    end
    
elseif ~show_as_scale && strcmp(Model,'TF')
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
            tname = strrep(tname,'GA-','G. Depth ');
            tname = strrep(tname,'GB-','G. Pos. ');
            tname = strrep(tname,'GC-','G. Width ');
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
if show_as_percent_error && ~strcmpi(Col_Property,'loss')
    % Make sure the source map has an odd number of values
    map = flipud(max(min(cbrewer('div','RdBu',11,'spline'),1),0));
    
    % Center and expand color map
    map = center_colormap(map,Col_data_cut,0);
    
    colormap(figh,map)
else
    colormap(figh,flipud(hot))
end

x = 0;
P = gobjects(N,1);

T = tiledlayout(figh,2,ceil(N/2), 'Padding', 'loose', 'TileSpacing', 'tight'); 
for idx = 1:N
    nexttile;

    % do scatter plot
    X = search_pnts_cut{:,idx};
    
    P(idx) = scatter(X,Y_data_cut,10,Col_data_cut,'filled','SizeData',30,...
        'MarkerEdgeColor','k');
    cur_ax = P(idx).Parent;

    xlabel(cur_ax,pnames{idx},'Interpreter','latex','FontSize',fs);
    if Y_log
        set(cur_ax,'yscale','log');
    end

    % some axis properties
    set(cur_ax,'XGrid','On','YGrid','On','GridLineStyle','-','Layer','Top',...
        'FontSize',fs,'TickLabelInterpreter','latex')
        
    if Col_log
        set(cur_ax,'ColorScale','log');
    end
end
drawnow

% Add colorbar
sig = 1;
mincb = round(min(Col_data_cut),sig,'Significant');
maxcb = round(max(Col_data_cut),sig,'Significant');
while mincb == maxcb
    sig = sig+1;
    mincb = round(min(Col_data_cut),sig,'Significant');
    maxcb = round(max(Col_data_cut),sig,'Significant');
end

if show_as_percent_error && ~strcmpi(Col_Property,'loss')
    abs_max = max([abs(mincb) abs(maxcb)]);
    mincb = -abs_max;
    maxcb = abs_max;
end

cb = colorbar(P(1).Parent,'Position',[0.935 0.11 0.02 0.815],...
    'TickLabelInterpreter','latex','FontSize',fs,'Limits',[mincb maxcb]); % [left, bottom, width, height]

cb.Title.String = Col_Label;
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = fs+5;

% Add title
title(T,[Salt ' ' Model ' Model ' Trial ': Parameters vs Loss Function'],...
    'Interpreter','Latex','FontSize',fs+16);
xlabel(T,'Parameter','Interpreter','Latex','FontSize',fs+5)
ylabel(T,Y_Label,'Interpreter','Latex','FontSize',fs+5)
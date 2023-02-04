clear; %#ok<*UNRCH>
%% Data options
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; %  'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Theory = 'BH';
ModelID = 'MG';
BestOnly = true;
SelectOnly = [5 3 1 5]; %[2 4 5 1]; %[5 3 1 5];
Reps = [1:5];
savefile = false; % switch to save the final plots to file
saveloc = 'C:\Users\Hayden\Documents\Patey_Lab\Amgen_Presentation_Images';
DM_Multiplier = 1e5;

%% Plot options
fs = 28; % font size
markers = {'o' 's' '^' 'v' 'd' '>' '<' 'p' 'h' 'x'};
show_as_percent_error = true; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)
include_av_line = false;
wwbuffer = 0.075;
ttbuffer = 0.08;
bbbuffer = 0.15;

plot_LE = true;
plot_RLE = true;
plot_a = false;
plot_c = false;
plot_ac = false;
plot_volume = true;
plot_density = false;
plot_loss = false;
plot_finite_T_data = true;

ylim_LE = [-2 15];
ylim_RLE = [-2 20];
ylim_a = [];
ylim_c = [];
ylim_ac = [];
ylim_volume = [-15 15];
ylim_density = [];

ylim_MP = [-55 5];
ylim_dH = [-4 12];
ylim_dV = [-5 20];
ylim_LiqV = [-15 10];
ylim_SolV = [-15 10];


% ylim_LE = [-5 20];
% ylim_RLE = [0 80];
% ylim_a = [];
% ylim_c = [];
% ylim_ac = [];
% ylim_volume = [0 170];
% ylim_density = [];
% 
% ylim_MP = [];
% ylim_dH = [0 100];
% ylim_dV = [-120 120];
% ylim_LiqV = [0 100];
% ylim_SolV = [0 100];

%% Script begins
cm3_per_Ang3 = 1e-24; % cubic cm/cubic angstrom
N_A = 6.02214076e23; % formula units/mol
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';

Structures   = {'Rocksalt' 'Wurtzite'};
MinPlotTypes = {'LE' 'RLE' 'a' 'c' 'ac' 'V' 'Density'};
MinPlotIdxes = [plot_LE plot_RLE plot_a plot_c plot_ac plot_volume plot_density];
MinPlotYLims = {ylim_LE ylim_RLE ylim_a ylim_c ylim_ac ylim_volume ylim_density};
MinPlotTypes = MinPlotTypes(MinPlotIdxes);
MinPlotYLims = MinPlotYLims(MinPlotIdxes);
FiniteTTypes = {'MP' 'dH' 'dV' 'LiqV' 'SolV'};
FiniteTYLims = {ylim_MP ylim_dH ylim_dV ylim_LiqV ylim_SolV};
% MP_Titles = {'$T_{m} - T_{m}^{*}$' ...
%     '$\Delta H^{\circ}_{\mathrm{fus}} - \Delta H^{\circ *}_{\mathrm{fus}}$' ...
%     '$\Delta V_{\mathrm{fus}} - \Delta V_{\mathrm{fus}}^{*}$' ...
%     '$V_{\mathrm{liq}}(T_{m}^{*}) - V_{\mathrm{liq}}^{*}(T_{m}^{*})$' ...
%     '$V_{\mathrm{sol}}(T_{m}^{*}) - V_{\mathrm{sol}}^{*}(T_{m}^{*})$' ...
%     '$D_{\mathrm{liq,Li}^{+}} - D_{\mathrm{liq,Li}^{+}}^{*}$'};
MP_Titles = {'Error in $T_m$' ...
    '$\Delta H^{\circ}_{\mathrm{fus}}$ Error' ...
    '$\Delta V_{\mathrm{fus}}$ Error' ...
    '$V_{\mathrm{liq}}(T_{m})$ Error' ...
    '$V_{\mathrm{sol}}(T_{m})$ Error'};
FT_Loss_Names = {'MP' 'Fusion_Enthalpy' 'MP_Volume_Change' 'Liquid_MP_Volume' 'Solid_MP_Volume'};

if show_as_percent_error
    MinPlotTitles = {'##Struct## $E_{L}$ Err.' ...
        '$\Delta E_{L}$(##Struct##)' ...
        '$a - a^{*}$' ...
        '$c - c^{*}$' ...
        '$c/a - c^{*}/a^{*}$' ...
        '##Struct## $V(T = 0)$ Err.' ...
        '$\rho - \rho^{*}$'};
    MinPlotYLabels = {'[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
    
    FiniteTYLabels = {'[K]' '[\%]' '[\%]' '[\%]' '[\%]' '[\%]'};
else
    MinPlotTitles = {'##Struct## $E_{L}$ Err.' ...
        'E_{L}$(##Struct##) - $E_{L}$(rocksalt)' ...
        '$a - a^{*}$' ...
        '$c - c^{*}$' ...
        '$c/a - c^{*}/a^{*}$' ...
        '##Struct## $V(T = 0)$ Err.' ...
        '$\rho - \rho^{*}$'};
    MinPlotYLabels = {'[kJ mol$^{-1}$]' ...
                      '[kJ mol$^{-1}$]' ...
                      '[\AA]' ...
                      '[\AA]' ...
                      '[$(c/a)$]' ...
                      '[$\AA^{3}$]' ...
                      'g cm$^{-3}$]'};
                  
    FiniteTYLabels = {'[K]' '[kJ mol$^{-1}$]' '[$\AA^{3}$]' '[$\AA^{3}$]' '[$\AA^{3}$]' ['[$10^{-' num2str(log10(DM_Multiplier)) '}$ cm$^{2}$/s]']};
end
MinPlotYLabels = MinPlotYLabels(MinPlotIdxes);
MinPlotTitles = MinPlotTitles(MinPlotIdxes);

N_Salts = numel(Salts);
N_Markers = numel(markers);
N_Structures = numel(Structures);
N_UppPlot_Rows = sum([plot_loss plot_finite_T_data]);
N_MinPlot_Rows = sum(MinPlotIdxes);
N_FiniteTTypes = numel(FiniteTTypes);
N_Rows = N_UppPlot_Rows + N_MinPlot_Rows;

% Find reps of models
Models = cell(1,length(Salts));
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
Min_PlotData        = nan(N_MinPlot_Rows,N_Salts,N_Models,N_Structures);
Finite_T_PlotData   = nan(N_FiniteTTypes,N_Salts,N_Models);
Rel_Energies_DFT    = nan(N_Salts,N_Structures);
Total_loss          = nan(N_Salts,N_Models);

% Plot color scheme
Colours = cbrewer('qual','Set3',max(N_Salts,3));

% Load target DFT/Experimental data
DFT = Load_Best_DFT_Data;
Exp = Load_Experimental_Data;
for idx = 1:N_Salts
    % Correct DFT to experimental lattice energies
    Correction = Exp.(Salts{idx}).Rocksalt.E - DFT.(Salts{idx}).Rocksalt.Energy;
    for jdx = 1:N_Structures
        DFT.(Salts{idx}).(Structures{jdx}).Energy = DFT.(Salts{idx}).(Structures{jdx}).Energy + Correction;
    end
    
    % Update targets to experimental lattice parameters where available
    DFT.(Salts{idx}).Rocksalt.a = Exp.(Salts{idx}).Rocksalt.a_zero;
    DFT.(Salts{idx}).Rocksalt.b = Exp.(Salts{idx}).Rocksalt.b_zero;
    DFT.(Salts{idx}).Rocksalt.c = Exp.(Salts{idx}).Rocksalt.c_zero;
    DFT.(Salts{idx}).Rocksalt.V = Exp.(Salts{idx}).Rocksalt.V_zero;
    DFT.(Salts{idx}).Rocksalt.density = Exp.(Salts{idx}).Rocksalt.density_zero;
    if isfield(Exp.(Salts{idx}),'Wurtzite') && ~isnan(Exp.(Salts{idx}).Wurtzite.a)
        DFT.(Salts{idx}).Wurtzite.a = Exp.(Salts{idx}).Wurtzite.a_zero;
        DFT.(Salts{idx}).Wurtzite.b = Exp.(Salts{idx}).Wurtzite.b_zero;
        DFT.(Salts{idx}).Wurtzite.c = Exp.(Salts{idx}).Wurtzite.c_zero;
        DFT.(Salts{idx}).Wurtzite.V = Exp.(Salts{idx}).Wurtzite.V_zero;
        DFT.(Salts{idx}).Wurtzite.density = Exp.(Salts{idx}).Wurtzite.density_zero;
    end
    
    % Calculate target relative lattice energies
    if plot_RLE
        for jdx = 1:N_Structures
            % Relative Energies
            Target_gap = DFT.(Salts{idx}).(Structures{jdx}).Energy - DFT.(Salts{idx}).Rocksalt.Energy;
            Rel_Energies_DFT(idx,jdx) = Target_gap;
        end
    end
end

% Find the targets
data_found = false;
for idx = 1:N_Salts
    for jdx = 1:N_Models
        dat_file = fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' Models{jdx} '_data.mat']);
        if ~isfile(dat_file)
            continue
        else
            % Load data
            try
                data = load(dat_file).full_data;
                Bayesopt_model = data.Settings;
                Bayesopt_Loss_Options = Bayesopt_model.Loss_Options;
                data_found = true;
                break
            catch
                continue
            end
        end
    end
    if data_found
        break
    end
end
if ~data_found
    warning(['Unable to load targets for ' Theory ' Model ' ModelID])
    Bayesopt_Loss_Options = init_loss_options;
elseif plot_volume && ~plot_a
    for jdx = 1:N_Structures
        if Bayesopt_Loss_Options.(Structures{jdx}).V < sqrt(eps)
            Bayesopt_Loss_Options.(Structures{jdx}).V = Bayesopt_Loss_Options.(Structures{jdx}).a;
        end
    end
end

if isfield(Bayesopt_model,'Fix_Charge') && ~Bayesopt_model.Fix_Charge
    Fix_Charge = 'Parameter Charge';
else
    Fix_Charge = 'Fixed Charge';
end
if isfield(Bayesopt_model,'Additivity') && ~Bayesopt_model.Additivity
    Additivity = 'Additivity Off';
else
    Additivity = 'Additivity On';
end

% Minimization Data
if N_MinPlot_Rows
    Settings = Initialize_MD_Settings;
    for idx = 1:N_Salts
        Salt = Salts{idx};
        Settings.Salt = Salt;
        for iidx = 1:N_Models
            Model = Models{iidx};
            
            % Find the fully optimized model
            dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);
            
            if ~isfile(dat_file)
                disp(['Could not load results found for: ' Salt ', ' Theory ', Model ' Model '.']);
                continue
            end
            
            % Load data
            try
                data = load(dat_file).full_data;
%                 data.Settings.Model = Model;
%                 [data.Settings,ModelFound] = Load_Model_Params(data.Settings);
                
                Minimization_Data = data.Minimization_Data;
                
                if isfield(data,'secondary_result')
                    optimvals = nan(1,length(data.secondary_result));
                    for jdx = 1:length(data.secondary_result)
                        optimvals(jdx) = [data.secondary_result(jdx).optimValues.fval];
                    end
                    Total_loss(idx,iidx) = min(optimvals);
                elseif isfield(data,'bayesopt_results') && ~isempty(data.bayesopt_results)
                    %Total_loss(idx,iidx) = data.bayesopt_results.MinObjective;
                    Total_loss(idx,iidx) = min(data.bayesopt_results.ObjectiveTrace);
                else
                    Total_loss(idx,iidx) = nan;
                end
            catch
                disp(['Could not obtain crystal minimization data for: ' Salt ', ' Theory ', Model ' Model '.']);
                continue
            end

            % Grab available structures from data
            Min_Dat_Structs = cell(1,length(Minimization_Data));
            for jdx = 1:length(Minimization_Data)
                Min_Dat_Structs{jdx} = Minimization_Data{jdx}.Structure;
            end

            % Loop through structures
            for jdx = 1:N_Structures
                Structure = Structures{jdx};
                Settings.Structure = Structure;
                DefGeometry = Default_Crystal(Settings);

                md_idx = find(strcmpi(Min_Dat_Structs,Structure));
                rs_idx = find(strcmpi(Min_Dat_Structs,'Rocksalt'));
                if isempty(md_idx)
                    continue
                end

                for rowidx = 1:N_MinPlot_Rows
                    if show_as_percent_error
                        switch MinPlotTypes{rowidx}
                            case 'LE'
                                Min_PlotData(rowidx,idx,iidx,jdx) = (Minimization_Data{md_idx}.E - ...
                                    DFT.(Salt).(Structure).Energy)*100/ ...
                                    abs(DFT.(Salt).(Structure).Energy);
                            case 'RLE'
                                if md_idx ~= rs_idx
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Minimization_Data{md_idx}.E - ...
                                        Minimization_Data{rs_idx}.E;
                                end
                            case 'a'
                                Min_PlotData(rowidx,idx,iidx,jdx) = (Minimization_Data{md_idx}.a - ...
                                    DFT.(Salt).(Structure).a)*100/ ...
                                    DFT.(Salt).(Structure).a;
                            case 'c'
                                Min_PlotData(rowidx,idx,iidx,jdx) = (Minimization_Data{md_idx}.c - ...
                                    DFT.(Salt).(Structure).c)*100/ ...
                                    DFT.(Salt).(Structure).c;
                            case 'ac'
                                Target_ca = DFT.(Salt).(Structure).c / DFT.(Salt).(Structure).a;
                                Model_ca  = Minimization_Data{md_idx}.c / ...
                                             Minimization_Data{md_idx}.a;
                                Min_PlotData(rowidx,idx,iidx,jdx) = (Model_ca - Target_ca)*100 / Target_ca;
                            case {'V' 'Density'}
                                ABC_vec = DefGeometry.Transform.*...
                                    [Minimization_Data{md_idx}.a; ...
                                     Minimization_Data{md_idx}.b; ...
                                     Minimization_Data{md_idx}.c];
                                Volume = ( dot( cross(ABC_vec(1,:),ABC_vec(2,:)),ABC_vec(3,:) ) ) / DefGeometry.NF; % Angstrom^3 / Formula Unit
                                Density = Exp.(Salt).MM/(Volume*cm3_per_Ang3*N_A); % g/cm^3
                                if strcmp(MinPlotTypes{rowidx},'V')
                                    Min_PlotData(rowidx,idx,iidx,jdx) = (Volume - DFT.(Salt).(Structure).V)*100/ ...
                                        DFT.(Salt).(Structure).V;
                                else
                                    Min_PlotData(rowidx,idx,iidx,jdx) = (Density - DFT.(Salt).(Structure).density)*100/ ...
                                        DFT.(Salt).(Structure).density;
                                end
                        end
                    else
                        switch MinPlotTypes{rowidx}
                            case 'LE'
                                Min_PlotData(rowidx,idx,iidx,jdx) = Minimization_Data{md_idx}.E - ...
                                    DFT.(Salt).(Structure).Energy;
                            case 'RLE'
                                if md_idx ~= rs_idx
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Minimization_Data{md_idx}.E - ...
                                        Minimization_Data{rs_idx}.E;
                                end
                            case 'a'
                                Min_PlotData(rowidx,idx,iidx,jdx) =  Minimization_Data{md_idx}.a - ...
                                    DFT.(Salt).(Structure).a;
                            case 'c'
                                Min_PlotData(rowidx,idx,iidx,jdx) = Minimization_Data{md_idx}.c - ...
                                    DFT.(Salt).(Structure).c;
                            case 'ac'
                                Target_ca = DFT.(Salt).(Structure).c / DFT.(Salt).(Structure).a;
                                Model_ca  = Minimization_Data{md_idx}.c / ...
                                             Minimization_Data{md_idx}.a;
                                Min_PlotData(rowidx,idx,iidx,jdx) = Model_ca - Target_ca;
                            case {'V' 'Density'}
                                ABC_vec = DefGeometry.Transform.*...
                                    [Minimization_Data{md_idx}.a; ...
                                     Minimization_Data{md_idx}.b; ...
                                     Minimization_Data{md_idx}.c];
                                Volume = ( dot( cross(ABC_vec(1,:),ABC_vec(2,:)),ABC_vec(3,:) ) ) / DefGeometry.NF; % Angstrom^3 / Formula Unit
                                Density = Exp.(Salt).MM/(Volume*cm3_per_Ang3*N_A); % g/cm^3
                                if strcmp(MinPlotTypes{rowidx},'V')
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Volume - DFT.(Salt).(Structure).V;
                                else
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Density - DFT.(Salt).(Structure).density;
                                end
                        end
                    end
                end
            end
            
            % Finite T Data
            if plot_finite_T_data
                try
                    Finite_T_Data = data.Finite_T_Data;
                catch
                    Settings.Salt = Salt;
                    Finite_T_Data = Initialize_Finite_T_Data(Settings);
                    disp(['No finite T data for: ' Salt ', ' Theory ', Model ' Model '.'])
                end
                
                if ~isfield(Finite_T_Data,'Liquid_DM_MP')
                    Finite_T_Data.Liquid_DM_MP = nan;
                    Finite_T_Data.Exp_DM_MP = nan;
                end
                
                for ridx = 1:N_FiniteTTypes
                    if show_as_percent_error
                        switch FiniteTTypes{ridx}
                            case 'dH'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.Fusion_dH - Finite_T_Data.Exp_Fusion_dH)...
                                	/abs(Finite_T_Data.Exp_Fusion_dH);
                            case 'dV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.Fusion_dV - Finite_T_Data.Exp_Fusion_dV)...
                                    /abs(Finite_T_Data.Exp_Fusion_dV);
                            case 'LiqV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.Liquid_V_MP - Finite_T_Data.Exp_Liquid_V_MP)...
                                	/abs(Finite_T_Data.Exp_Liquid_V_MP);
                            case 'SolV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.Solid_V_MP - Finite_T_Data.Exp_Solid_V_MP)...
                                	/abs(Finite_T_Data.Exp_Solid_V_MP);
                            case 'MP'
%                                 Finite_T_PlotData(ridx,idx,iidx) = ...
%                                     100*(Finite_T_Data.MP - Finite_T_Data.Exp_MP)...
%                                 	/abs(Finite_T_Data.Exp_MP);
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.MP - Finite_T_Data.Exp_MP;
                            case 'DM'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.Liquid_DM_MP    - Finite_T_Data.Exp_DM_MP)...
                                    /abs(Finite_T_Data.Exp_DM_MP);
                        end
                    else
                        switch FiniteTTypes{ridx}
                            case 'dH'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.Fusion_dH - Finite_T_Data.Exp_Fusion_dH;
                            case 'dV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.Fusion_dV - Finite_T_Data.Exp_Fusion_dV;
                            case 'LiqV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.Liquid_V_MP - Finite_T_Data.Exp_Liquid_V_MP;
                            case 'SolV'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.Solid_V_MP - Finite_T_Data.Exp_Solid_V_MP;
                            case 'MP'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    Finite_T_Data.MP - Finite_T_Data.Exp_MP;
                            case 'DM'
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    (Finite_T_Data.Liquid_DM_MP - Finite_T_Data.Exp_DM_MP).*DM_Multiplier;
                        end
                    end
                end
            end
        end
    end
end

if BestOnly && numel(SelectOnly) ~= N_Salts
    N_Models = 1;
    
    Min_PlotData_BB        = nan(N_MinPlot_Rows,N_Salts,N_Structures);
    Finite_T_PlotData_BB   = nan(N_FiniteTTypes,N_Salts);
    
    [Total_loss,midx] = min(Total_loss,[],2);
    
    for jdx = 1:N_Salts
        Model_idx = midx(jdx);
        for idx = 1:N_MinPlot_Rows
            for kdx = 1:N_Structures
                Min_PlotData_BB(idx,jdx,kdx) = Min_PlotData(idx,jdx,Model_idx,kdx);
            end
        end
        for idx = 1:N_FiniteTTypes
            Finite_T_PlotData_BB(idx,jdx) = Finite_T_PlotData(idx,jdx,Model_idx);
        end
    end
elseif BestOnly && numel(SelectOnly) == N_Salts
    N_Models = 1;
    Min_PlotData_BB        = nan(N_MinPlot_Rows,N_Salts,N_Structures);
    Finite_T_PlotData_BB   = nan(N_FiniteTTypes,N_Salts);
    
    TL = nan(N_Salts,1);
    for jdx = 1:N_Salts
        Model_idx = SelectOnly(jdx);
        TL(jdx) = Total_loss(jdx,SelectOnly(jdx));
        for idx = 1:N_MinPlot_Rows
            for kdx = 1:N_Structures
                Min_PlotData_BB(idx,jdx,kdx) = Min_PlotData(idx,jdx,Model_idx,kdx);
            end
        end
        for idx = 1:N_FiniteTTypes
            Finite_T_PlotData_BB(idx,jdx) = Finite_T_PlotData(idx,jdx,Model_idx);
        end
    end
    Total_loss = TL;
end

%% Plotting
% Gen figure
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
switch Theory
    case 'JC'
        PubTheoryName = 'CLJ';
    case 'BH'
        PubTheoryName = 'CBH';
    case 'BF'
        PubTheoryName = 'Coulomb Wang-Buckingham';
    case 'TF'
        PubTheoryName = 'CBHM';
    case 'BD'
        PubTheoryName = 'CBH Modifification D';
    case 'BE'
        PubTheoryName = 'CBH Modifification E';
    case 'Mie'
        PubTheoryName = 'Coulomb Mie[n-6]';
end


%% Generate plotting panels
xx = 1:N_Salts;
row_height = 1/N_Rows;
MinPlot_Frac = N_MinPlot_Rows/N_Rows;
axobj_ft = gobjects(N_FiniteTTypes,1);

t = tiledlayout(figh,4,4);

%% start with the MP in a special box
axmp = axes(t);
axmp.Layout.Tile = 3;
axmp.Layout.TileSpan = [4 2];

for idx = 1:N_FiniteTTypes
    if idx == 1
        axobj_ft(idx) = axmp;
    else
        axobj_ft(idx) = nexttile(t); % [left bottom width height]
    end
    hold(axobj_ft(idx),'on');

    % Gather data
    plot_data = squeeze(Finite_T_PlotData_BB(idx,:));

    % plot data
    p = bar(axobj_ft(idx),1,plot_data,'FaceColor','flat','Visible','on','BarWidth',1,...
        'LineWidth',2);
    for jdx = 1:N_Salts
        p(jdx).CData = Colours(jdx,:);
    end

    % Add boxes for targets
    if isempty(FiniteTYLims{idx})
        ylim(axobj_ft(idx),'padded');
    else
        ylim(axobj_ft(idx),FiniteTYLims{idx});
    end
    xlim(axobj_ft(idx),'padded');
    ylmits = ylim;
    xlmits = xlim;
    if isfield(Bayesopt_Loss_Options,FT_Loss_Names{idx}) && Bayesopt_Loss_Options.(FT_Loss_Names{idx}) > sqrt(eps)
        rech = rectangle(axobj_ft(idx),'Position',[xlmits(1) ylmits(1) (xlmits(2)-xlmits(1)) (ylmits(2) - ylmits(1))],...
                'Curvature',0,'LineWidth',3,'EdgeColor','r');
    end

    % Set plot properties
    set(axobj_ft(idx),'box','on','TickLabelInterpreter','latex');
    set(axobj_ft(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
    xticklabels(axobj_ft(idx),[])
    axobj_ft(idx).XAxis.TickLength = [0,0];
    ylim(axobj_ft(idx),ylmits)
    xlim(axobj_ft(idx),xlmits)
    ylabel(axobj_ft(idx), FiniteTYLabels{idx},'Fontsize',fs-3,'Interpreter','latex')
    title(axobj_ft(idx),MP_Titles(idx),'Fontsize',fs,'Interpreter','latex')

    set(axobj_ft(idx), 'YTickMode', 'auto')
    set(axobj_ft(idx),'xminorgrid','off','yminorgrid','on')
end


% Create figure and axis
axobj_min = gobjects(N_Structures,N_MinPlot_Rows);

for jdx = 1:N_MinPlot_Rows
    cylim = zeros(N_Structures,2);
    for idx = 1:N_Structures
        
        if (strcmp(MinPlotTypes{jdx},'LE')  && strcmp(Structures{idx},'Wurtzite')) || ...
           (strcmp(MinPlotTypes{jdx},'RLE') && strcmp(Structures{idx},'Rocksalt'))
            continue
        end
        
        axobj_min(idx,jdx) = nexttile(t); % [left bottom width height]
        hold(axobj_min(idx,jdx),'on');

        % Gather data
        plot_data = squeeze(Min_PlotData_BB(jdx,:,idx));

        % plot data
        p = bar(axobj_min(idx,jdx),1,plot_data,'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        for kdx = 1:N_Salts
            p(kdx).CData = Colours(kdx,:);
        end

        % Set plot properties
        set(axobj_min(idx,jdx),'box','on','TickLabelInterpreter','latex');
        set(axobj_min(idx,jdx),'XMinorTick','off','YMinorTick','on','FontSize',fs);
        xticklabels(axobj_min(idx,jdx),[])
        axobj_min(idx,jdx).XAxis.TickLength = [0,0];
        
        if isempty(MinPlotYLims{jdx})
            ylim(axobj_min(idx,jdx),'padded');
        else
            ylim(axobj_min(idx,jdx),MinPlotYLims{jdx});
        end
        xlim(axobj_min(idx,jdx),'padded')
        xlmits = xlim;
        ylmits = ylim;

        if strcmp(MinPlotTypes{jdx},'RLE') && ~strcmp(Structures{idx},'Rocksalt')
            for ib = 1:numel(p)
                xx(ib) = p(ib).XData+p(ib).XOffset;
            end

            p = bar(axobj_min(idx,jdx),1,Rel_Energies_DFT(:,idx),'FaceColor',[0.5 0.5 0.5],'Visible','on','BarWidth',1,...
                'LineWidth',2,'EdgeColor','b','FaceAlpha',0.5,'LineStyle','--');
        end
        
        if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
            cylim(idx,:) = ylim(axobj_min(idx,jdx));
        end

        ylabel(axobj_min(idx,jdx), MinPlotYLabels{jdx},'Fontsize',fs,'Interpreter','latex')
        title(axobj_min(idx,jdx),strrep(MinPlotTitles(jdx),'##Struct##',Structures{idx}),'Fontsize',fs,'Interpreter','latex')
        ylmits = [min(cylim(:,1)) max(cylim(:,2))];
        set(axobj_min(idx,jdx), 'YTickMode', 'auto')
        set(axobj_min(idx,jdx),'xminorgrid','off','yminorgrid','on')

        % Add boxes for targets
        if isfield(Bayesopt_Loss_Options,Structures{idx}) && Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}) > sqrt(eps)
            rech = rectangle(axobj_min(idx,jdx),'Position',[xlmits(1) ylmits(1) (xlmits(2) - xlmits(1)) (ylmits(2) - ylmits(1))],...
                    'Curvature',0,'LineWidth',3,'EdgeColor','r');
        end
        xlim(axobj_min(idx,jdx),xlmits)
        ylim(axobj_min(idx,jdx),ylmits)
    end
end
                
h = bar(axobj_min(1,1),nan,nan(1,N_Salts),'FaceColor','flat','Visible','on','BarWidth',1,...
    'LineWidth',2);
for kdx = 1:N_Salts
    h(kdx).CData = Colours(kdx,:);
end
h(end+1) = bar(axobj_min(1,1),nan,nan,'FaceColor',[0.5 0.5 0.5],'Visible','on','BarWidth',1,...
    'LineWidth',2,'EdgeColor','b','FaceAlpha',0.5,'LineStyle','--');
legtxt = strcat({'\ \ '},Salts,{'\quad'});
legtxt{end+1}  = '\ \ DFT Target';
legh = legend(h,legtxt,'Position',[0.25 0 0.5 0.06],'Orientation','Horizontal',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts+1);






% if savefile
%     
%     %set(mainpanelh,'renderer','opengl')
%     %saveas(figh,filename,'pdf')
%     %print(mainpanelh,filename,'-bestfit')
% %    exportgraphics(mainpanelh,fullfile(saveloc,filename))
%     
% %     figh.PaperPositionMode='manual';
% %     figh.PaperUnits ='normalized';
% %     figh.PaperPosition = [0 0 1 1];
% %     exportgraphics(mainpanelh,filename)
%     %exportapp(figh,filename)
%     %filename = ['Model_Overview_' Theory '_Model_' ModelID '.eps'];
%     %hgexport(figh,fullfile(saveloc,filename))
%     %print(figh,fullfile(saveloc,filename),'-deps','-noui')
%     figh=gcf;
%     filename = ['Model_Overview_' Theory '_Model_' ModelID '.png'];
%     print(figh,fullfile(saveloc,filename),'-dpng','-r300','-noui')
% end
filename = ['Model_Overview_' Theory '_Model_' ModelID '.emf'];
exportgraphics(figh,fullfile(saveloc,filename),'BackgroundColor','none')

function Output = LegStructure(Structure)
    switch Structure
        case 'BetaBeO'
            Output = '$\beta$-BeO';
        case 'FiveFive'
            Output = '5-5';
        otherwise
            Output = Structure;
    end
end

function ylim = GetYlims(FiniteTType,ylims)
    FiniteTTypes

end

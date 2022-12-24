clear; %#ok<*UNRCH>
%% Data options
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; %  'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Theory = 'BF';
ModelID = 'MK';
BestOnly = false;
SelectOnly = [];
Reps = [1:5];
savefile = false; % switch to save the final plots to file
saveloc = 'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Thesis\Thesis_Draft\BO_Figures';
DM_Multiplier = 1e5;

%% Plot options
fs = 28; % font size
markers = {'o' 's' '^' 'v' 'd' '>' '<' 'p' 'h' 'x'};
show_as_percent_error = false; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)
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

%% Script begins
cm3_per_Ang3 = 1e-24; % cubic cm/cubic angstrom
N_A = 6.02214076e23; % formula units/mol
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';

Structures   = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
MinPlotTypes = {'LE' 'RLE' 'a' 'c' 'ac' 'V' 'Density'};
MinPlotIdxes = [plot_LE plot_RLE plot_a plot_c plot_ac plot_volume plot_density];
MinPlotTypes = MinPlotTypes(MinPlotIdxes);
FiniteTTypes = {'MP' 'dH' 'dV' 'LiqV' 'SolV' 'DM'};
% MP_Titles = {'$T_{m} - T_{m}^{*}$' ...
%     '$\Delta H^{\circ}_{\mathrm{fus}} - \Delta H^{\circ *}_{\mathrm{fus}}$' ...
%     '$\Delta V_{\mathrm{fus}} - \Delta V_{\mathrm{fus}}^{*}$' ...
%     '$V_{\mathrm{liq}}(T_{m}^{*}) - V_{\mathrm{liq}}^{*}(T_{m}^{*})$' ...
%     '$V_{\mathrm{sol}}(T_{m}^{*}) - V_{\mathrm{sol}}^{*}(T_{m}^{*})$' ...
%     '$D_{\mathrm{liq,Li}^{+}} - D_{\mathrm{liq,Li}^{+}}^{*}$'};
MP_Titles = {'$T_{m}$' ...
    '$\Delta H^{\circ}_{\mathrm{fus}}$' ...
    '$\Delta V_{\mathrm{fus}}$' ...
    '$V_{\mathrm{liq}}(T_{m}^{*})$' ...
    '$V_{\mathrm{sol}}(T_{m}^{*})$' ...
    '$D_{\mathrm{liq,Li}^{+}}(T_{m}^{*})$'};
FT_Loss_Names = {'MP' 'Fusion_Enthalpy' 'MP_Volume_Change' 'Liquid_MP_Volume' 'Solid_MP_Volume' 'Liquid_DM_MP'};

if show_as_percent_error
    MinPlotTitles = {'$\left( E_{L} - E_{L}^{*} \right) / \left| E_{L}^{*} \right|$' ...
        'Error in $E_{L}$: $E_{L}$[Structure] - $E_{L}$[Rocksalt]' ...
        '$\left( a - a^{*} \right) / \left| a^{*} \right|$' ...
        '$\left( c - c^{*} \right) / \left| c^{*} \right|$' ...
        '$\left( c/a - c^{*}/a^{*} \right) / \left| c^{*}/a^{*} \right|$' ...
        '$\left( V(T = 0) - V(T = 0)^{*} \right) / \left| V^{*} \right|$' ...
        '$\left( \rho - \rho^{*} \right) / \left| \rho^{*} \right|$'};
    MinPlotYLabels = {'[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
    
    FiniteTYLabels = {'[\%]' '[\%]' '[\%]' '[\%]' '[\%]' '[\%]'};
else
    MinPlotTitles = {'Error in $E_{L}$: $E_{L} - E_{L}^{*}$' ...
        '$\Delta E_{L} = E_{L}$(other structure) - $E_{L}$(rocksalt)' ...
        '$a - a^{*}$' ...
        '$c - c^{*}$' ...
        '$c/a - c^{*}/a^{*}$' ...
        'Error in $V(T = 0)$: $V(T = 0) - V^{*}(T = 0)$' ...
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
                                Finite_T_PlotData(ridx,idx,iidx) = ...
                                    100*(Finite_T_Data.MP - Finite_T_Data.Exp_MP)...
                                	/abs(Finite_T_Data.Exp_MP);
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
    'Name','','Visible','On','color','w');
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

if N_Salts > 1 && ~BestOnly
    % Generate plotting panels
    xx = 1:N_Salts;
    row_height = 1/N_Rows;
    MinPlot_Frac = N_MinPlot_Rows/N_Rows;
    if plot_loss
        loss_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height 1 row_height],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        FTdrop = row_height;
        axobj_loss = gobjects(N_Salts,1);

        ww0 = 1 - 0.04;
        bb0 = 0.05;
        ww1 = ww0/N_Salts - ww0/(N_Salts+0.5);
        ww2 = ww1/(max(N_Salts-1,1));
        ww = ww0/(N_Salts+0.5);
        hh = 1 - bb0 - 0.3;
        bb = bb0;
        for idx = 1:N_Salts
            ll = (0.995-ww0) + (ww0/N_Salts+ww2)*(idx-1);
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
                scatter(axobj_loss(idx),xpos,plot_data(mdx)./(10.^expon),200,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(idx,:),'linewidth',2,'marker',markers{mdx},'MarkerFaceAlpha',0.5);
            end

            % Set plot properties
            set(axobj_loss(idx),'box','on','TickLabelInterpreter','latex');
            set(axobj_loss(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
            %set(axobj_loss(idx),'YScale','log');
            xticks(axobj_loss(idx),1:N_Models)
            xticklabels(axobj_loss(idx),[])
            axobj_loss(idx).XAxis.TickLength = [0,0];
            ylim(axobj_loss(idx),[0 max(plot_data./(10.^expon))*1.1])
            if savefile
                sgtitle(loss_panel,'Minimized Objective Function $f\left(\mathbf{x}^{*}\right)$',...
                    'Fontsize',fs-2,'Interpreter','latex')
            else
                sgtitle(loss_panel,['Minimized Objective Function: ' PubTheoryName ' Model ' ModelID ' [' Fix_Charge ' / ' Additivity ']'],...
                    'Fontsize',fs-2,'Interpreter','latex')
            end
%             if idx == 1
%                 ylabel(axobj_loss(idx), '$f(\mathbf{x}^{*})$','Fontsize',fs,'Interpreter','latex')
%             end
            xlim(axobj_loss(idx),[0 N_Models+1])
            set(axobj_loss(idx), 'YTickMode', 'auto');
            set(axobj_loss(idx),'xminorgrid','off','yminorgrid','on');
            axobj_loss(idx).YAxis.Exponent = 0;
            ylabel(axobj_loss(idx), ['$f\left(\mathbf{x}^{*}\right)/10^{' num2str(expon) '}$'],'Fontsize',fs-4,'Interpreter','latex')
        end
    else
        FTdrop = 0;
    end

    if plot_finite_T_data
        finite_T_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height-FTdrop-0.0247 1 row_height+0.0227],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        
        axobj_ft = gobjects(N_FiniteTTypes,1);

        ww0 = 1 - 0.07;
        bb0 = 0.14;
        ww1 = ww0/N_FiniteTTypes - ww0/(N_FiniteTTypes+2.6);
        ww2 = ww1/(N_FiniteTTypes-1);
        ww = ww0/(N_FiniteTTypes+2.6);
        hh = 1 - bb0 - 0.25;
        bb = bb0;
        for idx = 1:N_FiniteTTypes
            ll = (0.995-ww0) + (ww0/N_FiniteTTypes+ww2)*(idx-1);
            axobj_ft(idx) = subplot('Position',[ll bb ww hh],...
                'parent',finite_T_panel); % [left bottom width height]
            hold(axobj_ft(idx),'on');

            % Gather data
            plot_data = squeeze(Finite_T_PlotData(idx,:,:));
            salt_ave = mean(plot_data,2,'omitnan'); % average the data across salts

            % Plot data
            if include_av_line
                plot(axobj_ft(idx),xx,salt_ave,'--k','linewidth',2);
            end
            for jdx = 1:N_Models
                yy = plot_data(:,jdx);
                scatter(axobj_ft(idx),xx,yy,200,Colours(xx,:),'filled','MarkerEdgeColor','k',...
                    'linewidth',1,'marker',markers{jdx},'MarkerFaceAlpha',0.5);
            end
            

            % Add boxes for targets
            ylim(axobj_ft(idx),'padded');
            ylmits = ylim;
            plot(axobj_ft(idx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)
            ylim(axobj_ft(idx),ylmits);
            if isfield(Bayesopt_Loss_Options,FT_Loss_Names{idx}) && Bayesopt_Loss_Options.(FT_Loss_Names{idx}) > sqrt(eps)
                rech = rectangle(axobj_ft(idx),'Position',[0.55 ylmits(1) N_Salts-0.1 (ylmits(2) - ylmits(1))],...
                        'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(FT_Loss_Names{idx}),4),'EdgeColor','r');
            end

            % Set plot properties
            set(axobj_ft(idx),'box','on','TickLabelInterpreter','latex');
            set(axobj_ft(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
            xticks(axobj_ft(idx),1:N_Salts)
            xticklabels(axobj_ft(idx),[])
            axobj_ft(idx).XAxis.TickLength = [0,0];
            ylim(axobj_ft(idx),ylmits)
            ylabel(axobj_ft(idx), FiniteTYLabels{idx},'Fontsize',fs-3,'Interpreter','latex')
            title(axobj_ft(idx),MP_Titles(idx),'Fontsize',fs,'Interpreter','latex')
            xlim(axobj_ft(idx),[0.5 N_Salts+0.5])
            set(axobj_ft(idx), 'YTickMode', 'auto')
            set(axobj_ft(idx),'xminorgrid','off','yminorgrid','on')
        end
    end

    if N_MinPlot_Rows > 0
        min_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 0 1 MinPlot_Frac],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]

        % Some widths and heights
        bbz = 1 - bbbuffer;
        hh = bbz/N_MinPlot_Rows - ttbuffer*bbz;
        % Create figure and axis
        axobj_min = gobjects(N_Structures,N_MinPlot_Rows);
        ij=1;

        for jdx = 1:N_MinPlot_Rows
            cylim = zeros(N_Structures,2);
            for idx = 1:N_Structures

                if strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE')
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures + 0.19*((1-wwbuffer)/N_Structures);
                    ww = (1/N_Structures)*(1-wwbuffer) - 0.22*(1/N_Structures);
                else
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures-0.005;
                    ww = (1/N_Structures)*(1-wwbuffer);
                end
                bb = 1 - (jdx)*bbz/N_MinPlot_Rows;

                axobj_min(idx,jdx) = subplot('Position',[ll bb ww hh],...
                    'parent',min_panel); % [left bottom width height]
                hold(axobj_min(idx,jdx),'on');

                % Gather data
                plot_data = squeeze(Min_PlotData(jdx,:,:,idx));
                salt_ave = mean(plot_data,2,'omitnan'); % average the data across salts

                % Plot data
                if include_av_line
                    plot(axobj_min(idx,jdx),xx,salt_ave,'--k','linewidth',2);
                end
                if strcmp(MinPlotTypes{jdx},'RLE') && ~strcmp(Structures{idx},'Rocksalt')
                    plot(axobj_min(idx,jdx),xx,Rel_Energies_DFT(:,idx),'-','linewidth',2,'color','b')
                    scatter(axobj_min(idx,jdx),xx,Rel_Energies_DFT(:,idx),200,Colours(xx,:),...
                        'filled','MarkerEdgeColor','b',...
                        'linewidth',2,'marker','o');
                end

                for kdx = 1:N_Models
                    yy = plot_data(:,kdx);
                    scatter(axobj_min(idx,jdx),xx,yy,200,Colours(xx,:),'filled','MarkerEdgeColor','k',...
                        'linewidth',1,'marker',markers{kdx},'MarkerFaceAlpha',0.5);
                end
                
                % Set plot properties
                set(axobj_min(idx,jdx),'box','on','TickLabelInterpreter','latex');
                set(axobj_min(idx,jdx),'XMinorTick','off','YMinorTick','on','FontSize',fs);
                xticks(axobj_min(idx,jdx),1:N_Salts)
                xticklabels(axobj_min(idx,jdx),[])
                axobj_min(idx,jdx).XAxis.TickLength = [0,0];
                ylim(axobj_min(idx,jdx),'padded');
                if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    cylim(idx,:) = ylim(axobj_min(idx,jdx));
                    plot(axobj_min(idx,jdx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)
                elseif (strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    csclylim = ylim(axobj_min(idx,jdx));
                    plot(axobj_min(idx,jdx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)
                    ylim(axobj_min(idx,jdx),csclylim);
                end
                plot(axobj_min(idx,jdx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)
                
                if idx == 1
                    ylabel(axobj_min(idx,jdx), MinPlotYLabels{jdx},'Fontsize',fs,'Interpreter','latex')
                end
                if jdx == N_MinPlot_Rows
                    xlabel(axobj_min(idx,jdx), LegStructure(Structures{idx}),'Fontsize',fs,'Interpreter','latex')
                end
                if idx == ceil(N_Structures/2)
                    title(axobj_min(idx,jdx),MinPlotTitles(jdx),'Fontsize',fs,'Interpreter','latex')
                    axobj_min(idx,jdx).TitleHorizontalAlignment = 'Left';
                    %axobj_min(idx,jdx).Title.Position = [1.073991031390134,49.12118639249388,0];
                end
                ij = ij+1;
            end
            ylmits = [min(cylim(:,1)) max(cylim(:,2))];
            for idx = 1:N_Structures
                if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    ylim(axobj_min(idx,jdx),ylmits);
                end
                set(axobj_min(idx,jdx), 'YTickMode', 'auto')
                set(axobj_min(idx,jdx),'xminorgrid','off','yminorgrid','on')
                if idx ~= 1 && ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    yticklabels(axobj_min(idx,jdx),[])
                end

                % Add boxes for targets
                if isfield(Bayesopt_Loss_Options,Structures{idx}) && Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}) > sqrt(eps)
                    rech = rectangle(axobj_min(idx,jdx),'Position',[0.55 ylmits(1) N_Salts-0.1 (ylmits(2) - ylmits(1))],...
                            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}),4),'EdgeColor','r');
                end
                xlim(axobj_min(idx,jdx),[0.5 N_Salts+0.5])
            end
        end

        h = gobjects(length(Salts)+1,1);
        for idx = 1:length(Salts)
            h(idx) = plot(axobj_min(1,1),nan, nan, 'o', 'MarkerSize', 16, 'DisplayName', [Salts{idx} '\quad'],...
                'MarkerFaceColor',Colours(idx,:),'linewidth',2,'Color','k');
        end
        h(end) = plot(axobj_min(1,1),nan, nan, 'o', 'MarkerSize', 16, 'DisplayName', 'DFT',...
            'MarkerFaceColor','w','linewidth',2,'Color','b');
        legh = legend(h,'Position',[0.36 0 0.33 0.06],'Orientation','Horizontal',...
            'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts+1);
    end
elseif N_Salts > 1 && BestOnly
    
    % Generate plotting panels
    xx = 1:N_Salts;
    row_height = 1/N_Rows;
    MinPlot_Frac = N_MinPlot_Rows/N_Rows;
    if plot_loss
        loss_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height 1 row_height],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        FTdrop = row_height;
        
        axobj_loss = subplot('Position',[0.1 0.05 0.8 0.65],...
            'parent',loss_panel); % [left bottom width height]
        hold(axobj_loss,'on');

        % plot data
        p = bar(axobj_loss,1,Total_loss,'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        for idx = 1:N_Salts
            p(idx).CData = Colours(idx,:);
        end

        % Set plot properties
        set(axobj_loss,'box','on','TickLabelInterpreter','latex');
        set(axobj_loss,'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
        xticks(axobj_loss,1)
        xticklabels(axobj_loss,[])
        axobj_loss.XAxis.TickLength = [0,0];
        %ylim(axobj_loss,'padded')
        set(axobj_loss,'YScale','log')
        ylim(axobj_loss,[10^(floor( log10(min(Total_loss)))) 10^(ceil( log10(max(Total_loss))))])
        
        if savefile
            sgtitle(loss_panel,'Minimized Objective Function $f\left(\mathbf{x}^{*}\right)$',...
                'Fontsize',fs,'Interpreter','latex')
        else
            sgtitle(loss_panel,['Minimized Objective Function: ' PubTheoryName ' Model ' ModelID ' [' Fix_Charge ' / ' Additivity ']'],...
                'Fontsize',fs,'Interpreter','latex')
        end
        ylabel(axobj_loss, '$f(\mathbf{x}^{*})$','Fontsize',fs,'Interpreter','latex')
        xlim(axobj_loss,'padded')
        set(axobj_loss, 'YTickMode', 'auto')
        set(axobj_loss,'xminorgrid','off','yminorgrid','on')
    else
        FTdrop = 0;
    end

    if plot_finite_T_data
        finite_T_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height-FTdrop 1 row_height],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        axobj_ft = gobjects(N_FiniteTTypes,1);

        ww0 = 1 - 0.065;
        bb0 = 0.14;
        ww1 = ww0/N_FiniteTTypes - ww0/(N_FiniteTTypes+2.3);
        ww2 = ww1/(N_FiniteTTypes-1);
        ww = ww0/(N_FiniteTTypes+2.3);
        hh = 1 - bb0 - 0.25;
        bb = bb0;
        for idx = 1:N_FiniteTTypes
            ll = (0.995-ww0) + (ww0/N_FiniteTTypes+ww2)*(idx-1);
            axobj_ft(idx) = subplot('Position',[ll bb ww hh],...
                'parent',finite_T_panel); % [left bottom width height]
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
            ylim(axobj_ft(idx),'padded');
            xlim(axobj_ft(idx),'padded');
            ylmits = ylim;
            xlmits = xlim;
            if isfield(Bayesopt_Loss_Options,FT_Loss_Names{idx}) && Bayesopt_Loss_Options.(FT_Loss_Names{idx}) > sqrt(eps)
                rech = rectangle(axobj_ft(idx),'Position',[xlmits(1) ylmits(1) (xlmits(2)-xlmits(1)) (ylmits(2) - ylmits(1))],...
                        'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(FT_Loss_Names{idx}),4),'EdgeColor','r');
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
    end
    
    
    if N_MinPlot_Rows > 0
        min_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 0 1 MinPlot_Frac],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]

        % Some widths and heights
        bbz = 1 - bbbuffer;
        hh = bbz/N_MinPlot_Rows - ttbuffer*bbz;
        % Create figure and axis
        axobj_min = gobjects(N_Structures,N_MinPlot_Rows);
        ij=1;

        for jdx = 1:N_MinPlot_Rows
            cylim = zeros(N_Structures,2);
            for idx = 1:N_Structures

                if strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE')
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures + 0.25*((1-wwbuffer)/N_Structures);
                    ww = (1/N_Structures)*(1-wwbuffer) - 0.25*(1/N_Structures) - 0.0025;
                else
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures - 0.005;
                    ww = (1/N_Structures)*(1-wwbuffer);
                end
                bb = 1 - (jdx)*bbz/N_MinPlot_Rows;

                axobj_min(idx,jdx) = subplot('Position',[ll bb ww hh],...
                    'parent',min_panel); % [left bottom width height]
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
                ylim(axobj_min(idx,jdx),'padded')
                xlim(axobj_min(idx,jdx),'padded')
                xlmits = xlim;

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
                if idx == 1
                    ylabel(axobj_min(idx,jdx), MinPlotYLabels{jdx},'Fontsize',fs,'Interpreter','latex')
                end
                if jdx == N_MinPlot_Rows
                    xlabel(axobj_min(idx,jdx), LegStructure(Structures{idx}),'Fontsize',fs,'Interpreter','latex')
                end
                if idx == ceil(N_Structures/2)
                    title(axobj_min(idx,jdx),MinPlotTitles(jdx),'Fontsize',fs,'Interpreter','latex')
                    axobj_min(idx,jdx).TitleHorizontalAlignment = 'Left';
                end
                ij = ij+1;
            end
            ylmits = [min(cylim(:,1)) max(cylim(:,2))];
            for idx = 1:N_Structures
                if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    ylim(axobj_min(idx,jdx),ylmits);
                end
                set(axobj_min(idx,jdx), 'YTickMode', 'auto')
                set(axobj_min(idx,jdx),'xminorgrid','off','yminorgrid','on')
                if idx ~= 1 && ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    yticklabels(axobj_min(idx,jdx),[])
                end

                % Add boxes for targets
                if isfield(Bayesopt_Loss_Options,Structures{idx}) && Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}) > sqrt(eps)
                    rech = rectangle(axobj_min(idx,jdx),'Position',[xlmits(1) ylmits(1) (xlmits(2) - xlmits(1)) (ylmits(2) - ylmits(1))],...
                            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}),4),'EdgeColor','r');
                end
                xlim(axobj_min(idx,jdx),xlmits)
            end
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
    legtxt{end+1}  = '\ \ DFT';
    legh = legend(h,legtxt,'Position',[0.36 0 0.33 0.06],'Orientation','Horizontal',...
        'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts+1);
else
    % Generate plotting panels
    xx = 1:N_Models;
    row_height = 1/N_Rows;
    MinPlot_Frac = N_MinPlot_Rows/N_Rows;
    if plot_loss
        loss_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height 1 row_height],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        FTdrop = row_height;
        
        axobj_loss = subplot('Position',[0.1 0.05 0.8 0.65],...
            'parent',loss_panel); % [left bottom width height]
        hold(axobj_loss,'on');

        % plot data
        p = bar(axobj_loss,1:N_Models,Total_loss,'FaceColor',Colours(1,:),'Visible','on','BarWidth',1,...
            'LineWidth',2);
        for mdx = 1:N_Models
            xpos = p.XEndPoints(mdx);
            scatter(axobj_loss,xpos,Total_loss(mdx),200,'MarkerEdgeColor','k',...
                'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx},'MarkerFaceAlpha',0.5);
        end

        % Set plot properties
        set(axobj_loss,'box','on','TickLabelInterpreter','latex');
        set(axobj_loss,'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
        xticks(axobj_loss,1:N_Models)
        xticklabels(axobj_loss,[])
        axobj_loss.XAxis.TickLength = [0,0];
        ylim(axobj_loss,'padded')
        if savefile
            sgtitle(loss_panel,'Minimized Objective Function $f\left(\mathbf{x}^{*}\right)$',...
                'Fontsize',fs,'Interpreter','latex')
        else
            sgtitle(loss_panel,['Minimized Objective Function: ' Salts{1} ' ' PubTheoryName ' Model ' ModelID ' [' Fix_Charge ' / ' Additivity ']'],...
                'Fontsize',fs,'Interpreter','latex')
        end
        ylabel(axobj_loss, '$f(\mathbf{x}^{*})$','Fontsize',fs,'Interpreter','latex')
        xlim(axobj_loss,[0 N_Models+1])
        set(axobj_loss, 'YTickMode', 'auto')
        set(axobj_loss,'xminorgrid','off','yminorgrid','on')
    else
        FTdrop = 0;
    end

    if plot_finite_T_data
        finite_T_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height-FTdrop 1 row_height],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]
        axobj_ft = gobjects(N_FiniteTTypes,1);

        ww0 = 1 - 0.05;
        bb0 = 0.14;
        ww1 = ww0/N_FiniteTTypes - ww0/(N_FiniteTTypes+2);
        ww2 = ww1/(N_FiniteTTypes-1);
        ww = ww0/(N_FiniteTTypes+2);
        hh = 1 - bb0 - 0.25;
        bb = bb0;
        for idx = 1:N_FiniteTTypes
            ll = (1-ww0) + (ww0/N_FiniteTTypes+ww2)*(idx-1);
            axobj_ft(idx) = subplot('Position',[ll bb ww hh],...
                'parent',finite_T_panel); % [left bottom width height]
            hold(axobj_ft(idx),'on');

            % Gather data
            plot_data = squeeze(Finite_T_PlotData(idx,:,:));
            
            % plot data
            p = bar(axobj_ft(idx),1:N_Models,plot_data,'FaceColor',Colours(1,:),'Visible','on','BarWidth',1,...
                'LineWidth',2);
            for mdx = 1:N_Models
                xpos = p.XEndPoints(mdx);
                scatter(axobj_ft(idx),xpos,plot_data(mdx),200,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx},'MarkerFaceAlpha',0.5);
            end

            % Add boxes for targets
            ylim(axobj_ft(idx),'padded');
            ylmits = ylim;
            if isfield(Bayesopt_Loss_Options,FT_Loss_Names{idx}) && Bayesopt_Loss_Options.(FT_Loss_Names{idx}) > sqrt(eps)
                rech = rectangle(axobj_ft(idx),'Position',[0.05 ylmits(1) N_Models+0.9 (ylmits(2) - ylmits(1))],...
                        'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(FT_Loss_Names{idx}),4),'EdgeColor','r');
            end

            % Set plot properties
            set(axobj_ft(idx),'box','on','TickLabelInterpreter','latex');
            set(axobj_ft(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
            xticks(axobj_ft(idx),1:N_Models)
            xticklabels(axobj_ft(idx),[])
            axobj_ft(idx).XAxis.TickLength = [0,0];
            ylim(axobj_ft(idx),ylmits)
            ylabel(axobj_ft(idx), FiniteTYLabels{idx},'Fontsize',fs-3,'Interpreter','latex')
            title(axobj_ft(idx),MP_Titles(idx),'Fontsize',fs,'Interpreter','latex')
            xlim(axobj_ft(idx),[0 N_Models+1])
            set(axobj_ft(idx), 'YTickMode', 'auto')
            set(axobj_ft(idx),'xminorgrid','off','yminorgrid','on')
        end
    end

    if N_MinPlot_Rows > 0
        min_panel = uipanel(figh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 0 1 MinPlot_Frac],...
            'AutoResizeChildren','off','BackgroundColor','w'); % [left bottom width height]

        % Some widths and heights
        bbbuffer = 0.08;
        bbz = 1 - bbbuffer;
        hh = bbz/N_MinPlot_Rows - ttbuffer*bbz;
        % Create figure and axis
        axobj_min = gobjects(N_Structures,N_MinPlot_Rows);
        ij=1;

        for jdx = 1:N_MinPlot_Rows
            cylim = zeros(N_Structures,2);
            for idx = 1:N_Structures

                if strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE')
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures + 0.25*((1-wwbuffer)/N_Structures);
                    ww = (1/N_Structures)*(1-wwbuffer) - 0.25*(1/N_Structures);
                else
                    ll = wwbuffer + (1-wwbuffer)*(idx-1)/N_Structures;
                    ww = (1/N_Structures)*(1-wwbuffer);
                end
                bb = 1 - (jdx)*bbz/N_MinPlot_Rows;

                axobj_min(idx,jdx) = subplot('Position',[ll bb ww hh],...
                    'parent',min_panel); % [left bottom width height]
                hold(axobj_min(idx,jdx),'on');

                % Gather data
                plot_data = squeeze(Min_PlotData(jdx,:,:,idx));
                
                % plot data
                p = bar(axobj_min(idx,jdx),1:N_Models,plot_data,'FaceColor',Colours(1,:),'Visible','on','BarWidth',1,...
                    'LineWidth',2);
                for mdx = 1:N_Models
                    xpos = p.XEndPoints(mdx);
                    scatter(axobj_min(idx,jdx),xpos,plot_data(mdx),200,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx},'MarkerFaceAlpha',0.5);
                end
                
                if strcmp(MinPlotTypes{jdx},'RLE') && ~strcmp(Structures{idx},'Rocksalt')
                    plot(axobj_min(idx,jdx),[0 N_Models+1],[Rel_Energies_DFT(:,idx) Rel_Energies_DFT(:,idx)],'-.b','linewidth',2)
                end

                % Set plot properties
                set(axobj_min(idx,jdx),'box','on','TickLabelInterpreter','latex');
                set(axobj_min(idx,jdx),'XMinorTick','off','YMinorTick','on','FontSize',fs);
                xticks(axobj_min(idx,jdx),1:N_Models)
                xticklabels(axobj_min(idx,jdx),[])
                axobj_min(idx,jdx).XAxis.TickLength = [0,0];
                ylim(axobj_min(idx,jdx),'padded')

                if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    cylim(idx,:) = ylim(axobj_min(idx,jdx));
                end
                if idx == 1
                    ylabel(axobj_min(idx,jdx), MinPlotYLabels{jdx},'Fontsize',fs,'Interpreter','latex')
                end
                if jdx == N_MinPlot_Rows
                    xlabel(axobj_min(idx,jdx), LegStructure(Structures{idx}),'Fontsize',fs,'Interpreter','latex')
                end
                if idx == ceil(N_Structures/2)
                    title(axobj_min(idx,jdx),MinPlotTitles(jdx),'Fontsize',fs,'Interpreter','latex')
                    axobj_min(idx,jdx).TitleHorizontalAlignment = 'Left';
                end
                ij = ij+1;
            end
            ylmits = [min(cylim(:,1)) max(cylim(:,2))];
            for idx = 1:N_Structures
                if ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    ylim(axobj_min(idx,jdx),ylmits);
                end
                set(axobj_min(idx,jdx), 'YTickMode', 'auto')
                set(axobj_min(idx,jdx),'xminorgrid','off','yminorgrid','on')
                if idx ~= 1 && ~(strcmp(Structures{idx},'CsCl') && strcmp(MinPlotTypes{jdx},'RLE'))
                    yticklabels(axobj_min(idx,jdx),[])
                end

                % Add boxes for targets
                if isfield(Bayesopt_Loss_Options,Structures{idx}) && Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}) > sqrt(eps)
                    rech = rectangle(axobj_min(idx,jdx),'Position',[0.05 ylmits(1) N_Models+0.9 (ylmits(2) - ylmits(1))],...
                            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}),4),'EdgeColor','r');
                end
                xlim(axobj_min(idx,jdx),[0 N_Models+1])
            end
        end
    end
end





if savefile
    
    %set(mainpanelh,'renderer','opengl')
    %saveas(figh,filename,'pdf')
    %print(mainpanelh,filename,'-bestfit')
%    exportgraphics(mainpanelh,fullfile(saveloc,filename))
    
%     figh.PaperPositionMode='manual';
%     figh.PaperUnits ='normalized';
%     figh.PaperPosition = [0 0 1 1];
%     exportgraphics(mainpanelh,filename)
    %exportapp(figh,filename)
    %filename = ['Model_Overview_' Theory '_Model_' ModelID '.eps'];
    %hgexport(figh,fullfile(saveloc,filename))
    %print(figh,fullfile(saveloc,filename),'-deps','-noui')
    figh=gcf;
    filename = ['Model_Overview_' Theory '_Model_' ModelID '.png'];
    print(figh,fullfile(saveloc,filename),'-dpng','-r300','-noui')
end


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

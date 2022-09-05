clear; %#ok<*UNRCH>
%% Data options
Salts = {'NaCl'};
Theory = 'JA';
ModelID = 'LA';
Reps = [1:5];
savefile = true; % switch to save the final plots to file
saveloc = 'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Thesis\Thesis_Draft\BO_Figures';

%% Plot options
fs = 22; % font size
markers = {'o' 's' '^' 'v' 'd' '>' '<' 'p' 'h' 'x'};
show_as_percent_error = false; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)
wwbuffer = 0.055;
ttbuffer = 0.08;
bbbuffer = 0.14;

plot_LE = true;
plot_RLE = true;
plot_a = true;
plot_c = true;
plot_ac = false;
plot_volume = false;
plot_density = false;
plot_loss = true;
plot_finite_T_data = true;

%% Script begins
cm3_per_Ang3 = 1e-24; % cubic cm/cubic angstrom
N_A = 6.02214076e23; % formula units/mol
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

Structures   = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
MinPlotTypes = {'LE' 'RLE' 'a' 'c' 'ac' 'V' 'Density'};
MinPlotIdxes = [plot_LE plot_RLE plot_a plot_c plot_ac plot_volume plot_density];
MinPlotTypes = MinPlotTypes(MinPlotIdxes);
FiniteTTypes = {'MP' 'dH' 'dV' 'LiqV' 'SolV' 'DM'};
MinPlotTitles = {'Error in $E_L$: $\left( E_{L} - E_{L}^{*} \right) / \left| E_{L}^{*} \right|$' ...
    '$E_{L}$[Structure] - $E_{L}$[Rocksalt]' ...
    'Error in $a$: $\left( a - a^{*} \right) / \left| a^{*} \right|$' ...
    'Error in $c$: $\left( c - c^{*} \right) / \left| c^{*} \right|$' ...
    'Error in $c/a$: $\left( c/a - c^{*}/a^{*} \right) / \left| c^{*}/a^{*} \right|$' ...
    'Error in $\rho$: $\left( \rho - \rho^{*} \right) / \left| \rho^{*} \right|$' ...
    'Error in $V$: $\left( V - V^{*} \right) / \left| V^{*} \right|$'};
MinPlotTitles = MinPlotTitles(MinPlotIdxes);
MP_Titles = {'MP Er.' '$\Delta H_{\mathrm{Fus.}}$ Er.' '$\Delta V_{\mathrm{Fus.}}$ Er.' ...
    '$V_{\mathrm{Liq.}}$($T_{m}$) Er.' '$V_{\mathrm{Sol.}}$($T_{m}$) Er.' 'D$_{\mathrm{Li}}$($T_{m}$) Er.'};
FT_Loss_Names = {'MP' 'Fusion_Enthalpy' 'MP_Volume_Change' 'Liquid_MP_Volume' 'Solid_MP_Volume' 'Liquid_DM_MP'};

if show_as_percent_error    
    MinPlotYLabels = {'[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
    
    FiniteTYLabels = {'[\%]' '[\%]' '[\%]' '[\%]' '[\%]' '[\%]'};
else
    MinPlotYLabels = {'[kJ mol$^{-1}$]' ...
                      '[kJ mol$^{-1}$]' ...
                      '[\AA]' ...
                      '[\AA]' ...
                      '[$(c/a)$]' ...
                      '[$\AA^{3}$]' ...
                      'g cm$^{-3}$]'};
                  
    FiniteTYLabels = {'[K]' '[kJ mol$^{-1}$]' '[$\AA^{3}$]' '[$\AA^{3}$]' '[$\AA^{3}$]' '[$10^{-5}$ cm$^{2}$/s]'};
end
MinPlotYLabels = MinPlotYLabels(MinPlotIdxes);

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
                Bayesopt_model = functions(data.bayesopt_results.ObjectiveFcn).workspace{1}.Model;
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
    error(['Unable to load targets for ' Theory ' Model ' ModelID])
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
                Minimization_Data = data.Minimization_Data;
                Total_loss(idx,iidx) = data.loss;
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
                                    Min_PlotData(rowidx,idx,iidx,jdx) = (Density - DFT.(Salt).(Structure).density)*100/ ...
                                        DFT.(Salt).(Structure).density;
                                else
                                    Min_PlotData(rowidx,idx,iidx,jdx) = (Volume - DFT.(Salt).(Structure).V)*100/ ...
                                        DFT.(Salt).(Structure).V;
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
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Density - DFT.(Salt).(Structure).density;
                                else
                                    Min_PlotData(rowidx,idx,iidx,jdx) = Volume - DFT.(Salt).(Structure).V;
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
                                    (Finite_T_Data.Liquid_DM_MP - Finite_T_Data.Exp_DM_MP).*1e5;
                        end
                    end
                end
            end
        end
    end
end

%% Plotting
% Gen figure
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
mainpanelh = uipanel(figh,'BorderType','none',...
    'Position',[0 0 1 1]);


if N_Salts > 1
    % Generate plotting panels
    xx = 1:N_Salts;
    row_height = 1/N_Rows;
    MinPlot_Frac = N_MinPlot_Rows/N_Rows;
    if plot_loss
        loss_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height 1 row_height]); % [left bottom width height]
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
            ll = (1-ww0) + (ww0/N_Salts+ww2)*(idx-1);
            axobj_loss(idx) = subplot('Position',[ll bb ww hh],...
                'parent',loss_panel); % [left bottom width height]
            hold(axobj_loss(idx),'on');


            % Gather data
            plot_data = squeeze(Total_loss(idx,:));

            p = bar(axobj_loss(idx),1:N_Models,plot_data,'FaceColor',Colours(idx,:),'Visible','on','BarWidth',1,...
                'LineWidth',2);
            for mdx = 1:N_Models
                xpos = p.XEndPoints(mdx);
                scatter(axobj_loss(idx),xpos,plot_data(mdx),100,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(idx,:),'linewidth',2,'marker',markers{mdx});
            end

            % Set plot properties
            set(axobj_loss(idx),'box','on','TickLabelInterpreter','latex');
            set(axobj_loss(idx),'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
            xticks(axobj_loss(idx),1:N_Models)
            xticklabels(axobj_loss(idx),[])
            axobj_loss(idx).XAxis.TickLength = [0,0];
            ylim(axobj_loss(idx),'padded')
            sgtitle(loss_panel,['Minimized Objective Function: ' Theory ' Model ' ModelID ' [' Fix_Charge ' / ' Additivity ']'],...
                'Fontsize',fs,'Interpreter','latex')
            if idx == 1
                ylabel(axobj_loss(idx), '$f(\mathbf{x})$','Fontsize',fs,'Interpreter','latex')
            end
            xlim(axobj_loss(idx),[0 N_Models+1])
            set(axobj_loss(idx), 'YTickMode', 'auto')
            set(axobj_loss(idx),'xminorgrid','off','yminorgrid','on')
        end
    else
        FTdrop = 0;
    end

    if plot_finite_T_data
        finite_T_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height-FTdrop 1 row_height]); % [left bottom width height]
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
            salt_ave = mean(plot_data,2,'omitnan'); % average the data across salts

            % Plot data
            plot(axobj_ft(idx),xx,salt_ave,'--k','linewidth',2);
            for jdx = 1:N_Models
                yy = plot_data(jdx,:);
                scatter(axobj_ft(idx),xx,yy,100,Colours(xx,:),'filled','MarkerEdgeColor','k',...
                    'linewidth',1,'marker',markers{jdx});
            end
            plot(axobj_ft(idx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)

            % Add boxes for targets
            ylim(axobj_ft(idx),'padded');
            ylmits = ylim;
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
        min_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 0 1 MinPlot_Frac]); % [left bottom width height]

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
                salt_ave = mean(plot_data,2,'omitnan'); % average the data across salts

                % Plot data
                plot(axobj_min(idx,jdx),xx,salt_ave,'--k','linewidth',2);
                if strcmp(MinPlotTypes{jdx},'RLE') && ~strcmp(Structures{idx},'Rocksalt')
                    plot(axobj_min(idx,jdx),xx,Rel_Energies_DFT(:,idx),'-','linewidth',2,'color','r')
                    scatter(axobj_min(idx,jdx),xx,Rel_Energies_DFT(:,idx),100,Colours(xx,:),...
                        'filled','MarkerEdgeColor','r',...
                        'linewidth',2,'marker','o');
                end

                for kdx = 1:N_Models
                    yy = plot_data(kdx,:);
                    scatter(axobj_min(idx,jdx),xx,yy,100,Colours(xx,:),'filled','MarkerEdgeColor','k',...
                        'linewidth',1,'marker',markers{kdx});
                end
                plot(axobj_min(idx,jdx),[0.5 N_Salts+0.5],[0 0],':k','linewidth',1)

                % Set plot properties
                set(axobj_min(idx,jdx),'box','on','TickLabelInterpreter','latex');
                set(axobj_min(idx,jdx),'XMinorTick','off','YMinorTick','on','FontSize',fs);
                xticks(axobj_min(idx,jdx),1:N_Salts)
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
                    xlabel(axobj_min(idx,jdx), Structures{idx},'Fontsize',fs,'Interpreter','latex')
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
                    rech = rectangle(axobj_min(idx,jdx),'Position',[0.55 ylmits(1) N_Salts-0.1 (ylmits(2) - ylmits(1))],...
                            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}),4),'EdgeColor','r');
                end
                xlim(axobj_min(idx,jdx),[0.5 N_Salts+0.5])
            end
        end

        h = gobjects(length(Salts)+1,1);
        for idx = 1:length(Salts)
            h(idx) = plot(axobj_min(1,1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', Salts{idx},...
                'MarkerFaceColor',Colours(idx,:),'linewidth',2,'Color','k');
        end
        h(end) = plot(axobj_min(1,1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'DFT',...
            'MarkerFaceColor','w','linewidth',2,'Color','r');
        legh = legend(h,'Position',[0.36 0 0.33 0.06],'Orientation','Horizontal',...
            'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts+1);
    end
else
    
    
    
    
    
    
    % Generate plotting panels
    xx = 1:N_Models;
    row_height = 1/N_Rows;
    MinPlot_Frac = N_MinPlot_Rows/N_Rows;
    if plot_loss
        loss_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height 1 row_height]); % [left bottom width height]
        FTdrop = row_height;
        
        axobj_loss = subplot('Position',[0.1 0.05 0.8 0.65],...
            'parent',loss_panel); % [left bottom width height]
        hold(axobj_loss,'on');

        % plot data
        p = bar(axobj_loss,1:N_Models,Total_loss,'FaceColor',Colours(1,:),'Visible','on','BarWidth',1,...
            'LineWidth',2);
        for mdx = 1:N_Models
            xpos = p.XEndPoints(mdx);
            scatter(axobj_loss,xpos,Total_loss(mdx),100,'MarkerEdgeColor','k',...
                'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx});
        end

        % Set plot properties
        set(axobj_loss,'box','on','TickLabelInterpreter','latex');
        set(axobj_loss,'XMinorTick','off','YMinorTick','on','FontSize',fs-3);
        xticks(axobj_loss,1:N_Models)
        xticklabels(axobj_loss,[])
        axobj_loss.XAxis.TickLength = [0,0];
        ylim(axobj_loss,'padded')
        sgtitle(loss_panel,['Minimized Objective Function: ' Theory ' Model ' ModelID ' [' Fix_Charge ' / ' Additivity ']'],...
            'Fontsize',fs,'Interpreter','latex')
        ylabel(axobj_loss, '$f(\mathbf{x})$','Fontsize',fs,'Interpreter','latex')
        xlim(axobj_loss,[0 N_Models+1])
        set(axobj_loss, 'YTickMode', 'auto')
        set(axobj_loss,'xminorgrid','off','yminorgrid','on')
    else
        FTdrop = 0;
    end

    if plot_finite_T_data
        finite_T_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 1-row_height-FTdrop 1 row_height]); % [left bottom width height]
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
                scatter(axobj_ft(idx),xpos,plot_data(mdx),100,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx});
            end

            % Add boxes for targets
            ylim(axobj_ft(idx),'padded');
            ylmits = ylim;
            if isfield(Bayesopt_Loss_Options,FT_Loss_Names{idx}) && Bayesopt_Loss_Options.(FT_Loss_Names{idx}) > sqrt(eps)
                rech = rectangle(axobj_ft(idx),'Position',[0.05 ylmits(1) N_Models+0.95 (ylmits(2) - ylmits(1))],...
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
        min_panel = uipanel(mainpanelh,'FontSize',fs,...
            'BorderType','none',...
            'Position',[0 0 1 MinPlot_Frac]); % [left bottom width height]

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
                    scatter(axobj_min(idx,jdx),xpos,plot_data(mdx),100,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',Colours(1,:),'linewidth',2,'marker',markers{mdx});
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
                    xlabel(axobj_min(idx,jdx), Structures{idx},'Fontsize',fs,'Interpreter','latex')
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
                    rech = rectangle(axobj_min(idx,jdx),'Position',[0.05 ylmits(1) N_Models+0.95 (ylmits(2) - ylmits(1))],...
                            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structures{idx}).(MinPlotTypes{jdx}),4),'EdgeColor','r');
                end
                xlim(axobj_min(idx,jdx),[0 N_Models+1])
            end
        end
    end
end

if savefile
    filename = ['Model_Overview_' Theory '_Model_' ModelID '.pdf'];
    %set(mainpanelh,'renderer','opengl')
    %saveas(figh,filename,'pdf')
    %print(mainpanelh,filename,'-bestfit')
    exportgraphics(mainpanelh,fullfile(saveloc,filename))
end

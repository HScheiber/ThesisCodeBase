% Script for plotting an error comparison for multiple models, to show how
% each performs.
%#ok<*UNRCH>

%% Analysis Parameters:

% Data options
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Molar_masses = [25.939 42.394 86.845 133.85];  % g/mol
Theory = 'JC';
Basenum = 'E';
Midnum = 'D';
savefile = true; % switch to save the final plots to file
filename = ['Target_Compare_' Theory '_' Basenum  Midnum '.png'];

% Plot options
Structures = {'Rocksalt' 'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
plot_LE = true;
plot_RLE = true;
plot_a = true;
plot_c = true;
plot_ac = false;
plot_volume = false;
plot_density = false;
plot_loss = true;
show_targets = true;
fs = 22; % font size

% Scale parameters
set_ylim = false; % auto sets y-limits when false
LE_ylim = [-5 10];
RLE_ylim = [-5 35];
a_ylim = [-6 10];
c_ylim = [-6 10];
ca_ylim = [-Inf Inf];
density_ylim = [-Inf Inf];
volume_ylim = [-8 20];
loss_ylog = true;

% Other parameters
show_as_percent_error = true; % Plot as percent error. If false, plot as the numerical error value (i.e. including units)
% Calculates error in lattice energy with respect to experimental lattice energies when available
Target_Experimental_Energies = true; 
% Calculates error in lattice parameters, volumes, and densities with respect to experimental lattice parameters when available
Target_Experimental_latpar = true;
Exp_latpart_zero = true; % When true, the experimental lattice param is at 0 K, when false, it is at 300 K
BestOnly = false; % If true, include only the lowest loss model for each salt
Title_txt = '';
cm3_per_Ang3 = 1e-24; % cubic cm/cubic angstrom
N_A = 6.02214076e23; % formula units/mol

shift_CsCl_RLE = false;

%% Script begins
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

% Find reps of models
files = {dir(fullfile(ML_results_dir,['*' Theory '_Model_' Basenum Midnum '*'])).name};
Models = unique(regexp(files,[Basenum Midnum '([0-9])+'],'match','once'));
nonempty_idx = ~cellfun(@isempty,Models);
Models = Models(nonempty_idx);

% Make sure the models are sorted
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end


% Assign a marker to each rep
markers = {'o' 's' '^' 'v' 'd' '>' '<' 'p' 'h' 'x'};
N_Models = length(Models);
N_Markers = length(markers);
if N_Markers < N_Models
    Nrep = ceil(N_Models/N_Markers);
    markers = repmat(markers,1,Nrep);
end

if plot_LE
    LE_Structures = Structures;
else
    LE_Structures = {};  
end
if plot_RLE
    RLE_Structures = Structures(~strcmpi(Structures,'Rocksalt'));
else
    RLE_Structures = {};
end
if plot_a
    a_Structures = Structures;
else
    a_Structures = {};
end
if plot_c
    c_Structures = setdiff(Structures,{'Rocksalt' 'CsCl' 'Sphalerite'},'stable');
else
    c_Structures = {};
end
if plot_ac
    ca_Structures = setdiff(Structures,{'Rocksalt' 'CsCl' 'Sphalerite'},'stable');
else
    ca_Structures = {};
end
if plot_density
    density_Structures = Structures;
else
    density_Structures = {};
end
if plot_volume
    volume_Structures = Structures;
else
    volume_Structures = {};
end

%% Pre-allocate data arrays
N_Salts = length(Salts); % number of models plus the DFT results
N_Structures = length(Structures);

Energies = nan(N_Salts,N_Models,N_Structures);
Rel_Energies = nan(N_Salts,N_Models,N_Structures);
Rel_Energies_Target = nan(N_Salts,N_Structures);
LatPars_a = nan(N_Salts,N_Models,N_Structures);
LatPars_c = nan(N_Salts,N_Models,N_Structures);
LatPars_ca = nan(N_Salts,N_Models,N_Structures);
Densities = nan(N_Salts,N_Models,N_Structures);
Volumes = nan(N_Salts,N_Models,N_Structures);
Total_loss = nan(N_Salts,N_Models);

% Plot color scheme
Colours = cbrewer('qual','Set3',N_Salts);
%Colours = [0 0 0; cbrewer('qual','Paired',N_Models)];

Geometry = struct;
for idx = 1:length(Structures)
    Structure = Structures{idx};
    Settings = Initialize_MD_Settings;
    Settings.Salt = 'LiF';
    Settings.Structure = Structure;
    Geometry.(Structure) = Default_Crystal(Settings);
end

% Load target data
DFT = Load_Best_DFT_Data;

if Target_Experimental_Energies || Target_Experimental_latpar
    EXP = Load_Experimental_Data;
    
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        
        if Target_Experimental_Energies
            Correction = EXP.(Salt).Rocksalt.E - DFT.(Salt).Rocksalt.Energy;
            for jdx = 1:length(Structures)
                DFT.(Salt).(Structures{jdx}).Energy = DFT.(Salt).(Structures{jdx}).Energy + Correction;
            end
        end

        if Target_Experimental_latpar && Exp_latpart_zero
            DFT.(Salt).Rocksalt.a = EXP.(Salt).Rocksalt.a_zero;
            DFT.(Salt).Rocksalt.b = EXP.(Salt).Rocksalt.b_zero;
            DFT.(Salt).Rocksalt.c = EXP.(Salt).Rocksalt.c_zero;
            DFT.(Salt).Rocksalt.V = EXP.(Salt).Rocksalt.V_zero;
            DFT.(Salt).Rocksalt.density = EXP.(Salt).Rocksalt.density_zero;
            if isfield(EXP.(Salt),'Wurtzite') && ~isnan(EXP.(Salt).Wurtzite.a)
                DFT.(Salt).Wurtzite.a = EXP.(Salt).Wurtzite.a_zero;
                DFT.(Salt).Wurtzite.b = EXP.(Salt).Wurtzite.b_zero;
                DFT.(Salt).Wurtzite.c = EXP.(Salt).Wurtzite.c_zero;
                DFT.(Salt).Wurtzite.V = EXP.(Salt).Wurtzite.V_zero;
                DFT.(Salt).Wurtzite.density = EXP.(Salt).Wurtzite.density_zero;
            end
        elseif Target_Experimental_latpar && ~Exp_latpart_zero
            DFT.(Salt).Rocksalt.a = EXP.(Salt).Rocksalt.a;
            DFT.(Salt).Rocksalt.b = EXP.(Salt).Rocksalt.b;
            DFT.(Salt).Rocksalt.c = EXP.(Salt).Rocksalt.c;
            DFT.(Salt).Rocksalt.V = EXP.(Salt).Rocksalt.V;
            DFT.(Salt).Rocksalt.density = EXP.(Salt).Rocksalt.density;
            if isfield(EXP.(Salt),'Wurtzite') && ~isnan(EXP.(Salt).Wurtzite.a)
                DFT.(Salt).Wurtzite.a = EXP.(Salt).Wurtzite.a;
                DFT.(Salt).Wurtzite.b = EXP.(Salt).Wurtzite.b;
                DFT.(Salt).Wurtzite.c = EXP.(Salt).Wurtzite.c;
                DFT.(Salt).Wurtzite.V = EXP.(Salt).Wurtzite.V;
                DFT.(Salt).Wurtzite.density = EXP.(Salt).Wurtzite.density;
            end
        end
    end
    clear('EXP')
end


% Calculate target relative lattice energies
if plot_RLE
    for idx = 1:N_Salts
        Salt = Salts{idx};
        for jdx = 1:N_Structures
            Structure = Structures{jdx};
            
            if shift_CsCl_RLE && strcmpi(Structure,'cscl')
                RLE_Shift = -50;
            else
                RLE_Shift = 0;
            end
            
            % Relative Energies
            if contained_in_cell(Structure,RLE_Structures)
                Target_gap = DFT.(Salt).(Structure).Energy - DFT.(Salt).Rocksalt.Energy;
                Rel_Energies_Target(idx,jdx) = Target_gap + RLE_Shift;
            end
        end
    end
end


% Find the model targets (ASSUMING THEY ARE ALL THE SAME!!)
if show_targets
    data_found = false;
    for idx = 1:N_Salts
        for jdx = 1:N_Models
            sres = dir([ML_results_dir filesep  Salts{idx} '_' Theory '_Model_' Models{jdx} '_bayesopt.mat']);

            if isempty(sres)
                continue
            else
                model_filename = fullfile(sres.folder,sres.name);

                % Load data
                try
                    data = load(model_filename);
                    Bayesopt_Loss_Options = functions(data.results.ObjectiveFcn).workspace{1}.Model.Loss_Options;
                    data_found = true;
                    break
                catch
                    %warning(['Unable to load Bayesopt data for ' Salts{1} ' ' Theory ' Model ' Models{1}])
                    continue
                end
            end
        end
        if data_found
            break
        end
    end
    if ~data_found
        warning(['Unable to load Bayesopt data for ' Theory ' Model ' Models{1}(1:2)])
    end
else
    Bayesopt_Loss_Options = init_loss_options;
end


% Loop over chosen models
for idx = 1:N_Salts
    Salt = Salts{idx};
    MM = Molar_masses(idx);
    for iidx = 1:N_Models
        Model = Models{iidx};

        % Find the fully optimized model
        sres = dir([ML_results_dir filesep '*' Salt '*' Theory '*' 'Model_' Model '_*fullopt.mat']);

        if length(sres) > 1
            warning(['Multiple results found for: ' Salt ', ' Theory ', Model ' Model '. Using first result.']);
            sres = sres(1);
        elseif isempty(sres)
            warning(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end

        model_filename = fullfile(sres.folder,sres.name);

        % Load data
        data = load(model_filename);
        try
            Minimization_Data = data.Minimization_Data;
            Loss = data.loss;
            clear('data')
        catch
            warning(['No saved crystal data for ' Salt ' ' Theory ' Model ' Model])
            continue
        end
        
        % Grab available structures from data
        Min_Dat_Structs = cell(1,length(Minimization_Data));
        for jdx = 1:length(Minimization_Data)
            Min_Dat_Structs{jdx} = Minimization_Data{jdx}.Structure;
        end

        % Loop through structures
        for jdx = 1:length(Structures)
            Structure = Structures{jdx};
            
            if shift_CsCl_RLE && strcmpi(Structure,'cscl')
                RLE_Shift = -50;
            else
                RLE_Shift = 0;
            end
            
            if ~any(strcmpi(Min_Dat_Structs,Structure))
                continue
            else
                md_idx = find(strcmpi(Min_Dat_Structs,Structure));
                md_idx = md_idx(1);
            end
            
            % Absolute energies
            if contained_in_cell(Structure,LE_Structures)
                if show_as_percent_error
                    Energies(idx,iidx,jdx) = (Minimization_Data{md_idx}.E - ...
                        DFT.(Salt).(Structure).Energy)*100/ ...
                        abs(DFT.(Salt).(Structure).Energy);
                else
                    Energies(idx,iidx,jdx) = Minimization_Data{md_idx}.E - ...
                        DFT.(Salt).(Structure).Energy; 
                end
                disp([Salt ' ' Model ' ' Structure ' Energy [kJ/mol]: ' num2str(Minimization_Data{md_idx}.E,'%.15f')])
            end

            % Relative Energies
            if contained_in_cell(Structure,RLE_Structures)
                
                Model_gap  = Minimization_Data{md_idx}.E - ...
                             Minimization_Data{strcmpi(Min_Dat_Structs,'Rocksalt')}.E;
                
                Rel_Energies(idx,iidx,jdx) = Model_gap + RLE_Shift;
            end

            % lattice param a
            if contained_in_cell(Structure,a_Structures)
                if show_as_percent_error
                    LatPars_a(idx,iidx,jdx) = (Minimization_Data{md_idx}.a - DFT.(Salt).(Structure).a)*100/ ...
                        DFT.(Salt).(Structure).a;
                else
                    LatPars_a(idx,iidx,jdx) = Minimization_Data{md_idx}.a - DFT.(Salt).(Structure).a;
                end
                disp([Salt ' ' Model ' ' Structure ' a [Angstrom]: ' num2str(Minimization_Data{md_idx}.a,'%.15f')])
            end

            % lattice param c
            if contained_in_cell(Structure,c_Structures)
                if show_as_percent_error
                    LatPars_c(idx,iidx,jdx) = (Minimization_Data{md_idx}.c - DFT.(Salt).(Structure).c)*100/ ...
                        DFT.(Salt).(Structure).c;
                else
                    LatPars_c(idx,iidx,jdx) = Minimization_Data{md_idx}.c - DFT.(Salt).(Structure).c;
                end
            end
            
            % c/a
            if contained_in_cell(Structure,ca_Structures)
                Target_ca = DFT.(Salt).(Structure).c / DFT.(Salt).(Structure).a;

                Model_ca  = Minimization_Data{md_idx}.c / ...
                             Minimization_Data{md_idx}.a;

                if show_as_percent_error
                    LatPars_ca(idx,iidx,jdx) = (Model_ca - Target_ca)*100 / Target_ca;
                else
                    LatPars_ca(idx,iidx,jdx) = Model_ca - Target_ca;
                end
            end
            
            
            % densities
            switch Structure
                case 'Rocksalt'
                    NF_per_Cell = 4;
                case 'Wurtzite'
                    NF_per_Cell = 2;
                case 'CsCl'
                    NF_per_Cell = 1;
                case 'FiveFive'
                    NF_per_Cell = 4;
                case 'BetaBeO'
                    NF_per_Cell = 4;
                case 'Sphalerite'
                    NF_per_Cell = 4;
                case {'NiAs' 'AntiNiAs'}
                    NF_per_Cell = 2;
            end
            
            ABC_vec = Geometry.(Structure).Transform.*[Minimization_Data{md_idx}.a; Minimization_Data{md_idx}.b; Minimization_Data{md_idx}.c];
            Volume = (dot(cross(ABC_vec(1,:),ABC_vec(2,:)),ABC_vec(3,:)))/NF_per_Cell; % Angstrom^3 / Formula Unit
            Density = MM/(Volume*cm3_per_Ang3*N_A); % g/cm^3
            
            if contained_in_cell(Structure,density_Structures)
                if show_as_percent_error
                    Densities(idx,iidx,jdx) = (Density - DFT.(Salt).(Structure).density)*100/ ...
                        DFT.(Salt).(Structure).density;
                else
                    Densities(idx,iidx,jdx) = Density - DFT.(Salt).(Structure).density;
                end
            end
            
            % Volumes
            if contained_in_cell(Structure,volume_Structures)
                if show_as_percent_error
                    Volumes(idx,iidx,jdx) = (Volume - DFT.(Salt).(Structure).V)*100/ ...
                        DFT.(Salt).(Structure).V;
                else
                    Volumes(idx,iidx,jdx) = Volume - DFT.(Salt).(Structure).V;
                end
                disp([Salt ' ' Model ' ' Structure ' V [Angstrom^3/Formula Unit]: ' num2str(Volume,'%.15f')])
            end
        end

        % total loss
        if plot_loss || BestOnly
            Total_loss(idx,iidx) = Loss;
        end
    end
end

if BestOnly
    [Total_loss,min_idx] = min(Total_loss,[],2);
    
    Energies_base = zeros(size(Energies,1),1,size(Energies,3));
    Rel_Energies_base = zeros(size(Rel_Energies,1),1,size(Rel_Energies,3));
    LatPars_a_base = zeros(size(LatPars_a,1),1,size(LatPars_a,3));
    LatPars_c_base = zeros(size(LatPars_c,1),1,size(LatPars_c,3));
    LatPars_ca_base = zeros(size(LatPars_ca,1),1,size(LatPars_ca,3));
    Densities_base = zeros(size(Densities,1),1,size(Densities,3));
    Volumes_base = zeros(size(Volumes,1),1,size(Volumes,3));
    
    for idx = 1:N_Salts
        Energies_base(idx,:,:) = Energies(idx,min_idx(idx),:);
        Rel_Energies_base(idx,:,:) = Rel_Energies(idx,min_idx(idx),:);
        LatPars_a_base(idx,:,:) = LatPars_a(idx,min_idx(idx),:);
        LatPars_c_base(idx,:,:) = LatPars_c(idx,min_idx(idx),:);
        LatPars_ca_base(idx,:,:) = LatPars_ca(idx,min_idx(idx),:);
        Densities_base(idx,:,:) = Densities(idx,min_idx(idx),:);
        Volumes_base(idx,:,:) = Volumes(idx,min_idx(idx),:);
    end
    
    Energies = Energies_base;
    Rel_Energies = Rel_Energies_base;
    LatPars_a = LatPars_a_base;
    LatPars_c = LatPars_c_base;
    LatPars_ca = LatPars_ca_base;
    Densities = Densities_base;
    Volumes = Volumes_base;
    
    N_Models = 1;
end

plot_switches = [plot_loss plot_LE plot_RLE plot_a plot_c plot_ac plot_density plot_volume];
N_compares = sum(plot_switches);
all_y = {Total_loss, Energies, Rel_Energies, LatPars_a, LatPars_c, LatPars_ca Densities Volumes};
all_ylim = {[],LE_ylim, RLE_ylim, a_ylim, c_ylim, ca_ylim density_ylim volume_ylim};
plot_types = {'loss' 'LE' 'RLE' 'a' 'c' 'ca' 'density' 'V'};

bar_x = 1:N_Structures;

if savefile
    spec_ther = '';
else
    spec_ther = [': ' Theory ' Model ' Basenum Midnum];
end

Titles = {['Minimized Objective Function' spec_ther] ...
    'Error in $E_L$: $\left( E_{L} - E_{L}^{*} \right) / \left| E_{L}^{*} \right|$' ...
    '$E_{L}$[Structure] - $E_{L}$[Rocksalt]' ...
    'Error in $a$: $\left( a - a^{*} \right) / \left| a^{*} \right|$' ...
    'Error in $c$: $\left( c - c^{*} \right) / \left| c^{*} \right|$' ...
    'Error in $c/a$: $\left( c/a - c^{*}/a^{*} \right) / \left| c^{*}/a^{*} \right|$' ...
    'Error in $\rho$: $\left( \rho - \rho^{*} \right) / \left| \rho^{*} \right|$' ...
    'Error in $V$: $\left( V - V^{*} \right) / \left| V^{*} \right|$'};
if show_as_percent_error
    ylabs = {'Loss' ...
        '[\%]' ...
        '[kJ mol$^{-1}$]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]' ...
        '[\%]'};
else
    ylabs = {'Loss' ...
        '[kJ mol$^{-1}$]' ...
        '[kJ mol$^{-1}$]' ...
        '[\AA]' ...
        '[\AA]' ...
        '$(c/a)$' ...
        'g cm$^{-3}$' ...
        '$\AA^{3}$ / formula unit'};
end

% Create figure and axis
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
t = tiledlayout(figh,N_compares,1,'TileSpacing','tight');

for idx = 1:length(all_y)
    bar_y = all_y{idx};
    plot_switch = plot_switches(idx);
    if ~plot_switch
        continue
    end
    plot_type = plot_types{idx};
    cTitle = Titles{idx};
    ylab = ylabs{idx};
    cYLIM = all_ylim{idx};
    
    axobj = nexttile(t);
    hold(axobj,'on');
    if strcmp(plot_type,'loss')
        p = bar(axobj,bar_y,'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        
        xpos = zeros(N_Salts,N_Models);
        for kdx = 1:N_Salts
            for mdx = 1:N_Models
                xpos(:,mdx) = p(mdx).XEndPoints';
                scatter(axobj,xpos(kdx,mdx),bar_y(kdx,mdx),100,'MarkerEdgeColor','k',...
                    'MarkerFaceColor',Colours(kdx,:),'linewidth',2,'marker',markers{mdx});
            end
        end
        
        for cdx = 1:length(p)
            p(cdx).CData = Colours;
        end
        
        %xticks(axobj,reshape(xpos',1,[]))
        %xticklabels(axobj,repmat(strrep(Models,[Basenum Midnum],''),1,N_Salts))
        
        xticks(axobj,1:N_Salts)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Salts+0.5])
        xticklabels(axobj,[])
        if loss_ylog
            set(axobj, 'YScale', 'log')
        end
        
        %ylim(axobj,[round(min(bar_y,[],'all')*0.1,1,'significant') round(max(bar_y,[],'all')*10, 1,'significant')])
    elseif strcmp(plot_type,'RLE')
        
        % Make model-averaged data for the line plot
        bar_av = reshape(mean(bar_y,2,'omitnan'),[],N_Structures);
        
        % Reshape the data
        bar_y = reshape(permute(bar_y,[2 1 3]),[],N_Structures);
        
        for jdx = 1:N_Structures
            
            c_bar_y = bar_y(:,jdx);
            c_bar_y_av = bar_av(:,jdx);
            c_bar_x = repelem(linspace(jdx-0.3,jdx+0.3,N_Salts),N_Models);
            c_bar_x_av = linspace(jdx-0.3,jdx+0.3,N_Salts);
            
            if strcmpi(Structures{jdx},'CsCl')
                mainax = axobj;
                axobj = axes(figh,'YAxisLocation','right',...
                    'YColor','blue','XColor', 'blue','box','on'); % [left bottom width height]
                hold(axobj,'on')
                cscl_ax = axobj;
                cscl_idx = jdx;
            end
            
            plot(axobj,c_bar_x_av,Rel_Energies_Target(:,jdx),'-','linewidth',2,'color','r')
            for kdx = 1:N_Salts
                scatter(axobj,c_bar_x_av(kdx),Rel_Energies_Target(kdx,jdx),100,'MarkerEdgeColor','r',...
                    'MarkerFaceColor',Colours(kdx,:),'linewidth',2)
            end
            
            plot(axobj,c_bar_x_av,c_bar_y_av,'--k','linewidth',2);
            sca_idx = 1;
            for kdx = 1:N_Salts
                
                xx = c_bar_x(sca_idx:sca_idx+N_Models-1);
                yy = c_bar_y(sca_idx:sca_idx+N_Models-1);
                for mdx = 1:N_Models
                    scatter(axobj,xx(mdx),yy(mdx),100,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',Colours(kdx,:),'linewidth',1,'marker',markers{mdx});
                    
                end
%                 scatter(axobj,c_bar_x(sca_idx:sca_idx+N_Models-1),...
%                     c_bar_y(sca_idx:sca_idx+N_Models-1),100,'MarkerEdgeColor','k',...
%                     'MarkerFaceColor',Colours(kdx,:),'linewidth',2)
                sca_idx = sca_idx + N_Models;
            end
            if strcmpi(Structures{jdx},'CsCl')
                axobj = mainax;
            end
        end
        
        %yline(axobj,0,':k','linewidth',1)
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        plot(axobj,[0.5 N_Structures+0.5],[0 0],':k','linewidth',1)
        xticks(axobj,1:N_Structures)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        xticklabels(axobj,[])
        if set_ylim
            ylim(axobj,cYLIM);
        else
            ylim(axobj,'padded')
        end
        cylim = ylim(axobj);
        
        % Add boxes for targets
        for jdx = 1:N_Structures
            Structure = Structures{jdx};
            if isfield(Bayesopt_Loss_Options,Structure) && Bayesopt_Loss_Options.(Structure).RLE > sqrt(eps)
                
                if strcmpi(Structures{jdx},'CsCl')
                    continue
                end
                
                % [x y w h]. The x and y elements define the coordinate for the lower left corner of the rectangle. 
                % The w and h elements define the dimensions of the rectangle.
                rech = rectangle(axobj,'Position',[jdx-0.4 cylim(1) 0.8 (cylim(2) - cylim(1))],...
                    'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structure).RLE,4),'EdgeColor','r');
            end
        end
        ylim(axobj, cylim);
    else
        
        % Make model-averaged data for the line plot
        bar_av = reshape(mean(bar_y,2,'omitnan'),[],N_Structures);
        
        % Reshape the data
        bar_y = reshape(permute(bar_y,[2 1 3]),[],N_Structures);
        
        for jdx = 1:N_Structures
            c_bar_y = bar_y(:,jdx);
            c_bar_y_av = bar_av(:,jdx);
            c_bar_x = repelem(linspace(jdx-0.3,jdx+0.3,N_Salts),N_Models);
            c_bar_x_av = linspace(jdx-0.3,jdx+0.3,N_Salts);
            
            hold on
            
            plot(axobj,c_bar_x_av,c_bar_y_av,'--k','linewidth',2);
            sca_idx = 1;
            for kdx = 1:N_Salts
                xx = c_bar_x(sca_idx:sca_idx+N_Models-1);
                yy = c_bar_y(sca_idx:sca_idx+N_Models-1);
                
                for mdx = 1:N_Models
                    scatter(axobj,xx(mdx),yy(mdx),100,'MarkerEdgeColor','k',...
                        'MarkerFaceColor',Colours(kdx,:),'linewidth',1,'marker',markers{mdx});
                end
                sca_idx = sca_idx + N_Models;
            end
        end
        
        yline(axobj,0,':k','linewidth',1)
        xticks(axobj,1:N_Structures)
        axobj.XAxis.TickLength = [0,0];
        set(axobj,'box','on','TickLabelInterpreter','latex');
        set(axobj,'XMinorTick','off','YMinorTick','on','FontSize',fs);
        set(axobj,'XLim',[0.5 N_Structures+0.5])
        xticklabels(axobj,[])
        if set_ylim
            ylim(axobj,cYLIM);
        else
            ylim(axobj,'padded');
        end
        cylim = ylim(axobj);
        
        % Add boxes for targets
        for jdx = 1:N_Structures
            Structure = Structures{jdx};
            if isfield(Bayesopt_Loss_Options,Structure) && isfield(Bayesopt_Loss_Options.(Structure),plot_type) && ...
                    Bayesopt_Loss_Options.(Structure).(plot_type) > sqrt(eps)
                % [x y w h]. The x and y elements define the coordinate for the lower left corner of the rectangle. 
                % The w and h elements define the dimensions of the rectangle.
                rech = rectangle(axobj,'Position',[jdx-0.4 cylim(1) 0.8 (cylim(2) - cylim(1))],...
                    'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.(Structure).(plot_type),4),'EdgeColor','r');
            end
        end
        ylim(axobj, cylim);
    end
    
    title(axobj, cTitle,'Fontsize',fs,'Interpreter','latex')
    grid(axobj,'minor');
    axobj.YGrid = 'on';
    axobj.XGrid = 'off';
    ylabel(axobj, ylab,'Fontsize',fs-1,'Interpreter','latex')
    %set(axobj,'children',flipud(get(axobj,'children')))
end

title(t, Title_txt,'Fontsize',fs+5,'Interpreter','latex')
xticklabels(axobj,strrep(strrep(Structures,'FiveFive','5-5'),'BetaBeO','$\beta$-BeO'))

h(1) = plot(t.Children(1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiF',...
    'MarkerFaceColor',Colours(1,:),'linewidth',2,'Color','k');
h(2) = plot(t.Children(1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiCl',...
    'MarkerFaceColor',Colours(2,:),'linewidth',2,'Color','k');
h(3) = plot(t.Children(1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiBr',...
    'MarkerFaceColor',Colours(3,:),'linewidth',2,'Color','k');
h(4) = plot(t.Children(1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'LiI',...
    'MarkerFaceColor',Colours(4,:),'linewidth',2,'Color','k');
h(5) = plot(t.Children(1),nan, nan, 'o', 'MarkerSize', 12, 'DisplayName', 'DFT',...
    'MarkerFaceColor','w','linewidth',2,'Color','r');


legh = legend(h,'Location','SouthOutside','Orientation','Horizontal',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', N_Salts+1);
%legh.Layout.Tile = 'East';
%titleh = title(legh,'\bf{\underline{Salts}}','Fontsize',fs,'Interpreter','latex')

%Set CsCl axis properties last
if contained_in_cell('CsCl',Structures)
    drawnow;
    ll = mainax.Position(1) + (cscl_idx-1)*(mainax.Position(3)/(N_Structures));
    bb = mainax.Position(2);
    ww = mainax.Position(3)/(N_Structures);
    hh = mainax.Position(4);
    
    ylim(cscl_ax,'padded')
    xticks(cscl_ax,[])
    xlim(cscl_ax,[cscl_idx-0.5 cscl_idx+0.5])
    set(cscl_ax,'TickLabelInterpreter','latex','YMinorTick','on','FontSize',fs-2,...
        'Position',[ll bb ww hh],'LineWidth',1);
    grid(cscl_ax,'minor');
    cscl_ax.YGrid = 'on';
    cscl_ax.XGrid = 'off';
    xticks(cscl_ax,1:N_Structures)
    xticklabels(cscl_ax,[])
    %cscl_ax.YRuler.TickLabelGapOffset = -5;
    
    cylim_cscl = ylim(cscl_ax);
    if isfield(Bayesopt_Loss_Options,'CsCl') && Bayesopt_Loss_Options.CsCl.RLE > sqrt(eps)
        % [x y w h]. The x and y elements define the coordinate for the lower left corner of the rectangle. 
        % The w and h elements define the dimensions of the rectangle.
        rech = rectangle(cscl_ax,'Position',[cscl_idx-0.4 cylim_cscl(1) 0.8 (cylim_cscl(2) - cylim_cscl(1))],...
            'Curvature',0,'LineWidth',min(3*Bayesopt_Loss_Options.CsCl.RLE,4),'EdgeColor','r');
    end
    ylim(cscl_ax,cylim_cscl);
end

if savefile
    set(figh,'renderer','opengl')
    exportgraphics(figh,fullfile(pwd,'figures',filename),'Resolution',300)
end
%clear;
global sldr_h results X Y idx_X idx_Y input_table p full_opt_point ...
    full_opt_eval bayes_opt_eval bayes_opt_point bayes_opt_est_eval bayes_opt_est_point show_error transforms_extra
% Analysis parameters
Salt = 'LiI';
Theory = 'BE';
Model = 'OA1'; % C5squaredexponential
Plot_X = 'epsilon_MM'; % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ' 
Plot_Y = 'epsilon_XX';   % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'

% 'r0_MM'    'r0_XX'    'epsilon_MM'    'epsilon_XX'    'gamma_MX'

grid_density = 50; % Parameter to set the grid density
fs = 24;
show_error = false;

Model_type = 'error'; % constraints, error, objective

% Load data
A = load(fullfile('C:\Users\Hayden\Documents\Patey_Lab\BO_Models',Salt,[Salt '_' Theory '_Model_' Model '_data.mat'])).full_data;
results = A.bayesopt_results;

bayes_opt_eval = results.MinObjective;
bayes_opt_point = table2array(results.XAtMinObjective);

bayes_opt_est_eval = nan;
bayes_opt_est_point = table2array(results.XAtMinEstimatedObjective);

% Load best point
full_opt_point = A.full_opt_point;
full_opt_eval = A.loss;

% Variable names
varnames = {results.VariableDescriptions.Name};
transforms = {results.VariableDescriptions.Transform};
N_vars = length(varnames);

% Get indexes of the options of interest
[~,idx_X] = contained_in_cell(Plot_X,varnames);
[~,idx_Y] = contained_in_cell(Plot_Y,varnames);

% Grab projection of best point, bayes best point, and bayes estimated best point
X_best = full_opt_point(idx_X);
Y_best = full_opt_point(idx_Y);
X_bayes = bayes_opt_point(idx_X);
Y_bayes = bayes_opt_point(idx_Y);
X_bayes_est = bayes_opt_est_point(idx_X);
Y_bayes_est = bayes_opt_est_point(idx_Y);

% Get limits of X and Y coordinate
lims_X = results.VariableDescriptions(idx_X).Range;
Transform_X = transforms{idx_X};
lims_Y = results.VariableDescriptions(idx_Y).Range;
Transform_Y = transforms{idx_Y};

% Generate grids of X and Y coordinate
if strcmp(Transform_X,'log')
    Grid_X = linspace(log(lims_X(1)),log(lims_X(2)),grid_density+1);
    Grid_X = exp(Grid_X);
else
    Grid_X = linspace(lims_X(1),lims_X(2),grid_density+1);
end
if strcmp(Transform_Y,'log')
    Grid_Y = linspace(log(lims_Y(1)),log(lims_Y(2)),grid_density+1);
    Grid_Y = exp(Grid_Y);
else
    Grid_Y = linspace(lims_Y(1),lims_Y(2),grid_density+1);
end

[X,Y] = meshgrid(Grid_X,Grid_Y);

% Flatten the grids for input
flat_X = reshape(X,[],1);
flat_Y = reshape(Y,[],1);

% Create input table
flat_input = ones(length(flat_X),N_vars);
flat_input(:,idx_X) = flat_X;
flat_input(:,idx_Y) = flat_Y;

% Set to slice of best point
idxs = 1:N_vars;
idxs([idx_X idx_Y]) = [];
flat_input(:,idxs) = repmat(full_opt_point(idxs),length(flat_X),1);

input_table = array2table(flat_input,'VariableNames',varnames);

% Sample the model at the grid points
switch lower(Model_type)
    case 'objective'
        [flat_Z,flat_C] = results.predictObjective(input_table);
    case 'error'
        [flat_Z,flat_C] = results.predictError(input_table);
    case 'constraints'
        [flat_Z,flat_C] = results.predictConstraints(input_table);
end

% Reshape the model outputs to match the grids
Z = reshape(flat_Z,size(X));

if show_error
    C = reshape(flat_C,size(X));
else
    C = Z;
end

minZ = min(Z,[],'all');
maxZ = max(Z,[],'all');

if maxZ - minZ > 100
    Z = real(log1p(Z));
    full_opt_eval = log1p(full_opt_eval);
    bayes_opt_eval = max(real(log1p(bayes_opt_eval)),0);
    bayes_opt_est_eval = max(real(log1p(bayes_opt_est_eval)),0);
end

% Create Figure
figh = figure('WindowState','Maximized');
colormap(figh,'turbo')
ax = axes(figh,'Position',[0.200 0.1100 0.700 0.8150]); % [left bottom width height]

% Add sliders
varnames_extra = varnames;
transforms_extra = transforms;
varnames_extra([idx_X idx_Y]) = [];
transforms_extra([idx_X idx_Y]) = [];
N_vars_extra = N_vars-2;

sldr_h = gobjects(N_vars_extra*2,1);
for idx = 1:N_vars_extra
    
    [~,idx_Q] = contained_in_cell(varnames_extra{idx},varnames);
    Transform_Q = transforms_extra{idx};
    lims_Q = results.VariableDescriptions(idx_Q).Range;
    x0 = full_opt_point(idxs(idx));
    
    if strcmp(Transform_Q,'log')
        lims_Q = log(lims_Q);
        x0 = log(x0);
        x0_str = num2str(exp(x0),'%.4f');
    else
        x0_str = num2str(x0,'%.4f');
    end
    
    
    textpos  = [0 0.7-0.1*idx    0.15 0.05]; % [left bottom width height]
    slidepos = [0 0.65-0.1*idx 0.15 0.05]; % [left bottom width height]
    str = [varnames_extra{idx} ' = ' x0_str];
    
    sldr_h(idx) = uicontrol(figh,'Style','Text','Units','Normalized','Position',textpos,...
        'String',str,'FontSize',fs-4,'Tag',varnames_extra{idx});
    sldr_h(idx+N_vars_extra) = uicontrol(figh,'Units','Normalized','Style','Slider','Position',slidepos,...
        'String',str,'Value',x0,'Max',lims_Q(2),'Min',lims_Q(1),'FontSize',fs,...
        'Callback',@UpdatePlot,'UserData',[idx idx+N_vars_extra],...
        'SliderStep',[0.0001 0.1]);
end

% Add Buttons
button_h = gobjects(3,1);
button_h(1) = uicontrol(figh,'Units','Normalized','Style','PushButton','Position',[0.02 0.90 0.13 0.05],...
    'String','Full Opt. Point','Value',0,'FontSize',fs-5,...
    'Callback',@UpdatePlot,'UserData',1);
button_h(2) = uicontrol(figh,'Units','Normalized','Style','PushButton','Position',[0.02 0.80 0.13 0.05],...
    'String','Bayes Opt. Point','Value',0,'FontSize',fs-5,...
    'Callback',@UpdatePlot,'UserData',2);
button_h(3) = uicontrol(figh,'Units','Normalized','Style','PushButton','Position',[0.02 0.70 0.13 0.05],...
    'String','Est. Bayes Opt. Point','Value',0,'FontSize',fs-5,...
    'Callback',@UpdatePlot,'UserData',3);

% Plot
hold(ax,'on')
p = surf(ax,X,Y,Z,C,'EdgeColor','none','FaceAlpha',0.50);

% Plot the best point.
h(1) = scatter3(ax,X_best,Y_best,full_opt_eval,100,full_opt_eval,'filled',...
    'MarkerEdgeColor','k','MarkerFaceColor','g');
h(2) = scatter3(ax,X_bayes,Y_bayes,bayes_opt_eval,100,bayes_opt_eval,'filled',...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
h(3) = scatter3(ax,X_bayes_est,Y_bayes_est,bayes_opt_est_eval,100,bayes_opt_est_eval,'filled',...
    'MarkerEdgeColor','k','MarkerFaceColor','m');

xlabel(ax,param_name_map(Plot_X,Salt),'Interpreter','latex','FontSize',fs);
ylabel(ax,param_name_map(Plot_Y,Salt),'Interpreter','latex','FontSize',fs);
zlabel(ax,'Log($f(\bf{x}) + 1$)','Interpreter','latex','FontSize',fs);

set(ax,'XGrid','On','YGrid','On','GridLineStyle','-','Layer','Top',...
    'TickLength',[0 0],'FontSize',fs,'TickLabelInterpreter','latex')

title(ax,[Salt ' ' Theory ' Model ' Model ': Bayesian Surrogate Model'],...
    'Interpreter','Latex','FontSize',fs+16);

legend(ax,h,{'Full Opt.','Bayes. Opt.','Est. Bayes. Opt.'},...
    'Interpreter','Latex','FontSize',fs)

xlim(ax,lims_X)
ylim(ax,lims_Y)

if strcmp(Transform_X,'log')
    set(ax, 'XScale', 'log')
end
if strcmp(Transform_Y,'log')
    set(ax, 'YScale', 'log')
end


function output = param_name_map(p_name,Salt)
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    switch lower(p_name)
        case 'sdmm'
            output = ['$S_D$[' Metal '-' Metal ']'];
        case 'sdxx'
            output = ['$S_D$[' Halide '-' Halide ']'];
        case 'srmm'
            output = ['$S_R$[' Metal '-' Metal ']'];
        case 'srxx'
            output = ['$S_R$[' Halide '-' Halide ']'];
        case 'snmm'
            output = ['$S_N$[' Metal '-' Metal ']'];
        case 'snxx'
            output = ['$S_N$[' Halide '-' Halide ']'];
        case 'snmx'
            output = ['$S_N$[' Metal '-' Halide ']'];
        case 'sq'
            output = '$S_Q$';
        case 'sigma_mm'
            output = ['$\sigma$[' Metal '-' Metal ']'];
        case 'sigma_xx'
            output = ['$\sigma$[' Halide '-' Halide ']'];
        case 'sigma_mx'
            output = ['$\sigma$[' Metal '-' Halide ']'];
        case 'epsilon_mm'
            output = ['$\varepsilon$[' Metal '-' Metal ']'];
        case 'epsilon_xx'
            output = ['$\varepsilon$[' Halide '-' Halide ']'];
        case 'epsilon_mx'
            output = ['$\varepsilon$[' Metal '-' Halide ']'];
        case 'r0_xx'
            output = ['$r_{0}$[' Halide '-' Halide ']'];
        case 'r0_mm'
            output = ['$r_{0}$[' Metal '-' Metal ']'];
        case 'r0_mx'
            output = ['$r_{0}$[' Metal '-' Halide ']'];
        case 'gamma_mx'
            output = ['\gamma$[' Metal '-' Halide ']'];
        case 'gamma_mm'
            output = ['\gamma$[' Metal '-' Metal ']'];
        case 'gamma_xx'
            output = ['\gamma$[' Halide '-' Halide ']'];
        otherwise
            output = p_name;
    end
end

function UpdatePlot(src,~)
    global sldr_h results X input_table p full_opt_point idx_X idx_Y...
        bayes_opt_point bayes_opt_est_point show_error transforms_extra
    
    if strcmp(src.Style,'pushbutton')
        
        N_vars =  length(input_table.Properties.VariableNames);
        fidx = 1:N_vars;
        fidx([idx_X idx_Y]) = [];
        N_vars = N_vars-2;
        
        switch src.UserData
            case 1 % Full Opt. Point
                
                % Update sliders
                for idx = 1:N_vars
                    jdx = fidx(idx);
                    varname = input_table.Properties.VariableNames{jdx};
                    Transform_Q = results.VariableDescriptions(jdx).Transform;
                    Val = full_opt_point(jdx);
                    str = [varname ' = ' num2str(Val,'%.4f')];
                    
                    if strcmp(Transform_Q,'log')
                        Val = log(Val);
                    end
                    sldr_h(N_vars+idx).Value = Val;
                    set(sldr_h(idx),'String',str)
                end
                
                % Prepare input table
                input_table{:,fidx} = repmat(full_opt_point(fidx),height(input_table),1);
                
            case 2 % Bayes Opt. Point
                
                % Update sliders
                for idx = 1:N_vars
                    jdx = fidx(idx);
                    varname = input_table.Properties.VariableNames{jdx};
                    Transform_Q = results.VariableDescriptions(jdx).Transform;
                    Val = bayes_opt_point(jdx);
                    str = [varname ' = ' num2str(Val,'%.4f')];
                    
                    if strcmp(Transform_Q,'log')
                        Val = log(Val);
                    end
                    
                    sldr_h(N_vars+idx).Value = Val;
                    set(sldr_h(idx),'String',str)
                end
                
                % Prepare input table
                input_table{:,fidx} = repmat(bayes_opt_point(fidx),height(input_table),1);
                
            case 3 % Est. Bayes Opt. Point
                
                % Update sliders
                for idx = 1:N_vars
                    jdx = fidx(idx);
                    varname = input_table.Properties.VariableNames{jdx};
                    Transform_Q = results.VariableDescriptions(jdx).Transform;
                    Val = bayes_opt_est_point(jdx);
                    str = [varname ' = ' num2str(Val,'%.4f')];
                    
                    if strcmp(Transform_Q,'log')
                        Val = log(Val);
                    end
                    
                    sldr_h(N_vars+idx).Value = Val;
                    set(sldr_h(idx),'String',str)
                end
                
                % Prepare input table
                input_table{:,fidx} = repmat(bayes_opt_est_point(fidx),height(input_table),1);
        end
        
    else
        t_idx = src.UserData(1); % text index
        s_idx = src.UserData(2); % slider index
        Scaling = sldr_h(t_idx).Tag; % Scaling selected
        Val = sldr_h(s_idx).Value; % Current value
        
        Transform_Q = transforms_extra{t_idx};
        if strcmp(Transform_Q,'log')
            Val = exp(Val);
        end
        
        % Update the string
        str = [Scaling ' = ' num2str(Val,'%.4f')];
        set(sldr_h(t_idx),'String',str)
        
        % Prepare input table
        input_table.(Scaling) = repmat(Val,height(input_table),1);
    end
    
    % Sample the model at the grid points
    [flat_Z,flat_C] = results.predictObjective(input_table);
    
    % Reshape the model outputs to match the grids
    Z = reshape(flat_Z,size(X));
    C = reshape(flat_C,size(X));

    minZ = min(Z,[],'all');
    maxZ = max(Z,[],'all');

    if maxZ - minZ > 100
        Z = max(real(log1p(Z)),0);
    end
    
    % Reload the plot
    p.ZData = Z;
    if show_error
        p.CData = C;
    else
        p.CData = Z;
    end
end

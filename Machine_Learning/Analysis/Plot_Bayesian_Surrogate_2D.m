clear;
global sldr_h results X Y idx_X idx_Y input_table p full_opt_point ...
    full_opt_eval bayes_opt_eval bayes_opt_point bayes_opt_est_eval bayes_opt_est_point show_error
% Analysis parameters
Salt = 'LiI';
Theory = 'JC';
Model = 'C5'; % C5squaredexponential
Plot_X = 'SDXX'; % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
Plot_Y = 'SRXX';   % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
grid_density = 50; % Parameter to set the grid density
fs = 24;
show_error = false;


% Load data
A = load(['C:\Users\Hayden\Documents\Patey_Lab\BO_Models\' Salt '_' Theory '_Model_' Model '_bayesopt.mat']);
results = A.results;
bayes_opt_eval = results.MinObjective;
bayes_opt_point = table2array(results.XAtMinObjective);

bayes_opt_est_eval = results.MinEstimatedObjective;
bayes_opt_est_point = table2array(results.XAtMinEstimatedObjective);

% Load best point
A = load(['C:\Users\Hayden\Documents\Patey_Lab\BO_Models\' Salt '_' Theory '_Model_' Model '_fullopt.mat']);
full_opt_point = A.full_opt_point;
full_opt_eval = A.loss;

% Variable names
varnames = {results.VariableDescriptions.Name};
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
lims_Y = results.VariableDescriptions(idx_Y).Range;

% Generate grids of X and Y coordinate
gridsize_X = diff(lims_X)/grid_density;
gridsize_Y = diff(lims_Y)/grid_density;
Grid_X = lims_X(1):gridsize_X:lims_X(2);
Grid_Y = lims_Y(1):gridsize_Y:lims_Y(2);

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
[flat_Z,flat_C] = results.predictObjective(input_table);

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
varnames_extra([idx_X idx_Y]) = [];
N_vars_extra = N_vars-2;

sldr_h = gobjects(N_vars_extra*2,1);
for idx = 1:N_vars_extra
    
    x0 = full_opt_point(idxs(idx));
    x0_str = num2str(x0,'%.4f');
    
    textpos  = [0 0.7-0.1*idx    0.15 0.05]; % [left bottom width height]
    slidepos = [0 0.65-0.1*idx 0.15 0.05]; % [left bottom width height]
    str = [varnames_extra{idx} ' = ' x0_str];
    
    [~,idx_Q] = contained_in_cell(varnames_extra{idx},varnames);
    lims_Q = results.VariableDescriptions(idx_Q).Range;
    
    sldr_h(idx) = uicontrol(figh,'Style','Text','Units','Normalized','Position',textpos,...
        'String',str,'FontSize',fs,'Tag',varnames_extra{idx});
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

function output = param_name_map(p_name,Salt)
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    switch p_name
        case 'SDMM'
            output = ['$S_D$[' Metal '-' Metal ']'];
        case 'SDXX'
            output = ['$S_D$[' Halide '-' Halide ']'];
        case 'SRMM'
            output = ['$S_R$[' Metal '-' Metal ']'];
        case 'SRXX'
            output = ['$S_R$[' Halide '-' Halide ']'];
        case 'SNMM'
            output = ['$S_N$[' Metal '-' Metal ']'];
        case 'SNXX'
            output = ['$S_N$[' Halide '-' Halide ']'];
        case 'SNMX'
            output = ['$S_N$[' Metal '-' Halide ']'];
        case 'SQ'
            output = '$S_Q$';
        otherwise
            output = p_name;
    end
end

function UpdatePlot(src,~)
    global sldr_h results X input_table p full_opt_point idx_X idx_Y...
        bayes_opt_point bayes_opt_est_point show_error
    
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
                    Val = full_opt_point(jdx);
                    
                    sldr_h(N_vars+idx).Value = Val;
                    
                    str = [varname ' = ' num2str(Val,'%.4f')];
                    set(sldr_h(idx),'String',str)
                end
                
                % Prepare input table
                input_table{:,fidx} = repmat(full_opt_point(fidx),height(input_table),1);
                
            case 2 % Bayes Opt. Point
                
                % Update sliders
                for idx = 1:N_vars
                    jdx = fidx(idx);
                    varname = input_table.Properties.VariableNames{jdx};
                    Val = bayes_opt_point(jdx);
                    
                    sldr_h(N_vars+idx).Value = Val;
                    
                    str = [varname ' = ' num2str(Val,'%.4f')];
                    set(sldr_h(idx),'String',str)
                end
                
                % Prepare input table
                input_table{:,fidx} = repmat(bayes_opt_point(fidx),height(input_table),1);
                
            case 3 % Est. Bayes Opt. Point
                
                % Update sliders
                for idx = 1:N_vars
                    jdx = fidx(idx);
                    varname = input_table.Properties.VariableNames{jdx};
                    Val = bayes_opt_est_point(jdx);
                    
                    sldr_h(N_vars+idx).Value = Val;
                    
                    str = [varname ' = ' num2str(Val,'%.4f')];
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

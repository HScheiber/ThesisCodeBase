% Analysis parameters
Salt = 'LiBr';
Theory = 'JC';
Model = 'CAa';
Plot_X = 'SRMM'; % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
Plot_Y = 'SRXX';   % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
Plot_Z = 'SDMM';   % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
grid_density = 50; % Parameter to set the grid density
fs = 24;

% Load data
A = load(['C:\Users\Hayden\Documents\Patey_Lab\BO_Models\' Salt '_' Theory '_Model_' Model '_bayesopt.mat']);
results = A.results;

% Variable names
varnames = {results.VariableDescriptions.Name};
N_vars = length(varnames);

% Get indexes of the options of interest
[~,idx_X] = contained_in_cell(Plot_X,varnames);
[~,idx_Y] = contained_in_cell(Plot_Y,varnames);
[~,idx_Z] = contained_in_cell(Plot_Z,varnames);

% Get limits of coordinates
lims_X = results.VariableDescriptions(idx_X).Range;
lims_Y = results.VariableDescriptions(idx_Y).Range;
lims_Z = results.VariableDescriptions(idx_Z).Range;

% Generate grids of X and Y coordinate
gridsize_X = diff(lims_X)/grid_density;
gridsize_Y = diff(lims_Y)/grid_density;
gridsize_Z = diff(lims_Z)/grid_density;
Grid_X = lims_X(1):gridsize_X:lims_X(2);
Grid_Y = lims_Y(1):gridsize_Y:lims_Y(2);
Grid_Z = lims_Z(1):gridsize_Z:lims_Z(2);

[X,Y,Z] = meshgrid(Grid_X,Grid_Y,Grid_Z);

% Flatten the grids for input
flat_X = reshape(X,[],1);
flat_Y = reshape(Y,[],1);
flat_Z = reshape(Z,[],1);

% Create input table
flat_input = ones(length(flat_X),N_vars);
flat_input(:,idx_X) = flat_X;
flat_input(:,idx_Y) = flat_Y;
flat_input(:,idx_Z) = flat_Z;
input_table = array2table(flat_input,'VariableNames',varnames);

% Sample the model at the grid points
flat_C = results.predictObjective(input_table);

% Reshape the model outputs to match the grids
C = reshape(flat_C,size(X));

minC = min(C,[],'all');
maxC = max(C,[],'all');

if maxC - minC > 100
    C = real(log1p(C));
end

% Plot
figh = figure('WindowState','Maximized');
%colormap(figh,'hot')
ax = axes(figh);

max_iso = max(C,[],'all');
min_iso = 0;
isovalues = 0.1:0.1:0.9;
N_iso = length(isovalues);
alphas = linspace(1,0.1,N_iso);
%colors = cbrewer('seq','YlOrRd',N_iso,'spline');
colors = hot(N_iso);

for idx = 1:N_iso
    iso = isovalues(idx);
    alpha = alphas(idx);
    col = colors(idx,:);

    isovalue = iso*(max_iso - min_iso) + min_iso;
    surf = isosurface(X,Y,Z,C,isovalue);

    p = patch(ax,surf);
    set(p,'FaceColor',col,'EdgeColor','none','FaceAlpha',alpha)
end


xlabel(ax,param_name_map(Plot_X,Salt),'Interpreter','latex','FontSize',fs);
ylabel(ax,param_name_map(Plot_Y,Salt),'Interpreter','latex','FontSize',fs);
zlabel(ax,param_name_map(Plot_Z,Salt),'Interpreter','latex','FontSize',fs);

set(ax,'XGrid','On','YGrid','On','ZGrid','On','GridLineStyle','-','Layer','Top',...
    'FontSize',fs,'TickLabelInterpreter','latex')

title(ax,[Salt ' ' Theory ' Model ' Model ': Bayesian Surrogate Model'],...
    'Interpreter','Latex','FontSize',fs+16);


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
        case 'SQ'
            output = '$S_Q$';
    end
end

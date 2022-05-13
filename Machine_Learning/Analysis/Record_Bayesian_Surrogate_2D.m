clear;

% Analysis parameters
Salt = 'LiCl';
Theory = 'JC';
Model = 'D1'; % C5squaredexponential
Plot_X = 'SRMM'; % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
Plot_Y = 'SRXX';   % options: 'SDMM' 'SDXX' 'SRMM' 'SRXX' 'SQ'
video_filename = [Salt '_' Theory '_Model_' Model '.mp4'];
grid_density = 200; % Parameter to set the grid density
step_density = 200; % Sets the resolution in time
fps = 24; % Frames per second
fs = 24; % font size
show_error = true;
axpos = [0.10 0.100 0.800 0.850];  % [left bottom width height]
initial_rot_angle = 75;
height_angle = 25;
Z_limit = [-0.2 2];
savevideo = true;
bg_color = hex2rgb('#010940');

% Load data
A = load(['C:\Users\Hayden\Documents\Patey_Lab\BO_Models\' Salt '_' Theory '_Model_' Model '_bayesopt.mat']);
results = A.results;
bayes_opt_eval = results.MinObjective;
bayes_opt_point = table2array(results.XAtMinObjective);

bayes_opt_next_point = table2array(results.NextPoint);

% Variable names
varnames = {results.VariableDescriptions.Name};
N_vars = length(varnames);

% Get indexes of the options of interest
[~,idx_X] = contained_in_cell(Plot_X,varnames);
[~,idx_Y] = contained_in_cell(Plot_Y,varnames);

% Grab projection of best point, and next point
X_bayes = bayes_opt_point(idx_X);
Y_bayes = bayes_opt_point(idx_Y);

X_next = bayes_opt_next_point(idx_X);
Y_next = bayes_opt_next_point(idx_Y);
[Z_next,C_next] = results.predictObjective(array2table(bayes_opt_next_point,'VariableNames',varnames));

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
N_flat_X = length(flat_X);

% Create input table
flat_input = ones(N_flat_X,N_vars);
flat_input(:,idx_X) = flat_X;
flat_input(:,idx_Y) = flat_Y;

% Set to slice of best point
idxs = 1:N_vars;
idxs([idx_X idx_Y]) = [];
flat_input(:,idxs) = repmat(bayes_opt_point(idxs),N_flat_X,1);

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
    bayes_opt_eval = max(real(log1p(bayes_opt_eval)),0);
    transform_fx = true;
else
    transform_fx = false;
end

% Create Figure
figh = figure('WindowState','Maximized','Color',bg_color,'ToolBar','none','MenuBar','none');
colormap(figh,'turbo')
ax = axes(figh,'Position',axpos,'ZGrid','On','YGrid','On','XGrid','On',...
    'Color',[0.4667 0.7725 0.9765],'XColor',[1 1 1], 'YColor',[1 1 1], 'ZColor',[1 1 1],...
    'GridColor',[0 0 0],'ZLim',Z_limit);
hold(ax,'on')
%pbaspect(ax,[1 1 1])

% Plot the X-Y plane
% Generate grids of X and Y coordinate
gridsize_X = diff(lims_X)/10;
gridsize_Y = diff(lims_Y)/10;
Grid_X = lims_X(1):gridsize_X:lims_X(2);
Grid_Y = lims_Y(1):gridsize_Y:lims_Y(2);

[Xa,Ya] = meshgrid(Grid_X,Grid_Y);

Z_zero = zeros(size(Xa, 1)); % Generate z data
surf(ax,Xa,Ya,Z_zero,'FaceColor','none','EdgeAlpha',0.75) % Plot the surface

% Add colorbar
hCB = colorbar(ax,'Position',[0.94 0.15 0.025 0.75],... % [left bottom width height]
    'TickLabelInterpreter','latex','FontSize',fs,'Color',[1 1 1]);
set(hCB.Title,{'String','Interpreter','FontSize','Color'},{'$\sigma(\mathbf{x})$','latex',fs,[1 1 1]})
caxis(ax,[0 max(C,[],'all')*1])

% variable names of non-plotted variables
varnames_extra = varnames;
varnames_extra([idx_X idx_Y]) = [];
N_vars_extra = N_vars-2;

% Plot
hold(ax,'on')
p = surf(ax,X,Y,Z,C,'EdgeColor','none','FaceAlpha',1);

% Plot the best point.
h_scatter = scatter3(ax,X_bayes,Y_bayes,bayes_opt_eval,100,bayes_opt_eval,'filled',...
    'MarkerEdgeColor','k','MarkerFaceColor','r');

xlabel(ax,param_name_map(Plot_X,Salt),'Interpreter','latex','FontSize',fs);
ylabel(ax,param_name_map(Plot_Y,Salt),'Interpreter','latex','FontSize',fs);
zlabel(ax,'Log($f(\bf{x}) + 1$)','Interpreter','latex','FontSize',fs);

set(ax,'XGrid','On','YGrid','On','GridLineStyle','-','Layer','Top',...
    'TickLength',[0 0],'FontSize',fs,'TickLabelInterpreter','latex');%,...
    %'Zlim',[-Inf Inf]);

arrow_start = [mean(lims_X) mean(lims_Y) maxZ*0.9];
min_point = [X_bayes Y_bayes bayes_opt_eval];
diff_arrow = min_point - arrow_start;

arrow_end = min_point - 0.05.*diff_arrow;
text_point = arrow_start - 0.05.*diff_arrow;

texth = text(text_point(1),text_point(2),text_point(3),'Best Known Point',...
    'Interpreter','latex','FontSize',fs,'HorizontalAlignment', 'center');

titletxt = '';
extra_point = bayes_opt_point(idxs);
for kdx = 1:N_vars_extra
    titletxt = [titletxt param_name_map(varnames_extra{kdx},Salt) ': ' num2str(extra_point(kdx),'%.4f') '\quad '];
end

title_h = title(ax,titletxt,'Interpreter','Latex','FontSize',fs,'Color',[1 1 1]);

xlim(ax,lims_X)
ylim(ax,lims_Y)


%% Full rotation at best point
% Calculate movement of frames
az = linspace(initial_rot_angle,initial_rot_angle+360,step_density);
el = [linspace(height_angle,height_angle+20,step_density/2) linspace(height_angle+20,height_angle,step_density/2)];
F1 = struct('cdata', cell(1,step_density), 'colormap', cell(1,step_density));

drawnow;

% Start off with a full rotation around the best point
for kdx = 1:step_density
    view(ax,az(kdx),el(kdx))
    harrow = arrow(arrow_start,arrow_end,30,50);
    F1(kdx) = getframe(figh);
    delete(harrow);
end

% Remove the best point plot and arrows
delete(h_scatter);
delete(texth);

%% Move through search space

% Now move towards the next point in a straight line while simultaneously rotating
diff_vector = bayes_opt_next_point(idxs) - bayes_opt_point(idxs);
tp = linspace(0,1,step_density);
maps = bayes_opt_point(idxs) + kron(tp',diff_vector);

%az = linspace(0,360,step_density);
%el = [linspace(10,20,step_density/2) linspace(20,10,(step_density)/2)];
F2 = struct('cdata', cell(1,step_density), 'colormap', cell(1,step_density));
for kdx = 1:step_density
    
    new_point = maps(kdx,:);
    flat_input(:,idxs) = repmat(new_point,N_flat_X,1);

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

    if transform_fx
        Z = real(log1p(Z));
    end
    
    % Update title
    titletxt = '';
    for mdx = 1:N_vars_extra
        titletxt = [titletxt param_name_map(varnames_extra{mdx},Salt) ': ' num2str(new_point(mdx),'%.4f') '\quad '];
    end
    title_h.String = titletxt;
    
    % update plot
    p.ZData = Z;
    p.CData = C;
    %ax.ZLim = [-Inf Inf];
    
    % update view and arrow
    %view(ax,az(kdx),45)
    
    % Record frame
    F2(kdx) = getframe(figh);
end

%% Rotate around next point

% Plot the next point.
h_scatter = scatter3(ax,X_next,Y_next,Z_next,100,Z_next,'filled',...
    'MarkerEdgeColor','k','MarkerFaceColor','r');

% Label next point
arrow_start = [mean(lims_X) mean(lims_Y) max(Z,[],'all')];
min_point = [X_next Y_next Z_next];
diff_arrow = min_point - arrow_start;
arrow_end = min_point - 0.05.*diff_arrow;
text_point = arrow_start - 0.05.*diff_arrow;
texth = text(text_point(1),text_point(2),text_point(3),'Next Search Point',...
    'Interpreter','latex','FontSize',fs,'HorizontalAlignment', 'center');

% Calculate movement of frames
az = linspace(initial_rot_angle,initial_rot_angle+360,step_density);
el = [linspace(height_angle,height_angle+20,step_density/2) linspace(height_angle+20,height_angle,step_density/2)];
F3 = struct('cdata', cell(1,step_density), 'colormap', cell(1,step_density));

% Full rotation around the next point
for kdx = 1:step_density
    view(ax,az(kdx),el(kdx))
    harrow = arrow(arrow_start,arrow_end,30,50);
    F3(kdx) = getframe(figh);
    if kdx ~= step_density
        delete(harrow);
    end
end

% Marge F1, F2, F3 together
F = [F1 F2 F3];
N_F = length(F);

%% Write video to file
if savevideo
    % Rotate the plot and record images
    writerObj = VideoWriter(video_filename,'MPEG-4');
    writerObj.FrameRate = fps;   % set the seconds per image

    open(writerObj);
    f = waitbar(0,'Rendering...');
    for i = 1:N_F
        % convert the image to a frame
        try
            writeVideo(writerObj, F(i));
        catch
        end
        waitbar(i/N_F,f,'Rendering...');
    end
    close(writerObj);
    close(f);
end


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

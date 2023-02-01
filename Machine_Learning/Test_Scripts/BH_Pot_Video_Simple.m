PlotType = 'LJ';
fps = 60; % Frames per second
fs = 32; % font size
lw = 2; % line width
Trail_length = 40;
Density = 50;
savevideo = true;
video_filename = 'BH_Pot.mp4';

epsilon = 1;
r0 = 1;

gammas =  [linspace(7,15,50*Density)];

rmin = 0;
rmax = 4;
ylims = [-3.5 3.5];

% Initialize frames
N_frame = numel(gammas);
F = struct('cdata', cell(1,N_frame), 'colormap', cell(1,N_frame));

% Generate figure
%bg_color = hex2rgb('#010940');
bg_color = [1 1 1];
figh = figure('WindowState','Maximized','Color',bg_color,'ToolBar','none','MenuBar','none');
axpos = [0.10 0.120 0.880 0.850];  % [left bottom width height]
ax = axes(figh,'Position',axpos,'YGrid','On','XGrid','On',...
    'Color',[1 1 1],'XColor',[0 0 0], 'YColor',[0 0 0],... % [0.4667 0.7725 0.9765]
    'GridColor',[0 0 0],'TickLabelInterpreter','latex','FontSize',fs,...
    'TickLength',[0 0],'Ylim',ylims,'Xlim',[rmin rmax],'TitleHorizontalAlignment','left',...
    'box','on');
hold(ax,'on')
r = rmin:0.001:rmax;
y_set = nan(Trail_length,numel(r));
plot_objs = gobjects(Trail_length,1);
Colours = flipud(cbrewer('seq','Blues',Trail_length,'spline'));
title_obj = title(ax,['$\gamma_{ij} = 0.0000$' newline '$\varepsilon_{ij} = -1$' newline '$\sigma_{ij} = 1$'],...
    'FontSize',fs,'Interpreter','latex','Color',[0 0 0],...
    'Position',[2 1.5 0]);
ylabel(ax,'$u_{ij}(r)$','FontSize',fs,'Interpreter','latex','Color',[0 0 0])
xlabel(ax,'$r$','FontSize',fs,'Interpreter','latex','Color',[0 0 0])

% Initialize plots
for idx = Trail_length:-1:1
    plot_objs(idx) = plot(ax,r,y_set(idx,:),...
        '-','LineWidth',lw,'Color',Colours(idx,:));
end
yline(ax,0,':k',"LineWidth",1)
lin_objy = plot(ax,[rmin r0],[epsilon epsilon],'--k',"LineWidth",2);
lin_objx = plot(ax,[r0 r0],[ylims(1) epsilon],'--k',"LineWidth",2);
scat_obj = scatter(ax,r0,epsilon,200,'r','LineWidth',2);
txt_obj = text(r0,epsilon+diff(ylims)*0.05,'$(\sigma_{ij},-\varepsilon_{ij})$','fontsize',fs,'interpreter','latex');

for idx = 1:numel(gammas)
    gamma = gammas(idx);
    
    if gamma < 6 
        epsil = -epsilon;
    else
        epsil = epsilon;
    end
    lin_objy.YData = [-epsil -epsil];
    lin_objx.YData = [ylims(1) -epsil];
    scat_obj.YData = -epsil;
    if gamma < 7
        txt_obj.Position = [r0,-epsil+diff(ylims)*0.05, 0];
    else
        txt_obj.Position = [r0*1.01,-epsil-diff(ylims)*0.03, 0];
    end
    
    
    k = 1/(gamma - 6);
    
    % Gen plot data
    ffun = @(r) 6*epsil*k*exp(gamma*(1 - r./r0)) - epsil*gamma*k*((r0./r).^6);
%     fpfun = @(r) -6*(gamma/r0)*epsil*k*exp(gamma*(1 - r./r0)) + 6*epsil*gamma*k*(r0^6)./(r.^7);
%     if gamma < 7
%         rt = fzero(fpfun,r0+10);
%     else
%         rt = fzero(fpfun,r0*0.6);
%     end
% 
%     fmin = ffun(rt);
%     fpmin = fpfun(rt);
    
    % Calculate potential
    f = ffun(r);
    
    % Derivative of potential
%     fprime = fpfun(r);
    % Update text
    title_obj.String = ['$\gamma_{ij} = ' num2str(gamma,'%.4f') '$' newline ...
        '$\varepsilon_{ij} = ' num2str(epsil,'%.1f') '$' newline '$\sigma_{ij} = 1.0$'];

    % Update plot Y-data
    for jdx = Trail_length:-1:2
        plot_objs(jdx).YData = plot_objs(jdx-1).YData;
    end
    plot_objs(1).YData = f;
    F(idx) = getframe(figh);
end





if savevideo
    N_F = numel(F);
    
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
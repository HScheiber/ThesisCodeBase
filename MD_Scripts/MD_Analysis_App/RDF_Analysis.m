% RDF_Analysis
function RDF_Analysis
%% Main Figure and Options
figtitle = 'RDF Analysis';
app.home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
app.icons = fullfile(app.home,'analysis','MD_Scripts','RDF_App','icons');
app.figh = figure('WindowState','maximized','NumberTitle','off',...
	'Name',figtitle,'Visible','on','Units','normalized');
app.FS = 24; % Font size
app.LW = 2; % Line width
app.Normalize_RDF = false; % When true, all RDF are normalized to their maximum = 1
app.IntegratePeakOn = false; % Switch that tells the script when the integrate peak button is 'on'

%% Initialize Axes Grid Properties
app.Num_LR = 3; % Number of panels in grid left to right
app.Num_UD = 3; % Number of panels in grid top to bottom
app.Xdat = cell(app.Num_UD,app.Num_LR);
app.Ydat = cell(app.Num_UD,app.Num_LR);
app.Left_Border = 0.18;
app.Right_Border = 0.01;
app.Top_Border = 0.08;
app.Bottom_Border = 0.1;
app.Text_Width = 0.06; % Room required for vertical text
app.Panel_Width = (1 - app.Left_Border - app.Right_Border)/app.Num_LR;
app.Panel_Height = (1 - app.Top_Border - app.Bottom_Border)/app.Num_UD;
app.RDF_xlim = [0 1.6];

%% Option panel
app.OptionsTitle = {'Salt' 'Thermostat'; 'Structure' 'Barostat'; 'Model' 'Sim. Anneal'};
app.OptionsMenu = {{'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'} {'nose-hoover T=250 C=0.2' 'v-rescale T=250 C=0.2' }; ...
    {'Rocksalt' 'Wurtzite' 'FiveFive'} {'Berendsen P=1.0 C=0.2' 'Parrinello-Rahman P=1.0 C=0.2' }; ...
    {'JC (SPC/E)' 'JC (TIP4PEW)' 'JC (TIP3P)' 'TF'} {'MD' 'Cell Optimized' 'Full Optimized'}};
BG_Num_LR = 2; % Number of panels in grid left to right
BG_Num_UD = 3; % Number of panels in grid top to bottom
BG_Left_Border = 0.01;
BG_Right_Border = 0.01;
BG_Top_Border = 0.00;
BG_Bottom_Border = 0.2;
BG_Panel_Width = (1 - BG_Left_Border - BG_Right_Border)/BG_Num_LR;
BG_Panel_Height = (1 - BG_Top_Border - BG_Bottom_Border)/BG_Num_UD;

%% Generate Buttons
for idx=1:app.Num_UD

    % Add button group
    app.bg(idx) = uibuttongroup(app.figh,'Position',[0 1-app.Top_Border-(idx)*app.Panel_Height ...
        app.Left_Border-app.Text_Width app.Panel_Height]); % [left bottom width height]

    % Add grid buttons to button group
    for ii=1:BG_Num_UD
        for jj=1:BG_Num_LR
            app.BGLabel{idx}(ii,jj) = uicontrol(app.bg(idx),'Style','text','String',...
                app.OptionsTitle{ii,jj},'FontSize',app.FS-8,'Units','Normalized','Position',...
                [BG_Left_Border+(jj-1)*BG_Panel_Width,1-BG_Top_Border-(ii)*BG_Panel_Height+app.Panel_Height*0.65,BG_Panel_Width,BG_Panel_Height*0.35]); % [left bottom width height]
            app.BGBox{idx}(ii,jj) = uicontrol(app.bg(idx),'Style','popupmenu','String',...
                app.OptionsMenu{ii,jj},'FontSize',app.FS-8,'Units','Normalized',...
                'Position',[BG_Left_Border+(jj-1)*BG_Panel_Width,1-BG_Top_Border-(ii)*BG_Panel_Height,BG_Panel_Width,BG_Panel_Height*0.65],...
                'Min',1,'Max',1,'Callback',{@LoadPlotsRDF});
        end
    end
    
    % Add time slider
    app.BGLabel{idx}(BG_Num_UD+1,1) = uicontrol(app.bg(idx),'Style','text','String',...
        'Time (0 ns)','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 BG_Bottom_Border*0.5 1-BG_Left_Border-BG_Right_Border BG_Bottom_Border*0.5]); % [left bottom width height]
    app.BGBox{idx}(BG_Num_UD+1,1) = uicontrol(app.bg(idx),'Style','slider','String',...
        '','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 0 1-BG_Left_Border-BG_Right_Border BG_Bottom_Border*0.5],...
        'Callback',{@LoadPlotsRDF},'SliderStep', [0.1 0.1], ...
        'Min', 0, 'Max', 10, 'Value',0,'Callback', @Round_slider);
end

% Make invisible background axis
app.InvisAxis = axes(app.figh,'Position',[0 0 1 1],'FontSize',app.FS,...
    'TickLabelInterpreter','latex','Box','off','XTick',[],'YTick',[],...
    'Color','none');

% Add x and y label placeholders
app.Xlabel = text(app.InvisAxis,'String','','FontSize',app.FS,...
    'HorizontalAlignment','center','Units','Normalized','Rotation',0,...
    'Interpreter','latex','Position',[0.5844 0.0347],...
    'VerticalAlignment','middle');

app.Ylabel = text(app.InvisAxis,'String','','FontSize',app.FS,...
    'HorizontalAlignment','center','Units','Normalized','Rotation',90,...
    'Interpreter','latex','Position',[0.1308 0.5166],...
    'VerticalAlignment','middle');

%% Turn off some toolbar options
set(0,'Showhidden','on')
ch = findall(app.figh,'tag','FigureToolBar');
UT = get(ch,'children');
set(UT(1:5),'Visible','Off','Separator','Off')
set(UT(1),'Separator','On')
set(UT(8),'Visible','Off','Separator','Off')

tbh = findall(app.figh,'Type','uitoolbar');

%% Add crosshair togglebutton to toolbar:
% 
% % Read an image
% [cdata,map] = imread(fullfile(app.home,'analysis','AppData','crosshair.gif'));
% 
% % Convert white pixels into a transparent background
% map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;
% 
% % Convert into 3D RGB-space
% cdata = ind2rgb(cdata,map);
% 
% % Set the button icon
% tbh = findall(app.figh,'Type','uitoolbar');
% app.CrosshairButton = uitoggletool(tbh,'CData',cdata,...
%     'Separator','off','Enable','off',...
%     'TooltipString','Data Crosshair','OffCallback',@CrosshairOff,...
%     'BusyAction','Cancel','tag','Exploration.Crosshair');

% Update crosshair button
%set(CrosshairButton,'OnCallback',{@crosshairfcn,X_Axis,Y_Axis,plot_extra,fs});
%CrosshairButton.Enable = 'on';
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Add integrate peak pushbutton to toolbar:
%-------------------------------------------------------------------------%
% Read an image
[cdata,map] = imread(fullfile(app.icons,'Integration.gif'));

% Convert white pixels into a transparent background
map(find(map(:,1)+map(:,2)+map(:,3)==3)) = NaN; %#ok<FNDSB>

% Convert into 3D RGB-space
cdata = ind2rgb(cdata,map);

% Set the button icon
inb = uipushtool(tbh,'CData',cdata,...
    'Separator','on','ClickedCallback',@IntegratePeak,...
    'TooltipString','Integrate Peak(s)',...
    'BusyAction','Cancel','tag','Integrate.Peaks');
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%% Add clear integrations pushbutton to toolbar:
%-------------------------------------------------------------------------%
% Read an image
[cdata,map] = imread(fullfile(app.icons,'ClearIntegration.gif'));

% Convert white pixels into a transparent background
map(find(map(:,1)+map(:,2)+map(:,3)==3)) = NaN; %#ok<FNDSB>

% Convert into 3D RGB-space
cdata = ind2rgb(cdata,map);

% Set the button icon
clrnt = uipushtool(tbh,'CData',cdata,...
    'Separator','off','ClickedCallback',@ClearIntegration,...
    'TooltipString','Clear Integrations',...
    'BusyAction','Cancel','tag','Integrate.Clear');
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%% Add normalize integration peaks pushbutton to toolbar:
%-------------------------------------------------------------------------%
% Read an image
[cdata,map] = imread(fullfile(app.icons,'Normalize.gif'));

% Convert white pixels into a transparent background
map(find(map(:,1)+map(:,2)+map(:,3)==3)) = NaN; %#ok<FNDSB>

% Convert into 3D RGB-space
cdata = ind2rgb(cdata,map);

% Set the button icon
inp = uipushtool(tbh,'CData',cdata,...
    'Separator','off','ClickedCallback',@Normalize,...
    'TooltipString','Normalize Integrations',...
    'BusyAction','Cancel','tag','Integrate.Normalize');
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% Add Analysis Type
app.AnalysisLabel = uicontrol(app.figh,'Style','text','String',...
    'Analysis: ','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0,0.9449,0.0650,0.0409]);
app.AnalysisBox = uicontrol(app.figh,'Style','popupmenu','String',...
    {'RDF' 'Energy vs Time' 'Temperature vs Time' 'Pressure vs Time' 'Volume vs Time'},...
    'FontSize',app.FS-5,'Units','Normalized','Position',[0.0625,0.946,0.1016,0.0415],...
    'Min',1,'Max',1,'Callback',{@LoadPlotsRDF});

%% Load Data
LoadPlotsRDF([],[]);

%% Set up empty stuff for integrate peaks function
% Set up cursor text
app.hText = uicontrol(app.figh,'Style','text','Units','Normalized',...
    'Position',[0.500 0.94 0.16 0.032],'backg',get(app.figh,'Color'),...
    'Fontsize',app.FS-5,'Fontweight','bold','Visible','off','Tag','Integrate.hText');

% Set up empty cursor lines
app.hCurv1 = line(nan(1,1),nan(1,1),'Color','red','Parent',app.InvisAxis,...
    'Visible','off','Tag','Integrate.hCurv');
app.hCurv2 = line(nan(1,1),nan(1,1),'Color','red','Parent',app.InvisAxis,...
    'Visible','off','Tag','Integrate.hCurv');

app.IntFirstStage = true;
app.Currenti = nan(1,1);
app.Currentj = nan(1,1);

%% Load Plots function
function LoadPlotsRDF(~,~)
    ps_per_ns = 1000;
    
    %% Grab plot settings
    Plot_Type = app.AnalysisBox.String{app.AnalysisBox.Value};
    Salt = cell(1,app.Num_UD);
    Structure = cell(1,app.Num_UD);
    Model = cell(1,app.Num_UD);
    Thermostat = cell(1,app.Num_UD);
    Barostat = cell(1,app.Num_UD);
    Sim_Anneal = cell(1,app.Num_UD);
    Time = cell(1,app.Num_UD); 
    for i=1:app.Num_UD % loop through the setting sections
        Salt{i} = app.OptionsMenu{1,1}{app.BGBox{i}(1,1).Value};

        Structure{i} = app.OptionsMenu{2,1}{app.BGBox{i}(2,1).Value};

        Model{i} = app.OptionsMenu{3,1}{app.BGBox{i}(3,1).Value};
        Model{i} = strrep(Model{i},' (SPC/E)','');
        Model{i} = strrep(Model{i},' (TIP3P)','3P');
        Model{i} = strrep(Model{i},' (TIP4PEW)','4P');

        Thermostat{i} = app.OptionsMenu{1,2}{app.BGBox{i}(1,2).Value};
        Thermostat{i} = strrep(strrep(strrep(strrep(Thermostat{i},'-','_'),' ','_'),'=',''),'.','P');

        Barostat{i} = app.OptionsMenu{2,2}{app.BGBox{i}(2,2).Value};
        Barostat{i} = strrep(strrep(strrep(strrep(Barostat{i},'-','_'),' ','_'),'=',''),'.','P');

        Sim_Anneal{i} = app.OptionsMenu{3,2}{app.BGBox{i}(3,2).Value};
        Sim_Anneal{i} = strrep(strrep(strrep(strrep(Sim_Anneal{i},'-','_'),' ','_'),'=',''),'.','P');
        
        Time{i} = app.BGBox{i}(4,1).Value; % Time in ns
    end
    
    % Clear old axes if they exist
    ClearAxes
    
    if strcmpi(Plot_Type,'RDF')
        app.Xlabel.String = 'r (nm)';
        if app.Normalize_RDF
            app.Ylabel.String = 'Normalized Radial Dist. Function';
        else
            app.Ylabel.String = 'Radial Dist. Function';
        end
        
        app.Num_LR = 3; % Number of panels in grid left to right
        app = GenerateAxes(app);
        
        DOI = cell(app.Num_UD,app.Num_LR);
        Data_Type = {'MH' 'MM' 'HH'};
        Colours = {'k' 'r' 'b'};
        
        for i=1:app.Num_UD
            [Metal,Halide] = Separate_Metal_Halide(Salt{i});
            AxisTitle = {[Metal '-' Halide] [Metal '-' Metal] [Halide '-' Halide]};
            
            % Update time
            app.BGLabel{i}(4,1).String = ...
                regexprep(app.BGLabel{1}(4,1).String,'\([0-9|\.]+ ns\)',['(' num2str(Time{i}) ' ns)']);
            
            % Load data
            LH = load(fullfile(app.home,'data','PURE_SALT_MD',[Salt{i} '_' Structure{i} '_MD_Data.mat']));
            
            for j=1:app.Num_LR
            
                % Grab data of interest
                DOI{i,j} = LH.Data.(Salt{i}).(Structure{i}).(Model{i}). ...
                    (Thermostat{i}).(Barostat{i}).(Sim_Anneal{i}).RDF.(Data_Type{j});
                Times = [DOI{i,j}{:,1}];
                
                Index = abs(Times - Time{i}*ps_per_ns) < 1e-4;
                if sum(Index) == 1
                    app.Xdat{i,j} = DOI{i,j}{Index,2}(:,1);
                    app.Ydat{i,j} = DOI{i,j}{Index,2}(:,2);
                    
                    if app.Normalize_RDF
                        Ydat = app.Ydat{i,j}./max(app.Ydat{i,j});
                        ylimits = [0 1.1];
                        yticks = 0:0.25:1;
                    else
                        Ydat = app.Ydat{i,j};
                        if Time{i} == 0
                            ylimits = [0 115];
                            yticks = 0:25:110;
                        else
                            ylimits = [0 25];
                            yticks = 0:5:20;
                        end
                    end
                    
                    plot(app.axh(i,j),app.Xdat{i,j},Ydat,'Color',Colours{j},...
                        'LineWidth',app.LW,'LineStyle','-',...
                        'Visible','on');
                    
                    set(app.axh(i,j),'ylim',ylimits,'xlim',app.RDF_xlim,...
                        'FontSize',app.FS,...
                        'TickLabelInterpreter','latex','Box','on',...
                        'XMinorTick','on','YMinorTick','on')
                    
                    % Add in text label
                    text(app.axh(i,j),'String',AxisTitle{j},'FontSize',app.FS-3,...
                        'HorizontalAlignment','center','Units','Normalized','Rotation',0,...
                        'Interpreter','latex','Position',[0.5 0.9],...
                        'VerticalAlignment','middle');
                    
                    if i == app.Num_UD
                        set(app.axh(i,j),'xtick',0:0.5:3)
                    else
                        set(app.axh(i,j),'xtick',[])
                    end
                    
                    if j == 1
                        set(app.axh(i,j),'ytick',yticks)
                    else
                        set(app.axh(i,j),'ytick',[])
                    end
                end
            end
        end
    else
        app.Num_LR = 1; % Number of panels in grid left to right
        app = GenerateAxes(app);
        
        DOI = cell(app.Num_UD,app.Num_LR);
        app.Xlabel.String = 'Time (ps)';
        
        if strcmpi(Plot_Type,'Energy vs Time')
            Index = 2;
            app.Ylabel.String = 'Energy (kJ/mol)';
        elseif strcmpi(Plot_Type,'Temperature vs Time')
            Index = 3;
            app.Ylabel.String =  'Temperature (K)';
        elseif strcmpi(Plot_Type,'Pressure vs Time')
            Index = 4;    
            app.Ylabel.String =  'Pressure (Bar)';
        elseif strcmpi(Plot_Type,'Volume vs Time')
            Index = 5;
            app.Ylabel.String =  'Volume (nm$^{3}$)';
        end
        
        for i=1:app.Num_UD
            
            % Load data
            LH = load(fullfile(app.home,'data','PURE_SALT_MD',[Salt{i} '_' Structure{i} '_MD_Data.mat']));
            
            % Grab data of interest
            DOI{i} = LH.Data.(Salt{i}).(Structure{i}).(Model{i}). ...
                (Thermostat{i}).(Barostat{i}).(Sim_Anneal{i}).Energy(:,[1 Index]);
            
            % Plot data
            app.Xdat{i,1} = DOI{i}(:,1);
            app.Ydat{i,1} = DOI{i}(:,2);

            plot(app.axh(i,1),app.Xdat{i,1},app.Ydat{i,1},'Color','k',...
                'LineWidth',app.LW,'LineStyle','-',...
                'Visible','on');
            set(app.axh(i,1),'ylim',[-Inf Inf],'xlim',[-Inf Inf],...
                'FontSize',app.FS,...
                'TickLabelInterpreter','latex','Box','on',...
                'XMinorTick','on','YMinorTick','on')

            if i ~= app.Num_UD
                set(app.axh(i,1),'xtick',[])
            end
        end
    end
end

%% Function to generate grid of axes
function app = GenerateAxes(app)
    app.Panel_Width = (1 - app.Left_Border - app.Right_Border)/app.Num_LR;
    app.Panel_Height = (1 - app.Top_Border - app.Bottom_Border)/app.Num_UD;
    for i=1:app.Num_UD
        for j=1:app.Num_LR
            % Generate axis
            ax_pos = [app.Left_Border+(j-1)*app.Panel_Width ...
                1-app.Top_Border-(i)*app.Panel_Height app.Panel_Width app.Panel_Height]; % [left bottom width height]
            app.axh(i,j) = axes(app.figh,'Position',ax_pos,'FontSize',app.FS,...
                'TickLabelInterpreter','latex','Box','on','XMinorTick','on',...
                'YMinorTick','on','XTickLabel','','YTickLabel','');
        end
    end
end

%% Function to remove all axes
function ClearAxes
    app.Baseline = cell(0,0);
    app.DetrendLine = cell(0,0);
    app.Atxth = cell(0,0);
    app.Area = [];
    if isfield(app,'axh')
        for i=1:app.Num_UD
            for j=1:app.Num_LR
                delete(app.axh(i,j));
            end
        end
    end
end

%% Prevents slider from moving between steps
function Round_slider(hObject,~)
    newval = hObject.Value;	% get value from the slider
    newval = round(newval);	% round off this value
    set(hObject, 'Value', newval);	% set slider position to rounded off value
    LoadPlotsRDF([],[]);
end

function IntegratePeak(varargin)
    
    % If already enabled, turn off
    if app.IntegratePeakOn
        app.figh = EnableUIButtons(app.figh);
        app.IntegratePeakOn = false;
               
        set(app.hText,'Visible','off')
        set(app.hCurv1,'Parent',app.InvisAxis,'Visible','off');
        set(app.hCurv2,'Parent',app.InvisAxis,'Visible','off');
        set(app.figh,'WindowButtonMotionFcn','')
    else
        % Disable all buttons while active
        app.figh = DisableUIButtons(app.figh);
        app.IntegratePeakOn = true;
        app.IntFirstStage = true;
        
        set(app.figh,'WindowButtonMotionFcn',@IntegrateMouseMoveFcn)
    end
end

function Normalize(varargin)
    % This function serves to renormalize the measured peak areas
    N = length(app.Area);
    

        %CurrentAreas(i) = app.Area{i};
        %delete(app.Baseline{i})
        %delete(app.DetrendLine{i})
        %delete(app.Atxth{i})
    
    MinArea = min(app.Area);
    NewAreas = app.Area./MinArea;
    
    % Modify text
    for i=1:N
        app.Atxth{i}.String = num2str(NewAreas(i),'%4.3f');
    end

end

function ClearIntegration(varargin)
    % Clear peak integrations
    for i=1:length(app.DetrendLine)
        delete(app.Baseline{i})
        delete(app.DetrendLine{i})
        delete(app.Atxth{i})
    end
    app.Baseline = cell(0,0);
    app.DetrendLine = cell(0,0);
    app.Atxth = cell(0,0);
end

function IntegrateMouseMoveFcn(~,~)
    
    % Grab the mouse position
    MousePosition = get(app.figh,'CurrentPoint');
    MouseX = MousePosition(1);
    MouseY = MousePosition(2);
    AxisFound = false;
    
    if app.IntFirstStage
        % Loop through the axes to find if the mouse is over any of them
        for i=1:app.Num_UD
            for j=1:app.Num_LR
                
                % Position of axis within figure
                AxisPos = app.axh(i,j).Position; % [left bottom width height]

                % Create the x and y lower (L) and upper (U) bounds
                AxisXLower = AxisPos(1);
                AxisXUpper = (AxisPos(1)+AxisPos(3));
                AxisYLower = AxisPos(2);
                AxisYUpper = (AxisPos(2)+AxisPos(4));

                % Check to see if the mouse is over the axis
                if MouseX >= AxisXLower && MouseX <= AxisXUpper && ...
                        MouseY >= AxisYLower && MouseY <= AxisYUpper

                    AxisFound = true;
                    app.figh.CurrentAxes = app.axh(i,j);
                    app.Currenti = i;
                    app.Currentj = j;

                    % Get mouse location within current axis
                    Axispt = get(app.figh.CurrentAxes,'CurrentPoint');

                    % Create custom invisible cursor for use while within axes
                    set(app.figh,'PointerShapeCData',ones(16, 16)*nan,...
                        'WindowButtonDownFcn', @ClickFcn,'Pointer','custom');

                    % Get Axes limits
                    yl = ylim(app.figh.CurrentAxes);

                    % Update cursor line position
                    set(app.hCurv1,'YData',yl,'XData',[Axispt(1), Axispt(1)],...
                        'Visible','on','Parent',app.figh.CurrentAxes);

                    set(app.hText,'String', sprintf('r = %f nm',...
                        Axispt(1)),'Visible','on');
                    break
                end
            end
            if AxisFound
                break
            end
        end
        
    % Second step, the chosen axis is fixed
    else
        % Position of axis within figure
        AxisPos = app.figh.CurrentAxes.Position; % [left bottom width height]

        % Create the x and y lower (L) and upper (U) bounds
        AxisXLower = AxisPos(1);
        AxisXUpper = (AxisPos(1)+AxisPos(3));
        AxisYLower = AxisPos(2);
        AxisYUpper = (AxisPos(2)+AxisPos(4));
        
        % Check to see if the mouse is over the current axis
        if MouseX >= AxisXLower && MouseX <= AxisXUpper && ...
                MouseY >= AxisYLower && MouseY <= AxisYUpper

            AxisFound = true;

            % Get mouse location within current axis
            Axispt = get(app.figh.CurrentAxes,'CurrentPoint');

            % Create custom invisible cursor for use while within axes
            set(app.figh,'PointerShapeCData',ones(16, 16)*nan,...
                'WindowButtonDownFcn', @ClickFcn,'Pointer','custom');

            % Get Axes limits
            yl = ylim(app.figh.CurrentAxes);

            % Update cursor line position
            set(app.hCurv2,'YData',yl,'XData',[Axispt(1), Axispt(1)],...
                'Visible','on','Parent',app.figh.CurrentAxes);

            set(app.hText,'String', sprintf('r = %f nm',...
                Axispt(1)),'Visible','on');
        end
    end
    
    if ~AxisFound
        % If cursor not over an axis, return cursor to normal
        set(app.figh,'WindowButtonDownFcn','','Pointer','arrow');
        set(app.hCurv1,'Parent',app.InvisAxis,'Visible','off');
        set(app.hText,'Visible','off');
    end
end


function ClickFcn(varargin)
    
    % Get mouse location and axis info
    outpt = get(app.figh.CurrentAxes, 'CurrentPoint');
    yl = ylim(app.figh.CurrentAxes);
    xl = xlim(app.figh.CurrentAxes);
    if outpt(1) >= xl(1) && outpt(1) <= xl(2) && outpt(1,2) >= yl(1) && outpt(1,2) <= yl(2)
        if app.IntFirstStage
            app.X1_Integrate = outpt(1);
        else % Second stage
            app.X2_Integrate = outpt(1);
            
            % Determine which x is lower and which is upper bound, if x1 = x2, area = 0
            xlower = min([app.X1_Integrate app.X2_Integrate]);
            xupper = max([app.X1_Integrate app.X2_Integrate]);
            
            % No area selected then abort
            if xlower == xupper
                return
            end
            
            % Obtain data from selected section, correct for baseline
            Selectedx = app.Xdat{app.Currenti,app.Currentj};
            Selectedy = app.Ydat{app.Currenti,app.Currentj};

            Selectedy(Selectedx>xupper | Selectedx<xlower)= [];
            Selectedx(Selectedx>xupper | Selectedx<xlower)= [];
            
            % Calculate area under curve
            DetrendY = cumtrapz(Selectedx,Selectedy);
            app.Area(end+1) = DetrendY(end);
            
            % Plot the baseline
            hold(app.figh.CurrentAxes,'on')
            app.Baseline{end+1} = plot(app.figh.CurrentAxes,[Selectedx(1) Selectedx(end)],...
                [0 0],'Color','red');
            
            % Output a scaled detrend graph over the integrated peak
            app.DetrendLine{end+1} = plot(app.figh.CurrentAxes,Selectedx,DetrendY,'Color','green');
            
            % Determine x position for area under curve text
            xpos = (Selectedx(1)+Selectedx(end))/2;

            % Determine y position for area under curve text
            ypos = app.Area(end)*1.05;

            % Output the area under the curve on the plot
            Areatxt = num2str(app.Area(end),3);
            app.Atxth{end+1} = text(xpos,ypos,Areatxt,'HorizontalAlignment','center','Color','red');
            
            % Reset stuff
            set(app.hText,'Visible','off')
            set(app.hCurv1,'Parent',app.InvisAxis,'Visible','off');
            set(app.hCurv2,'Parent',app.InvisAxis,'Visible','off');
            set(app.figh,'windowbuttonmotionfcn','',...
                'WindowButtonDownFcn','',...
                'Pointer','arrow');
            app.IntegratePeakOn = false;
            EnableUIButtons(app.figh);
            return            
        end
    end
    
    % Set to second stage
    if app.IntFirstStage
        app.IntFirstStage = false;
    end
end

end
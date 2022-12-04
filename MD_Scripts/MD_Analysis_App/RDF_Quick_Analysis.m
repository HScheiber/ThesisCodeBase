% RDF_Analysis

function RDF_Quick_Analysis
warning('off','MATLAB:hg:uicontrol:StringMustBeNonEmpty')
apploc = fullfile(fileparts(mfilename('fullpath')),'AppData','savestate.mat');
global app
app = InitializeApp;
if isfile(apploc)
    % Load previous app state from file
    DeleteAllCallback
    app = LoadAppState(app,apploc);
    LoadPlotState;
    RefactorPlots; % Reload plots
else
    % Load app state and Generate Axis
    DeleteAllCallback;
end

% Attempt to clear any previous temp folders associated with program
try
    rmdir(app.tempfolder,'s')
catch
    % do nothing
end

% Make temporary folder to hold some temp files
mkdir(app.tempfolder);

% get parallel pool number of cores from local machine
app.PrefCores = feature('numcores');    

function app = InitializeApp
    % This function initializes all the app GUI parts
    %% Main Figure and Options
    figtitle = 'RDF Quick Analysis';
    if ispc
        app.datdir = 'D:\Molten_Salts_MD'; %'C:\Users\Hayden\Documents\Patey_Lab\Results';
        app.home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
        app.vmd_loc = 'C:\"Program Files (x86)"\VMD\vmd.exe';
        app.FS = 24; % Font size
        app.gmx = 'gmx_d';
        app.wsl = 'wsl source ~/.bashrc; ';
        app.pipe = ' ^| ';
    else
%         app.datdir = '/media/hayden/Windows/Users/Hayden/Documents/Patey_Lab/Results';
        app.datdir = '/home/user/project/Molten_Salts_MD';
%         app.home = '/media/hayden/Windows/Users/Hayden/Documents/Patey_Lab/ThesisCodeBase';
        app.home = '/home/user/ThesisCodeBase';
        app.vmd_loc = '/usr/local/bin/vmd';
        app.FS = 20; % Font size
        app.gmx = 'gmx_d';
        app.wsl = 'source ~/.bashrc; ';
        app.pipe = ' | ';
    end
    app.AppData = fullfile(fileparts(mfilename('fullpath')),'AppData');
    app.figh = figure('WindowState','maximized','NumberTitle','off',...
        'Name',figtitle,'Visible','on','Units','normalized',...
        'CloseRequestFcn',@CloseAppCallback);
    
    app.auxfig = figure('NumberTitle','off',...
        'Name','Trajectory Selector','Visible','off','Units','normalized',...
        'CloseRequestFcn',@CloseTrajSelCallback,...
        'MenuBar','none','toolbar','none','CloseRequestFcn',@CloseTrajSelCallback);
    
    app.tempfolder = tempname([tempdir filesep 'RDF_App']);
    app.LW = 3; % Line width
    app.Normalize_RDF = false; % When true, all RDF are normalized to their maximum = 1
    app.IntegratePeakOn = false; % Switch that tells the script when the integrate peak button is 'on'
    app.ColScheme = {'qual' 'Set2'}; % 2D Colour scheme
    app.Colours_3D = flipud(cbrewer('seq','YlOrRd',200,'PCHIP'));
    app.timestep = 0.001; % picoseconds
    app.HistBinQ = 0.01; % Histogram bin width for Ql distributions
    app.HistBinW = 0.01; % Histogram bin width for Wl distributions
    app.l = 2:2:10; % series of l to select for Wl and Ql calculations

    %% Option panel
    app.DataTypeMenu = {'MD' 'Cell Optimized' 'Full Optimized'};
    app.SaltMenu = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl' 'CsCl'};
    app.StructureMenu = {'Rocksalt' 'Wurtzite' 'BetaBeO' 'CsCl' 'FiveFive' 'NiAs' 'Sphalerite'};
    app.ModelMenu = {'JC (SPC/E)' 'JC (TIP4P-EW)' 'JC (TIP3P)' 'TF'};

    % Positions
    app.BG_Position = [0,0,0.139583,1];

    %% Generate Aux UI Interface
    app.TJS.bg_time  = uibuttongroup(app.auxfig,'Position',[0 0 1/3 0.5],'Visible','On'); % [left bottom width height]
    app.TJS.bg_type  = uibuttongroup(app.auxfig,'Position',[0 0.5 1/3 0.5],'Visible','On','Tag','type'); % [left bottom width height]
    app.TJS.bg_order = uibuttongroup(app.auxfig,'Position',[1/3 0.5 1/3 0.5],'Visible','On','Tag','order'); % [left bottom width height]
    app.TJS.bg_neighbour = uibuttongroup(app.auxfig,'Position',[2/3 0.5 1/3 0.5],'Visible','On','Tag','neighbour'); % [left bottom width height]
    app.TJS.bg_buttons = uibuttongroup(app.auxfig,'Position',[1/3 0 2/3 0.15],'Visible','On'); % [left bottom width height]
    app.TJS.bg_filter = uibuttongroup(app.auxfig,'Position',[1/3 0.15 2/3 0.35],'Visible','On','Tag','filter'); % [left bottom width height]
    
    % Primary Buttons
    app.TJS.Select_Viewer_text = uicontrol(app.TJS.bg_buttons,'Style','popupmenu','String',...
        {'OVITO' 'VMD' 'PyMOL'},'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 0 1/4 1/2]); % [left bottom width height]
    app.TJS.Select_Viewer = uicontrol(app.TJS.bg_buttons,'Style','text','String',...
        'Viewer:','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 1/2 1/4 1/2]); % [left bottom width height]
    app.TJS.View_Button = uicontrol(app.TJS.bg_buttons,'Style','PushButton','String',...
        'View','FontSize',app.FS-5,'Units','Normalized',...
        'Position',[1/4 0 1/4 1],'Callback',@ViewTrajectoryCallback); % [left bottom width height]
    app.TJS.Save_Button = uicontrol(app.TJS.bg_buttons,'Style','PushButton','String',...
        'Save','FontSize',app.FS-5,'Units','Normalized',...
        'Position',[2/4 0 1/4 1],'Callback',@SaveTrajSelCallback); % [left bottom width height]
    app.TJS.Cancel_Button = uicontrol(app.TJS.bg_buttons,'Style','PushButton','String',...
        'Cancel','FontSize',app.FS-5,'Units','Normalized',...
        'Position',[3/4 0 1/4 1],'Callback',@CloseTrajSelCallback); % [left bottom width height]
    
    % Time stuff
    app.TJS.InitTime_Text = uicontrol(app.TJS.bg_time,'Style','text','String',...
        'Initial Time (ns)','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 5/6 1 1/6]); % [left bottom width height]
    app.TJS.InitTime_Edit = uicontrol(app.TJS.bg_time,'Style','edit','String',...
        '','FontSize',app.FS-10,'Units','Normalized',...
        'Position',[1/3 4/6 1/3 1/6],... % [left bottom width height]
        'Callback',@InitialTimeCallback);
    app.TJS.FinTime_Text = uicontrol(app.TJS.bg_time,'Style','text','String',...
        'Final Time (ns)','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 3/6 1 1/6]); % [left bottom width height]
    app.TJS.FinTime_Edit = uicontrol(app.TJS.bg_time,'Style','edit','String',...
        '','FontSize',app.FS-10,'Units','Normalized',...
        'Position',[1/3 2/6 1/3 1/6],... % [left bottom width height]
        'Callback',@FinalTimeCallback);
    app.TJS.FrameTime_Text = uicontrol(app.TJS.bg_time,'Style','text','String',...
        'Time Per Frame (ps)','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 1/6 1 1/6]); % [left bottom width height]
    app.TJS.FrameTime_Edit = uicontrol(app.TJS.bg_time,'Style','edit','String',...
        '','FontSize',app.FS-10,'Units','Normalized',...
        'Position',[1/3 0/6 1/3 1/6],... % [left bottom width height]
        'Callback',@TimeStepsCallback);
    
    % Order module
    app.TJS.Order_Text = uicontrol(app.TJS.bg_order,'Style','text','String',...
        'Order Parameter: ','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 6/7 7/10 1/7]); % [left bottom width height]
    app.TJS.Order_Menu = uicontrol(app.TJS.bg_order,'Style','popupmenu','String',...
        {'Ql' 'Wl'},'FontSize',app.FS-11,'Units','Normalized',...
        'Position',[7/10 6/7 3/10 1/7],'Callback',@OrderParameterChanged); % [left bottom width height]
    app.TJS.Order_TextL = uicontrol(app.TJS.bg_order,'Style','text','String',...
        'L = ','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.0258 0.7072 0.2000 0.1429]); % [left bottom width height]
    app.TJS.Order_MenuL = uicontrol(app.TJS.bg_order,'Style','popupmenu','String',...
        {'2' '4' '6' '8' '10'},'FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.2 5/7 0.15 1/7]); % [left bottom width height]
    app.TJS.Order_TextRange = uicontrol(app.TJS.bg_order,'Style','text','String',...
        [char(8804) ' Ql ' char(8804)] ,'FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.5659 0.7072 0.3000 0.1429]); % [left bottom width height]
    app.TJS.Order_RangeMin = uicontrol(app.TJS.bg_order,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.4424 0.7339 0.1627 0.1227],'Callback',@OrderParameterRangeChanged); % [left bottom width height]
    app.TJS.Order_RangeMax = uicontrol(app.TJS.bg_order,'Style','edit','String',...
        '1','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.8274 0.7339 0.1627 0.1251],'Callback',@OrderParameterRangeChanged); % [left bottom width height]
    app.TJS.Order_Button = uicontrol(app.TJS.bg_order,'Style','PushButton','String',...
        'Add Filter','FontSize',app.FS-11,'Units','Normalized','tag','Type',...
        'Position',[1/7 0 2/3 1/7],'Callback',@AddFilterCallback); % [left bottom width height]
    
    % Ql and Wl cluster selection options
    app.TJS.bg_order2 = uibuttongroup(app.TJS.bg_order,'Position',[0 2/7 1 3/7],... % [left bottom width height]
        'BorderType','none');
    app.TJS.Order_QlRadius = uicontrol(app.TJS.bg_order2,'Style','radiobutton','String',...
        ['Search Rad. (' char(197) ')'],'FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 2/3 0.7 1/3],... % [left bottom width height]
        'Value',1,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_QlNeighbours = uicontrol(app.TJS.bg_order2,'Style','radiobutton','String',...
        'N. Neighbours','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 1/3 0.7 1/3],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_FirstShell = uicontrol(app.TJS.bg_order2,'Style','checkbox','String',...
        '1st Shell','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 0 0.5 1/3],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_SecondShell = uicontrol(app.TJS.bg_order2,'Style','checkbox','String',...
        '2nd Shell','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.5 0 0.5 1/3],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_QlEditMin = uicontrol(app.TJS.bg_order2,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.5888 0.6658 0.1433 0.3012],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_QlEdit = uicontrol(app.TJS.bg_order2,'Style','edit','String',...
        '3','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.7527 0.6658 0.1433 0.3012],... % [left bottom width height]
        'Callback',@QlEditCallbackTraj);
    app.TJS.Order_NeighbEdit = uicontrol(app.TJS.bg_order2,'Style','edit','String',...
        '6','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.5888 0.3519 0.1433 0.3012],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallbackTraj);
    app.TJS.Order_AllPoints = uicontrol(app.TJS.bg_order,'Style','radiobutton','String',...
        'Each t','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 1/7 1/3 1/7],... % [left bottom width height]
        'Value',0);
    app.TJS.Order_OnePoint = uicontrol(app.TJS.bg_order,'Style','radiobutton','String',...
        'Time (ps) = ','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.2664 0.1429 0.4147 0.1429],... % [left bottom width height]
        'Value',0);
    app.TJS.Order_OnePointEdit = uicontrol(app.TJS.bg_order,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.6532 0.1487 0.3369 0.1227],... % [left bottom width height]
        'Value',0,'Callback',@OrderTimeEditCallback);
    
    % Type module
    app.TJS.Type_Text = uicontrol(app.TJS.bg_type,'Style','text','String',...
        'Select by Type','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 2/3 1 1/6]); % [left bottom width height]
    app.TJS.Type_Menu = uicontrol(app.TJS.bg_type,'Style','popupmenu','String',...
        {'All' 'Metal','Halide'},'FontSize',app.FS-10,'Units','Normalized',...
        'Position',[1/3 0.5 1/3 1/6]); % [left bottom width height]
    app.TJS.Type_Button = uicontrol(app.TJS.bg_type,'Style','PushButton','String',...
        'Add Filter','FontSize',app.FS-10,'Units','Normalized','tag','Type',...
        'Position',[1/6 0 2/3 1/7],'Callback',@AddFilterCallback); % [left bottom width height]
    
    % Neighbour module
    app.TJS.Neighbour_Text = uicontrol(app.TJS.bg_neighbour,'Style','text','String',...
        'Number of Neighbours (N):','FontSize',app.FS-10,'Units','Normalized','Position',...
        [0 5/6 1 1/6]); % [left bottom width height]
    app.TJS.Neighbour_TextRange = uicontrol(app.TJS.bg_neighbour,'Style','text','String',...
        [char(8804) ' N ' char(8804)] ,'FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.5659-0.2 0.7072 0.3000 0.1429]); % [left bottom width height]
    app.TJS.Neighbour_RangeMin = uicontrol(app.TJS.bg_neighbour,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.4424-0.2 0.7339 0.1627 0.1227],'Callback',@NeighbourRangeChanged); % [left bottom width height]
    app.TJS.Neighbour_RangeMax = uicontrol(app.TJS.bg_neighbour,'Style','edit','String',...
        '6','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.8274-0.2 0.7339 0.1627 0.1251],'Callback',@NeighbourRangeChanged); % [left bottom width height]
    app.TJS.bg_neighbour2 = uibuttongroup(app.TJS.bg_neighbour,'Position',[0 2/6 1 2/6],... % [left bottom width height]
        'BorderType','none');
    app.TJS.Neighbour_Search_Radius = uicontrol(app.TJS.bg_neighbour2,'Style','radiobutton','String',...
        ['Search Rad. (' char(197) ')'],'FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 1/2 0.7 1/2],... % [left bottom width height]
        'Value',1,'Callback',@NNEditCallbackTraj);
    app.TJS.Neighbour_FirstShell = uicontrol(app.TJS.bg_neighbour2,'Style','checkbox','String',...
        '1st Shell','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 0 0.5 1/2],... % [left bottom width height]
        'Value',0,'Callback',@NNEditCallbackTraj);
    app.TJS.Neighbour_SecondShell = uicontrol(app.TJS.bg_neighbour2,'Style','checkbox','String',...
        '2nd Shell','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.5 0 0.5 1/2],... % [left bottom width height]
        'Value',0,'Callback',@NNEditCallbackTraj);
    app.TJS.Neighbour_Search_Radius_Min = uicontrol(app.TJS.bg_neighbour2,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.5888 0.57 0.1433 0.3012],... % [left bottom width height]
        'Value',0,'Callback',@NNEditCallbackTraj);
    app.TJS.Neighbour_Search_Radius_Max = uicontrol(app.TJS.bg_neighbour2,'Style','edit','String',...
        '3','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.7527 0.57 0.1433 0.3012],... % [left bottom width height]
        'Callback',@NNEditCallbackTraj);
    app.TJS.Neighbour_AllPoints = uicontrol(app.TJS.bg_neighbour,'Style','radiobutton','String',...
        'Each t','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0 1/6 1/3 1/6],... % [left bottom width height]
        'Value',0);
    app.TJS.Neighbour_OnePoint = uicontrol(app.TJS.bg_neighbour,'Style','radiobutton','String',...
        'Time (ps) = ','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.2664 1/6 0.4147 0.1429],... % [left bottom width height]
        'Value',0);
    app.TJS.Neighbour_OnePointEdit = uicontrol(app.TJS.bg_neighbour,'Style','edit','String',...
        '0','FontSize',app.FS-11,'Units','Normalized',...
        'Position',[0.6532 0.3422-1/6 0.3369 0.1227],... % [left bottom width height]
        'Value',0,'Callback',@OrderTimeEditCallback);    
    app.TJS.Neighbour_Button = uicontrol(app.TJS.bg_neighbour,'Style','PushButton','String',...
        'Add Filter','FontSize',app.FS-10,'Units','Normalized','tag','Type',...
        'Position',[1/6 0 2/3 1/7],'Callback',@AddFilterCallback); % [left bottom width height]
    
    % Filter module
    app.TJS.filter_notes = uicontrol(app.TJS.bg_filter,'Style','edit','String',...
        '','FontSize',app.FS-12,'Units','Normalized',...
        'Position',[0.25 0 0.75 6/7],'Enable','off',... % [left bottom width height]
        'Min',1,'Max',1000,'HorizontalAlignment','Left');
    app.TJS.filter_title = uicontrol(app.TJS.bg_filter,'Style','text','String',...
        'Current Filters','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.25 6/7 0.75 1/7]); % [left bottom width height]
    app.TJS.filter_clear = uicontrol(app.TJS.bg_filter,'Style','pushbutton','String',...
        'Clear Filters','FontSize',app.FS-11,'Units','Normalized','Position',...
        [0.025 3/7 0.20 1/7],'Callback',@ClearFiltersCallback); % [left bottom width height]
    app.TJS.current_filters = {};
    
    %% Generate Main UI Interface
    app.noteplusminus = uicontrol(app.figh,'Style','PushButton','String',...
        'Notes','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1391 0.9747 0.0419 0.0262],... % [left bottom width height]
        'Callback',@AddNoteCallback);
    
    app.Notepad = uicontrol(app.figh,'Style','edit','String',...
        '','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1393 0.6263 0.1354 0.3485],'Visible','Off',...
        'Min',1,'Max',1000,'HorizontalAlignment','Left','Callback',@EditNoteCallBack);
    
    % Button group
    app.bg = uibuttongroup(app.figh,'Position',app.BG_Position,'Visible','On'); % [left bottom width height]
    
    dv = 1/14; % Rows of buttons scaling within button group
    
    % Data Type
    app.BGLabel(1) = uicontrol(app.bg,'Style','text','String',...
        'Data Type','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623 1 0.0342]); % [left bottom width height]
    app.DataTypeh = uicontrol(app.bg,'Style','popupmenu','String',...
        app.DataTypeMenu,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1 0.9294 0.8 0.0342],... % [left bottom width height]
        'Min',1,'Max',1,'Callback',@DataTypeCallback);

    % Salt
    app.BGLabel(2) = uicontrol(app.bg,'Style','text','String',...
        'Salt','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-dv 1 0.0342]); % [left bottom width height]
    app.Salth = uicontrol(app.bg,'Style','popupmenu','String',...
        app.SaltMenu,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1 0.9294-dv 0.8 0.0342],... % [left bottom width height]
        'Min',1,'Max',1,'Callback',@SaltCallback);

    % Structure
    app.BGLabel(3) = uicontrol(app.bg,'Style','text','String',...
        'Structure','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-2*dv 1 0.0342],'Visible','Off'); % [left bottom width height]
    app.Structureh = uicontrol(app.bg,'Style','popupmenu','String',...
        app.StructureMenu,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1 0.9294-2*dv 0.8 0.0342],... % [left bottom width height]
        'Min',1,'Max',1,'Visible','Off');

    % Model
    % Find models associated with selected salt
    files = dir([app.datdir filesep app.SaltMenu{1}]);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    % remove '.' and '..'
    dfolders = {subFolders(~ismember({subFolders(:).name},{'.','..'})).name};

    app.Modelhtxt = uicontrol(app.bg,'Style','text','String',...
        'Model:','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-3*dv 0.4 0.0342]); % [left bottom width height]
    app.Modelh = uicontrol(app.bg,'Style','popupmenu','String',...
        dfolders,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.1 0.9294-3*dv 0.8 0.0342]); % [left bottom width height]
    
    app.Modelfilter = uibuttongroup(app.bg,'Position',[0.4 0.965-3*dv 0.6 0.0342],...
        'BorderType','none','Visible','On');
    app.Modelfilter_b1 = uicontrol(app.Modelfilter,'Style',...
                      'radiobutton',...
                      'String','New',...
                      'Units','Normalized',...
                      'Value',0,'FontSize',app.FS-8,...
                      'Position',[0 0 0.5 1],...
                      'Callback',@DataTypeCallback);
    app.Modelfilter_b2 = uicontrol(app.Modelfilter,'Style','radiobutton',...
                      'String','All',...
                      'Units','Normalized',...
                      'Value',1,'FontSize',app.FS-8,...
                      'Position',[0.5 0 0.5 1],...
                      'Callback',@DataTypeCallback);
    

    % Initial Time edit box
    app.BGLabel(5) = uicontrol(app.bg,'Style','text','String',...
        'Initial Time (ns)','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-4*dv 1 0.0342]); % [left bottom width height]
    app.InitialTimeh = uicontrol(app.bg,'Style','edit','String',...
        '0','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.3 0.9294-4*dv 0.4 0.0342],...
        'Callback',...
        @InitialTimeCallback);

    % Final Time edit box
    app.BGLabel(6) = uicontrol(app.bg,'Style','text','String',...
        'Final Time (ns)','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-5*dv 1 0.0342]); % [left bottom width height]
    app.FinalTimeh = uicontrol(app.bg,'Style','edit','String',...
        '30','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.3 0.9294-5*dv 0.4 0.0342],...
        'Callback',@FinalTimeCallback);

    % Add Time per Frame
    app.BGLabel(7) = uicontrol(app.bg,'Style','text','String',...
        'Time Per Frame (ps):','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-6*dv 1 0.0342]); % [left bottom width height]
    app.TimePerFrameh = uicontrol(app.bg,'Style','edit','String',...
        '1','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.3 0.9294-6*dv 0.4 0.0342],...
        'Callback',@TimeStepsCallback);

    % Plot Type 
    app.Options = {'' 'RDF' 'Potential Energy' 'Total Energy' 'Time' 'Kinetic Energy' ...
        'Temperature' 'Pressure' 'Volume' 'Density' 'Enthalpy' 'CRDF' 'Average Ql' ...
        'Average Wl' 'Ql Distribution' 'Wl Distribution' 'MSD' 'Nearest N.'};
    app.BGLabel(8) = uicontrol(app.bg,'Style','text','String',...
        'X-Axis','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0 0.9623-7*dv 0.5 0.0342]); % [left bottom width height]
    app.BGLabel(9) = uicontrol(app.bg,'Style','text','String',...
        'Y-Axis','FontSize',app.FS-8,'Units','Normalized','Position',...
        [0.5 0.9623-7*dv 0.5 0.0342]); % [left bottom width height]
    app.PlotTypeXh = uicontrol(app.bg,'Style','popupmenu','String',...
        app.Options,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0.9294-7*dv 0.48 0.0342],... % [left bottom width height]
        'Value',1,'Callback',@PlotTypeXCallback);
    app.PlotTypeYh = uicontrol(app.bg,'Style','popupmenu','String',...
        app.Options,'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.9294-7*dv 0.48 0.0342],... % [left bottom width height]
        'Value',1,'Callback',@PlotTypeYCallback);

    % Third axis
    app.ThirdAxis = uicontrol(app.bg,'Style','checkbox','String',...
        '3rd t-axis','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0.9294-7.48*dv 0.48 0.0342],... % [left bottom width height]
        'Value',0,'Callback',@ThirdAxisCallback,'Visible','off');
    app.ThirdAxisLabelA = uicontrol(app.bg,'Style','text','String',...
        'dT=','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.9294-7.48*dv 0.1472 0.0342],'Visible','off'); % [left bottom width height]
    app.ThirdAxisLabelB = uicontrol(app.bg,'Style','text','String',...
        '(ps)','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.8314 0.9294-7.48*dv 0.1786 0.0342],'Visible','off'); % [left bottom width height]
    app.ThirdAxisText = uicontrol(app.bg,'Style','edit','String',...
        '0','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.6572 0.9294-7.48*dv 0.1856 0.0342],... % [left bottom width height]
        'Value',0,'Callback',@ThirdAxisTextCallback,'Visible','off');
    app.TimeWindow = uicontrol(app.bg,'Style','text','String',...
        'Time Window (ps)','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 0.9294-8*dv 0.6686 0.0342],... % [left bottom width height]
        'Visible','off');
    app.TimeWindowText = uicontrol(app.bg,'Style','edit','String',...
        '','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.6875 0.933-8*dv 0.2992 0.026],... % [left bottom width height]
        'Value',0,'Callback',@TimeWindowCallback,'Visible','off');

    % Ql and Wl cluster selection options
    app.Qlbg = uibuttongroup(app.bg,'Position',[0 0.8994-9.2*dv 1.0 0.1155],... % [left bottom width height]
        'Visible','off','BorderType','none');
    app.QlRadius = uicontrol(app.Qlbg,'Style','radiobutton','String',...
        ['Search Rad. (' char(197) ')'],'FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 3/4 0.7 1/4],... % [left bottom width height]
        'Value',1,'Callback',@QlEditCallback);
    app.QlNeighbours = uicontrol(app.Qlbg,'Style','radiobutton','String',...
        'N. Neighbours','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 2/4 0.7 1/4],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallback);
    app.FirstShell = uicontrol(app.Qlbg,'Style','checkbox','String',...
        '1st Shell','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0 1/4 0.5 1/4],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallback);
    app.SecondShell = uicontrol(app.Qlbg,'Style','checkbox','String',...
        '2nd Shell','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.5 1/4 0.5 1/4],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallback);
    app.QlEditMin = uicontrol(app.Qlbg,'Style','edit','String',...
        '0','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.69 0.75 0.15 0.24],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallback);
    app.QlEdit = uicontrol(app.Qlbg,'Style','edit','String',...
        '3','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.84 0.75 0.15 0.24],... % [left bottom width height]
        'Callback',@QlEditCallback);
    app.NeighbEdit = uicontrol(app.Qlbg,'Style','edit','String',...
        '6','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.765 0.50 0.15 0.24],... % [left bottom width height]
        'Value',0,'Callback',@QlEditCallback);
    

    app.compareQlButton1 = uicontrol(app.Qlbg,'Style','pushbutton','String',...
        'Compare 1','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0 0.48 1/4],... % [left bottom width height]
        'Callback',@CompareQlCallback,'BackgroundColor',[0.9400    0.9400    0.9400]);
    app.compareQlButton2 = uicontrol(app.Qlbg,'Style','pushbutton','String',...
        'Compare 2','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0 0.48 1/4],... % [left bottom width height]
        'Callback',@CompareQlCallback,'BackgroundColor',[0.9400    0.9400    0.9400]);
    
    app.compareQlButton1b = uicontrol(app.figh,'Style','pushbutton','String',...
        'Compare 1','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.0014 0.2423 0.0670 0.0289],... % [left bottom width height]
        'Callback',@CompareQlCallback,'Visible','Off','BackgroundColor',[0.9400 0.9400 0.9400]);
    app.compareQlButton2b = uicontrol(app.figh,'Style','pushbutton','String',...
        'Compare 2','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.0712 0.2423 0.0670 0.0289],... % [left bottom width height]
        'Callback',@CompareQlCallback,'Visible','Off','BackgroundColor',[0.9400 0.9400 0.9400]);

    % Generate Plot Button
    app.GeneratePlot = uicontrol(app.bg,'Style','PushButton','String',...
        'Create Plot','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0.8994-10*dv 0.48 0.0342*1.5],... % [left bottom width height]
        'Callback',@GenPlotCallback);

    % Delete All Button
    app.ClearAxish = uicontrol(app.bg,'Style','PushButton','String',...
        'Delete All','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.8994-10*dv 0.48 0.0257],... % [left bottom width height]
        'Callback',@DeleteAllCallback);
    
    % Delete Current Button
    app.ClearAxish = uicontrol(app.bg,'Style','PushButton','String',...
        'Delete Plot','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.8994-9.64*dv 0.48 0.0257],... % [left bottom width height]
        'Callback',@DeletePlotCallback);

    % Save Data Button
    app.SaveData = uicontrol(app.bg,'Style','PushButton','String',...
        'Save Data','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0.9144-11*dv 0.48 0.0342*1.5],... % [left bottom width height]
        'Callback',@SaveDataCallback);

    % Load Data Button
    app.LoadData = uicontrol(app.bg,'Style','PushButton','String',...
        'Load Data','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.9144-11*dv 0.48 0.0342*1.5],... % [left bottom width height]
        'Callback',@LoadDataCallback);

    % View and Slice trajectory
    app.SliceTrajectoryh = uicontrol(app.bg,'Style','PushButton','String',...
        'View Traj.','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.01 0.9294-12*dv 0.48 0.0342*1.5],... % [left bottom width height]
        'Callback',@SliceTrajectoryCallback);
    
    % Load trajectory
    app.LoadTrajectoryh = uicontrol(app.bg,'Style','PushButton','String',...
        'Load Traj.','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.51 0.9294-12*dv 0.48 0.0342*1.5],... % [left bottom width height]
        'Callback',@LoadTrajectoryCallback);

    % View controls for Ql/Wl distributions
    app.TimeChangeText = uicontrol(app.bg,'Style','text','String',...
        'Time: 0 ns','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.0 0.9294-12.5*dv 0.6 0.030],... % [left bottom width height]
        'Visible','off');
    app.TimeChangeForward = uicontrol(app.bg,'Style','pushbutton','String',...
        '>>','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.8 0.9294-12.5*dv+0.0035 0.2 0.025],... % [left bottom width height]
        'Value',0,'Callback',@PlotChangeCallback,'Visible','off','tag','tforward');
    app.TimeChangeBackward = uicontrol(app.bg,'Style','pushbutton','String',...
        '<<','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.6 0.9294-12.5*dv+0.0035 0.2 0.025],... % [left bottom width height]
        'Value',0,'Callback',@PlotChangeCallback,'Visible','off','tag','tbackward');

    app.PlotChangeText = uicontrol(app.bg,'Style','text','String',...
        'Plot: 0 of 0','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.0 0.9294-13*dv 0.6 0.030],... % [left bottom width height]
        'Visible','on');
    app.PlotChangeForward = uicontrol(app.bg,'Style','pushbutton','String',...
        '>>','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.8 0.9294-13*dv+0.0035 0.2 0.025],... % [left bottom width height]
        'Value',0,'Callback',@PlotChangeCallback,'Visible','on',...
        'tag','pforward');
    app.PlotChangeBackward = uicontrol(app.bg,'Style','pushbutton','String',...
        '<<','FontSize',app.FS-8,'Units','Normalized',...
        'Position',[0.6 0.9294-13*dv+0.0035 0.2 0.025],... % [left bottom width height]
        'Value',0,'Callback',@PlotChangeCallback,'Visible','on',...
        'tag','pbackward');

    % Make invisible background axis
    app.InvisAxis = axes(app.figh,'Position',[0 0 1 1],'FontSize',app.FS,...
        'TickLabelInterpreter','latex','Box','off','XTick',[],'YTick',[],...
        'Color','none','Visible','off');

    %% Turn off some toolbar options
    set(0,'Showhidden','on')
    ch = findall(app.figh,'tag','FigureToolBar');
    UT = get(ch,'children');
    %set(UT(1:5),'Visible','Off','Separator','Off')
    set(UT(1),'Separator','On')
    set(UT(8),'Visible','Off','Separator','Off')
end

function GenPlotCallback(varargin) % Create Plots function
    %% Grab plot settings
    Data_Type = app.DataTypeh.Value; % {'MD' 'Cell Opt' 'Full Opt'}
    Salt = app.Salth.String{app.Salth.Value};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    Structure = app.Structureh.String{app.Structureh.Value}; % Only used for non-MD calculation
    if ~isempty(app.Modelh.String)
        Model = app.Modelh.String{app.Modelh.Value}; % {'JC (SPC/E)' 'JC (TIP4P-EW)' 'JC (TIP3P)' 'TF'};
    else
        warndlg('A model must be selected.')
        return
    end

    X_Axis = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    Y_Axis = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    app.Axh.ZLabel.String = '';
    
    Distr_Selx = ~isempty(regexp(X_Axis,'RDF','once')) || ~isempty(regexp(X_Axis,'[Q|W]l','once')) ...
        || ~isempty(regexp(X_Axis,'MSD','once')) || ~isempty(regexp(X_Axis,'Near','once'));
    RDF_sel = ~isempty(regexp(X_Axis,'RDF','once')) || ~isempty(regexp(Y_Axis,'RDF','once'));
    MSD_sel = ~isempty(regexp(X_Axis,'MSD','once')) || ~isempty(regexp(Y_Axis,'MSD','once'));
    NN_sel = ~isempty(regexp(X_Axis,'Near','once')) || ~isempty(regexp(Y_Axis,'Near','once'));
    Ql_sel = ~isempty(regexp(X_Axis,'Ql','once')) || ~isempty(regexp(Y_Axis,'Ql','once'));
    Wl_sel = ~isempty(regexp(X_Axis,'Wl','once')) || ~isempty(regexp(Y_Axis,'Wl','once'));
    
    if Data_Type == 1 % MD
        % Ensemble type
        Ensemble_c = regexp(Model,'_(.{3}$)','tokens','once');
        Ensemble = Ensemble_c{1};
        
        % Directory containing results of interest
        DatDir = [app.datdir filesep Salt filesep Model];
        
        if str2double(app.InitialTimeh.String) > str2double(app.FinalTimeh.String)
            warndlg('The Initial Time must be less than the final time.');
            return
        end
        
        Min_Time = str2double(app.InitialTimeh.String)*1000; % min time in ps
        Min_Time_Str = num2str(Min_Time); % min time in ps
        Max_Time = str2double(app.FinalTimeh.String)*1000; % max time in ps
        Max_Time_Str = num2str(Max_Time); % max time in ps

        if NN_sel % Nearest neighbour module selected
            
            TimePerFrame = app.TimePerFrameh.String; % Amount of time to skip between selected trajectory frames
            ThirdAxisSel = logical(app.ThirdAxis.Value); % is the third time-axis selected?
            ThirdAxisDT = str2double(app.ThirdAxisText.String); % Time step of third axis (ps)
            TimeWindow = str2double(app.TimeWindowText.String)/2; % Time window of third axis points
            
            First_Shell_Switch = logical(app.FirstShell.Value);
            Second_Shell_Switch = logical(app.SecondShell.Value);
            SelRadius = logical(app.QlRadius.Value); % Ql/Wl search space: radius or number of neighbours?
            
            if SelRadius
                leg_add = [' ' app.QlEditMin.String '-' app.QlEdit.String ' \AA'];
            elseif First_Shell_Switch && Second_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ \& 2$^{\textrm{nd}}$ Shell';
            elseif First_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ Shell';
            elseif Second_Shell_Switch
                leg_add = ' 2$^{\textrm{nd}}$ Shell';
            end
            
            if Distr_Selx
                Ref_Type = Y_Axis;
            else
                Ref_Type = X_Axis;
            end
            
            if ThirdAxisSel
                Time_Points = Min_Time:ThirdAxisDT:Max_Time;
                Start_Points = max(Time_Points - TimeWindow,0);
                End_Points = Time_Points + TimeWindow;
            else
                Start_Points = Min_Time;
                End_Points = Max_Time;
                Mid_Point = mean([Start_Points End_Points])/1000; % ns
            end
            
            % Find trajectory file
            Traj_info = dir([DatDir filesep '*.trr']);
            Tpr_info = dir([DatDir filesep '*.tpr']);
            Mdp_info = dir([DatDir filesep '*.mdp']);
            Gro_info = dir([DatDir filesep '*.gro']);
            
            if isempty(Traj_info) || isempty(Tpr_info) || isempty(Mdp_info) % no trajectory file found
                warndlg('Missing required file, one of: *.trr; *.tpr; or *.mdp')
                return
            end
            if length(Mdp_info) > 1
                Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
            end
            if length(Gro_info) > 1
                [~,gidx] = min(strlength({Gro_info.name})); % grab shortest name
                Gro_info = Gro_info(gidx);
            end
            
            Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
            Mdp_text = fileread(Mdp_file);
            
            % Check if PBC on
            pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));
            
            Traj_file = [DatDir filesep Traj_info.name]; % primary traj file
            Gro_file = [DatDir filesep Gro_info.name]; % gro molecule file
            
            % Calculate Ql with python script
            set(app.figh, 'pointer', 'watch')
            set(app.auxfig, 'pointer', 'watch')
            drawnow;
            try
                Results = py.LiXStructureDetector.Calculate_NN(Traj_file,Gro_file,Metal,Halide,...
                    Ref_Type,Start_Points,End_Points,str2double(TimePerFrame),First_Shell_Switch,Second_Shell_Switch,...
                    str2double(app.QlEditMin.String),str2double(app.QlEdit.String),pbc_on);
            catch e
                warndlg(['Error Executing Command:' newline e.message]);
                return
            end
            set(app.figh, 'pointer', 'arrow')
            set(app.auxfig, 'pointer', 'arrow')

            % pre allocate arrays
            N_SP = size(Results,2);
            NN_Hist_Data = cell(N_SP,1);
            
            NN_Edges = -0.5:1:20.5;
            NN_Centers = 0:1:20;
            for idx = 1:N_SP
                
                % Data from current time window
                Results_T = int64(Results{idx});
                
                % Calculate Ql histogram distribution for each value of l
                NN_Hist_Data{idx} = histcounts(Results_T(~isnan(Results_T)),NN_Edges)./(length(Results_T(~isnan(Results_T)))); %normalize count by number of reference atoms
            end
            
            if ThirdAxisSel
                
                y = Time_Points/1000; % ns
                [X_NN,Y_NN] = meshgrid(NN_Centers,y);

                Z_NN = cell2mat(NN_Hist_Data);
                
                [DX,DY] = size(Z_NN);
                NN_Data = zeros(3,DX,DY);
                NN_Data(1,:,:) = X_NN;
                NN_Data(2,:,:) = Y_NN;
                NN_Data(3,:,:) = Z_NN;
                
                app.NN{end+1} = NN_Data;
                app.Time{end+1} = y;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan;
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','-') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.NumPlots = app.NumPlots + 1;
                app.PlotPoint = app.NumPlots;
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 0;
                app.Notes{end+1} = 'Nearest Neighbours vs Time plot';
            else
                app.NN{end+1} = [NN_Centers' NN_Hist_Data{1}'];
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 5;
                if isnan(app.NN_type_plot) % First NN plot
                    app.NumPlots = app.NumPlots + 1;
                    app.PlotPoint = app.NumPlots;
                    app.NN_type_plot = app.PlotPoint;
                    app.Notes{end+1} = 'Nearest Neighbours Distribution';
                else
                    app.PlotPoint = app.NN_type_plot;
                end
            end

        elseif MSD_sel % MSD selected
            % Set up matlab parallel features
            if ~isempty(gcp('nocreate'))
                Cur_Pool = gcp;
                Cur_Workers = Cur_Pool.NumWorkers;
                if Cur_Workers ~= app.PrefCores
                    delete(Cur_Pool);
                    set(app.figh, 'pointer', 'watch')
                    set(app.auxfig, 'pointer', 'watch')
                    drawnow;
                    parpool(app.PrefCores);
                    set(app.figh, 'pointer', 'arrow')
                    set(app.auxfig, 'pointer', 'arrow')
                end
            else
                set(app.figh, 'pointer', 'watch')
                set(app.auxfig, 'pointer', 'watch')
                drawnow;
                parpool(app.PrefCores);
                set(app.figh, 'pointer', 'arrow')
                set(app.auxfig, 'pointer', 'arrow')
            end
            
            TimePerFrame = app.TimePerFrameh.String; % Amount of time to skip between selected trajectory frames
            ThirdAxisSel = logical(app.ThirdAxis.Value); % is the third time-axis selected?
            ThirdAxisDT = str2double(app.ThirdAxisText.String); % Time step of third axis (ps)
            TimeWindow = str2double(app.TimeWindowText.String)/2; % Time window of third axis points
            
            if Distr_Selx
                Ref_Type = Y_Axis;
            else
                Ref_Type = X_Axis;
            end
            
            switch Ref_Type
                case Metal
                    ref = ['echo 2' app.pipe];
                case Halide
                    ref = ['echo 3' app.pipe];
                case 'All'
                    ref = ['echo 0' app.pipe];
            end
            
            if ThirdAxisSel
                Time_Points = Min_Time:ThirdAxisDT:Max_Time;
                Start_Points = max(Time_Points - TimeWindow,0);
                End_Points = Time_Points + TimeWindow;
            else
                Start_Points = Min_Time;
                End_Points = Max_Time;
                Mid_Point = mean([Start_Points End_Points])/1000; % ns
            end
        
            % Find trajectory file
            Traj_info = dir([DatDir filesep '*.trr']);
            Tpr_info = dir([DatDir filesep '*.tpr']);
            Mdp_info = dir([DatDir filesep '*.mdp']);
            
            if isempty(Traj_info) || isempty(Tpr_info) || isempty(Mdp_info) % no trajectory file found
                warndlg('Missing required file, one of: *.trr; *.tpr; or *.mdp')
                return
            end
            if length(Mdp_info) > 1
                Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
            end
            Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
            Mdp_text = fileread(Mdp_file);
            
            % Check if PBC on
            pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));
            if pbc_on
                % Is a cluster
                if ~isempty(regexp(Model,'_N[0-9]+','once'))
                    pbc_mod = ' -pbc nojump'; % pbc % needed eg FavSphal_N500
                    ech = ['echo "0"' app.pipe];
                % Is not a cluster
                else
                    pbc_mod = ''; % pbc
                    ech = ['echo "0"' app.pipe];
                end
            % PBC off
            else
                pbc_mod = ' -box 5 5 5'; % place into a box
                ech = ['echo "0"' app.pipe];
            end
            
            Traj_file = [DatDir filesep Traj_info.name]; % primary traj file
            Tpr_file = [DatDir filesep Tpr_info.name]; % traj input file
            
            N_SP = length(Start_Points);
            
            % Create a set of parallel futures and a waitbar
            f(1:N_SP) = parallel.FevalFuture;
            wb = PoolWaitbar(N_SP,'Calculating MSD(s)');
            for idx = 1:N_SP
            
                Start_Point = num2str(Start_Points(idx)); % ps
                End_Point = num2str(End_Points(idx)); % ps
                
                MSD_Out_file = [app.tempfolder filesep 'tempRDF-' Start_Point '.xvg'];
                Trajectory_Out_file = [app.tempfolder filesep 'tempTraj-' Start_Point '.xtc'];
                
                f(idx) = parfeval(@calc_MSD,3,wb,app.wsl,ech,app.gmx,Traj_file,...
                    Trajectory_Out_file,Tpr_file,Start_Point,End_Point,...
                    TimePerFrame,pbc_mod,pbc_on,MSD_Out_file,ref);
                
%                 [MSD,DiffCoeff,DiffCoeff_e] = calc_MSD(wb,app.wsl,ech,app.gmx,Traj_file,Trajectory_Out_file,Tpr_file,Start_Point,End_Point,...
%                     TimePerFrame,pbc_mod,pbc_on,MSD_Out_file,ref);
                
            end
            wait(f)
            delete(wb)
            
            % Collect data
            [MSD_Data,DiffCoeff,DiffCoeff_e] = fetchOutputs(f,'UniformOutput',false);            
            
            if ThirdAxisSel
                
                xsz = min(cellfun(@(x) size(x,1), MSD_Data));
                
                x_MSD = MSD_Data{1}(1:xsz,1); % nm
                
                y = Time_Points/1000; % ns
                [X_MSD,Y_MSD] = meshgrid(x_MSD,y);
                
                Z_MSD = nan(size(X_MSD)); % compose Z shape

                for id = 1:N_SP % Add in data
                    Z_MSD(id,:) = MSD_Data{id}(1:xsz,2)';
                end
                
                [DX,DY] = size(Z_MSD);
                MSD_Data = zeros(3,DX,DY);
                
                MSD_Data(1,:,:) = X_MSD;
                MSD_Data(2,:,:) = Y_MSD;
                MSD_Data(3,:,:) = Z_MSD;
                
                app.Time{end+1} = y;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = MSD_Data;
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','-') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]'];
                
                % Increment stuff
                app.NumPlots = app.NumPlots + 1;
                app.PlotPoint = app.NumPlots;
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 0;
                app.Notes{end+1} = 'MSD vs Time plot';
            else
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = MSD_Data{1};
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps] D = ' ...
                    num2str(DiffCoeff{1}) ' +/-' num2str(DiffCoeff_e{1}) ' nm$^{2}$ ns$^{-1}$'];
                
                % Increment stuff
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 3;
                
                if isnan(app.msd_type_plot) % First MSD plot
                    app.NumPlots = app.NumPlots + 1;
                    app.PlotPoint = app.NumPlots;
                    app.msd_type_plot = app.PlotPoint;
                    app.Notes{end+1} = 'MSD plot';
                else
                    app.PlotPoint = app.msd_type_plot;
                end
            end
            
        elseif RDF_sel % RDF/CRDF selected

            TimePerFrame = app.TimePerFrameh.String; % Amount of time to skip between selected trajectory frames
            ThirdAxisSel = logical(app.ThirdAxis.Value); % is the third time-axis selected?
            ThirdAxisDT = str2double(app.ThirdAxisText.String); % Time step of third axis (ps)
            TimeWindow = str2double(app.TimeWindowText.String)/2; % Time window of third axis points
        
            if Distr_Selx
                Ref_Type = Y_Axis;
            else
                Ref_Type = X_Axis;
            end
        
            if ThirdAxisSel
                Time_Points = Min_Time:ThirdAxisDT:Max_Time;
                Start_Points = max(Time_Points - TimeWindow,0);
                End_Points = Time_Points + TimeWindow;
            else
                Start_Points = Min_Time;
                End_Points = Max_Time;
                Mid_Point = mean([Start_Points End_Points])/1000; % ns
            end
        
            % Find trajectory file
            Traj_info = dir([DatDir filesep '*.trr']);
            Gro_info = dir([DatDir filesep '*.tpr']);
            Mdp_info = dir([DatDir filesep '*.mdp']);
            
            if isempty(Traj_info) || isempty(Gro_info) || isempty(Mdp_info) % no trajectory file found
                warndlg('Missing required file, one of: *.trr; *.tpr; or *.mdp')
                return
            end
            if length(Mdp_info) > 1
                Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
            end
            if length(Gro_info) > 1
                [~,gidx] = min(strlength({Gro_info.name})); % grab shortest name
                Gro_info = Gro_info(gidx);
            end
            
            Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
            Mdp_text = fileread(Mdp_file);
            
            % Check if PBC on
            pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));
            Traj_file = [DatDir filesep Traj_info.name]; % primary traj file
            Gro_file = [DatDir filesep Gro_info.name]; % gro molecule file
            
            % Calculate RDF with python script (using Freud)
            set(app.figh, 'pointer', 'watch')
            set(app.auxfig, 'pointer', 'watch')
            drawnow;
            try
                Results = py.LiXStructureDetector.Calculate_RDF(Traj_file,Gro_file,...
                    Metal,Halide,Ref_Type,Start_Points,End_Points,...
                    str2double(TimePerFrame),pbc_on);
            catch e
                warndlg(['Error Executing Command:' newline e.message]);
                return
            end
            RDF_Data = double(Results{1});
            CRDF_Data = double(Results{2});
            set(app.figh, 'pointer', 'arrow')
            set(app.auxfig, 'pointer', 'arrow')
            
            if ThirdAxisSel
                
                x_CRDF = squeeze(CRDF_Data(1,1,:)); % nm
                x_RDF = squeeze(RDF_Data(1,1,:)); % nm
                
                y = Time_Points/1000; % ns
                [X_CRDF,Y_CRDF] = meshgrid(x_CRDF,y);
                [X_RDF,Y_RDF] = meshgrid(x_RDF,y);
                
                Z_CRDF = squeeze(CRDF_Data(:,2,:)); % compose Z shape
                Z_RDF = squeeze(RDF_Data(:,2,:)); % compose Z shape
                
                [DX,DY] = size(Z_RDF);
                RDF_Data = zeros(3,DX,DY);
                [DX,DY] = size(Z_CRDF);
                CRDF_Data = zeros(3,DX,DY);
                
                RDF_Data(1,:,:) = X_RDF;
                RDF_Data(2,:,:) = Y_RDF;
                RDF_Data(3,:,:) = Z_RDF;
                CRDF_Data(1,:,:) = X_CRDF;
                CRDF_Data(2,:,:) = Y_CRDF;
                CRDF_Data(3,:,:) = Z_CRDF;
                
                app.Time{end+1} = y;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan;
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = RDF_Data;
                app.CRDF{end+1} = CRDF_Data;
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','-') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]'];
                
                % Increment stuff
                app.NumPlots = app.NumPlots + 1;
                app.PlotPoint = app.NumPlots;
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 0;
                app.Notes{end+1} = 'RDF vs Time plot';
            else
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = squeeze(RDF_Data)';
                app.CRDF{end+1} = squeeze(CRDF_Data)';
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]'];
                
                % Increment stuff
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 2;
                if isnan(app.rdf_type_plot) % First RDF plot
                    app.NumPlots = app.NumPlots + 1;
                    app.PlotPoint = app.NumPlots;
                    app.rdf_type_plot = app.PlotPoint;
                    app.Notes{end+1} = 'RDF plot';
                else
                    app.PlotPoint = app.rdf_type_plot;
                end
            end
            
        elseif Wl_sel

            TimePerFrame = app.TimePerFrameh.String; % Amount of time to skip between selected trajectory frames
            ThirdAxisSel = logical(app.ThirdAxis.Value); % is the third time-axis selected?
            ThirdAxisDT = str2double(app.ThirdAxisText.String); % Time step of third axis (ps)
            TimeWindow = str2double(app.TimeWindowText.String)/2; % Time window of third axis points
            
            First_Shell_Switch = logical(app.FirstShell.Value);
            Second_Shell_Switch = logical(app.SecondShell.Value);
            SelRadius = logical(app.QlRadius.Value); % Ql/Wl search space: radius or number of neighbours?
            
            if SelRadius
                leg_add = [' ' app.QlEditMin.String '-' app.QlEdit.String ' \AA'];
            elseif First_Shell_Switch && Second_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ \& 2$^{\textrm{nd}}$ Shell';
            elseif First_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ Shell';
            elseif Second_Shell_Switch
                leg_add = ' 2$^{\textrm{nd}}$ Shell';
            else
                leg_add = [' ' app.QlEdit.String ' Near N.'];
            end
            
            if Distr_Selx
                Ref_Type = Y_Axis;
                Calc_type = X_Axis;
            else
                Ref_Type = X_Axis;
                Calc_type = Y_Axis;
            end
            
            if ThirdAxisSel
                Time_Points = Min_Time:ThirdAxisDT:Max_Time;
                Start_Points = max(Time_Points - TimeWindow,0);
                End_Points = Time_Points + TimeWindow;
            else
                Start_Points = Min_Time;
                End_Points = Max_Time;
                Mid_Point = mean([Start_Points End_Points])/1000; % ns
            end
            
            % Find trajectory file
            Traj_info = dir([DatDir filesep '*.trr']);
            Tpr_info = dir([DatDir filesep '*.tpr']);
            Mdp_info = dir([DatDir filesep '*.mdp']);
            Gro_info = dir([DatDir filesep '*.gro']);
            
            if isempty(Traj_info) || isempty(Tpr_info) || isempty(Mdp_info) % no trajectory file found
                warndlg('Missing required file, one of: *.trr; *.tpr; or *.mdp')
                return
            end
            if length(Mdp_info) > 1
                Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
            end
            if length(Gro_info) > 1
                [~,gidx] = min(strlength({Gro_info.name})); % grab shortest name
                Gro_info = Gro_info(gidx);
            end
            
            Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
            Mdp_text = fileread(Mdp_file);
            
            % Check if PBC on
            pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));
            
            Traj_file = [DatDir filesep Traj_info.name]; % primary traj file
            Gro_file = [DatDir filesep Gro_info.name]; % gro molecule file
            
            % Calculate Wl and Ql with python script
            set(app.figh, 'pointer', 'watch')
            set(app.auxfig, 'pointer', 'watch')
            drawnow;
            try
                Results = py.LiXStructureDetector.Calculate_Wl(Traj_file,Gro_file,...
                    Metal,Halide,Ref_Type,Start_Points,End_Points,...
                    str2double(TimePerFrame),First_Shell_Switch,Second_Shell_Switch,...
                    SelRadius,str2double(app.NeighbEdit.String),...
                    str2double(app.QlEditMin.String),str2double(app.QlEdit.String),pbc_on,app.l);
            catch e
                warndlg(['Error Executing Command:' newline e.message]);
                return
            end
            ResultsWl = Results{1};
            ResultsQl = Results{2};
            set(app.figh, 'pointer', 'arrow')
            set(app.auxfig, 'pointer', 'arrow')

            % pre allocate arrays
            N_SP = size(ResultsWl,2);
            Wl_Ave_Data = cell(N_SP,1);
            Wl_Hist_Data = cell(N_SP,1);
            Ql_Ave_Data = cell(N_SP,1);
            Ql_Hist_Data = cell(N_SP,1);
            
            Ql_Edges = 0:app.HistBinQ:1;
            Ql_Centers = Ql_Edges(1:end-1)+app.HistBinQ/2;
            Wl_Edges = -0.5:app.HistBinW:0.5;
            Wl_Centers = Wl_Edges(1:end-1)+app.HistBinW/2;
            for idx = 1:N_SP
                
                % Data from current time window
                ResultsWl_T = double(ResultsWl{idx});
                ResultsQl_T = double(ResultsQl{idx});
                
                % Average values of Wl and Ql
                Wl_Ave_Data{idx} = mean(ResultsWl_T,2,'omitnan');
                Ql_Ave_Data{idx} = mean(ResultsQl_T,2,'omitnan');
                
                % Calculate Ql histogram distribution for each value of l
                Wl_Hist = zeros(length(app.l),length(Wl_Centers));
                Ql_Hist = zeros(length(app.l),length(Ql_Centers));
                for lidx = 1:length(app.l)
                    Wl_Hist(lidx,:) = histcounts(ResultsWl_T(lidx,:),Wl_Edges)./(app.HistBinW*length(ResultsWl_T(1,:))); %normalize count by number of reference atoms
                    Ql_Hist(lidx,:) = histcounts(ResultsQl_T(lidx,:),Ql_Edges)./(app.HistBinQ*length(ResultsQl_T(1,:))); %normalize count by number of reference atoms
                end
                
                Wl_Hist_Data{idx} = Wl_Hist;
                Ql_Hist_Data{idx} = Ql_Hist;
            end
            
            [XQlD,YQlD] = meshgrid(Ql_Centers,app.l);
            [XWlD,YWlD] = meshgrid(Wl_Centers,app.l);
            
            if ThirdAxisSel
                
                x_ql = app.l;
                y = Time_Points/1000; % ns
                [X_QL,Y_QL] = meshgrid(x_ql,y);

                Z_QL = nan(size(X_QL));
                Z_WL = nan(size(X_QL));

                for id = 1:length(Time_Points) % Add in data
                    Z_QL(id,:) = Ql_Ave_Data{id};
                    Z_WL(id,:) = Wl_Ave_Data{id};
                end

                [DXQl,DYQl] = size(Z_QL);
                Ql_Ave_Data = zeros(3,DXQl,DYQl);
                [DXWl,DYWl] = size(Z_WL);
                Wl_Ave_Data = zeros(3,DXWl,DYWl);

                Ql_Ave_Data(1,:,:) = X_QL;
                Ql_Ave_Data(2,:,:) = Y_QL;
                Ql_Ave_Data(3,:,:) = Z_QL;
                Wl_Ave_Data(1,:,:) = X_QL;
                Wl_Ave_Data(2,:,:) = Y_QL;
                Wl_Ave_Data(3,:,:) = Z_WL;

                [DXQl,DYQl] = size(Ql_Hist_Data{1});
                [DXWl,DYWl] = size(Wl_Hist_Data{1});
                Ql_Hist_Out = cell(N_SP,2);
                Wl_Hist_Out = cell(N_SP,2);
                for t = 1:N_SP
                    Ql_Hist_Out{t,2} = y(t);
                    Wl_Hist_Out{t,2} = y(t);

                    Ql_Hist_Out{t,1} = zeros(DXQl,DYQl,3);
                    Ql_Hist_Out{t,1}(:,:,1) = XQlD;
                    Ql_Hist_Out{t,1}(:,:,2) = YQlD;
                    Ql_Hist_Out{t,1}(:,:,3) = Ql_Hist_Data{t};

                    Wl_Hist_Out{t,1} = zeros(DXWl,DYWl,3);
                    Wl_Hist_Out{t,1}(:,:,1) = XWlD;
                    Wl_Hist_Out{t,1}(:,:,2) = YWlD;
                    Wl_Hist_Out{t,1}(:,:,3) = Wl_Hist_Data{t};
                end
                
                app.Ql{end+1} = Ql_Ave_Data;
                app.Wl{end+1} = Wl_Ave_Data;
                app.Ql_Dist{end+1} = Ql_Hist_Out;
                app.Wl_Dist{end+1} = Wl_Hist_Out;
                                
                app.Time{end+1} = y;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','-') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.PlotPoint   = app.PlotPoint   + 1;
                app.NumDataSets = app.NumDataSets + 1;
                app.NumPlots    = app.NumPlots    + 1;
                app.data_type(end+1) = 0;
                app.Notes{end+1} = 'Ql/Wl Average/Distribution vs Time plot';
                
            else
                Ql_Hist_Out = cell(1,2);
                Wl_Hist_Out = cell(1,2);
                Ql_Hist_Out{1,2} = Mid_Point;
                Wl_Hist_Out{1,2} = Mid_Point;

                [DXQl,DYQl] = size(Ql_Hist_Data{1});
                Ql_Hist_Out{1,1} = zeros(DXQl,DYQl,3);
                Ql_Hist_Out{1,1}(:,:,1) = XQlD;
                Ql_Hist_Out{1,1}(:,:,2) = YQlD;
                Ql_Hist_Out{1,1}(:,:,3) = Ql_Hist_Data{1};

                [DXQl,DYQl] = size(Wl_Hist_Data{1});
                Wl_Hist_Out{1,1} = zeros(DXQl,DYQl,3);
                Wl_Hist_Out{1,1}(:,:,1) = XWlD;
                Wl_Hist_Out{1,1}(:,:,2) = YWlD;
                Wl_Hist_Out{1,1}(:,:,3) = Wl_Hist_Data{1};

                app.Ql{end+1} = [app.l' Ql_Ave_Data{1}];
                app.Wl{end+1} = [app.l' Wl_Ave_Data{1}];
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan;
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                app.Ql{end+1} = nan;
                app.Wl{end+1} = nan;
                app.Ql_Dist{end+1} = Ql_Hist_Out;
                app.Wl_Dist{end+1} = Wl_Hist_Out;
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan;
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.NumDataSets = app.NumDataSets + 2;
                app.data_type(end+1) = 4;
                app.data_type(end+1) = 0;
                
                switch Calc_type
                    case 'Average Wl'
                        if isnan(app.AvQl_type_plot) % First Ave Wl plot
                            app.NumPlots  = app.NumPlots + 2;
                            app.PlotPoint = app.NumPlots - 1;
                            app.AvQl_type_plot = app.PlotPoint;
                            app.Notes{end+1} = 'Ql/Wl Average plots';
                            app.Notes{end+1} = 'Ql/Wl Distribution plot';
                        else
                            app.NumPlots  = app.NumPlots + 1;
                            app.PlotPoint = app.AvQl_type_plot;
                            app.Notes{end+1} = 'Ql/Wl Distribution plot';
                        end
                    case 'Wl Distribution'
                        if isnan(app.AvQl_type_plot) % First Ave Wl plot
                            app.NumPlots  = app.NumPlots + 2;
                            app.PlotPoint = app.NumPlots;
                            app.AvQl_type_plot = app.PlotPoint - 1;
                            app.Notes{end+1} = 'Ql/Wl Average plots';
                            app.Notes{end+1} = 'Ql/Wl Distribution plot';
                        else
                            app.NumPlots  = app.NumPlots + 1;
                            app.PlotPoint = app.NumPlots;
                            app.Notes{end+1} = 'Ql/Wl Distribution plot';
                        end
                end
            end
            
        elseif Ql_sel % Ql selected
            
            TimePerFrame = app.TimePerFrameh.String; % Amount of time to skip between selected trajectory frames
            ThirdAxisSel = logical(app.ThirdAxis.Value); % is the third time-axis selected?
            ThirdAxisDT = str2double(app.ThirdAxisText.String); % Time step of third axis (ps)
            TimeWindow = str2double(app.TimeWindowText.String)/2; % Time window of third axis points
            
            First_Shell_Switch = logical(app.FirstShell.Value);
            Second_Shell_Switch = logical(app.SecondShell.Value);
            SelRadius = logical(app.QlRadius.Value); % Ql/Wl search space: radius or number of neighbours?
            
            if SelRadius
                leg_add = [' ' app.QlEditMin.String '-' app.QlEdit.String ' \AA'];
            elseif First_Shell_Switch && Second_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ \& 2$^{\textrm{nd}}$ Shell';
            elseif First_Shell_Switch
                leg_add = ' 1$^{\textrm{st}}$ Shell';
            elseif Second_Shell_Switch
                leg_add = ' 2$^{\textrm{nd}}$ Shell';
            else
                leg_add = [' ' app.NeighbEdit.String ' Near N.'];
            end
            
            if Distr_Selx
                Ref_Type = Y_Axis;
                Calc_type = X_Axis;
            else
                Ref_Type = X_Axis;
                Calc_type = Y_Axis;
            end
            
            if ThirdAxisSel
                Time_Points = Min_Time:ThirdAxisDT:Max_Time;
                Start_Points = max(Time_Points - TimeWindow,0);
                End_Points = Time_Points + TimeWindow;
            else
                Start_Points = Min_Time;
                End_Points = Max_Time;
                Mid_Point = mean([Start_Points End_Points])/1000; % ns
            end
            
            % Find trajectory file
            Traj_info = dir([DatDir filesep '*.trr']);
            Tpr_info = dir([DatDir filesep '*.tpr']);
            Mdp_info = dir([DatDir filesep '*.mdp']);
            Gro_info = dir([DatDir filesep '*.gro']);
            
            if isempty(Traj_info) || isempty(Tpr_info) || isempty(Mdp_info) % no trajectory file found
                warndlg('Missing required file, one of: *.trr; *.tpr; or *.mdp')
                return
            end
            if length(Mdp_info) > 1
                Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
            end
            if length(Gro_info) > 1
                [~,gidx] = min(strlength({Gro_info.name})); % grab shortest name
                Gro_info = Gro_info(gidx);
            end
            
            Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
            Mdp_text = fileread(Mdp_file);
            
            % Check if PBC on
            pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));
            
            Traj_file = [DatDir filesep Traj_info.name]; % primary traj file
            Gro_file = [DatDir filesep Gro_info.name]; % gro molecule file
            
            % Calculate Ql with python script
            set(app.figh, 'pointer', 'watch')
            set(app.auxfig, 'pointer', 'watch')
            drawnow;
            try
                Results = py.LiXStructureDetector.Calculate_Ql(Traj_file,Gro_file,...
                    Metal,Halide,Ref_Type,Start_Points,End_Points,...
                    str2double(TimePerFrame),First_Shell_Switch,Second_Shell_Switch,...
                    SelRadius,str2double(app.NeighbEdit.String),...
                    str2double(app.QlEditMin.String),str2double(app.QlEdit.String),pbc_on,app.l);
            catch e
                warndlg(['Error Executing Command:' newline e.message]);
                return
            end
            set(app.figh, 'pointer', 'arrow')
            set(app.auxfig, 'pointer', 'arrow')

            % pre allocate arrays
            N_SP = size(Results,2);
            Ql_Ave_Data = cell(N_SP,1);
            Ql_Hist_Data = cell(N_SP,1);
            
            Ql_Edges = 0:app.HistBinQ:1;
            Ql_Centers = Ql_Edges(1:end-1)+app.HistBinQ/2;
            for idx = 1:N_SP
                
                % Data from current time window
                Results_T = double(Results{idx});
                
                % Average values of Ql
                Ql_Ave_Data{idx} = mean(Results_T,2,'omitnan');
                
                % Calculate Ql histogram distribution for each value of l
                Ql_Hist = zeros(length(app.l),length(Ql_Centers));
                for lidx = 1:length(app.l)
                    Results_Tl = Results_T(lidx,:);
                    Ql_Hist(lidx,:) = histcounts(Results_Tl(~isnan(Results_Tl)),Ql_Edges)./(app.HistBinQ.*length(Results_Tl(~isnan(Results_Tl)))); %normalize count by number of reference atoms
                end
                
                Ql_Hist_Data{idx} = Ql_Hist;
            end
            
            [XQlD,YQlD] = meshgrid(Ql_Centers,app.l);
            
            if ThirdAxisSel
                
                x_ql = app.l;
                y = Time_Points/1000; % ns
                [X_QL,Y_QL] = meshgrid(x_ql,y);

                Z_QL = nan(size(X_QL));

                for id = 1:length(Time_Points) % Add in data
                    Z_QL(id,:) = Ql_Ave_Data{id};
                end

                [DXQl,DYQl] = size(Z_QL);
                Ql_Ave_Data = zeros(3,DXQl,DYQl);

                Ql_Ave_Data(1,:,:) = X_QL;
                Ql_Ave_Data(2,:,:) = Y_QL;
                Ql_Ave_Data(3,:,:) = Z_QL;

                [DXQl,DYQl] = size(Ql_Hist_Data{1});
                Ql_Hist_Out = cell(N_SP,2);
                for t = 1:N_SP
                    Ql_Hist_Out{t,2} = y(t);

                    Ql_Hist_Out{t,1} = zeros(DXQl,DYQl,3);
                    Ql_Hist_Out{t,1}(:,:,1) = XQlD;
                    Ql_Hist_Out{t,1}(:,:,2) = YQlD;
                    Ql_Hist_Out{t,1}(:,:,3) = Ql_Hist_Data{t};
                end

                app.Ql{end+1} = Ql_Ave_Data;
                app.Ql_Dist{end+1} = Ql_Hist_Out;
                
                app.Wl{end+1} = nan(1,2);
                app.Wl_Dist{end+1} = nan;            
                app.Time{end+1} = y;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','-') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.NumPlots = app.NumPlots + 1;
                app.PlotPoint = app.NumPlots;
                app.NumDataSets = app.NumDataSets + 1;
                app.data_type(end+1) = 0;
                app.Notes{end+1} = 'Ql Average vs Time plots';
            else
                Ql_Hist_Out = cell(1,2);
                Ql_Hist_Out{1,2} = Mid_Point;

                [DXQl,DYQl] = size(Ql_Hist_Data{1});
                Ql_Hist_Out{1,1} = zeros(DXQl,DYQl,3);
                Ql_Hist_Out{1,1}(:,:,1) = XQlD;
                Ql_Hist_Out{1,1}(:,:,2) = YQlD;
                Ql_Hist_Out{1,1}(:,:,3) = Ql_Hist_Data{1};


                app.Ql{end+1} = [app.l' Ql_Ave_Data{1}];
                app.Ql_Dist{end+1} = nan;
                app.Wl{end+1} = nan(1,2);
                app.Wl_Dist{end+1} = nan;
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                app.Ql{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = Ql_Hist_Out;
                app.Wl{end+1} = nan(1,2);
                app.Wl_Dist{end+1} = nan;
                app.Time{end+1} = Mid_Point;
                app.PotE{end+1} = nan;
                app.KinE{end+1} = nan;
                app.TotE{end+1} = nan;
                app.Temp{end+1} = nan;
                app.Press{end+1} = nan;
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' ' Ref_Type ...
                    ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]' leg_add];
                
                % Increment stuff
                app.NumDataSets = app.NumDataSets + 2;
                app.data_type(end+1) = 4;
                app.data_type(end+1) = 0;
                
                switch Calc_type
                    case 'Average Ql'
                        if isnan(app.AvQl_type_plot) % First Ave Ql plot
                            app.NumPlots  = app.NumPlots + 2;
                            app.PlotPoint = app.NumPlots - 1;
                            app.AvQl_type_plot = app.PlotPoint;
                            app.Notes{end+1} = 'Ql Average plots';
                            app.Notes{end+1} = 'Ql Distribution plot';
                        else
                            app.NumPlots  = app.NumPlots + 1;
                            app.PlotPoint = app.AvQl_type_plot;
                            app.Notes{end+1} = 'Ql Distribution plot';
                        end
                    case 'Ql Distribution'
                        if isnan(app.AvQl_type_plot) % First Ave Ql plot
                            app.NumPlots  = app.NumPlots + 2;
                            app.PlotPoint = app.NumPlots;
                            app.AvQl_type_plot = app.NumPlots - 1;
                            app.Notes{end+1} = 'Ql Average plots';
                            app.Notes{end+1} = 'Ql Distribution plot';
                        else
                            app.NumPlots  = app.NumPlots + 1;
                            app.PlotPoint = app.NumPlots;
                            app.Notes{end+1} = 'Ql Distribution plot';
                        end
                end
            end
        else % Something other than MSD/RDF/Ql/Wl selected
            
            % Find energy file
            Ener_info = dir([DatDir filesep '*.edr']);
            Tpr_info = dir([DatDir filesep '*.tpr']);
            %Mat_info = dir([DatDir filesep Model '*.mat']);
            Gro_info = dir([DatDir filesep Model '.gro']);
            
            if isempty(Ener_info) || isempty(Tpr_info) || isempty(Gro_info) % missing input files
                warning('Missing input file(s).')
                return
            end
            
            Ener_file = [DatDir filesep Ener_info.name]; % traj file
            Tpr_file = [DatDir filesep Tpr_info.name]; % tpr file
            %Mat_file = [DatDir filesep Mat_info.name]; % mat file
            Gro_file = [DatDir filesep Gro_info.name]; % gro file
            Out_file = [app.tempfolder filesep 'tempEnergy.xvg'];
            
            % Check energy options
            gmx_command = [app.wsl 'echo 0' app.pipe app.gmx ' energy -f ' windows2unix(Ener_file)...
                ' -s ' windows2unix(Tpr_file)];
            [err,outpt] = system(gmx_command);
            if err ~= 1
                warndlg('Problem with energy check.')
                return
            end
            en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
            en_opts = en_opts{1};
            En_set = char(regexp(en_opts,'([0-9]{1,2})  Potential','tokens','once'));
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Kinetic-En.','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Total-Energy','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Temperature','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Pressure','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Density','tokens','once'))];
            En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Enthalpy','tokens','once'))];
            En_set = [En_set ' 0'];
            En_set = regexprep(En_set,' +',' ');
    
            % Grab data from results
            gmx_command = [app.wsl 'echo ' En_set app.pipe app.gmx ' energy -f ' windows2unix(Ener_file)...
                ' -o ' windows2unix(Out_file) ' -s ' windows2unix(Tpr_file) ...
                ' -b ' Min_Time_Str ' -e ' Max_Time_Str];
            
            set(app.figh, 'pointer', 'watch')
            set(app.auxfig, 'pointer', 'watch')
            drawnow;
            err = system(gmx_command,'-echo');
            set(app.figh, 'pointer', 'arrow')
            set(app.auxfig, 'pointer', 'arrow')
            if err ~= 0
                warndlg('Failed to collect data.')
                return
            end
            
            NSS = load_gro_file(Gro_file);
            NF = NSS.N_atoms/2; % number of ion pairs in supercell
            
            Data = import_xvg(Out_file); % Gather RDF
            delete(Out_file) % remove temp output file
            
            if strcmpi(Ensemble,'NVT')
                app.Time{end+1} = Data(:,1)./1000; % ns
                app.PotE{end+1} = Data(:,2)./NF;
                app.KinE{end+1} = Data(:,3)./NF;
                app.TotE{end+1} = Data(:,4)./NF;
                app.Temp{end+1} = Data(:,5);
                app.Press{end+1} = Data(:,6); % Bar
                app.Vol{end+1} = nan;
                app.Dens{end+1} = nan;
                app.Enth{end+1} = nan;
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]'];
            elseif strcmpi(Ensemble,'NPT')
                app.Time{end+1} = Data(:,1)./1000; % ns
                app.PotE{end+1} = Data(:,2)./NF;
                app.KinE{end+1} = Data(:,3)./NF;
                app.TotE{end+1} = Data(:,4)./NF;
                app.Temp{end+1} = Data(:,5);
                app.Press{end+1} = Data(:,6); % Bar
                app.Vol{end+1} = Data(:,7).*(1e-7^3).*(6.0221409e+23)./NF; % cm^3/mol
                app.Dens{end+1} = Data(:,8); % kg^3/m^3
                app.Enth{end+1} = Data(:,9)./NF; % kJ/mol
                app.MSD{end+1} = nan(1,2);
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.NN{end+1} = nan(1,2);
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                app.Ql_Dist{end+1} = nan;
                app.Wl_Dist{end+1} = nan;
                app.legtxt{end+1} = [Salt ' ' strrep(Model,'_','\_') ' [' num2str(Min_Time_Str) '-' num2str(Max_Time_Str) ' ps]'];
            end
            
            % Increment stuff
            app.NumDataSets = app.NumDataSets + 1;
            app.data_type(end+1) = 1;

            if isnan(app.energy_type_plot) % First energy type plot
                app.NumPlots = app.NumPlots + 1;
                app.PlotPoint = app.NumPlots;
                app.energy_type_plot = app.PlotPoint;
                app.Notes{end+1} = 'Energy vs Time plots';
            else
                app.PlotPoint = app.energy_type_plot;
            end
        end
        
    else % Cell/Full opt RDF/Wl/Ql
        
        Model = strrep(Model,' (SPC/E)','');
        Model = strrep(Model,' (TIP3P)','3P');
        Model = strrep(Model,' (TIP4P-EW)','4P');
        
        % Load a basic Gromacs Object and update its parameters
        Input = GromObject;
        Input.Structure = Structure;
        Input.Salt = Salt;
        Input.Model = Model;
        Input.Data_Type = Data_Type-1;
        
        [RDF,CRDF] = GenRDF(app,Input,0.02,15);
        
        NeighbourSelected = logical(app.QlNeighbours.Value);
        FirstShell = logical(app.FirstShell.Value);
        SecondShell = logical(app.SecondShell.Value);
        MX_Peak = RDF.R(RDF.MX ~= 0);
        XX_Peak = RDF.R(RDF.XX ~= 0);
        
        % Collapse into shells (tolerance = 0.25 Ang)
        MX_Shells = uniquetol(MX_Peak,0.25,'DataScale',1);
        XX_Shells = uniquetol(XX_Peak,0.25,'DataScale',1);
        Minima_X = [mean([MX_Shells(1) XX_Shells(1)]) mean([XX_Shells(1) MX_Shells(2)])];
        
        if FirstShell && SecondShell
            lims = [0 Minima_X(2)];
            leg_add = [' 0-' num2str(Minima_X(2)) ' \AA'];
        elseif FirstShell
            lims = [0 Minima_X(1)];
            leg_add = [' 0-' num2str(Minima_X(1)) ' \AA'];
        elseif SecondShell
            lims = [Minima_X(1) Minima_X(2)];
            leg_add = [' ' num2str(Minima_X(1)) '-' num2str(Minima_X(2)) ' \AA'];
        elseif NeighbourSelected
            lims = str2double(app.NeighbEdit.String);
            leg_add = [' ' app.NeighbEdit.String ' Near N.'];
        else
            lims = [str2double(app.QlEditMin.String) str2double(app.QlEdit.String)];
            leg_add = [' ' app.QlEditMin.String '-' app.QlEdit.String ' \AA'];
        end
        
        [Ql,Wl] = GenQlWl(app,Input,NeighbourSelected,lims);
        
        if isempty(Ql)
            warndlg('Unable to construct theoretical 0K Ql.')
            return
        end
        
        app.Time{end+1} = nan;
        app.PotE{end+1} = nan;
        app.KinE{end+1} = nan;
        app.TotE{end+1} = nan;
        app.Temp{end+1} = nan;
        app.Press{end+1} = nan;
        app.Vol{end+1} = nan;
        app.Dens{end+1} = nan;
        app.Enth{end+1} = nan;
        app.MSD{end+1} = nan(1,2);
        app.NN{end+1} = nan(1,2);
        app.Ql_Dist{end+1} = nan;
        app.Wl_Dist{end+1} = nan;
        
        app.Time{end+1} = nan;
        app.PotE{end+1} = nan;
        app.KinE{end+1} = nan;
        app.TotE{end+1} = nan;
        app.Temp{end+1} = nan;
        app.Press{end+1} = nan;
        app.Vol{end+1} = nan;
        app.Dens{end+1} = nan;
        app.Enth{end+1} = nan;
        app.MSD{end+1} = nan(1,2);
        app.NN{end+1} = nan(1,2);
        app.Ql_Dist{end+1} = nan;
        app.Wl_Dist{end+1} = nan;
        
        switch Y_Axis
            case [Metal '-' Halide]
                app.RDF{end+1} = [RDF.R'./10 RDF.MX'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.MX'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.MX'];
                app.Wl{end+1} = [Wl.l' Wl.MX'];
            case [Metal '-' Metal]
                app.RDF{end+1} = [RDF.R'./10 RDF.MM'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.MM'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.MM'];
                app.Wl{end+1} = [Wl.l' Wl.MM'];
            case [Halide '-' Halide]
                app.RDF{end+1} = [RDF.R'./10 RDF.XX'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.XX'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.XX'];
                app.Wl{end+1} = [Wl.l' Wl.XX'];
            case [Metal '-All']
                app.RDF{end+1} = [RDF.R'./10 RDF.MA'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.MA'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.MA'];
                app.Wl{end+1} = [Wl.l' Wl.MA'];
            case [Halide '-All']
                app.RDF{end+1} = [RDF.R'./10 RDF.XA'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.XA'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.XA'];
                app.Wl{end+1} = [Wl.l' Wl.XA'];
            case 'All-All'
                app.RDF{end+1} = [RDF.R'./10 RDF.AA'];
                app.CRDF{end+1} = [CRDF.R'./10 CRDF.AA'];
                app.Ql{end+1} = nan(1,2);
                app.Wl{end+1} = nan(1,2);
                
                app.RDF{end+1} = nan(1,2);
                app.CRDF{end+1} = nan(1,2);
                app.Ql{end+1} = [Ql.l' Ql.AA'];
                app.Wl{end+1} = [Wl.l' Wl.AA'];
        end
        
        if Data_Type == 2
            Opttxt = 'Cell Opt.';
        elseif Data_Type == 3
            Opttxt = 'Full Opt.';
        end
        app.legtxt{end+1} = [Salt ' ' Structure ' ' Model ' ' Opttxt ' ' Y_Axis];
        app.legtxt{end+1} = [Salt ' ' Structure ' ' Model ' ' Opttxt ' ' Y_Axis leg_add];
        
        % Increment stuff
        app.NumDataSets = app.NumDataSets + 2;
        app.data_type(end+1) = 2;
        app.data_type(end+1) = 4;

        if isnan(app.rdf_type_plot) % First RDF type plot
            app.NumPlots  = app.NumPlots + 1;
            app.rdf_type_plot = app.NumPlots;
            app.Notes{end+1} = 'RDF plots';
        end
        if isnan(app.AvQl_type_plot) % First Ql type plot
            app.NumPlots  = app.NumPlots + 1;
            app.AvQl_type_plot = app.NumPlots;
            app.Notes{end+1} = 'Average Wl/Ql plots';
        end
        
        if strcmp(X_Axis,'RDF') || strcmp(X_Axis,'CRDF')
            app.PlotPoint = app.rdf_type_plot;
        elseif strcmp(X_Axis,'Average Wl') || strcmp(X_Axis,'Average Ql')
            app.PlotPoint = app.AvQl_type_plot;
        end
    end
    SavePlotState;
    RefactorPlots;
end

function LoadPlots
    ReloadAxis; % Clear the current plots
    
    % Load notes for the current plot point
    app.Notepad.String = app.Notes{app.PlotPoint};
    
    % Load the data type of the current plot
    uniqplot = unique_plot(app);
    dt_uniq = app.data_type(uniqplot);
    DatType = dt_uniq(app.PlotPoint);
    
    N = length(app.XDat);
    if N == 0
        RefreshGUI;
        return
    elseif ~iscell(app.XDat)
        if all(isnan(app.XDat),'all')
            RefreshGUI;
            return
        end
    end
	
    X_axis = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    Y_axis = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    
    if app.compareQlPlots % Special comparison Ql/Wl plot
        
        if contains(X_axis,{'Ql'}) || contains(Y_axis,{'Ql'})
            QorW = 'Q';
            QorWLab = 'Q';
        else
            QorW = 'W';
            QorWLab = '10^{2} \times W';
        end
        
        XDat1 = app.compareQl1.XDat(1,:);
        XDat2 = app.compareQl2.XDat(1,:);
        
        YDat1 = app.compareQl1.YDat(:,1);
        YDat2 = app.compareQl2.YDat(:,1);
        
        ZDat1 = app.compareQl1.ZDat{1,1};
        ZDat2 = app.compareQl2.ZDat{1,1};

        % Create 2 waterfall plots
        hold(app.Axh,'on')
        plot1 = waterfall(app.Axh,XDat1,YDat1,ZDat1);
        plot2 = waterfall(app.Axh,XDat2,YDat2,ZDat2);
        cmap1 = cbrewer('qual', 'Pastel1', length(app.l),'PCHIP');
        cmap2 = cbrewer('qual', 'Set1', length(app.l),'PCHIP');
        if length(app.l) < 3
            cmap1 = cmap1(1:length(app.l),:);
            cmap2 = cmap2(1:length(app.l),:);
        end
        set(plot1, 'FaceColor', 'flat');
        set(plot1, 'Linewidth', app.LW);
        set(plot1, 'FaceAlpha', 0.7);
        set(plot1, 'EdgeColor', 'k');
        set(plot1, 'FaceVertexCData', cmap1)
        set(plot2, 'FaceColor', 'flat');
        set(plot2, 'Linewidth', app.LW);
        set(plot2, 'FaceAlpha', 0.7);
        set(plot2, 'EdgeColor', 'k');
        set(plot2, 'FaceVertexCData', cmap2)
        
        set(app.Axh, 'Color', [1 1 1]*0.85)
        set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
        app.Axh.XLabel.Interpreter = 'latex';
        app.Axh.YLabel.Interpreter = 'latex';
        app.Axh.ZLabel.Interpreter = 'latex';
        app.Axh.XLabel.FontSize = app.FS+10;
        app.Axh.YLabel.FontSize = app.FS+10;
        app.Axh.ZLabel.FontSize = app.FS+10;
        app.Axh.XLabel.String = ['$' QorWLab '_{\ell}$ Values'];
        app.Axh.YLabel.String = '$\ell$';
        app.Axh.ZLabel.String = ['$P(' QorW '_{\ell})$'];
        set(app.Axh,'XTickMode', 'auto', 'XTickLabelMode', 'auto','XMinorTick','on')
        set(app.Axh,'YTickLabel',YDat1,'YTick',YDat1,'YMinorTick','off',...
            'YLim',[1 YDat1(end)+1]);
        set(app.Axh,'ZTickMode', 'auto', 'ZTickLabelMode', 'auto','ZMinorTick','on')
        view(app.Axh,app.plotview(1),app.plotview(2));
        grid(app.Axh,'on')
        zlim(app.Axh,'auto')

        if isprop(app.Legend,'Visible')
            app.Legend.Visible = 'off';
        end
        titletxt = ['\textit{Light}: ' app.Comparetxt1 newline '\textbf{Dark}: ' app.Comparetxt2];
        app.Title = title(app.Axh,titletxt,'Interpreter','latex');
        app.Title.Position = [mean(app.Axh.XLim) app.Axh.YLim(2) app.Axh.ZLim(2)*0.9];
        set(app.figh,'renderer','painters');
        
    elseif DatType == 0 % Single point type plots: 3D or 4D plots
        % Get plot title
        legtxt = app.legtxt(uniqplot);
        titletxt = legtxt(app.PlotPoint);
        
        if isprop(app.Legend,'Visible')
            app.Legend.Visible = 'off';
        end
        
        if contains(X_axis,{'Ql' 'Wl'}) || contains(Y_axis,{'Ql' 'Wl'})

            if contains(X_axis,{'Ql'}) || contains(Y_axis,{'Ql'})
                QorW = 'Q';
                QorWLab = 'Q';
                Ql = true;
            else
                QorW = 'W';
                QorWLab = 'W';
                Ql = false;
            end

            if contains(X_axis,{'Average'}) || contains(Y_axis,{'Average'}) % Average Q1/W1 plots
                surf(app.Axh,app.XDat,app.YDat,app.ZDat,'MeshStyle','column',...
                    'LineWidth',app.LW,'LineStyle','-');

                colormap(app.Axh,app.Colours_3D);
                set(app.Axh,'ZLim',[0 10])
                app.Title = title(app.Axh,titletxt,'Interpreter','latex');
                set(app.Axh, 'XTickMode', 'auto', 'XTickLabelMode', 'auto',...
                    'YTickMode','auto','YTickLabelMode','auto',...
                    'ZTickMode','auto','ZTickLabelMode','auto');
                view(app.Axh,app.plotview(1),app.plotview(2));
                grid(app.Axh,'on')
                set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
                app.Axh.XLabel.Interpreter = 'latex';
                app.Axh.YLabel.Interpreter = 'latex';
                app.Axh.ZLabel.Interpreter = 'latex';
                app.Axh.XLabel.FontSize = app.FS+10;
                app.Axh.YLabel.FontSize = app.FS+10;
                app.Axh.ZLabel.FontSize = app.FS+10;
                zlim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
                app.Axh.XLabel.String = '$\ell$';
                app.Axh.YLabel.String = 'Time [ns]';
                app.Axh.XTick = app.l;

                switch X_axis
                    case 'Average Ql'
                        app.Axh.ZLabel.String = '$\langle Q_{\ell} \rangle$';
                    case 'Average Wl'
                        app.Axh.ZLabel.String = '$\langle W_{\ell} \rangle$';
                end

                switch Y_axis
                    case 'Average Ql'
                        app.Axh.ZLabel.String = '$\langle Q_{\ell} \rangle$';
                    case 'Average Wl'
                        app.Axh.ZLabel.String = '$\langle W_{\ell} \rangle$';
                end        
                RefreshGUI;
            else % Q1/W1 distribution plots
                XDat = app.XDat(1,:);
                if isnan(XDat)
                    RefreshGUI;
                    return
                end
                YDat = app.YDat(:,1);
                ZDat = app.ZDat{1,1};

                if Ql
                    XTicklabs = [(XDat - app.HistBinQ/2) (XDat(end)+app.HistBinQ/2)];
                else
                    XTicklabs = [(XDat - app.HistBinW/2) (XDat(end)+app.HistBinW/2)];
                end

                % Create a waterfall plot
                app.barh = waterfall(app.Axh,XDat,YDat,ZDat);
                
                [cmap] = cbrewer('qual', 'Pastel1', length(app.l),'PCHIP');
                if length(app.l) < 3
                    cmap = cmap(1:length(app.l),:);
                end
                set(app.barh, 'FaceColor', 'flat');
                set(app.barh, 'Linewidth', app.LW);
                set(app.barh, 'FaceAlpha', 0.7);
                set(app.barh, 'EdgeColor', 'k');
                set(app.barh, 'FaceVertexCData', cmap)
                set(app.Axh, 'Color', [1 1 1]*0.85)
                set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
                app.Axh.XLabel.Interpreter = 'latex';
                app.Axh.YLabel.Interpreter = 'latex';
                app.Axh.ZLabel.Interpreter = 'latex';
                app.Axh.XLabel.FontSize = app.FS+10;
                app.Axh.YLabel.FontSize = app.FS+10;
                app.Axh.ZLabel.FontSize = app.FS+10;
                app.Axh.XLabel.String = ['$' QorWLab '_{\ell}$ Values'];
                app.Axh.YLabel.String = '$\ell$';
                app.Axh.ZLabel.String = ['$P(' QorW '_{\ell})$'];
                set(app.Axh,'XMinorTick','on',...
                    'XLim',[XTicklabs(1) XTicklabs(end)]);
                set(app.Axh,'YTickLabel',YDat,'YTick',YDat,'YMinorTick','off',...
                    'YLim',[1 YDat(end)+1]);
                
                zlim(app.Axh,[0 max(max(ZDat,[],'all'),20)])
                set(app.Axh,'ZMinorTick','off');
                view(app.Axh,app.plotview(1),app.plotview(2));
                grid(app.Axh,'on')

                if isprop(app.Legend,'Visible')
                    app.Legend.Visible = 'off';
                end
                app.Title = title(app.Axh,titletxt,'Interpreter','latex');
                RefreshGUI;
            end
        elseif contains(X_axis,{'RDF' 'MSD'}) || contains(Y_axis,{'RDF' 'MSD'})% 3D MSD/RDF/CRDF plots
            hold(app.Axh,'On')

            surf(app.Axh,app.XDat,app.YDat,app.ZDat,'MeshStyle','row',...
                'LineWidth',app.LW-2,'LineStyle','-');% RDF plots
            if isprop(app.Legend,'Visible')
                app.Legend.Visible = 'off';
            end

            colormap(app.Axh,app.Colours_3D);
            set(app.Axh, 'Color', [1 1 1]*0.85)
            app.Title = title(app.Axh,titletxt,'Interpreter','latex');
            set(app.Axh, 'XTickMode', 'auto', 'XTickLabelMode', 'auto',...
                'YTickMode','auto','YTickLabelMode','auto',...
                'ZTickMode','auto','ZTickLabelMode','auto');
            view(app.Axh,app.plotview(1),app.plotview(2));
            grid(app.Axh,'on')
            set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
            app.Axh.XLabel.Interpreter = 'latex';
            app.Axh.YLabel.Interpreter = 'latex';
            app.Axh.ZLabel.Interpreter = 'latex';
            app.Axh.XLabel.FontSize = app.FS+10;
            app.Axh.YLabel.FontSize = app.FS+10;
            app.Axh.ZLabel.FontSize = app.FS+10;

            % RDF plots
            if contains(X_axis,{'RDF'}) || contains(Y_axis,{'RDF'})
                app.Axh.XLabel.String = 'r [nm]';
                app.Axh.YLabel.String = 't [ns]';
                app.Axh.ZLabel.String = 'g(r)';
                set(app.Axh,'ZLim',[0 10])
                set(app.Axh,'XLim',[0 1.5])
            else
                app.Axh.XLabel.String = 't [ps]';
                app.Axh.YLabel.String = 'Traj. t [ns]';
                app.Axh.ZLabel.String = 'MSD (t)';
                set(app.Axh,'ZLim',[0 2])
                set(app.Axh,'XLim',[0 Inf])
                set(app.Axh,'TickLength',[0.01 0.025])
            end

            switch X_axis
                case 'CRDF'
                    zlim(app.Axh,[0 15])
                    set(app.Axh,'TickLength',[0.01 0.025])
                case 'RDF'
                    zlim(app.Axh,'auto')
                    set(app.Axh,'TickLength',[0.01 0.025])
            end

            switch Y_axis
                case 'CRDF'
                    zlim(app.Axh,[0 15])
                    set(app.Axh,'TickLength',[0.01 0.025])
                case 'RDF'
                    zlim(app.Axh,'auto')
                    set(app.Axh,'TickLength',[0.01 0.025])
            end
            RefreshGUI;
        elseif contains(X_axis,'Nearest N.') || contains(Y_axis,'Nearest N.')
            
            hold(app.Axh,'On')
            XDat = app.XDat(1,:);
            YDat = app.YDat(:,1);
            ZDat = app.ZDat;
            
            % Remove empty columns
            XDat( :, ~any(ZDat,1) ) = [];
            ZDat( :, ~any(ZDat,1) ) = [];
            
            % Anonymous function to scale Xdata
            scaler = @(vals)XDat(1) + ((vals-1) * (XDat(end) - XDat(1)) / (numel(XDat) - 1));

            % Create the bar plot
            app.barh = bar3(app.Axh,YDat,ZDat);
            
            set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
            app.Axh.XLabel.Interpreter = 'latex';
            app.Axh.YLabel.Interpreter = 'latex';
            app.Axh.ZLabel.Interpreter = 'latex';
            app.Axh.XLabel.FontSize = app.FS+10;
            app.Axh.YLabel.FontSize = app.FS+10;
            app.Axh.ZLabel.FontSize = app.FS+10;
            app.Axh.XLabel.String = 'Near Neighbours';
            app.Axh.YLabel.String = 'Time [ns]';
            app.Axh.ZLabel.String = 'Frac.';

            % Change the XData, make plot pretty
            bar_xdata = get(app.barh, 'XData');
            bar_xdata = cellfun(scaler, bar_xdata, 'uni', 0);
            set(app.barh, {'XData'}, bar_xdata);
            set(app.Axh,'XTick',XDat,'XTickLabel',XDat,'XMinorTick','off',...
                'XLim',[XDat(1)-0.5 XDat(end)+0.5]);
            set(app.Axh, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
            set(app.Axh, 'ZTickMode', 'auto', 'ZTickLabelMode', 'auto')
            set(app.Axh,'YMinorTick','on');
            set(app.Axh,'ZLim',[0 1],'ZMinorTick','on');
            yint = (YDat(2)-YDat(1))/2;
            set(app.Axh,'YLim',[YDat(1)-yint YDat(end)+yint],'ZMinorTick','on');
            view(app.Axh,app.plotview(1),app.plotview(2));
            grid(app.Axh,'on')
            colormap(app.Axh,app.Colours_3D);
            shading(app.Axh,'interp')
            for ii = 1:length(app.barh)
                 set(app.barh(ii),'Cdata',get(app.barh(ii),'Zdata'))
                 set(app.barh(ii),'EdgeColor','k')
            end
            set(app.Axh, 'Color', [1 1 1]*0.85)
            if isprop(app.Legend,'Visible')
                app.Legend.Visible = 'off';
            end
            RefreshGUI;
        end
        
    elseif any(DatType == [1 2 3]) % Any 2D plot other than average Wl, Ql, or NN
        legtxt = app.legtxt(app.data_type == DatType);
        
        % Check for missing data
        for ii = N:-1:1
            if all(isnan(app.XDat{ii})) || all(isnan(app.YDat{ii}))
                app.XDat(:,ii) = [];
                app.YDat(:,ii) = [];
                legtxt(:,ii) = [];
            end
        end
        if isempty(legtxt)
            RefreshGUI;
            return
        end
        
        ploth = gobjects(N,1);
        Colours_2D = cbrewer(app.ColScheme{1}, app.ColScheme{2}, max(N,3),'PCHIP');
        
        hold(app.Axh,'On')
        for ii = 1:N
            ploth(ii) = plot(app.Axh,app.XDat{ii},app.YDat{ii},'Color',Colours_2D(ii,:),...
                'LineWidth',app.LW,'LineStyle','-');
        end
        
        if isempty(ploth)
            return
        end
        app.Legend = legend(app.Axh,ploth,legtxt,'Interpreter','latex',...
            'Fontsize',app.FS,'Visible','on');
        set(app.Axh, 'XTickMode', 'auto', 'XTickLabelMode', 'auto',...
            'YTickMode','auto','YTickLabelMode','auto','XLim',[-Inf Inf],'YLim',[-Inf Inf]);
        
        switch X_axis
            case 'Potential Energy'
                app.Axh.XLabel.String = 'Potential Energy [kJ mol$^{-1}$]';
            case 'Total Energy'
                app.Axh.XLabel.String = 'Total Energy [kJ mol$^{-1}$]';
            case 'Time'
                app.Axh.XLabel.String = 'Time [ns]';
            case 'Kinetic Energy'
                app.Axh.XLabel.String = 'Kinetic Energy [kJ mol$^{-1}$]';
            case 'Temperature'
                app.Axh.XLabel.String = 'Temperature [K]';
            case 'Pressure'
                app.Axh.XLabel.String = 'Pressure [GPa]';
            case 'Volume'
                app.Axh.XLabel.String = 'Volume [cm$^{3}$ mol$^{-1}$]';
            case 'Density'
                app.Axh.XLabel.String = 'Density [kg m$^{-3}$]';
            case 'Enthalpy'
                app.Axh.XLabel.String = 'Enthalpy [kJ mol$^{-1}$]';
            case 'CRDF'
                app.Axh.XLabel.String = 'r [nm]';
                app.Axh.YLabel.String = 'Cumulat. g(r)';
                ylim(app.Axh,[0 15])
                set(app.Axh,'TickLength',[0.01 0.025])
            case 'RDF'
                app.Axh.XLabel.String = 'r [nm]';
                app.Axh.YLabel.String = 'g(r)';
                ylim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
            case 'MSD'
                app.Axh.XLabel.String = 't [ps]';
                app.Axh.YLabel.String = 'MSD(t) [nm]';
                ylim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
            otherwise
                ylim(app.Axh,'auto')
                xlim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
        end
        
        switch Y_axis
            case 'Potential Energy'
                app.Axh.YLabel.String = 'Potential Energy [kJ mol$^{-1}$]';
            case 'Total Energy'
                app.Axh.YLabel.String = 'Total Energy [kJ mol$^{-1}$]';
            case 'Time'
                app.Axh.YLabel.String = 'Time [ns]';
            case 'Kinetic Energy'
                app.Axh.YLabel.String = 'Kinetic Energy [kJ mol$^{-1}$]';
            case 'Temperature'
                app.Axh.YLabel.String = 'Temperature [K]';
            case 'Pressure'
                app.Axh.YLabel.String = 'Pressure [Bar]';
            case 'Volume'
                app.Axh.YLabel.String = 'Volume [cm$^{3}$ mol$^{-1}$]';
            case 'Density'
                app.Axh.YLabel.String = 'Density [kg m$^{-3}$]';
            case 'Enthalpy'
                app.Axh.YLabel.String = 'Enthalpy [kJ mol$^{-1}$]';
            case 'CRDF'
                app.Axh.XLabel.String = 'r [nm]';
                app.Axh.YLabel.String = 'Cumulat. g(r)';
                ylim(app.Axh,[0 15])
                set(app.Axh,'TickLength',[0.01 0.025])
            case 'RDF'
                app.Axh.XLabel.String = 'r [nm]';
                app.Axh.YLabel.String = 'g(r)';
                ylim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
            case 'MSD'
                app.Axh.XLabel.String = 't [ps]';
                app.Axh.YLabel.String = 'MSD(t) [nm]';
                ylim(app.Axh,'auto')
                set(app.Axh,'TickLength',[0.01 0.025])
        end
        
        set(app.Legend,'Location','Best')
        RefreshGUI;
    elseif DatType == 4 % Average Ql/Wl
        legtxt = app.legtxt(app.data_type == DatType);
        Colours_2D = cbrewer(app.ColScheme{1}, app.ColScheme{2}, max(N,3),'PCHIP');
        
        if contains(X_axis,{'Ql'}) || contains(Y_axis,{'Ql'})
            QorW = 'Q';
        else
            QorW = 'W';
        end
        
        hold(app.Axh,'on')
        bardata = nan(length(app.l),length(app.YDat));
        XDat = app.XDat{1};
        for ii = 1:N
            bardata(:,ii) = app.YDat{ii};
        end

        ploth = bar(app.Axh,XDat,bardata,'EdgeColor','k','LineWidth',app.LW);
        for ii = 1:length(ploth)
            ploth(ii).FaceColor = Colours_2D(ii,:);
        end

        view(app.Axh,2);
        set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
        app.Axh.XLabel.Interpreter = 'latex';
        app.Axh.YLabel.Interpreter = 'latex';
        app.Axh.XLabel.FontSize = app.FS+10;
        app.Axh.YLabel.FontSize = app.FS+10;
        app.Axh.XLabel.String = '$\ell$';
        app.Axh.YLabel.String = ['$\langle ' QorW '_{\ell} \rangle$'];
        app.Axh.ZLabel.String = '';
        XDat_lab = cellstr(   num2str( ((XDat(1)-1):1:(XDat(end)+1))' )  );
        XDat_lab(1:2:end) = {''};

        yticks(app.Axh,'auto')
        yticklabels(app.Axh,'auto')
        ylim(app.Axh,'auto')
        set(app.Axh,'XTick',(XDat(1)-1):1:(XDat(end)+1),'XTickLabel',XDat_lab,'XMinorTick','off',...
            'XLim',[XDat(1)-1 XDat(end)+1]);
        set(app.Axh,'YMinorTick','on');%,'YLim',zlimits);
        app.Legend = legend(app.Axh,ploth,legtxt,'Interpreter','latex',...
            'Fontsize',app.FS,'Visible','on');
        RefreshGUI;
    elseif DatType == 5 % Nearest Neighbour
        legtxt = app.legtxt(app.data_type == DatType);
        Colours_2D = cbrewer(app.ColScheme{1}, app.ColScheme{2}, max(N,3),'PCHIP');
        
        hold(app.Axh,'on')
        
        bardata = nan(length(app.XDat{1}),length(app.YDat));
        XDat = app.XDat{1};
        for ii = 1:N
            bardata(:,ii) = app.YDat{ii};
        end
        
        % Remove empty rows
        XDat(~any(bardata,2),:) = [];
        bardata(~any(bardata,2),:) = [];
        
        ploth = bar(app.Axh,XDat,bardata,'EdgeColor','k','LineWidth',app.LW);
        for ii = 1:length(ploth)
            ploth(ii).FaceColor = Colours_2D(ii,:);
        end

        view(app.Axh,2);
        set(app.Axh,'TickLabelInterpreter','latex','Box','on','FontSize',app.FS)
        app.Axh.XLabel.Interpreter = 'latex';
        app.Axh.YLabel.Interpreter = 'latex';
        app.Axh.XLabel.FontSize = app.FS+10;
        app.Axh.YLabel.FontSize = app.FS+10;
        app.Axh.XLabel.String = 'Nearest Neighbours';
        app.Axh.YLabel.String = 'Frac.';
        app.Axh.ZLabel.String = '';
        XDat_lab = cellstr(   num2str( ((XDat(1)):1:(XDat(end)))' )  );

        yticks(app.Axh,'auto')
        yticklabels(app.Axh,'auto')
        ylim(app.Axh,[0 1])
        set(app.Axh,'XTick',(XDat(1)):1:(XDat(end)),'XTickLabel',XDat_lab,'XMinorTick','off',...
            'XLim',[XDat(1)-0.5 XDat(end)+0.5]);
        set(app.Axh,'YMinorTick','on');%,'YLim',zlimits);
        app.Legend = legend(app.Axh,ploth,legtxt,'Interpreter','latex',...
            'Fontsize',app.FS,'Visible','on');
        RefreshGUI;
    end
end

function RefactorPlots
    app.TimeInUse = false;
    
    if app.NumPlots == 0
        % Nothing to plot
        return
    end
    
    % Update notes
    app.Notepad.String = app.Notes{app.PlotPoint};
    
    % Grab plot settings
    X_Axis = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    Y_Axis = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    
    if app.PlotPoint == app.energy_type_plot % 'Energy type' plots
        noplot = false;
        switch X_Axis
            case 'Potential Energy'
                objnameX = 'PotE';
            case 'Total Energy'
                objnameX = 'TotE';
            case 'Time'
                objnameX = 'Time';
            case 'Kinetic Energy'
                objnameX = 'KinE';
            case 'Temperature'
                objnameX = 'Temp';
            case 'Pressure'
                objnameX = 'Press';
            case 'Volume'
                objnameX = 'Vol';
            case 'Density'
                objnameX = 'Dens';
            case 'Enthalpy'
                objnameX = 'Enth';
            otherwise
                noplot = true;
        end
        
        switch Y_Axis
            case 'Potential Energy'
                objnameY = 'PotE';
            case 'Total Energy'
                objnameY = 'TotE';
            case 'Time'
                objnameY = 'Time';
            case 'Kinetic Energy'
                objnameY = 'KinE';
            case 'Temperature'
                objnameY = 'Temp';
            case 'Pressure'
                objnameY = 'Press';
            case 'Volume'
                objnameY = 'Vol';
            case 'Density'
                objnameY = 'Dens';
            case 'Enthalpy'
                objnameY = 'Enth';
            otherwise
                noplot = true;
        end
        
        if noplot
            app.Axh.XLabel.String = '';
            app.Axh.YLabel.String = '';
            app.XDat = nan;
            app.YDat = nan;
            app.ZDat = nan;
        else
            app.XDat = app.(objnameX)(app.data_type == 1);
            app.YDat = app.(objnameY)(app.data_type == 1);
            app.ZDat = nan;
        end
        
    elseif app.PlotPoint == app.rdf_type_plot % 'RDF type' plots
        
        if strcmpi(X_Axis,'RDF') || strcmpi(Y_Axis,'RDF')% RDF
            
            DOI = app.RDF(app.data_type == 2);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        elseif strcmpi(X_Axis,'CRDF') || strcmpi(Y_Axis,'CRDF')% CRDF
            DOI = app.CRDF(app.data_type == 2);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        else
            app.Axh.XLabel.String = '';
            app.Axh.YLabel.String = '';
            app.XDat = nan;
            app.YDat = nan;
            app.ZDat = nan;
        end
    elseif app.PlotPoint == app.NN_type_plot % Nearest neighbour plots
        if strcmpi(X_Axis,'Nearest N.') || strcmpi(Y_Axis,'Nearest N.')% Nearest Neighbour
            
            DOI = app.NN(app.data_type == 5);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        else
            app.Axh.XLabel.String = '';
            app.Axh.YLabel.String = '';
            app.XDat = nan;
            app.YDat = nan;
            app.ZDat = nan;
        end
    elseif app.PlotPoint == app.msd_type_plot % 'MSD type' plots
        
        if strcmpi(X_Axis,'MSD') || strcmpi(Y_Axis,'MSD')% RDF
            
            DOI = app.MSD(app.data_type == 3);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        else
            app.Axh.XLabel.String = '';
            app.Axh.YLabel.String = '';
            app.XDat = nan;
            app.YDat = nan;
            app.ZDat = nan;
        end
    elseif app.PlotPoint == app.AvQl_type_plot % 'Average Ql type' plots

        if strcmpi(X_Axis,'Average Ql') || strcmpi(Y_Axis,'Average Ql') % Ql
            
            DOI = app.Ql(app.data_type == 4);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        elseif strcmpi(X_Axis,'Average Wl') || strcmpi(Y_Axis,'Average Wl')
            DOI = app.Wl(app.data_type == 4);
            N = length(DOI);
            
            app.XDat = cell(1,N);
            app.YDat = cell(1,N);
            app.ZDat = nan;
            
            for idx = 1:N
                app.XDat{idx} = DOI{idx}(:,1);
                app.YDat{idx} = DOI{idx}(:,2);
            end
        else
            app.Axh.XLabel.String = '';
            app.Axh.YLabel.String = '';
            app.XDat = nan;
            app.YDat = nan;
            app.ZDat = nan;
        end
        
    else % Plots that cannot be stacked (3D plots, 4D Distributions)
        unique_idx = unique_plot(app);
        
        if strcmpi(X_Axis,'MSD') || strcmpi(Y_Axis,'MSD')
            Unique_MSD = app.MSD(unique_idx);
            Curr_Data  = Unique_MSD{app.PlotPoint};
            if size(Curr_Data,3) <= 1
                app.XDat   = nan;
                app.YDat   = nan;
                app.ZDat   = nan;
            else
                app.XDat   = squeeze(Curr_Data(1,:,:));
                app.YDat   = squeeze(Curr_Data(2,:,:));
                app.ZDat   = squeeze(Curr_Data(3,:,:));
            end

        elseif strcmpi(X_Axis,'RDF') || strcmpi(Y_Axis,'RDF')% RDF
            Unique_RDF = app.RDF(unique_idx);
            Curr_Data  = Unique_RDF{app.PlotPoint};
            if size(Curr_Data,3) <= 1
                app.XDat   = nan;
                app.YDat   = nan;
                app.ZDat   = nan;
            else
                app.XDat   = squeeze(Curr_Data(1,:,:));
                app.YDat   = squeeze(Curr_Data(2,:,:));
                app.ZDat   = squeeze(Curr_Data(3,:,:));
            end
        elseif strcmpi(X_Axis,'CRDF') || strcmpi(Y_Axis,'CRDF') % CRDF
            
            Unique_CRDF = app.CRDF(unique_idx);
            Curr_Data  = Unique_CRDF{app.PlotPoint};
            if size(Curr_Data,3) <= 1
                app.XDat   = nan;
                app.YDat   = nan;
                app.ZDat   = nan;
            else
                app.XDat   = squeeze(Curr_Data(1,:,:));
                app.YDat   = squeeze(Curr_Data(2,:,:));
                app.ZDat   = squeeze(Curr_Data(3,:,:));
            end
        elseif strcmpi(X_Axis,'Nearest N.') || strcmpi(Y_Axis,'Nearest N.')
            Unique_NN = app.NN(unique_idx);
            Curr_Data  = Unique_NN{app.PlotPoint};
            if size(Curr_Data,3) <= 1
                app.XDat   = nan;
                app.YDat   = nan;
                app.ZDat   = nan;
            else
                app.XDat   = squeeze(Curr_Data(1,:,:));
                app.YDat   = squeeze(Curr_Data(2,:,:));
                app.ZDat   = squeeze(Curr_Data(3,:,:));
            end
        elseif strcmpi(X_Axis,'Average Wl') || strcmpi(Y_Axis,'Average Wl')
            Unique_Wl = app.Wl(unique_idx);
            Curr_Data  = Unique_Wl{app.PlotPoint};
            if size(Curr_Data,3) <= 1
                app.XDat   = nan;
                app.YDat   = nan;
                app.ZDat   = nan;
            else
                app.XDat   = squeeze(Curr_Data(1,:,:));
                app.YDat   = squeeze(Curr_Data(2,:,:));
                app.ZDat   = squeeze(Curr_Data(3,:,:));
            end
            
        elseif strcmpi(X_Axis,'Ql Distribution') || strcmpi(Y_Axis,'Ql Distribution')
            Unique_Ql_Dist = app.Ql_Dist(unique_idx);
            if ~iscell(Unique_Ql_Dist{app.PlotPoint})
                app.XDat = nan;
                app.YDat = nan;
                app.ZDat = nan;
            else
                app.TimeInUse = true;
                try
                    Curr_Data  = Unique_Ql_Dist{app.PlotPoint}(app.TimePoint,:);
                catch
                    app.TimePoint = 1;
                    Curr_Data  = Unique_Ql_Dist{app.PlotPoint}(app.TimePoint,:);
                end

                app.XDat   = Curr_Data{1,1}(:,:,1);
                app.YDat   = Curr_Data{1,1}(:,:,2);

                ZDat = cell(1,2);
                ZDat{1,1} = Curr_Data{1,1}(:,:,3);
                ZDat{1,2} = Curr_Data{1,2};
                app.ZDat = ZDat;
            end
            
        elseif strcmpi(X_Axis,'Wl Distribution') || strcmpi(Y_Axis,'Wl Distribution')         
            Unique_Wl_Dist = app.Wl_Dist(unique_idx);
            if ~iscell(Unique_Wl_Dist{app.PlotPoint})
                app.XDat = nan;
                app.YDat = nan;
                app.ZDat = nan;
            else
                app.TimeInUse = true;
                Curr_Data  = Unique_Wl_Dist{app.PlotPoint}(app.TimePoint,:);

                app.XDat   = Curr_Data{1,1}(:,:,1);
                app.YDat   = Curr_Data{1,1}(:,:,2);

                ZDat = cell(1,2);
                ZDat{1,1} = Curr_Data{1,1}(:,:,3);
                ZDat{1,2} = Curr_Data{1,2};
                app.ZDat = ZDat;
            end
        end
    end
    LoadPlots;
end

function TimeStepsCallback(src,~)
    str=get(src,'String');
    if isempty(str2double(str))
        set(src,'string','1');
        warndlg('The time per used frame must be a positive number.');
    elseif str2double(str) < 0
        set(src,'string','1');
        warndlg('The time per used frame must be a positive number.');
    else
        return
    end
end

function PlotTypeXCallback(~,~)
    XStr = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    YStr = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    RDF_Sely = ~isempty(regexp(YStr,'RDF','once')) || ~isempty(regexp(YStr,'[Q|W]l','once')) || ~isempty(regexp(YStr,'Near','once'));
    RDF_Selx = ~isempty(regexp(XStr,'RDF','once')) || ~isempty(regexp(XStr,'[Q|W]l','once')) || ~isempty(regexp(XStr,'Near','once'));
    MSD_Sely = ~isempty(regexp(YStr,'MSD','once'));
    MSD_Selx = ~isempty(regexp(XStr,'MSD','once'));
    
    if RDF_Selx % [C]RDF or [Q|W]l selected on X axis
        Salt = app.Salth.String{app.Salth.Value};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        RDFOptions = {[Metal '-' Halide] [Metal '-' Metal] [Halide '-' Halide] ...
                [Metal '-All'] [Halide '-All'] 'All-All'};
        set(app.PlotTypeYh,'String',RDFOptions)
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        if isempty(idy)
            idy = 6;
        end
    elseif MSD_Selx
        Salt = app.Salth.String{app.Salth.Value};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        MSDOptions = {Metal Halide 'All'};
        set(app.PlotTypeYh,'String',MSDOptions)
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        if isempty(idy)
            idy = 3;
        end
    elseif ~RDF_Sely && ~MSD_Sely % [C]RDF nor [Q|W]l selected on neither axis
        set(app.PlotTypeYh,'String',app.Options)
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        if isempty(idy)
            idy = 1;
        end
    else
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        if isempty(idy)
            idy = 1;
        end
    end
    
    % Set Y axis to its previous setting if it exists
    app.PlotTypeYh.Value = idy;
    app.TimePoint = 1;
    
    % Disable/enable 3rd axis stuff
    ThirdAxisCallback;
    RefactorPlots;
end

function PlotTypeYCallback(~,~)
    XStr = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    YStr = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    Distr_Sely = ~isempty(regexp(YStr,'RDF','once')) || ...
        ~isempty(regexp(YStr,'[Q|W]l','once')) || ~isempty(regexp(YStr,'Near','once'));
    Distr_Selx = ~isempty(regexp(XStr,'RDF','once')) || ...
        ~isempty(regexp(XStr,'[Q|W]l','once')) || ~isempty(regexp(XStr,'Near','once'));
    MSD_Sely = ~isempty(regexp(YStr,'MSD','once'));
    MSD_Selx = ~isempty(regexp(XStr,'MSD','once'));
    
    if Distr_Sely % [C]RDF or [Q|W]l selected on Y axis
        Salt = app.Salth.String{app.Salth.Value};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        RDFOptions = {[Metal '-' Halide] [Metal '-' Metal] [Halide '-' Halide] ...
                [Metal '-All'] [Halide '-All'] 'All-All'};
        set(app.PlotTypeXh,'String',RDFOptions)
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        if isempty(idx)
            idx = 6;
        end
    elseif MSD_Sely % MSD selected on X axis
        Salt = app.Salth.String{app.Salth.Value};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        MSDOptions = {Metal Halide 'All'};
        set(app.PlotTypeXh,'String',MSDOptions)
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        if isempty(idx)
            idx = 3;
        end
    elseif ~Distr_Selx && ~MSD_Selx % [C]RDF nor [Q|W]l selected on neither axis
        set(app.PlotTypeXh,'String',app.Options)
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        if isempty(idx)
            idx = 1;
        end
    else
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        if isempty(idx)
            idx = 1;
        end
    end
    
    % Set X axis to its previous setting if it exists
    app.PlotTypeXh.Value = idx;
    app.TimePoint = 1;
    
    % Disable/enable 3rd axis GUI
    ThirdAxisCallback;
    RefactorPlots;
end

function DataTypeCallback(~,~)
    val = app.DataTypeh.Value;
    XStr = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    YStr = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    
    if val == 1 % MD selected
        set(app.Structureh,'Visible','Off')
        set(app.BGLabel(3),'Visible','Off')
        
        Salt = app.Salth.String{app.Salth.Value};
        
        % Search for models in results folder
        % Find models associated with selected salt
        filesi = dir([app.datdir filesep Salt]);
        % Get a logical vector that tells which is a directory.
        dirFlagsi = [filesi.isdir];
        % Extract only those that are directories.
        subFoldersi = filesi(dirFlagsi);
        
        % Filter by date if selected
        % This is TRUE if user selects only new data
        if logical(app.Modelfilter_b1.Value)
            current_date = today('datenum');
            keep_idx = current_date - [subFoldersi.datenum] < 7;
            subFoldersi = subFoldersi(keep_idx);
        end
        
        % remove '.' and '..'
        dfoldersi = {subFoldersi(~ismember({subFoldersi(:).name},{'.','..'})).name};
        if isempty(app.Modelh.String)
            cur_sel = [];
        else
            cur_sel = find(cellfun(@(x) strcmp(x,app.Modelh.String{app.Modelh.Value}),dfoldersi),1);
        end

        if isempty(cur_sel)
            set(app.Modelh,'String',dfoldersi,'Value',1)
        else
            set(app.Modelh,'String',dfoldersi,'Value',cur_sel)
        end
        set(app.InitialTimeh,'Visible','On')
        set(app.FinalTimeh,'Visible','On')
        set(app.BGLabel(5),'Visible','On')
        set(app.BGLabel(6),'Visible','On')
        set(app.BGLabel(7),'Visible','On')
        set(app.TimePerFrameh,'Visible','On')
        
        set(app.SliceTrajectoryh,'Visible','On')
        set(app.PlotTypeXh,'String',app.Options)
        
        % Set X axis and Y axis to their previous setting
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        app.PlotTypeXh.Value = idx;
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        app.PlotTypeYh.Value = idy;
    else % minimized selected
        set(app.Structureh,'Visible','On')
        set(app.BGLabel(3),'Visible','On')
        set(app.Modelh,'String',app.ModelMenu,'Value',1)
        set(app.InitialTimeh,'Visible','Off')
        set(app.FinalTimeh,'Visible','Off')
        set(app.BGLabel(5),'Visible','Off')
        set(app.BGLabel(6),'Visible','Off')
        set(app.BGLabel(7),'Visible','Off')
        set(app.TimePerFrameh,'Visible','Off')
        options = {'RDF' 'CRDF' 'Average Ql' 'Average Wl'};
        cur_sel = find(cellfun(@(x) strcmp(x,app.PlotTypeXh.String{app.PlotTypeXh.Value}),options),1);

        if isempty(cur_sel)
            set(app.PlotTypeXh,'String',options,'Value',1)
        else
            set(app.PlotTypeXh,'String',options,'Value',cur_sel)
        end
        
        Salt = app.Salth.String{app.Salth.Value};
        [Metal,Halide] = Separate_Metal_Halide(Salt);
        RDFOptions = {[Metal '-' Halide] [Metal '-' Metal] [Halide '-' Halide] ...
        	[Metal '-All'] [Halide '-All'] 'All-All'};
        cur_sel = find(cellfun(@(x) strcmp(x,YStr),RDFOptions),1);
        if isempty(cur_sel)
            set(app.PlotTypeYh,'Value',1,'String',RDFOptions)
        else
            set(app.PlotTypeYh,'Value',cur_sel,'String',RDFOptions)
        end
        set(app.SliceTrajectoryh,'Visible','Off')
        
        % Set X axis and Y axis to their previous setting
        idx = find(strcmp(app.PlotTypeXh.String,XStr),1);
        idy = find(strcmp(app.PlotTypeYh.String,YStr),1);
        if isempty(idx)
            idx = 1;
        end
        if isempty(idy)
            idy = 1;
        end
        app.PlotTypeXh.Value = idx;
        app.PlotTypeYh.Value = idy;
    end
    
    % Disable/enable 3rd axis GUI stuff
    ThirdAxisCallback;
end

function SaltCallback(varargin)
    Salt = app.Salth.String{app.Salth.Value};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    dt = get(app.DataTypeh,'Value');
    
    if dt == 1 % MD
        % Search for models in results folder
        % Find models associated with selected salt
        filesi = dir([app.datdir filesep Salt]);
        % Get a logical vector that tells which is a directory.
        dirFlagsi = [filesi.isdir];
        % Extract only those that are directories.
        subFoldersi = filesi(dirFlagsi);
        
        % Filter by date if selected
        % This is TRUE if user selects only new data
        if logical(app.Modelfilter_b1.Value)
            current_date = today('datenum');
            keep_idx = current_date - [subFoldersi.datenum] < 7;
            subFoldersi = subFoldersi(keep_idx);
        end
        
        % remove '.' and '..'
        dfoldersi = {subFoldersi(~ismember({subFoldersi(:).name},{'.','..'})).name};
        if app.Modelh.Value > length(dfoldersi)
            app.Modelh.Value = length(dfoldersi);
        elseif ~isempty(dfoldersi) && app.Modelh.Value == 0
            app.Modelh.Value = 1;
        end
        set(app.Modelh,'String',dfoldersi)
    else % Minimized
        set(app.Modelh,'String',app.ModelMenu,'Value',1)
    end
    
    set(app.TJS.Type_Menu,'String',{'All' Metal Halide})

    RDFOptions = {[Metal '-' Halide] [Metal '-' Metal] [Halide '-' Halide] ...
            [Metal '-All'] [Halide '-All'] 'All-All'};
    MSDOptions = {Metal Halide 'All'};
        
    XStr = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    YStr = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    Distr_Selx = ~isempty(regexp(XStr,'RDF','once')) || ~isempty(regexp(XStr,'[Q|W]l','once'));
    Distr_Sely = ~isempty(regexp(YStr,'RDF','once')) || ~isempty(regexp(YStr,'[Q|W]l','once'));
    MSD_Selx = ~isempty(regexp(XStr,'MSD','once'));
    MSD_Sely = ~isempty(regexp(YStr,'MSD','once'));
        
    if Distr_Selx
        set(app.PlotTypeYh,'String',RDFOptions)
    elseif Distr_Sely
        set(app.PlotTypeXh,'String',RDFOptions)
    elseif MSD_Selx
        set(app.PlotTypeYh,'String',MSDOptions)
    elseif MSD_Sely
        set(app.PlotTypeXh,'String',MSDOptions)
    end

end

function InitialTimeCallback(src,~)
    str=get(src,'String');
    if isempty(str2double(str))
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    elseif str2double(str) < 0
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    else
        return
    end
end

function FinalTimeCallback(src,~)
    str=get(src,'String');
    if isempty(str2double(str))
        set(src,'string','30');
        warndlg('The final time point must be a non-negative number.');
    elseif str2double(str) < 0
        set(src,'string','30');
        warndlg('The final time point must be a non-negative number.');
    else
        return
    end
end

function DeleteAllCallback(~,~)
    % Generate Empty data arrays
    app.XDat = {};
    app.YDat = {};
    app.ZDat = {};
    app.Time = {};
    app.PotE = {};
    app.KinE = {};
    app.TotE = {};
    app.Temp = {};
    app.Press = {};
    app.Vol = {};
    app.Dens = {};
    app.Enth = {};
    app.MSD = {};
    app.RDF = {};
    app.CRDF = {};
    app.NN = {};
    app.Ql = {};
    app.Wl = {};
    app.Ql_Dist = {};
    app.Wl_Dist = {};
    app.legtxt = {}; % Keeps track of the legend or title info associated with each data set
    app.Notes = {}; % Keeps track of the notes associated with each plot
    app.PlotState = struct();  % Keeps track of the state associated with each plot
    app.data_type = []; % Keeps track of whether each data set is an "energy type" data set
    app.energy_type_plot = nan; % Keeps track of the index for the "energy type" plot
    app.rdf_type_plot = nan; % Keeps track of the index for the "rdf type" plot
    app.msd_type_plot = nan; % Keeps track of the index for the "msd type" plot
    app.AvQl_type_plot = nan; % Keeps track of the index for the "average Ql type" plot
    app.NN_type_plot = nan; % Keeps track of the index for the "nearest neighbour type" plot
    app.plotview = [-37.5000 30]; % default 3d angle
    app.TimePoint = 1; % Index to keeps track of the current time point in Ql and Wl plots
    app.PlotPoint = 0; % Index to track of the current plot
    app.NumDataSets = 0; % Keeps track of the number of data sets
    app.NumPlots = 0; % Keeps track of the number of number of number of plots
    app.TimeInUse = false; % Is the time in use
    app.compareQlPlots = false; % Sets whether or not Ql plots are to be compared
    app.compareQl1.XDat = [];
    app.compareQl1.YDat = [];
    app.compareQl1.ZDat = [];
    app.compareQl2.XDat = [];
    app.compareQl2.YDat = [];
    app.compareQl2.ZDat = [];
    app.Comparetxt1 = '';
    app.Comparetxt2 = '';
    app.compareQlButton1.BackgroundColor = [0.9400    0.9400    0.9400];
    app.compareQlButton1b.BackgroundColor = [0.9400    0.9400    0.9400];
    app.compareQlButton2.BackgroundColor = [0.9400    0.9400    0.9400];
    app.compareQlButton2b.BackgroundColor = [0.9400    0.9400    0.9400];
    
    ReloadAxis;
    RefreshGUI;
end

function DeletePlotCallback(~,~)
    if app.NumPlots == 0
        return
    end
    unipl = unique_plot(app);
    unidt = app.data_type(unipl);
    current_dt = unidt(app.PlotPoint); % data type to be deleted
    
    if current_dt == 0
        % Get the index of the specific dataset to be deleted
        fidx = find(unipl);
        didx = fidx(app.PlotPoint);
        
        % Delete Data Set Data
        app.Time(didx) = [];
        app.PotE(didx) = [];
        app.KinE(didx) = [];
        app.TotE(didx) = [];
        app.Temp(didx) = [];
        app.Press(didx) = [];
        app.Vol(didx) = [];
        app.Dens(didx) = [];
        app.Enth(didx) = [];
        app.MSD(didx) = [];
        app.RDF(didx) = [];
        app.CRDF(didx) = [];
        app.NN(didx) = [];
        app.Ql(didx) = [];
        app.Wl(didx) = [];
        app.Ql_Dist(didx) = [];
        app.Wl_Dist(didx) = [];
        app.legtxt(didx) = [];
        app.data_type(didx) = [];
        
        % Delete Plot point data
        app.Notes(app.PlotPoint) = [];
        app.PlotState(app.PlotPoint) = [];  % Keeps track of the state associated with each plot
        
        % Rearrange pointers if required
        if app.PlotPoint < app.energy_type_plot
            app.energy_type_plot = app.energy_type_plot-1;
        end
        if app.PlotPoint < app.rdf_type_plot
            app.rdf_type_plot = app.rdf_type_plot-1;
        end
        if app.PlotPoint < app.msd_type_plot
            app.msd_type_plot = app.msd_type_plot-1;
        end
        if app.PlotPoint < app.AvQl_type_plot
            app.AvQl_type_plot = app.AvQl_type_plot-1;
        end
        if app.PlotPoint < app.NN_type_plot
            app.NN_type_plot = app.NN_type_plot-1;
        end
        
        app.NumPlots = app.NumPlots-1;
        app.NumDataSets = app.NumDataSets-1;
        if app.PlotPoint > app.NumPlots
            app.PlotPoint = app.PlotPoint-1;
        end
        
    else % In other cases, delete all data associated with the specific data type, and remove the plot
        
        % Get the indeces of the specific dataset to be deleted
        didx = app.data_type == current_dt;
        locs = find(didx);% Indexes 
        num_del = length(locs); % Number of data sets to be deleted
        
        % Delete Data Set Data
        app.Time(didx) = [];
        app.PotE(didx) = [];
        app.KinE(didx) = [];
        app.TotE(didx) = [];
        app.Temp(didx) = [];
        app.Press(didx) = [];
        app.Vol(didx) = [];
        app.Dens(didx) = [];
        app.Enth(didx) = [];
        app.MSD(didx) = [];
        app.RDF(didx) = [];
        app.CRDF(didx) = [];
        app.Ql(didx) = [];
        app.Wl(didx) = [];
        app.Ql_Dist(didx) = [];
        app.Wl_Dist(didx) = [];
        app.legtxt(didx) = [];
        app.data_type(didx) = [];
        
        % Delete Plot point data
        app.Notes(app.PlotPoint) = [];
        app.PlotState(app.PlotPoint) = [];  % Keeps track of the state associated with each plot
        
        % Rearrange pointers as required
        if app.PlotPoint < app.energy_type_plot
            app.energy_type_plot = app.energy_type_plot-1;
        end
        if app.PlotPoint < app.rdf_type_plot
            app.rdf_type_plot = app.rdf_type_plot-1;
        end
        if app.PlotPoint < app.msd_type_plot
            app.msd_type_plot = app.msd_type_plot-1;
        end
        if app.PlotPoint < app.AvQl_type_plot
            app.AvQl_type_plot = app.AvQl_type_plot-1;
        end
        app.NumPlots = app.NumPlots-1;
        app.NumDataSets = app.NumDataSets-num_del;
        if app.PlotPoint > app.NumPlots
            app.PlotPoint = app.PlotPoint-1;
        end
        
        % Reset the pointer of the particular data set
        if current_dt == 1
            app.energy_type_plot = nan;
        elseif current_dt == 2
            app.rdf_type_plot = nan;
        elseif current_dt == 3
            app.msd_type_plot = nan;
        elseif current_dt == 4
            app.AvQl_type_plot = nan;
        end
    end
    ReloadAxis;
    LoadPlotState;
    RefactorPlots;
    RefreshGUI;
end

function ViewTrajectoryCallback(~,~)
    
    %% Grab plot settings
    Salt = app.Salth.String{app.Salth.Value};
    Model = app.Modelh.String{app.Modelh.Value}; % {'JC (SPC/E)' 'JC (TIP4P-EW)' 'JC (TIP3P)' 'TF'};    
        
    % Directory containing results of interest
    DatDir = [app.datdir filesep Salt filesep Model];
    
    % Time slice options
    Min_Time = str2double(app.TJS.InitTime_Edit.String)*1000; % min time in ps
    Max_Time = str2double(app.TJS.FinTime_Edit.String)*1000; % max time in ps
    TimePerFrame = str2double(app.TJS.FrameTime_Edit.String); % time per frame (ps)

    % Select the external viewer
    Viewer = app.TJS.Select_Viewer_text.String{app.TJS.Select_Viewer_text.Value};
    
    % Grab the current filters
    Filters = app.TJS.current_filters;
    
    % Sanity check
    if Min_Time > Max_Time
        warndlg('The Initial Time must be less than or equal to the final time.');
        return
    elseif isnan(Min_Time) || isnan(Max_Time) || isnan(TimePerFrame)
        warndlg('Unable to interpret the selected time slice.');
        return
    end
    
    % Find trajectory and gro files
    Traj_info = dir([DatDir filesep '*.trr']);
    Gro_info = dir([DatDir filesep '*.gro']);
    Mdp_info = dir([DatDir filesep '*.mdp']);
    if isempty(Traj_info) || isempty(Gro_info) || isempty(isempty(Mdp_info)) % no trajectory/gro file found
        warndlg('Unable to find required data file(s).')
        return
    end
    
    if length(Mdp_info) > 1
        Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
    end
    Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
    Mdp_text = fileread(Mdp_file);
    pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));

    % Select the initial time step gro file
    [~,idx] = min(cellfun(@length,{Gro_info.name}));
    Gro_names = {Gro_info.name};
    Gro_name = Gro_names{idx};

    Traj_file =[DatDir filesep Traj_info.name]; % traj file
    Gro_file = [DatDir filesep Gro_name]; % gro file
    % Tack on a random length 10 string to the name
    str = char(97 + floor(26 .* rand(10,1)))';
    Out_file = [app.tempfolder filesep Salt '_' Model '_' str '.pdb'];
    
    set(app.figh, 'pointer', 'watch')
    set(app.auxfig, 'pointer', 'watch')
    drawnow;
    try
        py.LiXStructureDetector.Slice_Traj(Traj_file,Gro_file,Out_file,Min_Time,Max_Time,...
            TimePerFrame,Filters,pbc_on,Viewer)
    catch e
        warndlg(['Error Executing Command:' newline e.message]);
        return
    end
    set(app.figh, 'pointer', 'arrow')
    set(app.auxfig, 'pointer', 'arrow')
    
    % Generate view trajectory command
    if strcmpi(Viewer,'OVITO')
        traj_view_command = ['start "C:\Program Files\OVITO Basic\ovito.exe" "' Out_file '"'];
    elseif strcmpi(Viewer,'VMD')
        if ispc
            traj_view_command = ['start ' app.vmd_loc ' "' Out_file '" -args "pbc box"'];
        else
            traj_view_command = ['gnome-terminal -- ' app.vmd_loc ' "' Out_file '"'];
        end
    elseif strcmpi(Viewer,'PYMOL')
        traj_view_command = ['start "C:\ProgramData\Anaconda3\PyMOLWin.exe" "' Out_file '"'];
    end
    system(traj_view_command);
end
    
function CloseAppCallback(~,~)
    
    % Delete the auxiliary figure
    delete(app.auxfig);
    
    % Save app state
    if app.PlotPoint > 1
        app.Notes{app.PlotPoint} = app.Notepad.String;
    end
    SaveAppState;

    % Close vmd
    if ispc
        vmdexe = 'vmd.exe';
        commandLine = sprintf('tasklist /FI "IMAGENAME eq %s"', vmdexe);
        % Now execute that command line and accept the result into "result".
        [~,result] = system(commandLine);
        % Look for our program's name in the result variable.
        itIsRunning = contains(lower(result), lower(vmdexe));
        if itIsRunning
            system(['Taskkill /IM ' vmdexe ' /F'],'-echo');
        end
    else
        commandLine = 'pgrep vmd_LINUXA+';
        [~,result] = system(commandLine);
        if ~isempty(result)
            PID = strtrim(result);
            system(['kill ' PID],'-echo');
        end
    end
    
    % Delete temporary files
    if exist(app.tempfolder,'dir')
        try
            rmdir(app.tempfolder,'s');
        catch
            warning('Unable to delete temp folder.')
        end
    end
    delete(app.figh)
end

function ThirdAxisTextCallback(src,~)
    
    str=get(src,'String');
    if isempty(str2double(str))
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    elseif str2double(str) < 0
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    else
        return
    end
    
end

function TimeWindowCallback(src,~)
    str=get(src,'String');
    if isempty(str2double(str))
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    elseif str2double(str) < 0
        set(src,'string','0');
        warndlg('The intial time point must be a non-negative number.');
    else
        return
    end
    
end

function ThirdAxisCallback(~,~)
    XStr = app.PlotTypeXh.String{app.PlotTypeXh.Value};
    YStr = app.PlotTypeYh.String{app.PlotTypeYh.Value};
    DTStr = app.DataTypeh.String{app.DataTypeh.Value};
    Distr_Sel = ~isempty(regexp(XStr,'RDF','once')) || ~isempty(regexp(YStr,'RDF','once')) ...
        || ~isempty(regexp(XStr,'[Q|W]l','once')) || ~isempty(regexp(YStr,'[Q|W]l','once')) ...
        || ~isempty(regexp(XStr,'MSD','once')) || ~isempty(regexp(YStr,'MSD','once')) || ...
        ~isempty(regexp(XStr,'Near','once')) || ~isempty(regexp(YStr,'Near','once'));
    
    Ql_Sel = ~isempty(regexp(XStr,'[Q|W]l','once')) || ~isempty(regexp(YStr,'[Q|W]l','once'));
    NN_Sel = ~isempty(regexp(XStr,'Near','once')) || ~isempty(regexp(YStr,'Near','once'));
    
    if Ql_Sel
        set(app.Qlbg,'Visible','On');
        set(app.QlNeighbours,'Visible','On');
        set(app.NeighbEdit,'Visible','On');
        set(app.compareQlButton1,'Visible','On')
        set(app.compareQlButton2,'Visible','On')
    elseif NN_Sel
        set(app.Qlbg,'Visible','On');
        set(app.QlNeighbours,'Visible','Off');
        set(app.NeighbEdit,'Visible','Off');
        set(app.compareQlButton1,'Visible','Off')
        set(app.compareQlButton2,'Visible','Off')
    else
        set(app.Qlbg,'Visible','Off');
    end
    
    if strcmpi(DTStr,'MD') && Distr_Sel
        if app.ThirdAxis.Value == 0 % button not pushed
            set(app.ThirdAxis,'Visible','On')
            set(app.ThirdAxisLabelA,'Visible','Off');
            set(app.ThirdAxisLabelB,'Visible','Off');
            set(app.ThirdAxisText,'Visible','Off');
            set(app.TimeWindowText,'Visible','Off');
            set(app.TimeWindow,'Visible','Off');
        else % button pushed
            set(app.ThirdAxis,'Visible','On')
            set(app.ThirdAxisLabelA,'Visible','On');
            set(app.ThirdAxisLabelB,'Visible','On');
            set(app.ThirdAxisText,'Visible','On');
            set(app.TimeWindowText,'Visible','On');
            set(app.TimeWindow,'Visible','On');
        end
    else
        set(app.ThirdAxis,'Visible','Off')
        set(app.ThirdAxisLabelA,'Visible','Off');
        set(app.ThirdAxisLabelB,'Visible','Off');
        set(app.ThirdAxisText,'Visible','Off');
        set(app.TimeWindowText,'Visible','Off');
        set(app.TimeWindow,'Visible','Off');
    end
end

function QlEditCallback(src,~)
    
    First_sel = logical(app.FirstShell.Value);
    Second_sel = logical(app.SecondShell.Value);
    Radius_sel = logical(app.QlRadius.Value);
    Neighbours_sel = logical(app.QlNeighbours.Value);
    Search_Min = str2double(app.QlEditMin.String);
    Search_Max = str2double(app.QlEdit.String);
    Neighbours = str2double(app.NeighbEdit.String);
    
    if (First_sel || Second_sel) && ~strcmp(src.Style,'radiobutton')
        app.QlRadius.Value = 0;
        app.QlNeighbours.Value = 0;
    elseif ~First_sel && ~Second_sel && ~Radius_sel && ~Neighbours_sel
        app.QlRadius.Value = 1; % Default
    end
    
    if strcmp(src.Style,'radiobutton')
        app.FirstShell.Value = 0;
        app.SecondShell.Value = 0;
    end
    
    if isnan(Search_Min)
        set(app.QlEditMin,'string','0');
        Search_Min = 0;
    elseif Search_Min < 0
        set(app.QlEditMin,'string','0');
        Search_Min = 0;
    end

    if isnan(Search_Max)
        set(app.QlEdit,'string',num2str(Search_Min+1));
    elseif Search_Max <= Search_Min
        set(app.QlEdit,'string',num2str(Search_Min+1));
    end

    if isnan(Neighbours)
        set(app.NeighbEdit,'string','6');
    elseif Neighbours <= 0
        set(app.NeighbEdit,'string','6');
    elseif mod(Neighbours,1) ~= 0
        set(app.NeighbEdit,'string',num2str(round(Neighbours)));
    end
end

function QlEditCallbackTraj(src,~)
    
    First_sel = logical(app.TJS.Order_FirstShell.Value);
    Second_sel = logical(app.TJS.Order_SecondShell.Value);
    Radius_sel = logical(app.TJS.Order_QlRadius.Value);
    Neighbours_sel = logical(app.TJS.Order_QlNeighbours.Value);
    Search_Min = str2double(app.TJS.Order_QlEditMin.String);
    Search_Max = str2double(app.TJS.Order_QlEdit.String);
    Neighbours = str2double(app.TJS.Order_NeighbEdit.String);
    
    if (First_sel || Second_sel) && ~strcmp(src.Style,'radiobutton')
        app.TJS.Order_QlRadius.Value = 0;
        app.TJS.Order_QlNeighbours.Value = 0;
    elseif ~First_sel && ~Second_sel && ~Radius_sel && ~Neighbours_sel
        app.TJS.Order_QlRadius.Value = 1; % Default
    end
    
    if strcmp(src.Style,'radiobutton')
        app.TJS.Order_FirstShell.Value = 0;
        app.TJS.Order_SecondShell.Value = 0;
    end
    
    if isnan(Search_Min)
        set(app.TJS.Order_QlEditMin,'string','0');
        Search_Min = 0;
    elseif Search_Min < 0
        set(app.TJS.Order_QlEditMin,'string','0');
        Search_Min = 0;
    end

    if isnan(Search_Max)
        set(app.TJS.Order_QlEdit,'string',num2str(Search_Min+1));
    elseif Search_Max <= Search_Min
        set(app.TJS.Order_QlEdit,'string',num2str(Search_Min+1));
    end

    if isnan(Neighbours)
        set(app.TJS.Order_NeighbEdit,'string','6');
    elseif Neighbours <= 0
        set(app.TJS.Order_NeighbEdit,'string','6');
    elseif mod(Neighbours,1) ~= 0
        set(app.TJS.Order_NeighbEdit,'string',num2str(round(Neighbours)));
    end
end

function NNEditCallbackTraj(src,~)
    First_sel = logical(app.TJS.Neighbour_FirstShell.Value);
    Second_sel = logical(app.TJS.Neighbour_SecondShell.Value);
    Radius_sel = logical(app.TJS.Neighbour_Search_Radius.Value);
    Search_Min = str2double(app.TJS.Neighbour_Search_Radius_Min.String);
    Search_Max = str2double(app.TJS.Neighbour_Search_Radius_Max.String);
    
    if (First_sel || Second_sel) && ~strcmp(src.Style,'radiobutton')
        app.TJS.Neighbour_Search_Radius.Value = 0;
    elseif ~First_sel && ~Second_sel && ~Radius_sel
        app.TJS.Neighbour_Search_Radius.Value = 1; % Default
    end
    
    if strcmp(src.Style,'radiobutton')
        app.TJS.Neighbour_FirstShell.Value = 0;
        app.TJS.Neighbour_SecondShell.Value = 0;
    end
    
    if isnan(Search_Min)
        set(app.TJS.Neighbour_Search_Radius_Min,'string','0');
        Search_Min = 0;
    elseif Search_Min < 0
        set(app.TJS.Neighbour_Search_Radius_Min,'string','0');
        Search_Min = 0;
    end

    if isnan(Search_Max)
        set(app.TJS.Neighbour_Search_Radius_Max,'string',num2str(Search_Min+1));
    elseif Search_Max <= Search_Min
        set(app.TJS.Neighbour_Search_Radius_Max,'string',num2str(Search_Min+1));
    end
end

function PlotChangeCallback(src,~)
    % Save the current plot state
    SavePlotState;

    % Which button was pressed?
    switch src.Tag
        case 'pforward'
            if app.PlotPoint < app.NumPlots
                app.PlotPoint = app.PlotPoint + 1;
            else
                app.PlotPoint = app.NumPlots;
            end
            LoadPlotState;
        case 'pbackward'
            if app.PlotPoint > 1
                app.PlotPoint = app.PlotPoint - 1;
            else
                app.PlotPoint = 1;
            end
            LoadPlotState;
        case 'tforward'
            if app.TimeInUse
                pidx = unique_plot(app);
                time_uniq = app.Time(pidx);
                if app.TimePoint < length(time_uniq{app.PlotPoint})
                    app.TimePoint = app.TimePoint + 1;
                    [caz,cel] = view(app.Axh);
                    app.plotview = [caz cel];
                else
                    return
                end
            else
                app.TimePoint = 0;
            end
        case 'tbackward'
            if app.TimeInUse
                if app.TimePoint > 1
                    app.TimePoint = app.TimePoint - 1;
                    [caz,cel] = view(app.Axh);
                    app.plotview = [caz cel];
                else
                    return
                end
            else
                app.TimePoint = 0;
            end
    end
    RefactorPlots; 
end

function SaveDataCallback(~,~)
    if app.NumPlots == 0
        warndlg('Nothing to save!')
        return
    end
    fn = [regexprep(strrep(strrep(app.legtxt{app.PlotPoint},' ','-'), '\',''),'\$|\^|\{|\}|textrm','') '.dat'];
    
    loc = fullfile(fileparts(mfilename('fullpath')),'Results',fn);
    
    % Function to save the current app object to file
    [filename,PathName] = uiputfile(...
        {'*.dat','Data Analysis Files (*.dat)'},'Save Data File',loc);
    if filename == 0
        return
    end
    
    % Save all data
    SavePlotState
    Data.Time = app.Time;
    Data.PotE = app.PotE;
    Data.KinE = app.KinE;
    Data.TotE = app.TotE;
    Data.Temp = app.Temp;
    Data.Press = app.Press;
    Data.Vol = app.Vol;
    Data.Dens = app.Dens;
    Data.Enth = app.Enth;
    Data.MSD = app.MSD;
    Data.RDF = app.RDF;
    Data.CRDF = app.CRDF;
    Data.NN = app.NN;
    Data.Ql = app.Ql;
    Data.Wl = app.Wl;
    Data.Ql_Dist = app.Ql_Dist;
    Data.Wl_Dist = app.Wl_Dist;
    Data.legtxt = app.legtxt;
    Data.Notes = app.Notes;
    Data.PlotState = app.PlotState;
    
    % Save state of app
    Data.PlotPoint = app.PlotPoint;
    Data.TimePoint = app.TimePoint;
    Data.NumDataSets = app.NumDataSets;
    Data.NumPlots = app.NumPlots;
    Data.TimeInUse = app.TimeInUse;
    Data.data_type = app.data_type;
    Data.energy_type_plot = app.energy_type_plot;
    Data.rdf_type_plot = app.rdf_type_plot;
    Data.msd_type_plot = app.msd_type_plot;
    Data.AvQl_type_plot = app.AvQl_type_plot;
    Data.PlotChangeForward.Enable = app.PlotChangeForward.Enable;
    Data.PlotChangeBackward.Enable = app.PlotChangeBackward.Enable;
    Data.PlotChangeText.String = app.PlotChangeText.String;
    Data.TimeChangeBackward.Visible = app.TimeChangeBackward.Visible;
    Data.TimeChangeBackward.Enable = app.TimeChangeBackward.Enable;
    Data.TimeChangeForward.Visible = app.TimeChangeForward.Visible;
    Data.TimeChangeForward.Enable = app.TimeChangeForward.Enable;
    Data.TimeChangeText.String = app.TimeChangeText.String;
    Data.TimeChangeText.Visible = app.TimeChangeText.Visible;
    Data.SliceTrajectoryh.Visible = app.SliceTrajectoryh.Visible;
    Data.QlEdit.String = app.QlEdit.String;
    Data.QlEditMin.String = app.QlEditMin.String;
    Data.NeighbEdit.String = app.NeighbEdit.String;
    Data.SecondShell.Value = app.SecondShell.Value;
    Data.FirstShell.Value = app.FirstShell.Value;
    Data.QlNeighbours.Value = app.QlNeighbours.Value;
    Data.QlRadius.Value = app.QlRadius.Value;
    Data.Qlbg.Visible = app.Qlbg.Visible;
    Data.QlNeighbours.Visible = app.QlNeighbours.Visible;
    Data.NeighbEdit.Visible = app.NeighbEdit.Visible;
    Data.compareQlButton1.Visible = app.compareQlButton1.Visible;
    Data.compareQlButton2.Visible = app.compareQlButton2.Visible;
    Data.TimeWindowText.String = app.TimeWindowText.String;
    Data.TimeWindowText.Visible = app.TimeWindowText.Visible;
    Data.TimeWindow.Visible = app.TimeWindow.Visible;
    Data.ThirdAxisText.String = app.ThirdAxisText.String;
    Data.ThirdAxisText.Visible = app.ThirdAxisText.Visible;
    Data.ThirdAxisLabelA.Visible = app.ThirdAxisLabelA.Visible;
    Data.ThirdAxisLabelB.Visible = app.ThirdAxisLabelB.Visible;
    Data.ThirdAxis.Value = app.ThirdAxis.Value;
    Data.ThirdAxis.Visible = app.ThirdAxis.Visible;
    Data.PlotTypeYh.String = app.PlotTypeYh.String;
    Data.PlotTypeYh.Value = app.PlotTypeYh.Value;
    Data.PlotTypeXh.String = app.PlotTypeXh.String;
    Data.PlotTypeXh.Value = app.PlotTypeXh.Value;
    Data.TimePerFrameh.String = app.TimePerFrameh.String;
    Data.TimePerFrameh.Visible = app.TimePerFrameh.Visible;
    Data.FinalTimeh.String = app.FinalTimeh.String;
    Data.FinalTimeh.Visible = app.FinalTimeh.Visible;
    Data.InitialTimeh.String = app.InitialTimeh.String;
    Data.InitialTimeh.Visible = app.InitialTimeh.Visible;
    Data.Modelh.String = app.Modelh.String;
    Data.Modelh.Value = app.Modelh.Value;
    Data.Modelfilter_b1.Value = app.Modelfilter_b1.Value;
    Data.Modelfilter_b2.Value = app.Modelfilter_b2.Value;
    Data.Structureh.Value = app.Structureh.Value;
    Data.Structureh.Visible = app.Structureh.Visible;
    Data.Salth.Value = app.Salth.Value;
    Data.DataTypeh.Value = app.DataTypeh.Value;
    Data.BGLabel_3.Visible = app.BGLabel(3).Visible;
    Data.BGLabel_5.Visible = app.BGLabel(5).Visible;
    Data.BGLabel_6.Visible = app.BGLabel(6).Visible;
    Data.BGLabel_7.Visible = app.BGLabel(7).Visible;
    Data.PlotState = app.PlotState;
    try
        unid = unique_plot(app);
        dt = app.data_type(unid);
        if dt(app.PlotPoint) == 0
            [caz,cel] = view(app.Axh);
            Data.plotview = [caz cel];
        else
            Data.plotview = app.plotview;
        end
    catch
        Data.plotview = app.plotview;
    end
    
    fn = fullfile(PathName,filename);
    save(fn,'Data')
end

function LoadDataCallback(~,~)

    loc = fullfile(fileparts(mfilename('fullpath')),'Results');
    % Function to load data from file
    [filenames,PathName] = uigetfile(...
        {'*.dat','Data Analysis Files (*.dat)';'*.*',...
        'All Files (*.*)'},'Select RDF Quick Analysis Data File(s)',...
        loc,'MultiSelect','on');
    
    if ~iscell(filenames)
        filenames = {filenames};
    end
    
    for ifn = 1:length(filenames)
        filename = filenames{ifn};
        
        if filename == 0
            continue
        end

        fn = fullfile(PathName,filename);
        try
            load(fn,'Data','-mat')
            Data.PlotState = Update_PlotState_Version(app,Data.PlotState);
            Data = Update_PlotState_Version(app,Data);
        catch
            warndlg(['Error loading selected data file: ' filename '. Data not received.'])
            continue
        end

        if isempty(Data.Time) % no data in file
            warndlg(['Error loading selected data file: ' filename '. No data found in file.'])
            continue
        else
            % Load data set data to app
            for idx = 1:Data.NumDataSets
                app.Time{end+1} = Data.Time{idx};
                app.PotE{end+1} = Data.PotE{idx};
                app.KinE{end+1} = Data.KinE{idx};
                app.TotE{end+1} = Data.TotE{idx};
                app.Temp{end+1} = Data.Temp{idx};
                app.Press{end+1} = Data.Press{idx};
                app.Vol{end+1} = Data.Vol{idx};
                app.Dens{end+1} = Data.Dens{idx};
                app.Enth{end+1} = Data.Enth{idx};
                app.MSD{end+1} = Data.MSD{idx};
                app.RDF{end+1} = Data.RDF{idx};
                app.CRDF{end+1} = Data.CRDF{idx};
                try
                    app.NN{end+1} = Data.NN{idx};
                catch
                    app.NN{end+1} = nan(1,2);
                end
                app.Ql{end+1} = Data.Ql{idx};
                app.Wl{end+1} = Data.Wl{idx};
                app.Ql_Dist{end+1} = Data.Ql_Dist{idx};
                app.Wl_Dist{end+1} = Data.Wl_Dist{idx};
                app.legtxt{end+1} = Data.legtxt{idx};
                app.data_type(end+1) = Data.data_type(idx);
            end

            if ~isfield(Data,'NN_type_plot') % backwards compatibility
                Data.NN_type_plot = nan;
            end
            
            for idx = 1:Data.NumPlots

                % energy type plot selected
                if idx == Data.energy_type_plot

                    if isnan(app.energy_type_plot)
                        % Assign a new energy_type plot if one does not exist
                        app.NumPlots = app.NumPlots + 1;
                        app.energy_type_plot = app.NumPlots;
                        app.Notes{end+1} = Data.Notes{idx}; % assign saved notes to this new plot

                        if isempty(fieldnames(app.PlotState)) || isempty(app.PlotState)
                            app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        else
                            app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        end

                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NumPlots;
                        end
                        
                    else
                        % Append saved energy_type plot notes to existing notes
                        app.Notes{app.energy_type_plot} = [app.Notes{app.energy_type_plot} Data.Notes{idx}];
                        if Data.PlotPoint == idx
                            app.PlotPoint = app.energy_type_plot;
                        end
                    end

                % RDF type plot selected
                elseif idx == Data.rdf_type_plot

                    if isnan(app.rdf_type_plot)
                        % Assign a new rdf_type plot if one does not exist
                        app.NumPlots = app.NumPlots + 1;
                        app.rdf_type_plot = app.NumPlots;
                        app.Notes{end+1} = Data.Notes{idx}; % assign saved notes to this new plot

                        if isempty(fieldnames(app.PlotState))
                            app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        else
                            app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        end

                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NumPlots;
                        end
                    else
                        % Append saved rdf_type plot notes to existing notes
                        app.Notes{app.rdf_type_plot} = [app.Notes{app.rdf_type_plot} Data.Notes{idx}];
                        if Data.PlotPoint == idx
                            app.PlotPoint = app.rdf_type_plot;
                        end
                    end

                % MSD type plot selected    
                elseif idx == Data.msd_type_plot

                    if isnan(app.msd_type_plot)

                        % Assign a new msd_type plot if one does not exist
                        app.NumPlots = app.NumPlots + 1;
                        app.msd_type_plot = app.NumPlots;
                        app.Notes{end+1} = Data.Notes{idx}; % assign saved notes to this new plot

                        if isempty(fieldnames(app.PlotState))
                            app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        else
                            app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        end

                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NumPlots;
                        end
                    else
                        % Append saved msd_type plot notes to existing notes
                        app.Notes{app.msd_type_plot} = [app.Notes{app.msd_type_plot} Data.Notes{idx}];
                        if Data.PlotPoint == idx
                            app.PlotPoint = app.msd_type_plot;
                        end
                    end
                % NN type plot selected   
                elseif idx == Data.NN_type_plot

                    if isnan(app.NN_type_plot)

                        % Assign a new AvQl_type plot if one does not exist
                        app.NumPlots = app.NumPlots + 1;
                        app.NN_type_plot = app.NumPlots;
                        app.Notes{end+1} = Data.Notes{idx}; % assign saved notes to this new plot

                        if isempty(fieldnames(app.PlotState))
                            app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        else
                            app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        end

                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NumPlots;
                        end
                    else
                        % Append saved NN_type plot notes to existing notes
                        app.Notes{app.NN_type_plot} = [app.Notes{app.NN_type_plot} Data.Notes{idx}];
                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NN_type_plot;
                        end
                    end
                % Average Ql type plot selected   
                elseif idx == Data.AvQl_type_plot

                    if isnan(app.AvQl_type_plot)

                        % Assign a new AvQl_type plot if one does not exist
                        app.NumPlots = app.NumPlots + 1;
                        app.AvQl_type_plot = app.NumPlots;
                        app.Notes{end+1} = Data.Notes{idx}; % assign saved notes to this new plot

                        if isempty(fieldnames(app.PlotState))
                            app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        else
                            app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                        end

                        if Data.PlotPoint == idx
                            app.PlotPoint = app.NumPlots;
                        end
                    else
                        % Append saved AvQl_type plot notes to existing notes
                        app.Notes{app.AvQl_type_plot} = [app.Notes{app.AvQl_type_plot} Data.Notes{idx}];
                        if Data.PlotPoint == idx
                            app.PlotPoint = app.AvQl_type_plot;
                        end
                    end

                % Other type of plot
                else
                    app.NumPlots = app.NumPlots + 1;
                    app.Notes{end+1} = Data.Notes{idx};
                    if isempty(fieldnames(app.PlotState)) || isempty(app.PlotState)
                        app.PlotState = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                    else
                    	app.PlotState(end+1) = Data.PlotState(idx); % assign the saved state of the plot to this new plot
                    end
                    if Data.PlotPoint == idx
                        app.PlotPoint = app.NumPlots;
                    end
                end
            end

        end
        app.TimePoint = Data.TimePoint;
        app.TimeInUse = Data.TimeInUse;
        app.NumDataSets = app.NumDataSets + Data.NumDataSets;
        app.PlotChangeForward.Enable = Data.PlotChangeForward.Enable;
        app.PlotChangeBackward.Enable = Data.PlotChangeBackward.Enable;
        app.PlotChangeText.String = Data.PlotChangeText.String;
        app.TimeChangeBackward.Visible = Data.TimeChangeBackward.Visible;
        app.TimeChangeBackward.Enable = Data.TimeChangeBackward.Enable;
        app.TimeChangeForward.Visible = Data.TimeChangeForward.Visible;
        app.TimeChangeForward.Enable = Data.TimeChangeForward.Enable;
        app.TimeChangeText.String = Data.TimeChangeText.String;
        app.TimeChangeText.Visible = Data.TimeChangeText.Visible;
        try
            app.SliceTrajectoryh.Visible = Data.SliceTrajectoryh.Visible;
        catch
            app.SliceTrajectoryh.Visible = 'On';
        end
        
        app.QlEdit.String = Data.QlEdit.String;
        try % Backwards compatiblity
            app.QlEditMin.String = Data.QlEditMin.String;
            app.NeighbEdit.String = Data.NeighbEdit.String;
            app.QlNeighbours.Visible = Data.QlNeighbours.Visible;
            app.NeighbEdit.Visible = Data.NeighbEdit.Visible;
            app.compareQlButton1.Visible = Data.compareQlButton1.Visible;
            app.compareQlButton2.Visible = Data.compareQlButton2.Visible;
        catch
            app.QlEditMin.String = '0';
            app.NeighbEdit.String = '6';
            app.QlNeighbours.Visible = 'On';
            app.NeighbEdit.Visible = 'On';
            app.compareQlButton1.Visible = 'On';
            app.compareQlButton2.Visible = 'On';
        end
        app.SecondShell.Value = Data.SecondShell.Value;
        app.FirstShell.Value = Data.FirstShell.Value;
        app.QlNeighbours.Value = Data.QlNeighbours.Value;
        app.QlRadius.Value = Data.QlRadius.Value;
        app.Qlbg.Visible = Data.Qlbg.Visible;
        app.TimeWindowText.String = Data.TimeWindowText.String;
        app.TimeWindowText.Visible = Data.TimeWindowText.Visible;
        app.TimeWindow.Visible = Data.TimeWindow.Visible;
        app.ThirdAxisText.String = Data.ThirdAxisText.String;
        app.ThirdAxisText.Visible = Data.ThirdAxisText.Visible;
        app.ThirdAxisLabelA.Visible = Data.ThirdAxisLabelA.Visible;
        app.ThirdAxisLabelB.Visible = Data.ThirdAxisLabelB.Visible;
        app.ThirdAxis.Value = Data.ThirdAxis.Value;
        app.ThirdAxis.Visible = Data.ThirdAxis.Visible;
        app.PlotTypeYh.String = Data.PlotTypeYh.String;
        app.PlotTypeYh.Value = Data.PlotTypeYh.Value;
        app.PlotTypeXh.String = Data.PlotTypeXh.String;
        app.PlotTypeXh.Value = Data.PlotTypeXh.Value;
        app.TimePerFrameh.String = Data.TimePerFrameh.String;
        app.TimePerFrameh.Visible = Data.TimePerFrameh.Visible;
        app.FinalTimeh.String = Data.FinalTimeh.String;
        app.FinalTimeh.Visible = Data.FinalTimeh.Visible;
        app.InitialTimeh.String = Data.InitialTimeh.String;
        app.InitialTimeh.Visible = Data.InitialTimeh.Visible;
        app.Modelh.String = Data.Modelh.String;
        app.Modelh.Value = Data.Modelh.Value;
        app.Modelfilter_b1.Value = Data.Modelfilter_b1.Value;
        app.Modelfilter_b2.Value = Data.Modelfilter_b2.Value;
        app.Structureh.Value = Data.Structureh.Value;
        app.Structureh.Visible = Data.Structureh.Visible;
        app.Salth.Value = Data.Salth.Value;
        app.DataTypeh.Value = Data.DataTypeh.Value;
        app.BGLabel(3).Visible = Data.BGLabel_3.Visible;
        app.BGLabel(5).Visible = Data.BGLabel_5.Visible;
        app.BGLabel(6).Visible = Data.BGLabel_6.Visible;
        app.BGLabel(7).Visible = Data.BGLabel_7.Visible;
        app.plotview = Data.plotview;
    end
    LoadPlotState;
    RefactorPlots; % Load data into plots
end

function SaveAppState
    SavePlotState
    filename = fullfile(app.AppData,'savestate.mat');
    savestate = struct();
	savestate.PlotPoint = app.PlotPoint;
	savestate.TimePoint = app.TimePoint;
    savestate.TimeInUse = app.TimeInUse;
	savestate.NumDataSets = app.NumDataSets;
    savestate.data_type = app.data_type;
    savestate.energy_type_plot = app.energy_type_plot;
    savestate.rdf_type_plot = app.rdf_type_plot;
    savestate.msd_type_plot = app.msd_type_plot;
    savestate.AvQl_type_plot = app.AvQl_type_plot;
    savestate.NumPlots = app.NumPlots;
	savestate.XDat = app.XDat;
	savestate.YDat = app.YDat;
	savestate.ZDat = app.ZDat;
	savestate.Time = app.Time;
	savestate.PotE = app.PotE;
	savestate.KinE = app.KinE;
	savestate.TotE = app.TotE;
	savestate.Temp = app.Temp;
	savestate.Press = app.Press;
	savestate.Vol = app.Vol;
	savestate.Dens = app.Dens;
	savestate.Enth = app.Enth;
    savestate.MSD = app.MSD;
	savestate.RDF = app.RDF;
	savestate.CRDF = app.CRDF;
    savestate.NN = app.NN;
	savestate.Ql = app.Ql;
	savestate.Wl = app.Wl;
	savestate.Ql_Dist = app.Ql_Dist;
	savestate.Wl_Dist = app.Wl_Dist;
    savestate.legtxt = app.legtxt;
    savestate.Notes = app.Notes;
	savestate.PlotChangeForward.Enable = app.PlotChangeForward.Enable;
	savestate.PlotChangeBackward.Enable = app.PlotChangeBackward.Enable;
	savestate.PlotChangeText.String = app.PlotChangeText.String;
	savestate.TimeChangeBackward.Visible = app.TimeChangeBackward.Visible;
	savestate.TimeChangeBackward.Enable = app.TimeChangeBackward.Enable;
	savestate.TimeChangeForward.Visible = app.TimeChangeForward.Visible;
	savestate.TimeChangeForward.Enable = app.TimeChangeForward.Enable;
	savestate.TimeChangeText.String = app.TimeChangeText.String;
	savestate.TimeChangeText.Visible = app.TimeChangeText.Visible;
	savestate.SliceTrajectoryh.Visible = app.SliceTrajectoryh.Visible;
	savestate.QlEdit.String = app.QlEdit.String;
    savestate.QlEditMin.String = app.QlEditMin.String;
    savestate.NeighbEdit.String = app.NeighbEdit.String;
	savestate.SecondShell.Value = app.SecondShell.Value;
	savestate.FirstShell.Value = app.FirstShell.Value;
	savestate.QlNeighbours.Value = app.QlNeighbours.Value;
	savestate.QlRadius.Value = app.QlRadius.Value;
	savestate.Qlbg.Visible = app.Qlbg.Visible;
	savestate.QlNeighbours.Visible = app.QlNeighbours.Visible;
	savestate.NeighbEdit.Visible = app.NeighbEdit.Visible;
	savestate.compareQlButton1.Visible = app.compareQlButton1.Visible;
	savestate.compareQlButton2.Visible = app.compareQlButton2.Visible;
	savestate.TimeWindowText.String = app.TimeWindowText.String;
	savestate.TimeWindowText.Visible = app.TimeWindowText.Visible;
	savestate.TimeWindow.Visible = app.TimeWindow.Visible;
	savestate.ThirdAxisText.String = app.ThirdAxisText.String;
	savestate.ThirdAxisText.Visible = app.ThirdAxisText.Visible;
	savestate.ThirdAxisLabelA.Visible = app.ThirdAxisLabelA.Visible;
	savestate.ThirdAxisLabelB.Visible = app.ThirdAxisLabelB.Visible;
	savestate.ThirdAxis.Value = app.ThirdAxis.Value;
	savestate.ThirdAxis.Visible = app.ThirdAxis.Visible;
	savestate.PlotTypeYh.String = app.PlotTypeYh.String;
	savestate.PlotTypeYh.Value = app.PlotTypeYh.Value;
	savestate.PlotTypeXh.String = app.PlotTypeXh.String;
	savestate.PlotTypeXh.Value = app.PlotTypeXh.Value;
	savestate.TimePerFrameh.String = app.TimePerFrameh.String;
	savestate.TimePerFrameh.Visible = app.TimePerFrameh.Visible;
	savestate.FinalTimeh.String = app.FinalTimeh.String;
	savestate.FinalTimeh.Visible = app.FinalTimeh.Visible;
	savestate.InitialTimeh.String = app.InitialTimeh.String;
	savestate.InitialTimeh.Visible = app.InitialTimeh.Visible;
	savestate.Modelh.String = app.Modelh.String;
	savestate.Modelh.Value = app.Modelh.Value;
    savestate.Modelfilter_b1.Value = app.Modelfilter_b1.Value;
    savestate.Modelfilter_b2.Value = app.Modelfilter_b2.Value;
	savestate.Structureh.Value = app.Structureh.Value;
	savestate.Structureh.Visible = app.Structureh.Visible;
	savestate.Salth.Value = app.Salth.Value;
	savestate.DataTypeh.Value = app.DataTypeh.Value;
	savestate.BGLabel_3.Visible = app.BGLabel(3).Visible;
	savestate.BGLabel_5.Visible = app.BGLabel(5).Visible;
	savestate.BGLabel_6.Visible = app.BGLabel(6).Visible;
	savestate.BGLabel_7.Visible = app.BGLabel(7).Visible;
    savestate.PlotState = app.PlotState;
    try
        unid = unique_plot(app);
        dt = app.data_type(unid);
        if dt(app.PlotPoint) == 0
            [caz,cel] = view(app.Axh);
            savestate.plotview = [caz cel];
        else
            savestate.plotview = app.plotview;
        end
    catch
        savestate.plotview = app.plotview;
    end
    
    save(filename,'savestate');
end

function app = LoadAppState(app,filename)
    savedata = open(filename);
    savestate = savedata.savestate;
    app.PlotPoint = savestate.PlotPoint;
    app.TimePoint = savestate.TimePoint;
    app.TimeInUse = savestate.TimeInUse;
    app.NumDataSets = savestate.NumDataSets;
    app.data_type = savestate.data_type;
    app.energy_type_plot = savestate.energy_type_plot;
    app.rdf_type_plot = savestate.rdf_type_plot;
    app.msd_type_plot = savestate.msd_type_plot;
    app.AvQl_type_plot = savestate.AvQl_type_plot;
    app.NumPlots = savestate.NumPlots;
    app.XDat = savestate.XDat;
    app.YDat = savestate.YDat;
    app.ZDat = savestate.ZDat;
    app.Time = savestate.Time;
    app.PotE = savestate.PotE;
    app.KinE = savestate.KinE;
    app.TotE = savestate.TotE;
    app.Temp = savestate.Temp;
    app.Press = savestate.Press;
    app.Vol = savestate.Vol;
    app.Dens = savestate.Dens;
    app.Enth = savestate.Enth;
    app.MSD = savestate.MSD;
    app.RDF = savestate.RDF;
    app.CRDF = savestate.CRDF;
    app.NN = savestate.NN;
    app.Ql = savestate.Ql;
    app.Wl = savestate.Wl;
    app.Ql_Dist = savestate.Ql_Dist;
    app.Wl_Dist = savestate.Wl_Dist;
    app.legtxt = savestate.legtxt;
    app.Notes = savestate.Notes;
    app.PlotChangeForward.Enable = savestate.PlotChangeForward.Enable;
    app.PlotChangeBackward.Enable = savestate.PlotChangeBackward.Enable;
    app.PlotChangeText.String = savestate.PlotChangeText.String;
    app.TimeChangeBackward.Visible = savestate.TimeChangeBackward.Visible;
    app.TimeChangeBackward.Enable = savestate.TimeChangeBackward.Enable;
    app.TimeChangeForward.Visible = savestate.TimeChangeForward.Visible;
    app.TimeChangeForward.Enable = savestate.TimeChangeForward.Enable;
    app.TimeChangeText.String = savestate.TimeChangeText.String;
    app.TimeChangeText.Visible = savestate.TimeChangeText.Visible;
    try
        app.SliceTrajectoryh.Visible = savestate.SliceTrajectoryh.Visible;
    catch
        app.SliceTrajectoryh.Visible = savestate.ViewTrajectoryh.Visible;
    end
    app.QlEdit.String = savestate.QlEdit.String;
    try % backwards compatible
        app.QlEditMin.String = savestate.QlEditMin.String;
        app.NeighbEdit.String = savestate.NeighbEdit.String;
        app.QlNeighbours.Visible = savestate.QlNeighbours.Visible;
        app.NeighbEdit.Visible = savestate.NeighbEdit.Visible;
        app.compareQlButton1.Visible = savestate.compareQlButton1.Visible;
        app.compareQlButton2.Visible = savestate.compareQlButton2.Visible;
    catch
        app.QlEditMin.String = '0';
        app.NeighbEdit.String = '6';
        app.QlNeighbours.Visible = 'On';
        app.NeighbEdit.Visible = 'On';
        app.compareQlButton1.Visible = 'On';
        app.compareQlButton2.Visible = 'On';
    end
    app.SecondShell.Value = savestate.SecondShell.Value;
    app.FirstShell.Value = savestate.FirstShell.Value;
    app.QlNeighbours.Value = savestate.QlNeighbours.Value;
    app.QlRadius.Value = savestate.QlRadius.Value;
    app.Qlbg.Visible = savestate.Qlbg.Visible;
    app.TimeWindowText.String = savestate.TimeWindowText.String;
    app.TimeWindowText.Visible = savestate.TimeWindowText.Visible;
    app.TimeWindow.Visible = savestate.TimeWindow.Visible;
    app.ThirdAxisText.String = savestate.ThirdAxisText.String;
    app.ThirdAxisText.Visible = savestate.ThirdAxisText.Visible;
    app.ThirdAxisLabelA.Visible = savestate.ThirdAxisLabelA.Visible;
    app.ThirdAxisLabelB.Visible = savestate.ThirdAxisLabelB.Visible;
    app.ThirdAxis.Value = savestate.ThirdAxis.Value;
    app.ThirdAxis.Visible = savestate.ThirdAxis.Visible;
    app.PlotTypeYh.String = savestate.PlotTypeYh.String;
    app.PlotTypeYh.Value = savestate.PlotTypeYh.Value;
    app.PlotTypeXh.String = savestate.PlotTypeXh.String;
    app.PlotTypeXh.Value = savestate.PlotTypeXh.Value;
    app.TimePerFrameh.String = savestate.TimePerFrameh.String;
    app.TimePerFrameh.Visible = savestate.TimePerFrameh.Visible;
    app.FinalTimeh.String = savestate.FinalTimeh.String;
    app.FinalTimeh.Visible = savestate.FinalTimeh.Visible;
    app.InitialTimeh.String = savestate.InitialTimeh.String;
    app.InitialTimeh.Visible = savestate.InitialTimeh.Visible;
    app.Modelh.String = savestate.Modelh.String;
    app.Modelh.Value = savestate.Modelh.Value;
    app.Modelfilter_b1.Value = savestate.Modelfilter_b1.Value;
    app.Modelfilter_b2.Value = savestate.Modelfilter_b2.Value;
    app.Structureh.Value = savestate.Structureh.Value;
    app.Structureh.Visible = savestate.Structureh.Visible;
    app.Salth.Value = savestate.Salth.Value;
    app.DataTypeh.Value = savestate.DataTypeh.Value;
    app.BGLabel(3).Visible = savestate.BGLabel_3.Visible;
    app.BGLabel(5).Visible = savestate.BGLabel_5.Visible;
    app.BGLabel(6).Visible = savestate.BGLabel_6.Visible;
    app.BGLabel(7).Visible = savestate.BGLabel_7.Visible;
    app.plotview = savestate.plotview;
    app.PlotState = savestate.PlotState;
end

function SavePlotState
    if app.NumPlots > 0
        app.PlotState(app.PlotPoint).PlotChangeForward.Enable = app.PlotChangeForward.Enable;
        app.PlotState(app.PlotPoint).PlotChangeBackward.Enable = app.PlotChangeBackward.Enable;
        app.PlotState(app.PlotPoint).PlotChangeText.String = app.PlotChangeText.String;
        app.PlotState(app.PlotPoint).TimeChangeBackward.Visible = app.TimeChangeBackward.Visible;
        app.PlotState(app.PlotPoint).TimeChangeBackward.Enable = app.TimeChangeBackward.Enable;
        app.PlotState(app.PlotPoint).TimeChangeForward.Visible = app.TimeChangeForward.Visible;
        app.PlotState(app.PlotPoint).TimeChangeForward.Enable = app.TimeChangeForward.Enable;
        app.PlotState(app.PlotPoint).TimeChangeText.String = app.TimeChangeText.String;
        app.PlotState(app.PlotPoint).TimeChangeText.Visible = app.TimeChangeText.Visible;
        app.PlotState(app.PlotPoint).SliceTrajectoryh.Visible = app.SliceTrajectoryh.Visible;
        app.PlotState(app.PlotPoint).QlEdit.String = app.QlEdit.String;
        app.PlotState(app.PlotPoint).QlEditMin.String = app.QlEditMin.String;
        app.PlotState(app.PlotPoint).NeighbEdit.String = app.NeighbEdit.String;
        app.PlotState(app.PlotPoint).SecondShell.Value = app.SecondShell.Value;
        app.PlotState(app.PlotPoint).FirstShell.Value = app.FirstShell.Value;
        app.PlotState(app.PlotPoint).QlNeighbours.Value = app.QlNeighbours.Value;
        app.PlotState(app.PlotPoint).QlRadius.Value = app.QlRadius.Value;
        app.PlotState(app.PlotPoint).Qlbg.Visible = app.Qlbg.Visible;
        app.PlotState(app.PlotPoint).QlNeighbours.Visible = app.QlNeighbours.Visible;
        app.PlotState(app.PlotPoint).NeighbEdit.Visible = app.NeighbEdit.Visible;
        app.PlotState(app.PlotPoint).compareQlButton1.Visible = app.compareQlButton1.Visible;
        app.PlotState(app.PlotPoint).compareQlButton2.Visible = app.compareQlButton2.Visible;
        app.PlotState(app.PlotPoint).TimeWindowText.String = app.TimeWindowText.String;
        app.PlotState(app.PlotPoint).TimeWindowText.Visible = app.TimeWindowText.Visible;
        app.PlotState(app.PlotPoint).TimeWindow.Visible = app.TimeWindow.Visible;
        app.PlotState(app.PlotPoint).ThirdAxisText.String = app.ThirdAxisText.String;
        app.PlotState(app.PlotPoint).ThirdAxisText.Visible = app.ThirdAxisText.Visible;
        app.PlotState(app.PlotPoint).ThirdAxisLabelA.Visible = app.ThirdAxisLabelA.Visible;
        app.PlotState(app.PlotPoint).ThirdAxisLabelB.Visible = app.ThirdAxisLabelB.Visible;
        app.PlotState(app.PlotPoint).ThirdAxis.Value = app.ThirdAxis.Value;
        app.PlotState(app.PlotPoint).ThirdAxis.Visible = app.ThirdAxis.Visible;
        app.PlotState(app.PlotPoint).PlotTypeYh.String = app.PlotTypeYh.String;
        app.PlotState(app.PlotPoint).PlotTypeYh.Value = app.PlotTypeYh.Value;
        app.PlotState(app.PlotPoint).PlotTypeXh.String = app.PlotTypeXh.String;
        app.PlotState(app.PlotPoint).PlotTypeXh.Value = app.PlotTypeXh.Value;
        app.PlotState(app.PlotPoint).TimePerFrameh.String = app.TimePerFrameh.String;
        app.PlotState(app.PlotPoint).TimePerFrameh.Visible = app.TimePerFrameh.Visible;
        app.PlotState(app.PlotPoint).FinalTimeh.String = app.FinalTimeh.String;
        app.PlotState(app.PlotPoint).FinalTimeh.Visible = app.FinalTimeh.Visible;
        app.PlotState(app.PlotPoint).InitialTimeh.String = app.InitialTimeh.String;
        app.PlotState(app.PlotPoint).InitialTimeh.Visible = app.InitialTimeh.Visible;
        app.PlotState(app.PlotPoint).Modelh.String = app.Modelh.String;
        app.PlotState(app.PlotPoint).Modelh.Value = app.Modelh.Value;
        app.PlotState(app.PlotPoint).Modelfilter_b1.Value = app.Modelfilter_b1.Value;
        app.PlotState(app.PlotPoint).Modelfilter_b2.Value = app.Modelfilter_b2.Value;
        app.PlotState(app.PlotPoint).Structureh.Value = app.Structureh.Value;
        app.PlotState(app.PlotPoint).Structureh.Visible = app.Structureh.Visible;
        app.PlotState(app.PlotPoint).Salth.Value = app.Salth.Value;
        app.PlotState(app.PlotPoint).DataTypeh.Value = app.DataTypeh.Value;
        app.PlotState(app.PlotPoint).BGLabel_3.Visible = app.BGLabel(3).Visible;
        app.PlotState(app.PlotPoint).BGLabel_5.Visible = app.BGLabel(5).Visible;
        app.PlotState(app.PlotPoint).BGLabel_6.Visible = app.BGLabel(6).Visible;
        app.PlotState(app.PlotPoint).BGLabel_7.Visible = app.BGLabel(7).Visible;
        [caz,cel] = view(app.Axh);
        app.PlotState(app.PlotPoint).plotview = [caz cel];
        app.PlotState(app.PlotPoint).TimePoint = app.TimePoint;
        app.PlotState(app.PlotPoint).Notepad.Visible = app.Notepad.Visible;
    end
end

function LoadPlotState
    if length(app.PlotState) >= app.PlotPoint && ~isempty(app.PlotState) && ~isempty(fieldnames(app.PlotState))
        if ~isempty(app.PlotState(app.PlotPoint).PlotChangeForward)
            app.PlotChangeForward.Enable = app.PlotState(app.PlotPoint).PlotChangeForward.Enable;
            app.PlotChangeBackward.Enable = app.PlotState(app.PlotPoint).PlotChangeBackward.Enable;
            app.PlotChangeText.String = app.PlotState(app.PlotPoint).PlotChangeText.String;
            app.TimeChangeBackward.Visible = app.PlotState(app.PlotPoint).TimeChangeBackward.Visible;
            app.TimeChangeBackward.Enable = app.PlotState(app.PlotPoint).TimeChangeBackward.Enable;
            app.TimeChangeForward.Visible = app.PlotState(app.PlotPoint).TimeChangeForward.Visible;
            app.TimeChangeForward.Enable = app.PlotState(app.PlotPoint).TimeChangeForward.Enable;
            app.TimeChangeText.String = app.PlotState(app.PlotPoint).TimeChangeText.String;
            app.TimeChangeText.Visible = app.PlotState(app.PlotPoint).TimeChangeText.Visible;
            try
                app.SliceTrajectoryh.Visible = app.PlotState(app.PlotPoint).SliceTrajectoryh.Visible;
            catch
                app.SliceTrajectoryh.Visible = 1;
            end
            app.QlEdit.String = app.PlotState(app.PlotPoint).QlEdit.String;
            app.SecondShell.Value = app.PlotState(app.PlotPoint).SecondShell.Value;
            app.FirstShell.Value = app.PlotState(app.PlotPoint).FirstShell.Value;
            app.QlNeighbours.Value = app.PlotState(app.PlotPoint).QlNeighbours.Value;
            app.QlRadius.Value = app.PlotState(app.PlotPoint).QlRadius.Value;
            app.Qlbg.Visible = app.PlotState(app.PlotPoint).Qlbg.Visible;
            try
                app.QlNeighbours.Visible = app.PlotState(app.PlotPoint).QlNeighbours.Visible;
                app.NeighbEdit.Visible = app.PlotState(app.PlotPoint).NeighbEdit.Visible;
                app.compareQlButton1.Visible = app.PlotState(app.PlotPoint).compareQlButton1.Visible;
                app.compareQlButton2.Visible = app.PlotState(app.PlotPoint).compareQlButton2.Visible;
            catch
                app.QlNeighbours.Visible = 'On';
                app.NeighbEdit.Visible = 'On';
                app.compareQlButton1.Visible = 'On';
                app.compareQlButton2.Visible = 'On';
            end
            app.TimeWindowText.String = app.PlotState(app.PlotPoint).TimeWindowText.String;
            app.TimeWindowText.Visible = app.PlotState(app.PlotPoint).TimeWindowText.Visible;
            app.TimeWindow.Visible = app.PlotState(app.PlotPoint).TimeWindow.Visible;
            app.ThirdAxisText.String = app.PlotState(app.PlotPoint).ThirdAxisText.String;
            app.ThirdAxisText.Visible = app.PlotState(app.PlotPoint).ThirdAxisText.Visible;
            app.ThirdAxisLabelA.Visible = app.PlotState(app.PlotPoint).ThirdAxisLabelA.Visible;
            app.ThirdAxisLabelB.Visible = app.PlotState(app.PlotPoint).ThirdAxisLabelB.Visible;
            app.ThirdAxis.Value = app.PlotState(app.PlotPoint).ThirdAxis.Value;
            app.ThirdAxis.Visible = app.PlotState(app.PlotPoint).ThirdAxis.Visible;
            app.PlotTypeYh.String = app.PlotState(app.PlotPoint).PlotTypeYh.String;
            app.PlotTypeYh.Value = app.PlotState(app.PlotPoint).PlotTypeYh.Value;
            app.PlotTypeXh.String = app.PlotState(app.PlotPoint).PlotTypeXh.String;
            app.PlotTypeXh.Value = app.PlotState(app.PlotPoint).PlotTypeXh.Value;
            app.TimePerFrameh.String = app.PlotState(app.PlotPoint).TimePerFrameh.String;
            app.TimePerFrameh.Visible = app.PlotState(app.PlotPoint).TimePerFrameh.Visible;
            app.FinalTimeh.String = app.PlotState(app.PlotPoint).FinalTimeh.String;
            app.FinalTimeh.Visible = app.PlotState(app.PlotPoint).FinalTimeh.Visible;
            app.InitialTimeh.String = app.PlotState(app.PlotPoint).InitialTimeh.String;
            app.InitialTimeh.Visible = app.PlotState(app.PlotPoint).InitialTimeh.Visible;
            app.Modelh.String = app.PlotState(app.PlotPoint).Modelh.String;
            app.Modelh.Value = app.PlotState(app.PlotPoint).Modelh.Value;
            app.Modelfilter_b1.Value = app.PlotState(app.PlotPoint).Modelfilter_b1.Value;
            app.Modelfilter_b2.Value = app.PlotState(app.PlotPoint).Modelfilter_b2.Value;
            app.Structureh.Value = app.PlotState(app.PlotPoint).Structureh.Value;
            app.Structureh.Visible = app.PlotState(app.PlotPoint).Structureh.Visible;
            app.Salth.Value = app.PlotState(app.PlotPoint).Salth.Value;
            app.DataTypeh.Value = app.PlotState(app.PlotPoint).DataTypeh.Value;
            app.BGLabel(3).Visible = app.PlotState(app.PlotPoint).BGLabel_3.Visible;
            app.BGLabel(5).Visible = app.PlotState(app.PlotPoint).BGLabel_5.Visible;
            app.BGLabel(6).Visible = app.PlotState(app.PlotPoint).BGLabel_6.Visible;
            app.BGLabel(7).Visible = app.PlotState(app.PlotPoint).BGLabel_7.Visible;
            app.plotview = app.PlotState(app.PlotPoint).plotview;
            app.TimePoint = app.PlotState(app.PlotPoint).TimePoint;
            app.Notepad.Visible = app.PlotState(app.PlotPoint).Notepad.Visible;
        end
    end
    DataTypeCallback;
end

function RefreshGUI
    if app.TimeInUse
        pidx = unique_plot(app);
        time_uniq = app.Time(pidx);
        
        app.TimeChangeText.String = ['Time: ' num2str(time_uniq{app.PlotPoint}(app.TimePoint)) ' ns'];
        app.TimeChangeText.Visible = 'On';
        app.TimeChangeForward.Visible = 'On';
        app.TimeChangeBackward.Visible = 'On';

        if app.TimePoint <= 1
            set(app.TimeChangeBackward,'Enable','Off')
        else
            set(app.TimeChangeBackward,'Enable','On')
        end

        if app.TimePoint >= length(time_uniq{app.PlotPoint})
            set(app.TimeChangeForward,'Enable','Off')
        else
            set(app.TimeChangeForward,'Enable','On')
        end
    else
        app.TimeChangeText.Visible = 'Off';
        app.TimeChangeForward.Visible = 'Off';
        app.TimeChangeBackward.Visible = 'Off';
    end

    app.PlotChangeText.String = ['Plot: ' num2str(app.PlotPoint) ' of ' num2str(app.NumPlots)];

    if app.PlotPoint >= app.NumPlots
        set(app.PlotChangeForward,'Enable','Off')
    else
        set(app.PlotChangeForward,'Enable','On')
    end
    if app.PlotPoint <= 1
        set(app.PlotChangeBackward,'Enable','Off')
    else
        set(app.PlotChangeBackward,'Enable','On')
    end
    
    if app.compareQlPlots
        set(app.bg,'Visible','Off')
        set(app.compareQlButton1b,'Visible','On')
        set(app.compareQlButton2b,'Visible','On')
    else
        set(app.bg,'Visible','On')
        set(app.compareQlButton1b,'Visible','Off')
        set(app.compareQlButton2b,'Visible','Off')
    end
end

function ReloadAxis
    
    if any(strcmp('Axh', fieldnames(app)))
        delete(app.Axh)
    end
    
    if any(strcmp('Legend', fieldnames(app)))
        delete(app.Legend)
    end
    
    if any(strcmp('Title', fieldnames(app)))
        delete(app.Title)
    end
    
    if any(strcmp('barh', fieldnames(app)))
    	delete(app.barh)
    end
    
    app.Legend = gobjects(1,1);
    app.Title = gobjects(1,1);
    app.barh = gobjects(1,1);
    app.Axh = axes('Parent',app.figh,'FontSize',app.FS,'TickLabelInterpreter','latex',...
        'Position',[0.1974+0.02 0.1088 0.7878-0.02 0.8352],'Box','on','XMinorTick','on',... % [left bottom width height]
        'YMinorTick','on','XTickLabel','','YTickLabel','');
    app.Axh.XLabel.Interpreter = 'latex';
    app.Axh.YLabel.Interpreter = 'latex';
    app.Axh.ZLabel.Interpreter = 'latex';
    app.Axh.XLabel.FontSize = app.FS+10;
    app.Axh.YLabel.FontSize = app.FS+10;
    app.Axh.ZLabel.FontSize = app.FS+10;
    app.Axh.XLabel.String = '';
    app.Axh.YLabel.String = '';
    app.Axh.ZLabel.String = '';
    app.Axh.XTick = [];
    app.Axh.YTick = [];
    app.Axh.ZTick = [];
    
    set(app.figh,'renderer','opengl');
end

function AddNoteCallback(~,~)
    
    if strcmpi(app.Notepad.Visible,'on')
        if app.PlotPoint < 1
            app.Notepad.Visible = 'Off';
        else
            app.Notes{app.PlotPoint} = app.Notepad.String;
            app.Notepad.Visible = 'Off';
        end
    else
        if app.PlotPoint < 1
            app.Notepad.Visible = 'Off';
        else
            app.Notepad.Visible = 'On';
        end
    end
    
end

function EditNoteCallBack(~,~)
    if app.PlotPoint > 1
        app.Notes{app.PlotPoint} = app.Notepad.String;
    end
end

function CompareQlCallback(src,~)
    Current_compare = ~isempty(app.compareQl1.XDat) && ~isempty(app.compareQl2.XDat);
    Button_red = (src.BackgroundColor(1) == 1);
    
    unipl = unique_plot(app);
    plots_legtxt = app.legtxt(unipl);
    
    switch src.String
        case 'Compare 1'
            if isempty(app.XDat)
                app.compareQl1.XDat = [];
                app.compareQl1.YDat = [];
                app.compareQl1.ZDat = [];
                app.Comparetxt1 = '';
                app.compareQlButton1.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton1b.BackgroundColor = [0.9400    0.9400    0.9400];
            elseif all(isnan(app.XDat),'all')
                app.compareQl1.XDat = [];
                app.compareQl1.YDat = [];
                app.compareQl1.ZDat = [];
                app.Comparetxt1 = '';
                app.compareQlButton1.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton1b.BackgroundColor = [0.9400    0.9400    0.9400];
            elseif Button_red
                app.compareQl1.XDat = [];
                app.compareQl1.YDat = [];
                app.compareQl1.ZDat = [];
                app.Comparetxt1 = '';
                app.compareQlButton1.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton1b.BackgroundColor = [0.9400    0.9400    0.9400];
            else
                app.compareQl1.XDat = app.XDat;
                app.compareQl1.YDat = app.YDat;
                app.compareQl1.ZDat = app.ZDat;
                app.Comparetxt1 = [plots_legtxt{app.PlotPoint} ' - ' num2str(app.ZDat{2}) ' ns'];
                app.compareQlButton1.BackgroundColor = [1 0 0];
                app.compareQlButton1b.BackgroundColor = [1 0 0];
            end
            
        case 'Compare 2'
            if isempty(app.XDat)
                app.compareQl2.XDat = [];
                app.compareQl2.YDat = [];
                app.compareQl2.ZDat = [];
                app.Comparetxt2 = '';
                app.compareQlButton2.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton2b.BackgroundColor = [0.9400    0.9400    0.9400];
            elseif all(isnan(app.XDat),'all')
                app.compareQl2.XDat = [];
                app.compareQl2.YDat = [];
                app.compareQl2.ZDat = [];
                app.Comparetxt2 = '';
                app.compareQlButton2.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton2b.BackgroundColor = [0.9400    0.9400    0.9400];
            elseif Button_red
                app.compareQl2.XDat = [];
                app.compareQl2.YDat = [];
                app.compareQl2.ZDat = [];
                app.Comparetxt2 = '';
                app.compareQlButton2.BackgroundColor = [0.9400    0.9400    0.9400];
                app.compareQlButton2b.BackgroundColor = [0.9400    0.9400    0.9400];
            else
                app.compareQl2.XDat = app.XDat;
                app.compareQl2.YDat = app.YDat;
                app.compareQl2.ZDat = app.ZDat;
                app.Comparetxt2 = [plots_legtxt{app.PlotPoint} ' - ' num2str(app.ZDat{2}) ' ns'];
                app.compareQlButton2.BackgroundColor = [1 0 0];
                app.compareQlButton2b.BackgroundColor = [1 0 0];
            end
    end
    
    if (app.compareQlButton1.BackgroundColor(1) == 1) && (app.compareQlButton2.BackgroundColor(1) == 1)
        app.compareQlPlots = true;
    else
        app.compareQlPlots = false;
    end
    
    RefreshGUI;
    if app.compareQlPlots || Current_compare
        RefactorPlots;
    end
    
end

function CloseTrajSelCallback(~,~)
    set(app.auxfig,'Visible','off');
end

function SliceTrajectoryCallback(~,~)
    app.TJS.InitTime_Edit.String = app.InitialTimeh.String;
    app.TJS.FinTime_Edit.String = app.FinalTimeh.String;
    app.TJS.FrameTime_Edit.String = app.TimePerFrameh.String;
    
    app.TJS.Order_FirstShell.Value = app.FirstShell.Value;
    app.TJS.Order_SecondShell.Value = app.SecondShell.Value;
    app.TJS.Order_QlRadius.Value = app.QlRadius.Value;
    app.TJS.Order_QlNeighbours.Value = app.QlNeighbours.Value;
    app.TJS.Order_QlEditMin.String = app.QlEditMin.String;
    app.TJS.Order_QlEdit.String = app.QlEdit.String;
    app.TJS.Order_NeighbEdit.String = app.NeighbEdit.String;
    
    app.TJS.Neighbour_FirstShell.Value = app.FirstShell.Value;
    app.TJS.Neighbour_SecondShell.Value = app.SecondShell.Value;
    if app.TJS.Neighbour_FirstShell.Value || app.TJS.Neighbour_SecondShell.Value
        app.TJS.Neighbour_Search_Radius.Value = 0;
    else
        app.TJS.Neighbour_Search_Radius.Value = 1;
    end
    app.TJS.Neighbour_Search_Radius_Min.String = app.QlEditMin.String;
    app.TJS.Neighbour_Search_Radius_Max.String = app.QlEdit.String;
    app.TJS.Neighbour_RangeMin.String = app.NeighbEdit.String;
    app.TJS.Neighbour_RangeMax.String = app.NeighbEdit.String;
    
    % Set Metal and Halide
    Salt = app.Salth.String{app.Salth.Value};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    set(app.TJS.Type_Menu,'String',{'All' Metal Halide});
    
    set(app.auxfig,'Visible','on','position',[0.275 0.3 0.45 0.4]); %[left bottom width height]
    figure(app.auxfig);
end

function OrderParameterChanged(~,~)
    if app.TJS.Order_Menu.Value == 1
        set(app.TJS.Order_TextRange,'String',[char(8804) ' Ql ' char(8804)])
    else
        set(app.TJS.Order_TextRange,'String',[char(8804) ' Wl ' char(8804)])
    end
end

function OrderParameterRangeChanged(~,~)
    MinStr = get(app.TJS.Order_RangeMin,'String');
    MaxStr = get(app.TJS.Order_RangeMax,'String');
    
    if isempty(str2double(MinStr)) || isnan(str2double(MinStr))
        set(app.TJS.Order_RangeMin,'string','0');
        warndlg('The minimum value must be a number.');
    end
    if isempty(str2double(MaxStr)) || isnan(str2double(MaxStr))
        set(app.TJS.Order_RangeMax,'string','1');
        warndlg('The maximum value must be a number.');
    end
    
    if str2double(MinStr) > str2double(MaxStr)
        set(app.TJS.Order_RangeMin,'string','0');
        set(app.TJS.Order_RangeMax,'string','1');
        warndlg('The minimum value should be less than or equal to the maximum value.');
    end
end

function OrderTimeEditCallback(src,~)
    Str = get(src,'String');
    if isempty(str2double(Str)) || isnan(str2double(Str)) || str2double(Str) < 0
        set(src,'string','0');
        warndlg('The time value must be a non-negative number.');        
    end
    
end

function NeighbourRangeChanged(~,~)
    
    MinStr = str2double(app.TJS.Neighbour_RangeMin.String);
    MaxStr = str2double(app.TJS.Neighbour_RangeMax.String);
    
    if isnan(MinStr)
        warndlg('The minimum value must be a number.');
        app.TJS.Neighbour_RangeMin.String = '0';
    elseif MinStr < 0
        app.TJS.Neighbour_RangeMin.String = '0';
    else
        app.TJS.Neighbour_RangeMin.String = round(MinStr);
    end
    
    if isnan(MaxStr)
        app.TJS.Neighbour_RangeMax.String = num2str(round(str2double(app.TJS.Neighbour_RangeMin.String)));
    elseif MaxStr < round(str2double(app.TJS.Neighbour_RangeMin.String))
        app.TJS.Neighbour_RangeMax.String = num2str(round(str2double(app.TJS.Neighbour_RangeMin.String)));
    else
        app.TJS.Neighbour_RangeMax.String = round(MaxStr);
    end
end

function ClearFiltersCallback(~,~)
    app.TJS.current_filters = {};
    app.TJS.filter_notes.String = '';
end

function AddFilterCallback(src,~)
    if strcmpi(src.Parent.Tag,'order')
        filter.Type = 'order';
        filter.OrderParam = app.TJS.Order_Menu.String{app.TJS.Order_Menu.Value};
        filter.L = str2double(app.TJS.Order_MenuL.String{app.TJS.Order_MenuL.Value});
        filter.Ql_Min = str2double(app.TJS.Order_RangeMin.String);
        filter.Ql_Max = str2double(app.TJS.Order_RangeMax.String);
        filter.First_Shell_Sel = logical(app.TJS.Order_FirstShell.Value);
        filter.Second_Shell_Sel = logical(app.TJS.Order_SecondShell.Value);
        filter.N_Neighbours_Sel = logical(app.TJS.Order_QlNeighbours.Value);
        filter.Search_Radius_Sel = logical(app.TJS.Order_QlRadius.Value);
        filter.Search_Radius_Max = str2double(app.TJS.Order_QlEdit.String);
        filter.Search_Radius_Min = str2double(app.TJS.Order_QlEditMin.String);
        filter.N_Neighbours = str2double(app.TJS.Order_NeighbEdit.String);
        filter.Time_Point_Sel = logical(app.TJS.Order_OnePoint.Value);
        filter.Time_Point = str2double(app.TJS.Order_OnePointEdit.String);
        app.TJS.current_filters{end+1} = filter;
        
        notes = strrep(filter.OrderParam,'l',app.TJS.Order_MenuL.String{app.TJS.Order_MenuL.Value});
        notes = [app.TJS.Order_RangeMin.String char(8804) ' ' notes ' ' char(8804) app.TJS.Order_RangeMax.String];
        
        if filter.First_Shell_Sel && filter.Second_Shell_Sel
            notes = [notes ' for 1st and 2nd shell bonds'];
        elseif filter.First_Shell_Sel
            notes = [notes ' for 1st shell bonds'];
        elseif filter.Second_Shell_Sel
            notes = [notes ' for 2nd shell bonds'];
        elseif filter.Search_Radius_Sel
            notes = [notes ' for ' app.TJS.Order_QlEditMin.String char(8804) ' Rij ' char(8804) app.TJS.Order_QlEdit.String];
        elseif filter.N_Neighbours_Sel
            notes = [notes ' for ' app.TJS.Order_NeighbEdit.String ' nearest neighbours'];
        end
        
        if filter.Time_Point_Sel
            notes = [notes ' at t = ' app.TJS.Order_OnePointEdit.String ' ps.'];
        else
            notes = [notes ' for each frame.'];
        end
    elseif strcmpi(src.Parent.Tag,'type')
        
        filter.Type = 'type';
        filter.Type_sel = app.TJS.Type_Menu.String{app.TJS.Type_Menu.Value};
        filter.Time_Point_Sel = true; % will apply at all time points
        app.TJS.current_filters{end+1} = filter;
        
        notes = ['Type: ' filter.Type_sel ' for each frame.'];
    elseif strcmpi(src.Parent.Tag,'neighbour')
        
        filter.Type = 'neighbour';
        filter.N_Neighbours_Min = str2double(app.TJS.Neighbour_RangeMin.String);
        filter.N_Neighbours_Max = str2double(app.TJS.Neighbour_RangeMax.String);
        filter.First_Shell_Sel = logical(app.TJS.Neighbour_FirstShell.Value);
        filter.Second_Shell_Sel = logical(app.TJS.Neighbour_SecondShell.Value);
        filter.Search_Radius_Sel = logical(app.TJS.Neighbour_Search_Radius.Value);
        filter.Search_Radius_Max = str2double(app.TJS.Neighbour_Search_Radius_Max.String);
        filter.Search_Radius_Min = str2double(app.TJS.Neighbour_Search_Radius_Min.String);
        filter.Time_Point_Sel = logical(app.TJS.Neighbour_OnePoint.Value);
        filter.Time_Point = str2double(app.TJS.Neighbour_OnePointEdit.String);
        app.TJS.current_filters{end+1} = filter;
        
        notes = [app.TJS.Neighbour_RangeMin.String char(8804) ' N ' char(8804) app.TJS.Neighbour_RangeMax.String];
        
        if filter.First_Shell_Sel && filter.Second_Shell_Sel
            notes = [notes ' for 1st and 2nd shell bonds'];
        elseif filter.First_Shell_Sel
            notes = [notes ' for 1st shell bonds'];
        elseif filter.Second_Shell_Sel
            notes = [notes ' for 2nd shell bonds'];
        elseif filter.Search_Radius_Sel
            notes = [notes ' for ' app.TJS.Neighbour_Search_Radius_Min.String char(8804) ' Rij ' char(8804) app.TJS.Neighbour_Search_Radius_Max.String];
        end
        
        if filter.Time_Point_Sel
            notes = [notes ' at t = ' app.TJS.Neighbour_OnePointEdit.String ' ps.'];
        else
            notes = [notes ' for each frame.'];
        end
    end
    
    if isempty(app.TJS.filter_notes.String)
        app.TJS.filter_notes.String = {notes};
    else
        app.TJS.filter_notes.String{end+1} = notes; %[app.TJS.filter_notes.String newline notes];
    end
    
end

function LoadTrajectoryCallback(~,~)
    
    loc = fullfile(fileparts(mfilename('fullpath')),'Results','');    
    
    % Function to save the current app object to file
    [filename,PathName] = uigetfile(...
        {'*.pdb','PDB Trajectory Files (*.pdb)'},'Load Trajectory File',loc);
    if filename == 0
        return
    end
    
    Out_file = fullfile(PathName,filename);
    
    % Select the external viewer
    Viewer = app.TJS.Select_Viewer_text.String{app.TJS.Select_Viewer_text.Value};
    
    % Generate view trajectory command
    if strcmpi(Viewer,'OVITO')
        traj_view_command = ['start "C:\Program Files\OVITO Basic\ovito.exe" "' Out_file '"'];
    elseif strcmpi(Viewer,'VMD')
        if ispc
            traj_view_command = ['start ' app.vmd_loc ' "' Out_file '" -args "pbc box"'];
        else
            traj_view_command = ['gnome-terminal -- ' app.vmd_loc ' "' Out_file '"'];
        end
    elseif strcmpi(Viewer,'PYMOL')
        traj_view_command = ['start "C:\ProgramData\Anaconda3\PyMOLWin.exe" "' Out_file '"'];
    end
    system(traj_view_command);
    
end

function SaveTrajSelCallback(~,~)
    
    % Grab currently selected trajectory
    Salt = app.Salth.String{app.Salth.Value};
    Model = app.Modelh.String{app.Modelh.Value}; % {'JC (SPC/E)' 'JC (TIP4P-EW)' 'JC (TIP3P)' 'TF'};
    
    % Time slice options
    Min_Time = str2double(app.TJS.InitTime_Edit.String)*1000; % min time in ps
    Max_Time = str2double(app.TJS.FinTime_Edit.String)*1000; % max time in ps
    TimePerFrame = str2double(app.TJS.FrameTime_Edit.String); % time per frame (ps)
    
    % Create a default filename
    fn = [Salt '-' Model '-[' num2str(Min_Time) '-' num2str(Max_Time) 'ps].pdb'];
    
    loc = fullfile(fileparts(mfilename('fullpath')),'Results',fn);
    
    % Sanity check
    if Min_Time > Max_Time
        warndlg('The Initial Time must be less than or equal to the final time.');
        return
    elseif isnan(Min_Time) || isnan(Max_Time) || isnan(TimePerFrame)
        warndlg('Unable to interpret the selected time slice.');
        return
    end
    
    % Function to save the current app object to file
    [filename,PathName] = uiputfile(...
        {'*.pdb','PDB Trajectory Files (*.pdb)'},'Save Trajectory File',loc);
    if filename == 0
        return
    end
    
    % Directory containing results of interest
    DatDir = [app.datdir filesep Salt filesep Model];

    % Select the external viewer
    Viewer = app.TJS.Select_Viewer_text.String{app.TJS.Select_Viewer_text.Value};
    
    % Grab the current filters
    Filters = app.TJS.current_filters;
    
    % Find trajectory and gro files
    Traj_info = dir([DatDir filesep '*.trr']);
    Gro_info = dir([DatDir filesep '*.gro']);
    Mdp_info = dir([DatDir filesep '*.mdp']);
    if isempty(Traj_info) || isempty(Gro_info) || isempty(isempty(Mdp_info)) % no trajectory/gro file found
        warndlg('Unable to find required data file(s).')
        return
    end
    
    if length(Mdp_info) > 1
        Mdp_info = Mdp_info(~contains({Mdp_info.name},'out'));
    end
    Mdp_file = [DatDir filesep Mdp_info.name]; % mdp file
    Mdp_text = fileread(Mdp_file);
    pbc_on = isempty(regexp(Mdp_text,'pbc += *no','once'));

    % Select the initial time step gro file
    [~,idx] = min(cellfun(@length,{Gro_info.name}));
    Gro_names = {Gro_info.name};
    Gro_name = Gro_names{idx};

    Traj_file =[DatDir filesep Traj_info.name]; % traj file
    Gro_file = [DatDir filesep Gro_name]; % gro file
    Out_file = fullfile(PathName,filename);
    
    set(app.figh, 'pointer', 'watch')
    set(app.auxfig, 'pointer', 'watch')
    drawnow;
    try
        py.LiXStructureDetector.Slice_Traj(Traj_file,Gro_file,Out_file,Min_Time,Max_Time,...
            TimePerFrame,Filters,pbc_on,Viewer)
    catch e
        warndlg(['Error Executing Command:' newline e.message]);
        return
    end
    set(app.figh, 'pointer', 'arrow')
    set(app.auxfig, 'pointer', 'arrow')
    
    % Generate view trajectory command
    if strcmpi(Viewer,'OVITO')
        traj_view_command = ['start "C:\Program Files\OVITO Basic\ovito.exe" "' Out_file '"'];
    elseif strcmpi(Viewer,'VMD')
        if ispc
            traj_view_command = ['start ' app.vmd_loc ' "' Out_file '" -args "pbc box"'];
        else
            traj_view_command = ['gnome-terminal -- ' app.vmd_loc ' "' Out_file '"'];
        end
    elseif strcmpi(Viewer,'PYMOL')
        traj_view_command = ['start "C:\ProgramData\Anaconda3\PyMOLWin.exe" "' Out_file '"'];
    end
    system(traj_view_command);
end

end
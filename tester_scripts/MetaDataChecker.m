%% Load metadata
global app
app.PlotIdx = 1;
app.SaltIdx = 1;
app.ModelIdx = 1;
app.MainDir = 'D:\Molten_Salts_MD';
app.Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl' 'CsCl'};
app.MetaDataFile = fullfile(app.MainDir,'Calculation_Metadata.csv');
app.MetaData = readtable(app.MetaDataFile);
app.tempfolder = tempname;
mkdir(app.tempfolder);
app.gmx = 'gmx_d';
app.wsl = 'wsl source ~/.bashrc; ';
app.pipe = ' ^| ';
app.ModelA = 'V1NN 5-1';
app.ModelB = 'V2NN 5-1';

app.n_traj = size(app.MetaData,1);
app.fs = 24;

% Aspect ratio of images
app.Aspect_Ratio = 1.4; % Width/Height
app.Width = 0.6;
app.Height = app.Width/app.Aspect_Ratio;

% Get info about available trajectories
n_traj = 0;
app.Calculations = cell(1,length(app.Salts));
for idx = 1:length(app.Salts)
    Salt = app.Salts{idx};
    sjobs = dir(fullfile(app.MainDir,Salt));
    app.Calculations{idx} = {sjobs(cellfun(@(x) ~strncmp( x,'.',1),{sjobs.name})).name};
    n_traj = n_traj + length(app.Calculations{idx});
end
app.n_traj = n_traj;

app.figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','on','Units','normalized');
app.axh = axes(app.figh,'Position',[0 0.5 1/3 app.Height],'FontSize',app.fs,...
    'TickLabelInterpreter','latex','Box','off','XTick',[],'YTick',[],...
    'Visible','On');  % [left bottom width height]
app.axhb = axes(app.figh,'Position',[1/3 1/2 1/3 app.Height],'FontSize',app.fs,...
    'TickLabelInterpreter','latex','Box','off','XTick',[],'YTick',[],...
    'Visible','On');  % [left bottom width height]
app.axh2 = axes(app.figh,'Position',[0.1 0.1 app.Width-0.1 0.38],'FontSize',app.fs,...
    'TickLabelInterpreter','latex','Box','off','XTick',[],'YTick',[],...
    'Visible','On');
hold(app.axh2,'off')

plottxt = ['Plot: 1 of ' num2str(app.n_traj) newline ...
     newline ...
     'Salt:' newline ...
     newline ...
     'Model:' newline ...
     newline ...
     'Start point:' newline ...
     'End point:' newline ...
     'Atoms (MetaData):' newline ...
     'Atoms (Gro file):' ];
app.CurrentPlotTxt = uicontrol(app.figh,'Style','text','String',...
    plottxt,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.67578125,0.634536082474227,0.4,0.3],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.SaltDropDown = uicontrol(app.figh,'Style','popupmenu','String',...
    app.Salts,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7148 0.8315 0.05 0.03],... % [left bottom width height]
    'Value',1,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelSalt');
app.ModelDropDown = uicontrol(app.figh,'Style','popupmenu','String',...
    app.Calculations{1},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7310 0.7582 0.27 0.03],... % [left bottom width height]
    'Value',1,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelModel');
app.PlotChangeForward = uicontrol(app.figh,'Style','pushbutton','String',...
    'Next Trajectory','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7 0.19 0.2 0.1],... % [left bottom width height]
    'Value',0,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','pforward');
app.PlotChangeReverse = uicontrol(app.figh,'Style','pushbutton','String',...
    'Previous Trajectory','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7 0.09 0.2 0.1],... % [left bottom width height]
    'Value',0,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','pbackward','Enable','off');
app.PlotChangeBar = uicontrol(app.figh,'Style','slider','String',...
    ['Trajectory: ' num2str(app.PlotIdx)],'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7 0.025 0.2 0.025],... % [left bottom width height]
    'Value',app.PlotIdx,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','pslider','Enable','on','Max',app.n_traj,'Min',1,...
    'SliderStep',[1/app.n_traj 10/app.n_traj]);

app.OpenTrajectoryTxt = uicontrol(app.figh,'Style','text','String',...
    'Open Trajectory','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7378 0.5402 0.1216 0.0412],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.TrajSelectA = uicontrol(app.figh,'Style','pushbutton','String',...
    app.ModelA,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.7 0.44 0.1 0.1],... % [left bottom width height]
    'Value',0,'Callback',@OpenTrajCallback,'Visible','on',...
    'tag','trajA');
app.TrajSelectB = uicontrol(app.figh,'Style','pushbutton','String',...
    app.ModelB,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.8 0.44 0.1 0.1],... % [left bottom width height]
    'Value',0,'Callback',@OpenTrajCallback,'Visible','on',...
    'tag','trajB');

app.TrajTypeA = uicontrol(app.figh,'Style','popupmenu','String',...
    {'V1NN 0-0' 'V1NN 5-1' 'V2NN 0-0' 'V2NN 0-1' 'V2NN 3-1' 'V2NN 5-1' 'V2NN 11-1' ...
    'V2NN 11-5' 'V2NN 21-5' 'V2NN 51-5' 'V2NN 51-10'},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.1185 0.9376 0.1 0.05],... % [left bottom width height]
    'Value',2,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','trajTypeA');
app.TrajTypeB = uicontrol(app.figh,'Style','popupmenu','String',...
    {'V1NN 0-0' 'V1NN 5-1' 'V2NN 0-0' 'V2NN 0-1' 'V2NN 3-1' 'V2NN 5-1' 'V2NN 11-1' ...
    'V2NN 11-5' 'V2NN 21-5' 'V2NN 51-5' 'V2NN 51-10'},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.4526 0.9376 0.1000 0.0500],... % [left bottom width height]
    'Value',6,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','trajTypeB');


app.PlotTypeXtxt = uicontrol(app.figh,'Style','text','String',...
    'X:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.70 0.388 0.0300 0.0500],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.PlotTypeX = uicontrol(app.figh,'Style','popupmenu','String',...
    {'Time' 'LJ-(SR)' 'Coulomb-(SR)' 'Coul.-recip.' 'Potential' 'Kinetic-En.' ...
    'Total-Energy' 'Conserved-En.' 'Temperature' 'Pressure' 'Box-X' ...
    'Box-Y' 'Box-Z' 'Volume' 'Density' 'pV' 'Enthalpy'},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.73 0.388 0.1500 0.0500],... % [left bottom width height]
    'Value',1,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','PlotTypeX');

app.PlotTypeY1txt = uicontrol(app.figh,'Style','text','String',...
    'Y1:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.70 0.338 0.03 0.050],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.PlotTypeY1 = uicontrol(app.figh,'Style','popupmenu','String',...
    {'LJ-(SR)' 'Coulomb-(SR)' 'Coul.-recip.' 'Potential' 'Kinetic-En.' ...
    'Total-Energy' 'Conserved-En.' 'Temperature' 'Pressure' 'Box-X' ...
    'Box-Y' 'Box-Z' 'Volume' 'Density' 'pV' 'Enthalpy'},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.73 0.338 0.15 0.05],... % [left bottom width height]
    'Value',4,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','PlotTypeY1');
app.PlotTypeY2txt = uicontrol(app.figh,'Style','text','String',...
    'Y2:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.70 0.288 0.03 0.050],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.PlotTypeY2 = uicontrol(app.figh,'Style','popupmenu','String',...
    {'LJ-(SR)' 'Coulomb-(SR)' 'Coul.-recip.' 'Potential' 'Kinetic-En.' ...
    'Total-Energy' 'Conserved-En.' 'Temperature' 'Pressure' 'Box-X' ...
    'Box-Y' 'Box-Z' 'Volume' 'Density' 'pV' 'Enthalpy'},'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.73 0.288 0.1500 0.0500],... % [left bottom width height]
    'Value',8,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','PlotTypeY2');

% Initialize plot
src.Tag = 'pslider';
PlotChangeCallback(src,[]);

function PlotChangeCallback(src,~)
    global app
    
    % Which button was pressed?
    switch src.Tag
        case 'pforward'
            app.PlotIdx = app.PlotIdx + 1;
            app = GetSaltModelIdx(app);
        case 'pbackward'
            app.PlotIdx = app.PlotIdx - 1;
            app = GetSaltModelIdx(app);
        case 'pslider'
            app.PlotIdx = round(app.PlotChangeBar.Value);
            app = GetSaltModelIdx(app);
        case 'SelSalt'
            if app.SaltDropDown.Value == app.SaltIdx
                return
            else
                app.SaltIdx = app.SaltDropDown.Value;
                app.ModelIdx = 1;
                app = GetPlotIdx(app);
            end
        case 'SelModel'
            if app.ModelDropDown.Value == app.ModelIdx
                return
            else
                app.ModelIdx = app.ModelDropDown.Value;
                app = GetPlotIdx(app);
            end
        case 'trajTypeA'
            app.TrajSelectA.String = app.TrajTypeA.String{app.TrajTypeA.Value};
            app.ModelA = app.TrajSelectA.String;
        case 'trajTypeB'
            app.TrajSelectB.String = app.TrajTypeB.String{app.TrajTypeB.Value};
            app.ModelB = app.TrajSelectB.String;
    end
    app.PlotChangeBar.Value = app.PlotIdx;
    app.SaltDropDown.Value = app.SaltIdx;
    app.ModelDropDown.String = app.Calculations{app.SaltIdx};
    app.ModelDropDown.Value = app.ModelIdx;
    
    if app.PlotIdx >= app.n_traj
        app.PlotChangeForward.Enable = 'off';
    elseif app.PlotIdx <= 1
        app.PlotChangeReverse.Enable = 'off';
    else
        app.PlotChangeForward.Enable = 'on';
        app.PlotChangeReverse.Enable = 'on';
    end
    
    Salt = app.Salts{app.SaltIdx};
    Model = app.Calculations{app.SaltIdx}{app.ModelIdx};
    
    % Find calculation in metadata file, if it exists    
    meta_idx = cellfun(@(x) strcmp(x,Salt),app.MetaData.Salt) & cellfun(@(x) strcmp(x,Model),app.MetaData.Model);
    if sum(meta_idx) > 0
        N_Atoms = app.MetaData{meta_idx,5};
        ns_0 = app.MetaData{meta_idx,6};
        ns_f = app.MetaData{meta_idx,7};
    else
        N_Atoms = nan;
        ns_0 = nan;
        ns_f = nan;
    end
    WorkDir = fullfile(app.MainDir,Salt,Model);
    
    % Update summary plots
    overview_fileA = fullfile(WorkDir,[strrep(app.ModelA,' ','-') '_' Salt '_' Model '.png']);
    
    try
        [X1,map1] = imread(overview_fileA);
        imshow(X1,map1,'Parent',app.axh,'border','tight','InitialMagnification','fit',...
            'XData', [0 1], 'YData', [0 1/app.Aspect_Ratio]);
    catch
        warndlg(['Missing file: ' overview_fileA])
    end
    
    overview_fileB = fullfile(WorkDir,[strrep(app.ModelB,' ','-') '_' Salt '_' Model '.png']);
    
    try
        [X1,map1b] = imread(overview_fileB);
        imshow(X1,map1b,'Parent',app.axhb,'border','tight','InitialMagnification','fit',...
            'XData', [0 1], 'YData', [0 1/app.Aspect_Ratio]);
    catch
        warndlg(['Missing file: ' overview_fileB])
    end
    
    %% Upadate energy overview
    % Find trajectory file
    Ener_file = fullfile(WorkDir,[Model '.edr']);
    Tpr_file = fullfile(WorkDir,[Model '.tpr']);
    Out_file = fullfile(app.tempfolder,'tempEnergy.xvg');
    Gro_file = fullfile(WorkDir,[Model '.gro']);

    if ~isfile(Ener_file) || ~isfile(Tpr_file) || ~isfile(Gro_file)
        warndlg('Missing required file, one of: *.edr; *.tpr; or *.gro')
        return
    end

    gmx_command = [app.wsl 'echo 0' app.pipe app.gmx ' energy -f ' windows2unix(Ener_file)...
        ' -s ' windows2unix(Tpr_file)];
    [err,outpt] = system(gmx_command);
    if err ~= 1
        warndlg('Problem with energy check.')
        return
    end
    
    Y1_sel = app.PlotTypeY1.String{app.PlotTypeY1.Value};
    Y2_sel = app.PlotTypeY2.String{app.PlotTypeY2.Value};
    X_sel = app.PlotTypeX.String{app.PlotTypeX.Value};
    
    en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
    en_opts = en_opts{1};
    Y1_num = char(regexp(en_opts,['([0-9]{1,2})  ' replace(Y1_sel,{'(' ')'},{'\(' '\)'})],'tokens','once'));
    Y2_num = char(regexp(en_opts,['([0-9]{1,2})  ' replace(Y2_sel,{'(' ')'},{'\(' '\)'})],'tokens','once'));
    X_num =  char(regexp(en_opts,['([0-9]{1,2})  ' replace(X_sel, {'(' ')'},{'\(' '\)'})],'tokens','once'));
    
    En_set = Y1_num;
    En_set = [En_set ' ' Y2_num];
    En_set = [En_set ' ' X_num];
    En_set = [En_set ' 0'];
    En_set = regexprep(En_set,' +',' ');
    
    % Grab data from results
    gmx_command = [app.wsl 'echo ' En_set app.pipe app.gmx ' energy -f ' windows2unix(Ener_file)...
        ' -o ' windows2unix(Out_file) ' -s ' windows2unix(Tpr_file) ...
        ' -b 0 -e 50000'];

    set(app.figh, 'pointer', 'watch')
    drawnow;
    err = system(gmx_command);
    set(app.figh, 'pointer', 'arrow')
    if err ~= 0
        warndlg('Failed to collect data.')
        return
    end
    
    NSS = load_gro_file(Gro_file);
    NF = NSS.N_atoms/2; % number of ion pairs in supercell
    Data = import_xvg(Out_file); % Gather energies
    delete(Out_file) % remove temp output file
	if isextensive(Y1_sel)
        NY1 = NF;
	else
        NY1 = 1;
	end
    if isextensive(Y2_sel)
        NY2 = NF;
    else
        NY2 = 1;
    end
    if isextensive(X_sel)
        NX = NF;
    else
        NX = 1;
    end
    
    if strcmp(Y1_sel,Y2_sel)
        Y1 = Data(:,2);
        Y2 =  Y1;
    else
        Y1 = Data(:,2);
        Y2 =  Data(:,3);
    end
    
    if strcmp(X_sel,'Time')
        X = Data(:,1)./1000; % ns
    else
        if strcmp(Y1_sel,X_sel)
            X = Y1;
        elseif strcmp(Y2_sel,X_sel)
            X = Y2;
        elseif strcmp(Y1_sel,Y2_sel)
            X = Data(:,3);
        else
            X = Data(:,4);
        end
    end
    
    nums = [str2double(Y1_num) str2double(Y2_num) str2double(X_num)];
    [~,sidx] = sort(nums);
    Data = [Y1 Y2 X];
    Data = Data(:,sidx);
    Y1 = Data(:,1)./NY1;
    Y2 = Data(:,2)./NY2;
    X = Data(:,3)./NX;
    
    yyaxis(app.axh2,'left')
    cla(app.axh2)
    %xline(app.axh2,ns_0,'--r','LineWidth',3)
    if ~isnan(ns_f) && strcmp(X_sel,'Time')
        xline(app.axh2,ns_f,'--r','LineWidth',3)
    end
    hold(app.axh2,'on')
    plot(app.axh2,X,Y1,'LineWidth',3);
    ylim(app.axh2,'auto')
    xlim(app.axh2,'auto')
    yticks(app.axh2,'auto')
	
    ylabel(app.axh2,GetAxisLabel(Y1_sel),'Interpreter','latex','FontSize',app.fs)
    hold(app.axh2,'off')
    yyaxis(app.axh2,'right')
    plot(app.axh2,X,Y2,'LineWidth',3);
    ylabel(app.axh2,GetAxisLabel(Y2_sel),'Interpreter','latex','FontSize',app.fs)
    ylim(app.axh2,'auto')
    yticks(app.axh2,'auto')
    
    
    xlabel(app.axh2,GetAxisLabel(X_sel),'Interpreter','latex','FontSize',app.fs)
    xticks(app.axh2,'auto')
    
    % Update figure text and title
    app.figh.Name = [Salt ' ' Model];
    
    plottxt = ['Plot: ' num2str(app.PlotIdx) ' of ' num2str(app.n_traj) newline ...
             newline ...
             'Salt:' newline ...
             newline ...
             'Model:' newline ...
             newline ...
             'Start point: '  num2str(ns_0) ' ns' newline ...
             'End point: '  num2str(ns_f) ' ns' newline ...
             'Atoms (MetaData): ' num2str(N_Atoms) newline ...
             'Atoms (Gro file): ' num2str(NSS.N_atoms)];
    app.CurrentPlotTxt.String = plottxt;
    
end


function OpenTrajCallback(src,~)
    global app
    
    Salt = app.Salts{app.SaltIdx};
    Model = app.Calculations{app.SaltIdx}{app.ModelIdx};
    WorkDir = fullfile(app.MainDir,Salt,Model);
    
    switch src.Tag
        case 'trajA'
        	traj_file = fullfile(WorkDir,[strrep(app.ModelA,' ','-') '_' Salt '_' Model '.xyz']);
            
        case 'trajB'
        	traj_file = fullfile(WorkDir,[strrep(app.ModelB,' ','-') '_' Salt '_' Model '.xyz']);
            
    end
    
    syscmd = ['start "C:\Program Files\OVITO Basic\ovito.exe" "' traj_file '"'];
    cpath = getenv('QT_PLUGIN_PATH');
    setenv('QT_PLUGIN_PATH','')
    system(syscmd);
    setenv(cpath)
end

function app = GetPlotIdx(app)
    n_plot = 0;
    for idx = 1:length(app.Salts)
        if idx == app.SaltIdx
            n_plot = n_plot + app.ModelIdx;
            app.PlotIdx = n_plot;
            return
        else
            n_plot = n_plot + length(app.Calculations{idx});
        end
    end
    
end

function app = GetSaltModelIdx(app)
    n_plot = 0;
    for idx = 1:length(app.Calculations)
        n_plot_prev = n_plot;
        n_plot = n_plot + length(app.Calculations{idx});
        
        if app.PlotIdx <= n_plot
            app.SaltIdx = idx;
            app.ModelIdx = app.PlotIdx - n_plot_prev;
            return
        end
    end
end


function label = GetAxisLabel(Property_type)
    switch Property_type
        case 'Time'
            label = 'Time [ns]';
        case 'LJ-(SR)'
            label = 'LJ Energy (SR) [kJ mol$^{-1}$]';
        case 'Coulomb-(SR)'
            label = 'Coulomb Energy (SR) [kJ mol$^{-1}$]';
        case 'Coul.-recip.'
            label = 'Coulomb Energy (Recip) [kJ mol$^{-1}$]';
        case 'Potential'
            label = 'Potential [kJ mol$^{-1}$]';
        case 'Kinetic-En.'
            label = 'Kinetic Energy [kJ mol$^{-1}$]';
        case 'Total-Energy'
            label = 'Total Energy [kJ mol$^{-1}$]';
        case 'Conserved-En.'
            label = 'Conserved Energy [kJ mol$^{-1}$]';
        case 'Temperature'
            label = 'Temperature [K]';
        case 'Pressure'
            label = 'Pressure [GPa]';
        case 'Box-X'
            label = 'Box Length X [nm]';
        case 'Box-Y'
            label = 'Box Length Y [nm]';
        case 'Box-Z'
            label = 'Box Length Z [nm]';
        case 'Volume'
            label = 'Volume [nm$^{3}$ molecule$^{-1}$]';
        case 'Density'
            label = 'Density [kg m$^{-3}$]';
        case 'pV'
            label = 'pV [kJ mol$^{-1}$]';
        case 'Enthalpy'
            label = 'Enthalpy [kJ mol$^{-1}$]';
    end
end

function extensive = isextensive(Property_type)

    switch Property_type
        case {'Time' 'Temperature' 'Pressure' 'Box-X' 'Box-Y' 'Box-Z' 'Density'}
            extensive = false;
        case {'LJ-(SR)' 'Coulomb-(SR)' 'Coul.-recip.' 'Potential' 'Kinetic-En.' ...
               'Total-Energy' 'Conserved-En.' 'Volume' 'pV' 'Enthalpy'}
            extensive = true;
    end
end
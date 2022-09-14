% @    xaxis  label "Relative position from center (nm)"
% @    yaxis  label "Density (kg m\S-3\N)"

% Density Profile Plotting
home = find_home;
%DataFile = fullfile(home,'data','Melting_Point_Data.mat');
DataFile = fullfile(home,'data','Melting_Point_Data.mat');

global app
app.Data = load(DataFile).Data;
app.SaltIdx = 1;
app.ModelIdx = 1;
app.CalcIdx = [];
app.Structure = 'Rocksalt';
app.Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
app.fs = 32;
app.lw = 3;

% Setup empty plot
app.figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
app.axh = axes(app.figh,'position',[0.1 0.1 0.7 0.8]);
hold(app.axh,'on')

app.Models = fieldnames(app.Data.(app.Salts{app.SaltIdx}));
app.Calcs = strrep(fieldnames(app.Data.(app.Salts{app.SaltIdx}).(app.Models{app.ModelIdx}).(app.Structure)),'MP_','');

% Setup UI
app.SaltTxt = uicontrol(app.figh,'Style','text','String',...
    'Salt:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.81,0.8,0.07,0.05],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.SaltDropDown = uicontrol(app.figh,'Style','popupmenu','String',...
    app.Salts,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.87 0.8 0.05 0.05],... % [left bottom width height]
    'Value',app.SaltIdx,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelSalt');

app.ModelTxt = uicontrol(app.figh,'Style','text','String',...
    'Model:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.81,0.7,0.07,0.05],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.ModelDropDown = uicontrol(app.figh,'Style','popupmenu','String',...
    app.Models,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.87 0.7 0.13 0.05],... % [left bottom width height]
    'Value',app.ModelIdx,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelModel');

app.CalcTxt = uicontrol(app.figh,'Style','text','String',...
    'Calc:','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.81,0.6,0.07,0.05],... % [left bottom width height]
    'Visible','on','HorizontalAlignment','left');
app.CalcDropDown = uicontrol(app.figh,'Style','list','String',...
    app.Calcs,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.87 0.4 0.13 0.25],... % [left bottom width height]
    'Value',app.CalcIdx,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelCalc','max',Inf,'min',0);

app.ClearButton = uicontrol(app.figh,'Style','pushbutton','String',...
    'Clear Plot Area','FontSize',app.fs,'Units','Normalized',...
    'Position',[0.825 0.26 0.15 0.1],... % [left bottom width height]
    'Value',0,'Callback',@ClearPlotCallBack,'Visible','on');

% Initialize plot
ClearPlotCallBack;
grid(app.axh,'minor');
grid(app.axh,'on');

function PlotChangeCallback(src,~)
    global app
    
    % Which button was pressed?
    switch src.Tag
        case 'SelSalt'
            if app.SaltDropDown.Value == app.SaltIdx
                return
            else
                app.SaltIdx = app.SaltDropDown.Value;
                app.ModelIdx = 1;
                app.CalcIdx = 1;
                
                % Update models and calcs
                app.Models = fieldnames(app.Data.(app.Salts{app.SaltIdx}));
                app.Calcs = strrep(fieldnames(app.Data.(app.Salts{app.SaltIdx}).(app.Models{app.ModelIdx}).(app.Structure)),'MP_','');
                
            end
        case 'SelModel'
            if app.ModelDropDown.Value == app.ModelIdx
                return
            else
                app.ModelIdx = app.ModelDropDown.Value;
                app.CalcIdx = [];
                app.Calcs = strrep(fieldnames(app.Data.(app.Salts{app.SaltIdx}).(app.Models{app.ModelIdx}).(app.Structure)),'MP_','');
            end
        case 'SelCalc'
            if isequal(app.CalcDropDown.Value,app.CalcIdx)
                return
            else
                app.CalcIdx = unique([app.CalcIdx app.CalcDropDown.Value]);
            end
        otherwise
            % Do nothing
    end
    app.SaltDropDown.Value = app.SaltIdx;
    app.ModelDropDown.String = app.Models;
    app.ModelDropDown.Value = app.ModelIdx;
    app.CalcDropDown.String = app.Calcs;
    app.CalcDropDown.Value = app.CalcIdx;
    
    % Add any new plots to the plot list
    Salt = app.Salts{app.SaltIdx};
    Model = app.Models{app.ModelIdx};
    for cidx = app.CalcIdx
        Calc = app.Calcs{cidx};
        LegText = [Salt ' ' strrep(Model,'_',' ') ' ' strrep(Calc,'_',' ')];
        
        % Check if already plotted, if not grab the plot information and save it
        if ~any(strcmp(app.LegTexts,LegText))
            app.LegTexts{end+1} = LegText;
            
            T_Trace = app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).T_Trace;
            T_Freeze_Trace = logical(app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).Freeze_Trace);
            T_Melt_Trace = logical(app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).Melt_Trace);
            Tm = mean(app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).dT);
            
            CData = repmat([1 1 1],length(T_Trace),1); % Default color is white
            CData(T_Melt_Trace',:) = repmat([1 0 0],sum(T_Melt_Trace),1); % Melted points are coloured red
            CData(T_Freeze_Trace',:) = repmat([0 0 1],sum(T_Freeze_Trace),1); % frozen points are coloured blue
            
            % "Fake" melting points are coloured light red
            Fake_MP_Melt_idx = T_Melt_Trace & app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).f_Trace < 0;
            CData(Fake_MP_Melt_idx,:) = repmat([1 0.75 0.75],sum(Fake_MP_Melt_idx),1);
            
            % "Fake" melting points are coloured light blue
            Fake_MP_Freeze_idx = T_Freeze_Trace & app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).f_Trace < 0;
            CData(Fake_MP_Freeze_idx,:) = repmat([0.75 0.75 1],sum(Fake_MP_Freeze_idx),1);
            
            app.T_Traces{end+1} = T_Trace;
            app.Tms{end+1} = Tm;
            app.CDatas{end+1} = CData;
            
        else
            continue
        end
    end
    
    % Reset the plot with new plot colours
    cla(app.axh)
    N_plots = length(app.LegTexts);
    colours = min(max(cbrewer('qual','Dark2',max(3,N_plots),'cubic'),0),1);
    plt_objs = gobjects(N_plots,1);
    X_max = 0;
    for pidx = 1:N_plots
        col = colours(pidx,:);
        X = 1:length(app.T_Traces{pidx});
        yline(app.axh,app.Tms{pidx},'-.','LineWidth',app.lw,'Color',col)
        plt_objs(pidx) = plot(app.axh,X,app.T_Traces{pidx},'-','LineWidth',app.lw,'Color',col);
        scatter(app.axh,X,app.T_Traces{pidx},150,app.CDatas{pidx},'filled','LineWidth',app.lw-1,'MarkerEdgeColor','k');
        X_max = max([X_max X]);
    end
    legend(app.axh,plt_objs,app.LegTexts,'Fontsize',app.fs,'Interpreter','latex',...
        'Location','southeast')
    title(app.axh,'Melting Point Search Temperature Trajectories','Fontsize',app.fs,'Interpreter','latex')
    set(app.axh,'XLim',[0 X_max+1],'FontSize',app.fs)
    ylim(app.axh,'auto')
    
    xlabel(app.axh, 'Algorithm Step','Fontsize',app.fs,'Interpreter','latex')
    ylabel(app.axh, '$T$ [K]','Fontsize',app.fs,'Interpreter','latex')
    set(app.axh,'XMinorTick','on','YMinorTick','on','FontSize',app.fs);
    set(app.axh,'box','on','TickLabelInterpreter','latex');
end

function ClearPlotCallBack(~,~)
    global app
    
    cla(app.axh)
    app.CalcDropDown.Value = [];
    app.CalcIdx = [];
    app.LegTexts = {};
    app.T_Traces = {};
    app.Tms      = {};
    app.CDatas   = {};
    src.Tag = '';
    PlotChangeCallback(src,[]);
end

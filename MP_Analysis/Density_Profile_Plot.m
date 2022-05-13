% @    xaxis  label "Relative position from center (nm)"
% @    yaxis  label "Density (kg m\S-3\N)"

% Density Profile Plotting
home = find_home;
DataFile = fullfile(home,'data','Melting_Point_Data.mat');

global app
app.Data = load(DataFile).Data;
app.SaltIdx = 1;
app.ModelIdx = 1;
app.CalcIdx = 1;
app.Structure = 'Rocksalt';
app.Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
app.fs = 24;

% Setup empty plot
app.figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
app.axh = axes(app.figh,'position',[0.1 0.1 0.7 0.8]);

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
app.CalcDropDown = uicontrol(app.figh,'Style','popupmenu','String',...
    app.Calcs,'FontSize',app.fs,'Units','Normalized',...
    'Position',[0.87 0.6 0.13 0.05],... % [left bottom width height]
    'Value',app.CalcIdx,'Callback',@PlotChangeCallback,'Visible','on',...
    'tag','SelCalc');

% Initialize plot
src.Tag = '';
PlotChangeCallback(src,[]);

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
                app.CalcIdx = 1;
                app.Calcs = strrep(fieldnames(app.Data.(app.Salts{app.SaltIdx}).(app.Models{app.ModelIdx}).(app.Structure)),'MP_','');
            end
        case 'SelCalc'
            if app.CalcDropDown.Value == app.CalcIdx
                return
            else
                app.CalcIdx = app.CalcDropDown.Value;
            end
        otherwise
            % Do nothing
    end
    app.SaltDropDown.Value = app.SaltIdx;
    app.ModelDropDown.String = app.Models;
    app.ModelDropDown.Value = app.ModelIdx;
    app.CalcDropDown.String = app.Calcs;
    app.CalcDropDown.Value = app.CalcIdx;
    
    
    Salt = app.Salts{app.SaltIdx};
    Model = app.Models{app.ModelIdx};
    Calc = app.Calcs{app.CalcIdx};
    
    ZProfile = app.Data.(Salt).(Model).(app.Structure).(['MP_' Calc]).DensityProfile;
    width = (ZProfile(end,1) - ZProfile(1,1))/2;

    plot(app.axh,ZProfile(:,1)+width,ZProfile(:,2),'-k','LineWidth',3);
    title(app.axh,['$\vec{Z}$-Dimension Density Profile: ' Salt ' ' Model ' ' ...
        strrep(strrep(Calc,'MP_',''),'_','\_')],'Fontsize',app.fs,'Interpreter','latex')
    set(app.axh,'XLim',[0 width*2],'FontSize',app.fs)
    
    xlabel(app.axh, 'Position along $\vec{Z}$ [nm]','Fontsize',app.fs,'Interpreter','latex')
    ylabel(app.axh, '$\rho$ [kg m$^{-3}$]','Fontsize',app.fs,'Interpreter','latex')
    grid(app.axh,'minor');
    set(app.axh,'XMinorTick','on','YMinorTick','on','FontSize',app.fs);
    set(app.axh,'box','on','TickLabelInterpreter','latex');
end
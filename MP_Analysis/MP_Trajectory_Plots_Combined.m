% @    xaxis  label "Relative position from center (nm)"
% @    yaxis  label "Density (kg m\S-3\N)"

% Density Profile Plotting
home = find_home;
DataFile = fullfile(home,'data','Melting_Point_Data.mat');

Data = load(DataFile).Data;
Salt = 'NaCl';
Structure = 'Rocksalt';
Model = 'JC';
Tag = 'Set2';
fs = 24;
linestyles = {'-',':','-.','--'};
mkrshapes = {'o' 's' 'd' '^' 'v' '>' '<' 'p' 'h'};

% Setup empty plot
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'position',[0.1 0.1 0.85 0.85]);
hold(axh,'on')

Calcs = strrep(fieldnames(Data.(Salt).(Model).(Structure)),'MP_','');
Calcs_matched = Calcs(cellfun(@(x) contains(x,Tag),Calcs));
N = length(Calcs_matched);
N_linestyle = length(linestyles);
N_markershape = length(mkrshapes);
colours = max(min(cbrewer('qual','Dark2',N,'spline'),1),0);
ph = gobjects(1,N);

switch Tag
    case 'Set1'
        Vars = regexp(Calcs_matched,'XY_([0-9]+)','tokens','once');
        [~,order_idx] = sort(str2double(vertcat(Vars{:})));
        Calcs_matched = Calcs_matched(order_idx);
        legtxts = strcat(regexprep(Calcs_matched,{'_Z_10' 'Set1_' '_'},{'' '' ' '}),' nm');
    case 'Set2'
        Vars = regexp(Calcs_matched,'Z_([0-9]+_*[0-9]*)','tokens','once');
        [~,order_idx] = sort(str2double( strrep(vertcat(Vars{:}),'_','.') ));
        Calcs_matched = Calcs_matched(order_idx);
        legtxts = strcat(regexprep(Calcs_matched,{'_XY_8_' 'Set2' 'Z_' '_'},{'' '' 'Z=' '.'}),' nm');
    case 'Set3'
        Vars = regexp(Calcs_matched,'R_([0-9]+)','tokens','once');
        [~,order_idx] = sort(str2double(vertcat(Vars{:})));
        Calcs_matched = Calcs_matched(order_idx);
        legtxts = strcat(regexprep(Calcs_matched,{'Set3_' '_'},{'' ' '}),' nm');
end

mkridx = 1;
linidx = 1;
for idx = 1:N
    if mkridx > N_markershape
        mkridx = 1;
    end
    if linidx > N_linestyle
        linidx = 1;
    end
    mkrshape = mkrshapes{mkridx};
    linestyle = linestyles{linidx};
    col = colours(idx,:);
    
    Calc = Calcs_matched{idx};
    T_Trace = Data.(Salt).(Model).(Structure).(['MP_' Calc]).T_Trace;
    T_Freeze_Trace = logical(Data.(Salt).(Model).(Structure).(['MP_' Calc]).Freeze_Trace);
    T_Melt_Trace = logical(Data.(Salt).(Model).(Structure).(['MP_' Calc]).Melt_Trace);
    Tm = mean(Data.(Salt).(Model).(Structure).(['MP_' Calc]).dT);


    X = 1:length(T_Trace);
    CData = repmat([1 1 1],length(X),1); % Default color is white
    CData(T_Melt_Trace',:) = repmat([1 0 0],sum(T_Melt_Trace),1); % Melted points are coloured red
    CData(T_Freeze_Trace',:) = repmat([0 0 1],sum(T_Freeze_Trace),1); % frozen points are coloured blue
    
    yline(axh,Tm,':','LineWidth',2,'Color',col)
    ph(idx) = plot(axh,X,T_Trace,'-.','LineWidth',2,'Color',col);
    scatter(axh,X,T_Trace,100,CData,'filled',mkrshape,'LineWidth',2,'MarkerEdgeColor',col);
    mkridx = mkridx + 1;
    linidx = linidx + 1;
end

switch Tag
    case 'Set1'
        legtxts = regexprep(Calcs_matched,{'_Z_10' 'Set1_' '_'},{'' '' ' '});
    case 'Set2'
    case 'Set3'
end

legend(axh,ph,legtxts)
title(axh,['Melting Point Search: Temperature Trajectory for ' Salt ' ' Model ' ' ...
    Tag],'Fontsize',fs,'Interpreter','latex')
set(axh,'XLim',[0 length(X)+1],'FontSize',fs)
ylim(axh,'auto')

xlabel(axh, 'Algorithm Step','Fontsize',fs,'Interpreter','latex')
ylabel(axh, '$T$ [K]','Fontsize',fs,'Interpreter','latex')
grid(axh,'minor');
set(axh,'XMinorTick','on','YMinorTick','on','FontSize',fs);
set(axh,'box','on','TickLabelInterpreter','latex');
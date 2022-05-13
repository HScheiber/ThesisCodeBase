% This function produces a crosshair that tracks the mouse

function crosshairfcn(~,~,x_axis,Y_Axis,plot_extra,fs)

if x_axis == "bl"
    xtext = 'BL';
    xunits = char(197); % Angstrom
elseif x_axis == "a"
    xtext = 'a';
    xunits = char(197); % Angstrom
elseif x_axis == "V"
    xtext = 'V';
    xunits = 'cm^{3}/mol'; % cubic cm per mol
elseif x_axis == "P"
    xtext = 'P';
    xunits = 'GPa'; % cubic cm per mol
elseif x_axis == "d"
    xtext = 'd';
    xunits = 'g/cm^{3}'; % cubic cm per mol
end

if strcmp(Y_Axis,'L')
    ytext = 'E';
    yunits = 'kJ/mol';
elseif strcmp(Y_Axis,'C')
    ytext = 'E';
    yunits = 'kJ/mol';
elseif strcmp(Y_Axis,'P')
    ytext = 'P';
    yunits = 'GPa';
end

ylim manual

%Disable data cursor and zoom, it doesn't interact well with this function
set(0,'Showhidden','on')
ch = findall(gcf,'tag','FigureToolBar');
UT = get(ch,'children');
UV(1) = findall(UT,'tag','Exploration.DataCursor');
UV(2) = findall(UT,'tag','Exploration.Brushing');
UV(3) = findall(UT,'tag','Exploration.Pan');
UV(4) = findall(UT,'tag','Exploration.ZoomOut');
UV(5) = findall(UT,'tag','Exploration.ZoomIn');
UV(6) = findall(UT,'tag','OptionsFunction');

set(UV(1:5),'Enable','Off','State','Off');
set(UV(6),'Enable','Off');
zoom on
zoom off

% Create custom invisible cursor for use while within axes
set(gcf,'PointerShapeCData',ones(16, 16)*nan);

% Set up cursor text
if plot_extra
    hText = annotation(gcf,'textbox','Units','Normalized',...
        'position',[0.0557291666666667 0.45 0.127083333333333 0.075],...
        'backg',get(gcf,'color'),'HorizontalAlignment','left','EdgeColor','none',...
        'fontsize',fs-3,'fontweight','bold','Visible','off','tag','crosshair');
else
    hText = annotation(gcf,'textbox','Units','Normalized',...
        'position',[0.0630208333333333 0.0176531671858775 0.113541666666667 0.0581516095534787],...
        'backg',get(gcf,'color'),'HorizontalAlignment','left','EdgeColor','none',...
        'fontsize',fs-3,'fontweight','bold','Visible','off','tag','crosshair');
end

% Get all axes
Ax = findall(gcf,'type','axes');
N = length(Ax);
hCurv = cell(1,N);
hCurh = cell(1,N);

for i = 1:N
    % Set up cursor lines
    hCurv{i} = line(xlim(Ax(i)), ylim(Ax(i)), ...
      'Color', 'red', 'Parent', Ax(i),'Visible','off','tag','crosshair');

    hCurh{i} = line(xlim(Ax(i)), ylim(Ax(i)), ...
      'Color', 'red', 'Parent', Ax(i),'Visible','off','tag','crosshair');
end

%Change the mouse movement callback
set(gcf,'windowbuttonmotionfcn', {@dragFcn,N,Ax,hCurv,hCurh,hText,xtext,xunits,ytext,yunits,x_axis,plot_extra,fs});

    function dragFcn(~,~,N,Ax,hCurv,hCurh,hText,xtext,xunits,ytext,yunits,x_axis,plot_extra,fs)
        for j = 1:N
            try
                % Get mouse location
                pt = get(Ax(j), 'CurrentPoint');
                % Update cursor text
                yl = ylim(Ax(j));
                xl = xlim(Ax(j));
                if pt(1) >= xl(1) && pt(1) <= xl(2) && pt(1,2) >= yl(1) && pt(1,2) <= yl(2)
                  set(hText,'String', sprintf([xtext ' = %3.3f ' xunits newline ytext ' = %6.3f ' yunits],...
                      pt(1), pt(1,2)),'Visible','on');
                  set(gcf,'Pointer','custom');
                  
                    % Update cursor line position
                    set(hCurv{j}, 'YData', ylim(Ax(j)));
                    set(hCurh{j}, 'XData', xlim(Ax(j))); 
                    set(hCurv{j}, 'XData', [pt(1), pt(1)],'Visible','on');
                    set(hCurh{j}, 'YData', [pt(1,2), pt(1,2)],'Visible','on');
                  
                  break
                else
                   set(hText, 'Visible', 'off')
                   set(hCurh{j}, 'Visible', 'off')
                   set(hCurv{j}, 'Visible', 'off')
                   set(gcf,'Pointer','arrow');
                end
            catch
                crosshairfcn([],[],x_axis,plot_extra,fs);
                return
            end
        end
    end
end
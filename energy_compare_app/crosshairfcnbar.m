% This function produces a crosshair that tracks the mouse

function crosshairfcnbar(~,~,x_axis,plot_extra)
if x_axis == "bl"
    xtext = 'BL';
else
    xtext = char(x_axis);
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
    hText = uicontrol('Style','text',...
        'Units','Normalized',...
        'position',[0.049 0.964 0.081291666666667 0.037383177570093],...
        'backg',get(gcf,'color'),'HorizontalAlignment','left',...
        'fontsize',12,'fontweight','bold','Visible','off','tag','crosshair');
else
    hText = uicontrol('Style','text',...
        'Units','Normalized',...
        'position',[0.047 0.963925233644859 0.081291666666667 0.037383177570093],...
        'backg',get(gcf,'color'),'HorizontalAlignment','left',...
        'fontsize',12,'fontweight','bold','Visible','off','tag','crosshair');
end

% Get all axes
Ax = findall(gcf,'type','axes');
N = length(Ax);
hCurv = cell(1,N);
hCurh = cell(1,N);
axlogical = false(1,N);

for i = 1:N
    % Set up cursor lines
    hCurv{i} = line(xlim(Ax(i)), ylim(Ax(i)), ...
      'Color', 'red', 'Parent', Ax(i),'Visible','off','tag','crosshair');
  
    hCurh{i} = line(xlim(Ax(i)), ylim(Ax(i)), ...
      'Color', 'red', 'Parent', Ax(i),'Visible','off','tag','crosshair');
  
    axlogical(i) = strcmp(Ax(i).Children(end).Type,'bar'); % true if bar plot
end

%Change the mouse movement callback
set(gcf,'windowbuttonmotionfcn', {@dragFcn,N,Ax,hCurv,hCurh,hText,xtext,x_axis,plot_extra,axlogical});

    function dragFcn(~,~,N,Ax,hCurv,hCurh,hText,xtext,x_axis,plot_extra,axlogical)
        for j = 1:N
            try
                % Get mouse location
                pt = get(Ax(j), 'CurrentPoint');
                % Update cursor line position
                set(hCurh{j}, 'XData', xlim(Ax(j))); 
                set(hCurh{j}, 'YData', [pt(1,2), pt(1,2)],'Visible','on');
                if ~axlogical(j)
                    set(hCurv{j}, 'YData', ylim(Ax(j)));
                    set(hCurv{j}, 'XData', [pt(1), pt(1)],'Visible','on');
                end
                % Update cursor text
                yl = ylim(Ax(j));
                xl = xlim(Ax(j));
                if pt(1) >= xl(1) && pt(1) <= xl(2) && pt(1,2) >= yl(1) && pt(1,2) <= yl(2)
                    
                    if axlogical(j)
                        set(hText,'String', sprintf([newline 'E = %6.3f kJ/mol'],pt(1,2)),'Visible','on');
                    else
                        set(hText,'String', sprintf([xtext ' = %3.3f ' char(197) newline 'E = %6.3f kJ/mol'],...
                          pt(1), pt(1,2)),'Visible','on');
                    end
                  
                    set(gcf,'Pointer','custom');
                    break
                else
                   set(hText, 'Visible', 'off')
                   set(hCurh{j}, 'Visible', 'off')
                   set(hCurv{j}, 'Visible', 'off')
                   set(gcf,'Pointer','arrow');
                end
            catch
                crosshairfcnbar([],[],x_axis,plot_extra);
                return
            end
        end
    end
end
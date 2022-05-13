function CrosshairOff(~,~)
% Return the mouse movement callback to nothing
set(gcf,'windowbuttonmotionfcn', '')

% Reset cursor
set(gcf,'Pointer','arrow');

% Find any remaining crosshair lines and disable them
fh = findall(gcf,'tag','crosshair');
set(fh,'visible','off');

% Enable the Data cursor and zoom
set(0,'Showhidden','on')
ch = findall(gcf,'tag','FigureToolBar');
UT = get(ch,'children');
UV(1) = findall(UT,'tag','Exploration.DataCursor');
UV(2) = findall(UT,'tag','Exploration.Brushing');
UV(3) = findall(UT,'tag','Exploration.Pan');
UV(4) = findall(UT,'tag','Exploration.ZoomOut');
UV(5) = findall(UT,'tag','Exploration.ZoomIn');
UV(6) = findall(UT,'tag','OptionsFunction');

set(UV,'Enable','On');

return
end
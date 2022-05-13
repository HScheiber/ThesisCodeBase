% This function simply turns off and disables all buttons on the user interface
% fh is the figure handle of the figure whose buttons should be disabled
function fh = DisableUIButtons(fh)
    if ~isvalid(fh)
        return
    end
    % Turn off any buttons that are currently on
    %InterfaceObj=findobj(fh,'State','on');
    %set(InterfaceObj,'State','off');

    % Disable any currently enabled buttons
    InterfaceObj = findobj(fh,'Enable','on','-not','Tag','Integrate.Peaks',...
        '-not','Tag','Integrate.hText','-not','Tag','Integrate.hCurv');    
    set(InterfaceObj,'Enable','off');
    zoom on
    zoom off
end

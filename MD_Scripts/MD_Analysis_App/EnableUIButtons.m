% This function simply turns on and enables all buttons on the user interface
% fh is the figure handle of the figure whose buttons should be enabled
function fh = EnableUIButtons(fh)
    if ~isvalid(fh)
        return
    end
    % Turn on any buttons that are currently off
    %InterfaceObj=findobj(fh,'State','off');
    %set(InterfaceObj,'State','on');

    % Enable any currently disabled buttons
    InterfaceObj=findobj(fh,'Enable','off');
    set(InterfaceObj,'Enable','on');
end
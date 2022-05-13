function varargout = clickableLegend(varargin)
% clickableLegend  Interactive legend for toggling or highlighting graphics
%
% clickableLegend is a wrapper around the LEGEND function that provides
% interactive display toggling or highlighting of lines or patches in a MATLAB
% plot. It enables you to,
% * Toggle (hide/show) a graphics object (or group) by clicking on an
%   associated text label in the legend
% * Highlight and identify a graphics object (or group) by clicking it. 
%
% Its usage is the same as the <a href="matlab: help legend">LEGEND</a> function with additional
% optional parameters. For further information please see the LEGEND documentation.
%
% ADDITIONAL ARGUMENTS specific to clickableLegend:
% These are passed in as parameter-value pairs
%
% * groups: A vector specifying the group membership for every line object.
%   The grouping parameter lets you have a group of lines represented and
%   controlled by one entry in the legend. This can be useful if you would
%   like to highlight or hide a group of lines by clicking on one legend entry.
%   CONSIDERATIONS:
%   # The number of unique entries in the group vector must match the
%     number of strings in the legend
%   # When passing a vector of line/patch handles as the first input
%     argument, the vector should contain all lines/patches whose group
%     membership is included in the group vector. ie. the length of the
%     vector of handles must match that of the length of the group vector.
% 
% * displayedLines: A vector of indices corresponding to the lines or groups that
%   should be displayed initially. This option is useful if there are too
%   many lines in the figure and you are only interested in looking at a
%   few at first.
%   Example: clickableLegend(..., 'displayedLines', [4 5 6])
%
% * plotOptions: A cell array of parameter value pairs that define
%   properties of the graphics objects in the legend. This can be useful,
%   for example, to ensure consistent marker sizes within the legend. 
% 
% Notes: 
% 1. If you save the figure and re-load it, the toggling functionality
% is not automatically re-enabled. To restore it, simply call clickableLegend
% with no arguments.
%
% 2. To prevent the axis from automatically scaling every time a line is
% turned on and off, issue the command: axis manual
%
% Example 1:
% z = peaks(100);
% plot(z(:,26:5:50))
% grid on;
% axis manual;
% clickableLegend({'Line1','Line2','Line3','Line4','Line5'}, 'Location', 'NorthWest');
%
% Example 2:
% f = plot([1:10;1:2:20]','x'); hold on;
% g = plot(sin([1:10;1:2:20]'),'r-');
% h = plot(11:20,rand(5,10)*5,'b:');
% clickableLegend([f;g;h], {'Line1','Line2','Line3'},...
%   'groups', [1 1 2 2 3 3 3 3 3], 'displayedLines', [2 3]);
%
% hgsave(gcf, 'testfig.fig');
% hgload testfig.fig
% clickableLegend
%
% See also legend, clickableLegend_examples

% Copyright 2009-2014 MathWorks, Inc.

% Extract any arguments for clickableLegend
[dispinds, groupmem, plotOptions, varargin] = ...
    extractOptionalArgs(varargin{:});

% Process group memberships
[groups, plotObj, varargin] = processGroups(groupmem, varargin{:});

% Process X and Y labels
[Xtext, Ytext, varargin] = ProcessText(varargin{:});

% Create legend
[varargout{1:nargout(@legend)}] = legend(varargin{:});

% Extract what is needed for the rest of the function and fix varargout
[~, objhan, plothan] = varargout{1:3}; 
% objhan: strings
% plothan: graphics objects
varargout = varargout(1:nargout);

if isempty(groupmem) % Default group membership
    groupmem = 1:length(plothan);
    plotObj = plothan;
    groups = groupmem;
end

if ~isempty(dispinds) % DisplayedLines parameter was specified
    hidden = true(1, length(plothan));
    dispinds(dispinds>length(plothan)) = [];
    hidden(dispinds) = false;
end

% Set the callbacks & plot options
for i = 1:length(plothan)
    set(objhan(i), 'HitTest', 'on', 'ButtonDownFcn',...
        @(varargin)togglevisibility(objhan(i),plotObj(groupmem==groups(i))),...
        'UserData', true);
    if ~isempty(dispinds) && hidden(i)
        togglevisibility(objhan(i), plotObj(groupmem==groups(i)));
    end
    set(plotObj(groupmem==groups(i)), 'HitTest', 'on', 'ButtonDownFcn', ...
        @(varargin)highlightObject(objhan(i),plothan(i),...
                   plotObj(groupmem==groups(i)),plotOptions,Xtext,Ytext),...
                    'UserData', false);
    if ~isempty(plotOptions)
        set(plothan(i), plotOptions{:}); 
    end
end

function togglevisibility(hObject, obj)
if get(hObject, 'UserData') % It is on, turn it off
    set(hObject,'UserData', false);
    set(obj,'Visible','off');
    set(findobj('Tag',obj.DisplayName),'Visible','off');
else
    set(hObject,'UserData', true);
    set(obj,'Visible','on');
    set(findobj('Tag',obj.DisplayName),'Visible','on');
    uistack(obj,'top')
    uistack(findobj('Tag',obj.DisplayName),'top')
end

function highlightObject(lTextObj, lMarkerObj, plotObj, plotOptions,Xtext,Ytext)
lw = get(plotObj,'LineWidth');
if ~iscell(lw), lw = {lw}; end
if strcmp(plotObj.Type,'scatter')
    ms = get(plotObj,'SizeData');
else
    ms = get(plotObj,'MarkerSize');
end
if ~iscell(ms), ms = {ms}; end

if ~get(plotObj(1), 'UserData') % It is not selected, highlight it
    XScale = diff(get(plotObj.Parent,'XLim'));
    YScale = diff(get(plotObj.Parent,'YLim'));
    fs = get(plotObj.Parent,'Fontsize');
    par = get(plotObj,'Parent');
    [Y,idx] = min(get(plotObj,'YData'));
    xdat = get(plotObj,'XData');
    X = xdat(idx);
    txt1 = [sprintf('%6.2f',X) Xtext ' , ' sprintf('%6.1f',Y) Ytext];
    txtObj = text(X-XScale*0.06,Y-YScale*0.025,txt1,'Interpreter','latex','FontSize',fs,...
        'Parent',par,'Visible','on','Fontweight','bold');
    txtObj.Tag = plotObj.DisplayName;
    
    %set(hObject, 'FontWeight', 'bold');
    set(lTextObj, 'EdgeColor', 'k');
    if strcmp(plotObj.Type,'scatter')
        set(plotObj, {'LineWidth', 'SizeData'}, [cellfun(@(x)x+0.499, lw, 'Uniformoutput', false) cellfun(@(x)x+20, ms, 'uniformoutput', false)]);
    elseif strcmp(plotObj.Type,'errorbar')
        set(plotObj, {'MarkerSize'}, cellfun(@(x)x+2, ms, 'uniformoutput', false))
    else
        set(plotObj, {'LineWidth', 'MarkerSize'}, [cellfun(@(x)x+2, lw, 'Uniformoutput', false) cellfun(@(x)x+2, ms, 'uniformoutput', false)]);
    end
    set(plotObj, 'UserData', true);
    uistack(plotObj,'top')
    uistack(txtObj,'top')
else
    %set(hObject, 'FontWeight', 'normal');
    set(lTextObj, 'EdgeColor', 'none');
    if strcmp(plotObj.Type,'scatter')
        set(plotObj, {'LineWidth', 'SizeData'}, [cellfun(@(x)x-0.499, lw, 'Uniformoutput', false) cellfun(@(x)x-20, ms, 'uniformoutput', false)]);
    elseif strcmp(plotObj.Type,'errorbar')
        set(plotObj, {'MarkerSize'}, cellfun(@(x)x-2, ms, 'uniformoutput', false))
    else
        set(plotObj, {'LineWidth', 'MarkerSize'}, [cellfun(@(x)x-2, lw, 'Uniformoutput', false) cellfun(@(x)x-2, ms, 'uniformoutput', false)]);
    end
    set(plotObj, 'UserData', false);
    delete(findobj('Tag',plotObj.DisplayName));
end
if ~isempty(plotOptions)
    set(lMarkerObj, plotOptions{:});
end

function [dispinds, groupmem, plotOpt, varargin] = extractOptionalArgs(varargin)
% Extract the displayedlines and/or groups arguments if specified

ind = find(strcmpi(varargin,'DisplayedLines'));
if ~isempty(ind)
    assert(ind<nargin, 'The DisplayedLines parameter value must be specified');
    dispinds = varargin{ind+1};
    varargin(ind:ind+1) = [];
else
    dispinds = [];
end

ind = find(strcmpi(varargin,'groups'));
if ~isempty(ind)
    assert(ind < nargin, 'The groups parameter value must be specified');
    groupmem = varargin{ind+1};
    varargin(ind:ind+1) = [];
else
    groupmem = [];
end

ind = find(strcmpi(varargin,'plotoptions'));
if ~isempty(ind)
    assert(ind < nargin, 'The plotOptions parameter value must be specified');
    plotOpt = varargin{ind+1};
    varargin(ind:ind+1) = [];
else
    plotOpt = {};
end

function [groups, obj, varargin] = processGroups(groupmem, varargin)
if isempty(groupmem)
    groups = []; obj = [];
    return;
end
if iscellstr(groupmem)
    groupmem = categorical(groupmem);
end
groups = unique(groupmem);
firstmem = zeros(size(groups));
if nargin > 1 && ishandle(varargin{1}(1))
    if strcmpi(get(varargin{1}(1),'Type'),'axes')
        hAxes = varargin{1}(1);
        obj = flipud([findobj(hAxes,'Type','line');findobj(hAxes,'Type','patch')]);
    else % It's a line/patch
        obj = varargin{1};
        [~,firstmem] = ismember(groups, groupmem);
        %for i = 1:length(groups)
        %    firstmem(i) = find(groupmem==groups(i),1);
        %end
        varargin{1} = obj(firstmem);
    end
else
    hAxes = gca;
    obj = flipud([findobj(hAxes,'Type','line');findobj(hAxes,'Type','patch')]);
end


function [Xtext, Ytext, varargin] = ProcessText(varargin)

ind = find(strcmpi(varargin,'X_Axis'));
if ~isempty(ind)
    assert(ind<nargin, 'The X_Axis parameter value must be specified');
    X_setting = varargin{ind+1};
    varargin(ind:ind+1) = [];
    
    if strcmp(X_setting,'a')
        Xtext = '\,\AA';
    elseif strcmp(X_setting,'bl')
        Xtext = '\,\AA';
    elseif strcmp(X_setting,'V')
        Xtext = '\,cm$^{3}$ mol$^{-1}$';
    elseif strcmp(X_setting,'P')
        Xtext = '\,GPa';
    elseif strcmp(X_setting,'d')
        Xtext = '\,g cm$^{-3}$';
    else
        error(['Unknown X Axis type: ' X_setting])
    end
    
else
    dispinds = [];
    Xtext = '\,\AA';
end

ind = find(strcmpi(varargin,'Y_Axis'));
if ~isempty(ind)
    assert(ind < nargin, 'The Y_Axis parameter value must be specified');
    Y_setting = varargin{ind+1};
    varargin(ind:ind+1) = [];
    if strcmp(Y_setting,'L')
        Ytext = '\,kJ mol$^{-1}$';
    elseif strcmp(Y_setting,'C')
        Ytext = '\,kJ mol$^{-1}$';
    elseif strcmp(Y_setting,'P')
        Ytext = '\,GPa';
    else
        error(['Unknown Y Axis type: ' Y_setting])
    end
else
    groupmem = [];
    Ytext = '\,kJ mol$^{-1}$';
end

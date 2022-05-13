function tf = isCharInCellstr(Thing, Cellstr)
%

%   Copyright 2016 The MathWorks, Inc.

if nargin > 0
    Thing = convertStringsToChars(Thing);
end

if nargin > 1
    if isstring(Cellstr)
        Cellstr = cellstr(Cellstr);
    end
end

tf = ischar(Thing) && ismember(Thing, Cellstr);
end
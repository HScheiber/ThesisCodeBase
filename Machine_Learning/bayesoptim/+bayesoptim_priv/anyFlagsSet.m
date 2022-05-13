function tf = anyFlagsSet(setflag)
%

%   Copyright 2016 The MathWorks, Inc.

c = struct2cell(setflag);
tf = any([c{:}]);
end
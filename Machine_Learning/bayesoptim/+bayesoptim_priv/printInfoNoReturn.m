function printInfoNoReturn(key, varargin)
%

%   Copyright 2016 The MathWorks, Inc.

fprintf('%s', bayesoptim.infoString(key, varargin{:}));
end
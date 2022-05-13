function s = infoString(key, varargin)
%

%   Copyright 2016 The MathWorks, Inc.

tag = ['stats:bayesoptim:bayesoptim:' key];
s = message(tag, varargin{:}).getString;
end

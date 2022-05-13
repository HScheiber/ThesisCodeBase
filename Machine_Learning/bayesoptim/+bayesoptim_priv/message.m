function m = message(key, varargin)
%

%   Copyright 2016 The MathWorks, Inc.

if nargin > 0
    key = convertStringsToChars(key);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

tag = ['stats:bayesoptim:bayesoptim:' key];
m = message(tag, varargin{:});
end

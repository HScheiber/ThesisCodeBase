function err(key, varargin)
%

%   Copyright 2016 The MathWorks, Inc.

msg = bayesoptim.message(key, varargin{:});
throwAsCaller(MException(msg.Identifier, getString(msg)));
end
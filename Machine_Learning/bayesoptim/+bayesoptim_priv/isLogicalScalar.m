function tf = isLogicalScalar(x)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isscalar(x) && (islogical(x) || isnumeric(x) && (x==0 || x==1));
end

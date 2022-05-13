function tf = isNonnegativeRealScalar(x)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isscalar(x) && isnumeric(x) && isreal(x) && x>=0;
end

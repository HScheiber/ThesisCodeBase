function tf = isNonnegativeInteger(x)
%

%   Copyright 2016 The MathWorks, Inc.

tf = bayesoptim.isLowerBoundedIntScalar(x, 1);
end
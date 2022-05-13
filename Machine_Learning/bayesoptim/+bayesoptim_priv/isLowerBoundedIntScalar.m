function tf = isLowerBoundedIntScalar(x, LB)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isscalar(x) && isnumeric(x) && isreal(x) && x==floor(x) && x>=LB;
end

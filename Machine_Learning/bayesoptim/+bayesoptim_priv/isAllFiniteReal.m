function tf = isAllFiniteReal(X)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isnumeric(X) && isreal(X) && all(isfinite(X(:)));
end

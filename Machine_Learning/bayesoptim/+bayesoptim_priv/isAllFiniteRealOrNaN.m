function tf = isAllFiniteRealOrNaN(X)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isnumeric(X) && isreal(X) && ~any(isinf(X(:)));
end

function TF = isInteger(X)
%

%   Copyright 2016-2017 The MathWorks, Inc.

TF = isnumeric(X) & isreal(X) & X==floor(X);
end
function tf = isPlusMinusOneOrNaN(x)
%

%   Copyright 2016 The MathWorks, Inc.

tf = isnumeric(x) && (isnan(x) || isequal(x,-1) || isequal(x,1));
end

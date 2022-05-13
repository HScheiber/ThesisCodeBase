function XGrid = NDGrid(LB, UB, NumDivs)
%

%   Copyright 2016 The MathWorks, Inc.

import bayesoptim.*
% Creates an N-Dimensional grid in the box defined by LB and UB
N = length(LB);
if N==0
    XGrid = [];
else
    NMinus1DGrid = NDGrid(LB(2:end), UB(2:end), NumDivs(2:end));
    XGrid = [];
    FirstLinspace = linspace(LB(1), UB(1), NumDivs(1));
    for i = 1:NumDivs(1)
        Block = [repmat(FirstLinspace(i), max(1, size(NMinus1DGrid,1)), 1), ...
            NMinus1DGrid];
        XGrid = [XGrid; Block];
    end
end
end
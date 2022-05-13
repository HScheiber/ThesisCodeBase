function EI = expectedImprovement(X, ObjectiveFcnGP, FBest)
% Returns a vector of EI scores for the rows of the design matrix X.
    
%   Copyright 2016 The MathWorks, Inc.

[PI, FSD, GammaX] = bayesoptim.probabilityOfImprovement(X, ObjectiveFcnGP, FBest, 0);
EI = FSD.*(GammaX.*PI + normpdf(GammaX, 0, 1));
end

function [PI, FSD, GammaX] = probabilityOfImprovement(X, ObjectiveFcnGP, FBest, Margin)
% Returns a vector of PI scores for the rows of the design matrix X.
% Returns the probability that the true function value improves upon FBest
% by at least 'Margin'.
    
%   Copyright 2016 The MathWorks, Inc.

[FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
FSD = sqrt(max(0, YSD.^2 - ObjectiveFcnGP.Sigma.^2));         	% Want SD of F, not Y. Need max to avoid complex sqrt.
GammaX = (FBest - Margin - FMean)./FSD;
PI = normcdf(GammaX, 0, 1);
end

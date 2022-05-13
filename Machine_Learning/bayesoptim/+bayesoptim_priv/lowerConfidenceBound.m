function LCB = lowerConfidenceBound(X, ObjectiveFcnGP, Kappa)
% Returns a vector of LCB scores for the rows of the design matrix X. The
% value returned is the negative of the lower confidence bound.
    
%   Copyright 2016 The MathWorks, Inc.

[FMean, YSD, ~] = predict(ObjectiveFcnGP, X);
LCB = Kappa*YSD - FMean;
end

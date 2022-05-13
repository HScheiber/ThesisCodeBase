function EIpC = expectedImprovementPerSecond(X, ObjectiveFcnGP, ObjectiveEvaluationTimeGP, FBest)
%

%   Copyright 2016 The MathWorks, Inc.

% Returns a vector of EIpC scores for the rows of the design matrix X. This
% function assumes that ObjectiveFcnGP exists but ObjectiveEvaluationTimeGP may be empty.
if isempty(ObjectiveEvaluationTimeGP)
    EIpC = bayesoptim.expectedImprovement(X, ObjectiveFcnGP, FBest);
else
    logObjectiveEvaluationTimePred = predict(ObjectiveEvaluationTimeGP, X);
    EIpC = bayesoptim.expectedImprovement(X, ObjectiveFcnGP, FBest)./exp(logObjectiveEvaluationTimePred);
end
end
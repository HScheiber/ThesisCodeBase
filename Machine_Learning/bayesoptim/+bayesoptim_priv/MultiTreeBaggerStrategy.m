classdef MultiTreeBaggerStrategy < bayesoptim.RFStrategy
    %MULTITREEBAGGERSTRATEGY This class contains methods in which a Multiple Random
    %Forests are used as estimators and cases specific to Random Forests

%   Copyright 2019 The MathWorks, Inc.

    methods
        function RF = iFitRobust(this, BayesOptModel, X, Y, PrivOptions)
            % Fit a Random Forests to the given data.
            minLeafSize = 3;
            n_features = size(X, 2);
            nPredictorsToSample = ceil((1/3) * n_features);
            nTrees = 50;
            if all(isnan(Y))
                bayesoptim.err('AllModelTargetsNaN');
            end
            
            % get indices of categorical variables
            idx = find(strcmp(PrivOptions.VarSpec.Types, 'categorical'));            
            if isempty(BayesOptModel)
                multiTreeBagger = bayesoptim.fitrmultirf(X, Y, nTrees, 'categorical', idx,...
                'MinLeafSize', minLeafSize, ...
                'NumPredictorsToSample', nPredictorsToSample);
                BayesOptModel = bayesoptim.BayesOptRFModel(multiTreeBagger);
            else
                newGroup = X(end, 1);
                BayesOptModel.Model = refit(BayesOptModel.Model, newGroup, X, Y, nTrees,...
                'categorical', idx,...
                'MinLeafSize', minLeafSize, ...
                'NumPredictorsToSample', nPredictorsToSample);
            end
            RF = BayesOptModel;
        end
    end
end


classdef SingleTreeBaggerStrategy < bayesoptim.RFStrategy
    %SingleTreeBaggerStrategy This class contains methods in which a single Random
    %Forests are used as estimators and cases specific to Random Forests

%   Copyright 2020 The MathWorks, Inc.
    
    methods
        function RF = iFitRobust(this, BayesOptModel, X, Y, PrivOptions)
            % Fit a Random Forests to the given data. 
            n_features = size(X, 2);
            nTrees = 50;
            minLeafSize = 3;
            nPredictorsToSample = ceil((1/3) * n_features);
            if all(isnan(Y))
                bayesoptim.err('AllModelTargetsNaN');
            end
            
            % get indices of categorical variables
            idx = find(strcmp(PrivOptions.VarSpec.Types, 'categorical'));
            
            % fit a random forest to the given data
            treeBagger = TreeBagger(nTrees, X, Y, 'Method', 'regression', ...
                'CategoricalPredictors', idx, 'MinLeafSize', minLeafSize, ...
                'NumPredictorsToSample', nPredictorsToSample);
            RF = bayesoptim.BayesOptRFModel(treeBagger);
        end
    end
end


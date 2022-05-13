classdef BayesOptRFModel < bayesoptim.BayesOptModel
    %BayesOptRFModel This class is a wrapper for Random Forests Model
    % This is created since we need two different behaviors depending on
    % Gaussian Process and Random Forests, and avoid extensive code changes
    % in BayesianOptimization Class

%   Copyright 2020 The MathWorks, Inc.
    
    properties

    end
    
    methods
        function obj = BayesOptRFModel(Model)
            %BAYESOPTRFMODEL Construct an instance of this class
            
            % The first parameter (Sigma) is set to zero since there is no
            % equivalent of Gaussian Process Sigma in Random Forests
            obj = obj@bayesoptim.BayesOptModel(0, Model);
        end
        
        function [mean, std, dummy] = predict(obj,inputArg)
            % Some methods in BayesianOptimization Class expects three
            % output values, the first two mean and std have equivalent values in
            % Gaussian Processes. However, Prediction intervals for true
            % response values have no equivalent in Random Forests. Thus,
            % zero is returned for the third return value.
            dummy = 0;  
            [mean, std] = predict(obj.Model, inputArg);
        end
        
        function shouldChoose = shouldChooseRandomPoint(obj, SigmaFTol)
            % Always return False
            shouldChoose = false;
        end
    end
end


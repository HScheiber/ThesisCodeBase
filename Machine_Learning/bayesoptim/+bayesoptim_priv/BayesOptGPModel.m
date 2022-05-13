classdef BayesOptGPModel < bayesoptim.BayesOptModel
    %BayesOptGPModel This class is a wrapper for Gaussian Process Model
    % This is created since we need two different behaviors depending on
    % Gaussian Process and Random Forests

%   Copyright 2020 The MathWorks, Inc.
    
    properties
        KernelInformation;
        Beta;
    end
    
    methods
        function obj = BayesOptGPModel(Sigma, Model, KernelInformation, Beta)
            % Initialize a Gaussian Process Model with Properties that will
            % be accessed by methods in BayesianOptimization class
            obj = obj@bayesoptim.BayesOptModel(Sigma, Model);
            obj.KernelInformation = KernelInformation;
            obj.Beta = Beta;
        end
        
        function [ypred, ysd, yint] = predict(obj,inputArg)
            % Some methods in BayesianOptimization Class expects three
            % output values, same output values of GP Regression 
            % ypred -> Predicted response values
            % ysd -> Estimated standard deviation of new response values
            % yint -> Prediction intervals for the true response values
            [ypred, ysd, yint] = predict(obj.Model, inputArg);
        end
        
        function shouldChoose = shouldChooseRandomPoint(obj, SigmaFTol)
            % Determine whether a random point should be chosen in next
            % trial
            if isa(obj.Model, 'bayesoptim.MultiGP')
                % Follow MultiTreeBagger's convention of always returning
                % false
                shouldChoose = false;
            else
                % Single GP model
                SigmaF = iGetSigmaF(obj.Model);
                shouldChoose = SigmaF/(SigmaF + obj.Sigma) < SigmaFTol;
            end
        end
    end
end

function s = iGetSigmaF(GPR)
    s = GPR.KernelInformation.KernelParameters(end);
end


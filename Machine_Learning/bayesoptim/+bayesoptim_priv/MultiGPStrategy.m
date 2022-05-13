classdef MultiGPStrategy < bayesoptim.GPStrategy
    % MULTIGPSTRATEGY This class contains methods in which a Multiple
    % Gaussian Process model is used as the estimator

%   Copyright 2020 The MathWorks, Inc.

    methods
        function BayesOptModel = iFitRobust(~, BayesOptModel, X, F, SigmaLowerBound, varargin)
            % Fit a Multi Gaussian Process model to the given data. 
            
            % default hyperparameters set in original SMAC paper
            
            if all(isnan(F))
                bayesoptim.err('AllModelTargetsNaN');
            end
            
            % Fit a MultiGP model         
            if isempty(BayesOptModel)
                multiGP = bayesoptim.fitrmultigp(X, F, SigmaLowerBound, varargin{:});
                GPModel = multiGP.Models;
                
                % There should be at least one (non-empty) cell in GPModel,
                % given that BayesOptModel is empty. Typically this is the
                % first model, however, if the first objective function was
                % a NaN, this could be a later cell value.
                idx = ~(cellfun('isempty',GPModel));
                if sum(idx) < 1 % Ensure that at least one model is non-empty 
                    BayesOptModel = [];
                else
                    idx = find(idx,1,'first');
                    GPModel = GPModel{idx};
                    BayesOptModel = bayesoptim.BayesOptGPModel(GPModel.Sigma, multiGP, ...
                        GPModel.KernelInformation, GPModel.Beta);
                end
            else
                newGroup = X(end, 1);
                BayesOptModel.Model = refit(BayesOptModel.Model, newGroup, X, F, SigmaLowerBound, varargin{:});
            end
        end
    end
end


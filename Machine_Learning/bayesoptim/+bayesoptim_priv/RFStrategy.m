classdef RFStrategy < bayesoptim.ModelStrategy
    %BAYESIANOPTIMIZATIONRF This class contains methods in which Random
    %Forests are used as estimators and cases specific to Random
    %Forests

%   Copyright 2020 The MathWorks, Inc.
    
    methods (Access = private)

        function ConstraintRF = fitConstraintModel(this, RF, ConstraintTrain, isDeterministic, useFitMethodNone, PrivOptions, XTrain)
            if all(isnan(ConstraintTrain))
                ConstraintRF = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitConstraint');
                end
            elseif useFitMethodNone && ~isempty(RF)
                ConstraintRF = refitWithFitMethodNone(RF,  Model, XTrain, ...
                    ConstraintTrain, PrivOptions);
            else
                ConstraintRF = iFitRobust(this, RF, XTrain, ConstraintTrain, ...
                    PrivOptions);
            end
        end
        
    end
    
    methods
        function ObjectiveFcnRF = fitObjectiveFcnModel(this, XTrain, FTrain, PrivOptions, Model, varargin)
            useFitMethodNone = varargin{1}{1};
            if all(isnan(FTrain))
                ObjectiveFcnRF = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitObjective');
                end
            elseif useFitMethodNone && ~isempty(this)
                ObjectiveFcnRF = refitWithFitMethodNone(this, Model, XTrain, FTrain, PrivOptions);
            else
                ObjectiveFcnRF = iFitRobust(this, Model, XTrain, FTrain, PrivOptions);
            end
        end
        
        function ObjectiveEvaluationTimeRF = fitObjectiveEvaluationTimeModel(this, XTrain, FTrain, PrivOptions, ObjectiveEvaluationTimeTrain, Model, varargin)
            useFitMethodNone = varargin{1}{1};
            % Only fit model to successful function evaluations
            FitIndices = ~isnan(ObjectiveEvaluationTimeTrain) & ~isnan(FTrain);
            if ~any(FitIndices)
                ObjectiveEvaluationTimeRF = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitObjectiveEvaluationTime');
                end
            elseif useFitMethodNone && ~isempty(this)
                ObjectiveEvaluationTimeRF = refitWithFitMethodNone(this, Model, ...
                                                                 XTrain(FitIndices,:), ...
                                                                 log(ObjectiveEvaluationTimeTrain(FitIndices)),...
                                                                 PrivOptions);
            else
                ObjectiveEvaluationTimeRF = iFitRobust(this, Model, XTrain(FitIndices, :), ...
                    log(ObjectiveEvaluationTimeTrain(FitIndices)), ...
                    PrivOptions);
            end
        end
        
        function ErrorRF = fitErrorModel(this,  XTrain, FTrain, PrivOptions, ErrorTrain, Model, varargin)
            useFitMethodNone = varargin{1}{1};
            
            if all(isnan(ErrorTrain))
                ErrorRF = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitError');
                end
            elseif all(ErrorTrain < 0)
                ErrorRF = [];
            elseif useFitMethodNone && ~isempty(this)
                ErrorRF = refitWithFitMethodNone(this, Model, XTrain, ...
                    ErrorTrain, PrivOptions);
            else
                ErrorRF = iFitRobust(this, Model, XTrain, ErrorTrain, PrivOptions);
            end
        end
        
        function ConstrtRFs = fitConstraintModels(this, XTrain, ConstraintTrain, ConstraintGPs, PrivOptions, varargin)
            useFitMethodNone = varargin{1}{1};
            ConstrtRFs = cell(1,PrivOptions.NumCoupledConstraints);
            if ~isempty(ConstraintTrain)
                for i = PrivOptions.NumCoupledConstraints:-1:1
                    if useFitMethodNone && ~isempty(ConstraintGPs{i})
                        ConstrtRFs{i} = fitConstraintModel(this, ConstraintGPs{i},...
                            ConstraintTrain(:,i), PrivOptions.AreCoupledConstraintsDeterministic(i), ...
                            useFitMethodNone, PrivOptions, XTrain);
                    else
                        ConstrtRFs{i} = fitConstraintModel(this, [],...
                            ConstraintTrain(:,i), PrivOptions.AreCoupledConstraintsDeterministic(i), ...
                            useFitMethodNone, PrivOptions, XTrain);
                    end
                end
            end
        end
        
        function RF = iFitRobust(this, BayesOptModel, X, Y, PrivOptions)
            % Fit a Random Forests to the given data. 
            n_features = size(X, 2);
            
            % default hyperparameters set in original SMAC paper
            nTrees = 10;
            minLeafSize = 10;
            nPredictorsToSample = ceil((5 / 6) * n_features);
            
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
            
%             if isempty(BayesOptModel)
%                 multiTreeBagger = fitrmultirf(X, Y, nTrees, 'cat', idx);
%                 BayesOptModel = bayesoptim.BayesOptRFModel(multiTreeBagger);
%             else
%                 newGroup = X(end, 1);
%                 BayesOptModel.Model = refit(BayesOptModel.Model, newGroup, X, Y, nTrees,'cat', idx);
%             end
%             RF = BayesOptModel;
        end

        function RF = refitWithFitMethodNone(this, RF, XTrain, FTrain, PrivOptions)
            % In Gaussian Process Strategy, iFitRobust and
            % refitWithFitMethodNone are different. However, this does not
            % apply to Random Forests.
            RF = iFitRobust(this, RF, XTrain, FTrain, PrivOptions);
        end
    end
end




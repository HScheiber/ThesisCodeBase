classdef RegressionMultiModelOptim < bayesoptim.MultiModelOptim
    %ClassificationMultiModelOptim This class performs bayesian
    %optimization across classification models and their hyperparameters

%   Copyright 2020 The MathWorks, Inc.
    
    properties(GetAccess=public,Constant)
        % Cell array that stores all possible classification learner names
        AllLearners = {'svm','tree','linear','ensemble','kernel','gp'};
        
        % Cell array that stores all possible fitrfunction names
        FitFunctions = {'fitrsvm','fitrtree','fitrlinear','fitrensemble',...
            'fitrkernel','fitrgp'};
    end
    
    methods
        % Constructor
        function this = RegressionMultiModelOptim(X,Y,multiModelOptimOptions)          
            this.X = X;
            this.Y = Y;
            this.Learners = multiModelOptimOptions.Learners;
            this.OptimizeHyperparameters = multiModelOptimOptions.OptimizeHyperparameters;
            this.OptimizationOptions = multiModelOptimOptions.OptimizationOptions;
            this.FitFunctionArgs = multiModelOptimOptions.RemainingArgs;
            this.UserPassedMaxObjectiveEvaluations = multiModelOptimOptions.UserPassedMaxObjectiveEvaluations;
            this.FitFunctionMap = containers.Map(this.AllLearners,this.FitFunctions);
            this.IrrelevantHyperParamsMap = [];
            
            [this.PredictorsMatrix,this.ResponseVector,~,dataSummary] = ...
                classreg.learning.regr.FullRegressionModel.prepareData(X, Y, this.FitFunctionArgs{:});
            
            this.DataSummary = dataSummary;
            updateLearners(this);
            checkObservationsInRows(this);
            checkForSparseData(this);
            checkAndUpdateValidationMethod(this);
        end
    end
    
    methods(Hidden)
        function loss = fitModel(this,parameterCombination)
        %FITMODEL fit the specified learner and return the loss value 
        %   LOSS = FITMODEL(THIS,PARAMETERCOMBINATION) fits the specified
        %   learner and its hyperparameters. It returns the loss value
        %   obtained by the fitting the learner to the input data.
            if size(this.Learners, 2) > 1
                model = char(parameterCombination.learner);
            else
                model = this.Learners{:};
            end

            % get the hyperparameters for the current algorithm
            optimizableHyperparams = getArgumentsTable(this,this.HyperparameterNames.(model), parameterCombination);
            nonOptimizableLearnerArgPairs = {};
            

            if isfield(this.NonOptimizableNVPs,model)
                nonOptimizableLearnerArgPairs = this.NonOptimizableNVPs.(model);
                if any(strcmp(model,{'svm'}))
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                elseif strcmp(model,'ensemble')
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.removeLearnrateForBaggedEnsemble(parameterCombination,nonOptimizableLearnerArgPairs);
                end
            end
            
            if strcmp(model,'ensemble')
                weakLearnerNonOptimizableNVPs = {};
                if isfield(this.WeakLearnerNonOptimizableNVPs,'tree')
                    weakLearnerNonOptimizableNVPs = this.WeakLearnerNonOptimizableNVPs.('tree');
                    [nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs] = ...
                        bayesoptim.MultiModelOptim.AddNumBinsForBoostedEnsemble(parameterCombination,nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs);
                end
                fitFunArgs = [{'Learners'},{templateTree(weakLearnerNonOptimizableNVPs{:})},this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            elseif strcmp(model,'linear')
                nonOptimizableLearnerArgPairs=AddSolverForLinear(this,optimizableHyperparams);
                fitFunArgs = [this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            else
                fitFunArgs = [this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            end

            % get objective function of current algorithm
            objFcn = getObjectiveFcn(this,this.BayesOptInfo.(model), fitFunArgs);
            loss = objFcn(optimizableHyperparams);
        end
    end
    
    methods(Access = {?tRegressionMultiModelOptim, ?bayesoptim.MultiModelOptim})
        
        function model = trainBestModel(this)        
        %TRAINBESTMODEL trains and returns the best model  
        %   MODEL = TRAINBESTMODEL(THIS) finds the best
        %   hyperparameters and trains the learner with best hyperparameter
        %   and returns the full classification model.
            
        
            [~, ~, MinObsLoc] = bestPoint(this.OptimizationResults, ...
                'Criterion', 'min-observed');
            
            if ~isfinite(MinObsLoc)
                error(message('stats:bayesoptim:bayesoptim:NoFeasibleResultFitauto')); 
            end
            
            optimizedHyperparameters = bestPoint(this.OptimizationResults,...
            'Criterion','min-visited-mean');
        
            if numel(this.Learners) > 1
                learner = char(optimizedHyperparameters.learner);
            else
                learner =this.Learners{:};
            end
            hyperparameters = optimizedHyperparameters.Properties.VariableNames;
            relevantHyperparameters = hyperparameters(startsWith(hyperparameters, learner));
            
            learnerArgPairs = {};
            weakLearnerArgPairs = {};
            
            for idx = 1 : size(relevantHyperparameters, 2)
                components = split(relevantHyperparameters{idx}, '_');
                hyperparameterName = components{end};

                % Extract weakLearner specific arguments if the model is
                % ensemble
                if ismember(relevantHyperparameters{idx}, {'ensemble_MinLeafSize', ...
                        'ensemble_MaxNumSplits', 'ensemble_SplitCriterion', 'ensemble_NumVariablesToSample'})
                    weakLearnerArgPairs = bayesoptim.MultiModelOptim.updateHyperparameterNVPs(...
                        weakLearnerArgPairs,hyperparameterName,...
                        relevantHyperparameters{idx},optimizedHyperparameters);
                end

                % Extract arguments specific to the classifier
                if ~strcmp(hyperparameterName, 'Coding') && ~strcmp(learner,'gp')...
                    && ~ismember(relevantHyperparameters{idx}, {'ensemble_MinLeafSize', ...
                        'ensemble_MaxNumSplits', 'ensemble_SplitCriterion', 'ensemble_NumVariablesToSample'})
                    learnerArgPairs = bayesoptim.MultiModelOptim.updateHyperparameterNVPs(learnerArgPairs,...
                        hyperparameterName,relevantHyperparameters{idx},...
                        optimizedHyperparameters);
                end
            end
            
            if strcmp(learner,'gp')
               % GP optimizes over KernelScale but does not support it as 
               % a Name-Value pair. Use BayesoptInfoGP's
               % updateArgsFromTable method to obtain args without
               % KernelScale
               % Remove "gp_" from each Name
               optimizedHyperparameters = optimizedHyperparameters(:,relevantHyperparameters);
               varNames = erase(optimizedHyperparameters.Properties.VariableNames,...
                   'gp_');
               optimizedHyperparameters.Properties.VariableNames = varNames;
               learnerArgPairs = this.BayesOptInfo.gp.updateArgsFromTable({},...
                   optimizedHyperparameters);
            end
            
            nonOptimizableLearnerArgPairs = {};
            
            if isfield(this.NonOptimizableNVPs,learner)
                nonOptimizableLearnerArgPairs = this.NonOptimizableNVPs.(learner);
            end
            
            weakLearnerNonOptimizableNVPs = {};

            
            if strcmp(learner,'ensemble')
                if isfield(this.WeakLearnerNonOptimizableNVPs,'tree')
                    weakLearnerNonOptimizableNVPs = this.WeakLearnerNonOptimizableNVPs.('tree');
                end
                nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.removeLearnrateForBaggedEnsemble(optimizedHyperparameters,nonOptimizableLearnerArgPairs);
                [nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs] = ...
                bayesoptim.MultiModelOptim.AddNumBinsForBoostedEnsemble(optimizedHyperparameters,nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs);
            end
           
            model = trainModels(this,learner,learnerArgPairs,nonOptimizableLearnerArgPairs,...
                weakLearnerArgPairs,weakLearnerNonOptimizableNVPs);
        end
        
        function updateSearchRangeToRemoveARDKernel(this)
            for i = 1:size(this.HyperparameterSearchRange,2)
                splits = split(this.HyperparameterSearchRange(i).Name,'_');
                
                if ~(numel(splits) == 2)
                    continue;
                end
                learner = splits{1};
                hyperparameter = splits{2};
                if strcmp(hyperparameter,'KernelFunction') && strcmp(learner,'gp')
                    KernelFunctionRange=this.HyperparameterSearchRange(i).Range;
                    this.HyperparameterSearchRange(i).Range(find(strncmp(KernelFunctionRange,'ard',3)))=[];
                    break;
                end
            end
        end
        
        function XTable = handleConditionalVariablesGPR(~,XTable)
        %HANDLECONDITIONALVARIABLESGPR  Gaussian Process conditional variable
        %function
        %   XTABLE = HANDLECONDITIONALVARIABLESGP(THIS,XTABLE)
        %   applies conditional constraints related to a KNN
        % Set KernelScale to NaN when KernelFunction is any ARD kernel.
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'gp_KernelFunction', 'gp_KernelScale'})
                ARDRows = ismember(XTable.gp_KernelFunction, {'ardexponential','ardmatern32',...
                    'ardmatern52','ardrationalquadratic','ardsquaredexponential'});
                XTable.gp_KernelScale(ARDRows) = NaN;
            end
        end
        
        function XTable = handleConditionalVariablesKernel(~,XTable)
            % Epsilon is irrelevant if Learner~='svm'
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'kernel_Learner', 'kernel_Epsilon'})
                XTable.kernel_Epsilon(XTable.kernel_Learner ~= 'svm') = NaN;
            % when observations <100000 kernel with learner svm is removed
            elseif ~bayesoptim.MultiModelOptim.hasVariables(XTable, {'kernel_Learner'})
                [r,~]=size(XTable);
                XTable.kernel_Epsilon = NaN(r,1);
            end
        end
        
        function updateLearners(this)
        %UPDATELEARNERS  update the learner types for optimization
        %   UPDATELEARNERS(THIS) checks for the learner options
        %   'auto','all','alllinear' and 'allnonlinear' and updates the
        %   learner types accordingly. If user passes specific learner
        %   types they are used directly without any change.
            learners = this.Learners{:};
            switch learners
                case 'all-linear'
                    this.LearnersIsAuto = false;
                    this.LearnersIsLinear = true;
                    this.Learners = {'linear'};
                case 'all-nonlinear'
                    this.LearnersIsAuto = false;
                    this.LearnersIsNoNLinear = true;
                    this.Learners = {'svm','tree','ensemble','kernel','gp'};
                case 'all'
                    this.LearnersIsAll = true;
                    this.LearnersIsAuto = false;
                    this.Learners = this.AllLearners;
                case 'auto'
                     dataProp = dataProperties(this);
                     updateLearnersForAuto(this,dataProp);
                otherwise
                    this.LearnersIsAuto = false;
            end
            checkIfLearnersSupportMixedData(this);
        end

        function model = trainModels(this,learner,LearnerArgPairs,nonOptimizableLearnerArgPairs,...
                weakLearnerArgPairs,weakLearnerNonOptimizableNVPs)
        %TRAINMODELS  trains a non ecoc model
        %   TRAINMODELS(THIS,LEARNER,LEARNERARGPAIRS,WEAKLEARNERARGPAIRS) 
        %   trains and returns a nonecoc model. The LEARNERARGPAIRS and
        %   THIS.FITFUNCTIONARGS are passed to the fit function. Only
        %   ensemble models use the additional WEAKLEARNERARGPAIRS

            switch learner
                case 'linear'
                    model = fitrlinear(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'svm'
                    LearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(LearnerArgPairs);
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                    model = fitrsvm(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'kernel'
                    model= fitrkernel(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'gp'
                    LearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(LearnerArgPairs);
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                    model = fitrgp(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'tree'
                    model = fitrtree(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'ensemble'
                    model = fitrensemble(this.X, this.Y, LearnerArgPairs{:}, nonOptimizableLearnerArgPairs{:},'Learners', ...
                        templateTree(weakLearnerArgPairs{:},weakLearnerNonOptimizableNVPs{:}),this.FitFunctionArgs{:});
            end
        end
        
        function dataProp =  dataProperties(this)
            dataProp = struct();
            [N,D] = size(this.PredictorsMatrix);
            dataProp.hasCategorical = ~isempty(this.DataSummary.CategoricalPredictors);
            isValueMissing = ismissing(this.PredictorsMatrix);
            fracObsWithMissingVals = sum(any(isValueMissing,2))/N;
            dataProp.hasManyMissing = fracObsWithMissingVals > 0.05;
            dataProp.isWide = D >= N;
            dataProp.isHighD = D > 100;
            dataProp.isBig = N > 50000;
        end
        
        function updateLearnersForAuto(this,dataProp)
            hasManyMissing =  dataProp.hasManyMissing;
            hasCategorical = dataProp.hasCategorical;
            isWide = dataProp.isWide;
            isBig = dataProp.isBig;
            
            if hasManyMissing
               surrogate = 'on';
               addTreeModel(this,surrogate);
               addBoostingModels(this,surrogate);
               addBaggingModels(this,isBig,surrogate);
            elseif hasCategorical
               surrogate = 'off';
               addTreeModel(this,surrogate);
               addBoostingModels(this,surrogate);
               addBaggingModels(this,isBig,surrogate);
               addSVMModel(this,'polynomial');
            elseif isWide                
                addLinearModel(this)
            else
                surrogate = 'off';
                addTreeModel(this,surrogate);
                addBoostingModels(this,surrogate);
                addBaggingModels(this,isBig,surrogate);
                addSVMModel(this,'rbf')
                addGPModel(this);
            end
            
        end
               
        function addBoostingModels(this,surrogate)
            numLearningCycles = 300;
            maxNumSplits = 10;
            addBoostingEnsembleModel(this,'LSBoost',numLearningCycles,maxNumSplits,surrogate);
        end
        
        function addBoostingEnsembleModel(this,method,numLearningCycles,maxNumSplits,surrogate)
            learnerName = 'ensemble';
            NoNOptimizableHyperparameters = [];
            weakLearnerName = 'tree';
            
            if strcmp(this.OptimizeHyperparameters,'auto')
                OptimizableHyperparameters = {'Method',method,...
                    'NumLearningCycles',numLearningCycles,'LearnRate',0.1};
                weakLearnerOptimizableHyperparameters = {'MaxNumSplits',maxNumSplits};
            else
                OptimizableHyperparameters = {'Method',method};
                weakLearnerOptimizableHyperparameters = {'MaxNumSplits',maxNumSplits};
            end
            
            weakLearnerNoNOptimizableHyperparameters = {'Surrogate',surrogate};
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        end
        
        function addBaggingEnsembleModel(this,numLearningCycles,surrogate,predictorSelection)
            learnerName = 'ensemble';
            OptimizableHyperparameters = {'Method','Bag',...
                'NumLearningCycles',numLearningCycles};
            NoNOptimizableHyperparameters = [];
            weakLearnerName = 'tree';
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = {'Surrogate',surrogate,...
                'PredictorSelection',predictorSelection};
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        end
        
        function addSVMModel(this,kernelFunction)
            learnerName = 'svm';
            OptimizableHyperparameters = {'KernelFunction',kernelFunction,'Standardize',true};            
            if strcmp(kernelFunction,'linear')
                OptimizableHyperparameters{end + 1}  = 'KernelScale';
                OptimizableHyperparameters{end + 1}  = 1;
            end
            
            NoNOptimizableHyperparameters = [];
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        
        end
        
        function addLinearModel(this)
            learnerName = 'linear';
            OptimizableHyperparameters = [];
            % removing the solver set to lbfgs 
            NoNOptimizableHyperparameters = {};
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters)
        end
        
        function addTreeModel(this,surrogate)
            learnerName = 'tree';
            OptimizableHyperparameters = [];
            NoNOptimizableHyperparameters = {'Surrogate',surrogate,...
                'MinParentSize',1,...
                'PredictorSelection','interaction-curvature'};
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,...
                weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        end
        
        function addGPModel(this)
            learnerName = 'gp';
            OptimizableHyperparameters = {'Standardize',true};
            NoNOptimizableHyperparameters = [];
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,...
                weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        end
    end
end
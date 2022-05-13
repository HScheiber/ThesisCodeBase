classdef (Abstract) MultiModelOptim < handle
    %MultiModelOptim This class performs bayesian
    %optimization across classification models and their hyperparameters

%   Copyright 2020 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess={?bayesoptim.MultiModelOptim,?tClassificationMultiModelOptim,...
            ?tRegressionMultiModelOptim})
        
        % Either a Table containing predictors and response, or a predictor
        % table or a predictor matrix
        X = [];
        
        % Either a char representing the response name if the input is a
        % table or a formula using predictors to obtain a response or a
        % response vector
        Y = [];
        
        % Either 'auto','all','all-linear' or 'all-nonlinear' or a cell array
        % of eligible learner names
        Learners = 'auto';
        
        % Scalar that stores number of classes in the response
        NumClasses = 0;
               
        % A map to store the learner and the corresponding fitc* function
        FitFunctionMap = [];
        
        % Either 'auto' or 'all', represents the set of hyperparameters to
        % be optimized
        OptimizeHyperparameters = 'auto'
        
        % Cell array that stores name value pairs to be passed to the fit
        % functions
        FitFunctionArgs = [];
        
        % Struct that the hyperparameter optimization options
        OptimizationOptions = [];
        
        % A map that stores the learners and hyperparameters relevant to it
        IrrelevantHyperParamsMap = [];
        
        % A cell array storing the BayesoptInfo objects for every learner
        % being optimized
        BayesOptInfo = [];
        
        % A cell array storing the hyperparameter names with the learner name
        % added to their names
        HyperparameterNames = [];
        
        % A array of OptimizableVariable objects that store the
        % hyperparameters to be optimized and their range
        HyperparameterSearchRange = [];
        
        % Matrix of predictors
        PredictorsMatrix = [];
        
        % True class labels used to train
        ResponseVector = [];
        
        % An object or table describing the results of hyperparameter optimization. 
        OptimizationResults = [];
        
        % A struct holding the summary of input data
        DataSummary = [];
        
        % A struct holding the non optimizable hyperparameter name value pairs for
        % every learner
        NonOptimizableNVPs = struct();        
         
        % A struct holding the non optimizable hyperparameter name value pairs for
        % every weakLearner
        WeakLearnerNonOptimizableNVPs  = struct();
        
        % A struct holding the hyperparameter search range when Learners is
        % 'auto'
        AutoSearchRange = struct();
        
        % A boolean indicating if the Learners is
        % 'all'
        LearnersIsAll = false;
        
        % A boolean indicating if the Learners is
        % 'auto'
        LearnersIsAuto = true;
        
        % A boolean indicating if the Learners is
        % 'linear'
        LearnersIsLinear = false;
        
        % A boolean indicating if the Learners is
        % 'nonlinear'
        LearnersIsNoNLinear = false;
        
        % A boolean to indicate if user passed MaxObjectiveEvaluations
        UserPassedMaxObjectiveEvaluations = false;
    end
    
    properties(Abstract,GetAccess=public,Constant)
        % Cell array that stores all possible classification learner names
        AllLearners
        
        % Cell array that stores all possible fitcfunction names
        FitFunctions
    end
    
    methods(Hidden)     
        function[optimizedModel,optimizationResults] = performOptimization(this)
            initializeLearnersSearchRange(this);
            runBayesopt(this);
            optimizedModel = compact(trainBestModel(this));
            optimizationResults = this.OptimizationResults;
        end
    end
    
    
    methods(Access = {?bayesoptim.MultiModelOptim,?tClassificationMultiModelOptim,...
            ?tRegressionMultiModelOptim})
        
        function initializeLearnersSearchRange(this)
        %INITIALIZELEARNERSSEARCHRANGE initialize the search range for all
        %learners
        %   INITIALIZELEARNERSSEARCHRANGE(THIS) 
        %   initializes the hyperparameter search range for all the learners and 
        %   generates BayesoptInfo objects for each learner.
            
            if this.LearnersIsAuto
                this.Learners = fieldnames(this.AutoSearchRange)';
            end
        
            % remove any duplicates passed by user
            this.Learners = unique(this.Learners);

            if size(this.Learners, 2) > 1
                this.HyperparameterSearchRange = optimizableVariable('learner',this.Learners, ...
                'Type', 'categorical');
            end

            % make it an empty struct
            this.HyperparameterNames = {};
            this.BayesOptInfo = {};
           
            % get optimizable variables for all learners
            for jdx = 1 : size(this.Learners, 2)
                if this.NumClasses > 2 && ismember(this.Learners{jdx}, {'kernel', 'svm', 'linear'})
                    % ECOC does not apply to regression. This is taken care of by NumClasses = 0
                    boInfo = classreg.learning.paramoptim.BayesoptInfo.makeBayesoptInfo('fitcecoc', ...
                             this.X,this.Y, {'Learners', this.Learners{jdx},this.FitFunctionArgs{:}});
                else
                    boInfo = classreg.learning.paramoptim.BayesoptInfo.makeBayesoptInfo(this.FitFunctionMap(this.Learners{jdx}), ...
                             this.X, this.Y, this.FitFunctionArgs); 
                end
                
                % add prefix to hyperparameter names to differentiate
                % hyperparameters with the same name from different
                % learner
                VariableDescriptions = getVariableDescriptions(boInfo, this.OptimizeHyperparameters);
                for k = 1 : size(VariableDescriptions, 1)
                    VariableDescriptions(k).Name = [this.Learners{jdx}, '_', VariableDescriptions(k).Name];
                end

                this.HyperparameterSearchRange = [this.HyperparameterSearchRange, VariableDescriptions'];

                % store the names of hyperparameters
                this.HyperparameterNames.(this.Learners{jdx}) = {VariableDescriptions.Name};
                
                % store bayesOptinfo of each learner
                this.BayesOptInfo.(this.Learners{jdx}) = boInfo;
            end
            
            if isempty(this.Learners)
                error(message('stats:bayesoptim:multimodeloptimoptions:NoLearner')); 
            end
            
            % remove ARD Kernel from the search range for Regression
            % Problem
            if(this.OptimizationOptions.IsClassregRegressionFunction)
                updateSearchRangeToRemoveARDKernel(this);
            end
            
            if this.LearnersIsAuto
                updateSearchRangeForAuto(this);
            end
            
            if this.LearnersIsAll
                updateSearchRangeForAll(this);
            end
            
            if this.LearnersIsLinear
                updateSearchRangeForAllLinear(this);
            end
            
            if this.LearnersIsNoNLinear
                updateSearchRangeForAllNoNLinear(this);
            end
            
            if ~this.UserPassedMaxObjectiveEvaluations
                this.OptimizationOptions.MaxObjectiveEvaluations = numel(this.Learners)*30;
            end
        end
        
        function runBayesopt(this)
        %RUNBAYESOPT perform Bayesian optimization
        %   RUNBAYESOPT(THIS) runs the Bayesian optimization across 
        %   all the learners and their hyperparameters
            
            fun = @(parameterCombination) fitModel(this,parameterCombination);
            condtionalVaraibleFunction = @(XTable) ConditionalVariableFunction(this,XTable);

            xConstraintFunction = @(XTable) XConstraintFunction(this,XTable);
            
            if this.OptimizationOptions.SaveIntermediateResults
                ouputFcn = @assignInBase;
            else
                ouputFcn = {};
            end
            
            if this.OptimizationOptions.ShowPlots
                plotFcn = @plotMinObjective;
            else
                plotFcn = {};
            end
            
            if this.OptimizationOptions.IsClassregRegressionFunction ==1
                this.OptimizationOptions.ModelType = 'MultiGaussianProcess';
                if numel(this.Learners) == 1
                    this.OptimizationOptions.ModelType = 'GaussianProcess';
                    fitautoMultipleLearners = false;
                    fitAutoSingleLearner = true;
                else
                    fitautoMultipleLearners = true;
                    fitAutoSingleLearner = false;
                end
            else
                this.OptimizationOptions.ModelType = 'MultiTreeBagger';
                if numel(this.Learners) == 1
                    this.OptimizationOptions.ModelType = 'SingleTreeBagger';
                    fitautoMultipleLearners = false;
                    fitAutoSingleLearner = true;
                else
                    fitautoMultipleLearners = true;
                    fitAutoSingleLearner = false;
                end
            end
            
            this.OptimizationResults = bayesopt(fun, this.HyperparameterSearchRange,...
                'AcquisitionFunctionName',  this.OptimizationOptions.AcquisitionFunctionName, ...
                'ModelType', this.OptimizationOptions.ModelType, ...
                'UseParallel', this.OptimizationOptions.UseParallel,...
                'MaxObjectiveEvaluations',this.OptimizationOptions.MaxObjectiveEvaluations, ...
                'MaxTime', this.OptimizationOptions.MaxTime, ...
                'ConditionalVariableFcn', condtionalVaraibleFunction,...
                'XConstraintFcn',xConstraintFunction,...
                'Verbose',this.OptimizationOptions.Verbose,...
                'OutputFcn', ouputFcn,...
                'PlotFcn',plotFcn,...
                'FitcautoMultipleLearners',fitautoMultipleLearners,...
                'FitcAutoSingleLearner',fitAutoSingleLearner,...
                'NumSeedPoints',15,...
                'IsClassregRegressionFunction',this.OptimizationOptions.IsClassregRegressionFunction);
        end
        
        function objFcn = getObjectiveFcn(this,BOInfo,fitFunArgs)   
        %GETOBJECTIVEFUNCTION returns the objective function for a fit function  
        %   OBJFCN = GETOBJECTIVEFUNCTION(THIS,BOINFO,FITFUNARGS) returns
        %   the objective function specified in BOINFO.FitFcn.
            objFcn = classreg.learning.paramoptim.createObjFcn(BOInfo, fitFunArgs, ...
                    this.X, this.Y, ...
                    this.OptimizationOptions.ValidationMethod, ...
                    this.OptimizationOptions.ValidationVal, ...
                    this.OptimizationOptions.Repartition, ...
                    this.OptimizationOptions.Verbose);
        end

        function args = getArgumentsTable(~,argNames, bayesoptParams)
            
            args = struct();
            for idx = 1 : size(argNames, 2)
                components = split(argNames{idx}, '_');
                argName = components{end};
                if ismember(argNames{idx}, bayesoptParams.Properties.VariableNames)
                    if ~ismissing(bayesoptParams.(argNames{idx}))
                        if iscategorical(bayesoptParams.(argNames{idx}))
                            args.(argName) = char(bayesoptParams.(argNames{idx}));
                        else
                            args.(argName) = bayesoptParams.(argNames{idx});
                        end
                    end
                end
            end
            args = struct2table(args);
        end

        function XTable = ConditionalVariableFunction(this,XTable)
        %CONDITIONALVARIABLEFUNCTION  conditionalVariableFunction for the entire
        %optimization process.
        %   XTABLE = CONDITIONALVARIABLEFUNCTION(THIS,XTABLE) invokes
        %   conditional variable functions for individual learners and
        %   updates valid hyperparameters for each iteration of the
        %   entire optimization process.

            learnerNames = fieldnames(this.HyperparameterNames);
            if isempty(this.IrrelevantHyperParamsMap)
                values = cellfun(@(x) this.HyperparameterNames.(x), fieldnames(this.HyperparameterNames), 'UniformOutput', 0);
                params = horzcat(values{:});
                this.IrrelevantHyperParamsMap = containers.Map;
                for idx = 1 : numel(learnerNames)
                    learner = learnerNames{idx};
                    irrelevantHyperParams = params(~ismember(params, this.HyperparameterNames.(learner)));
                    existingIrrelevantHyperParams = ismember(XTable.Properties.VariableNames, irrelevantHyperParams);
                    this.IrrelevantHyperParamsMap(learner) = existingIrrelevantHyperParams;
                end
            end

            % filter irrelevant hyperparameters that are not from each learner
            if numel(this.Learners) > 1
                if height(XTable) == 1
                        learner = char(XTable.learner);
                        XTable{1, this.IrrelevantHyperParamsMap(learner)} = missing;
                        XTable = handleConditionalVariables(this,XTable, learner);
                else
                    for idx = 1 : numel(learnerNames)
                        learner = learnerNames{idx};
                        XTable{XTable.learner == learner, this.IrrelevantHyperParamsMap(learner)} = missing;
                        XTable = handleConditionalVariables(this,XTable, learner);
                    end
                end
            else
                XTable = handleConditionalVariables(this,XTable, this.Learners);
            end
        end

        function XTable = handleConditionalVariables(this,XTable, learner)
        %HANDLECONDITIONALVARIABLE  invoke conditionalVariableFunction
        %based on the learner type.
        %   XTABLE = CONDITIONALVARIABLEFUNCTION(THIS,XTABLE,LEARNER) invokes
        %   conditional variable functions for individual learners based on
        %   the LEARNER type
            if strcmp(learner, 'ensemble')
                XTable = handleConditionalVariablesEnsemble(this,XTable);
            elseif strcmp(learner, 'nb')
                XTable = handleConditionalVariablesNaiveBayes(this,XTable);
            elseif strcmp(learner, 'discr')
                XTable = handleConditionalVariablesDiscr(this,XTable);
            elseif strcmp(learner, 'knn')
                XTable = handleConditionalVariablesKNN(this,XTable);
            elseif strcmp(learner, 'svm')
                XTable = handleConditionalVariablesSVM(this,XTable);
            elseif strcmp(learner, 'gp')
                XTable = handleConditionalVariablesGPR(this,XTable);
            elseif strcmp(learner, 'kernel')
                XTable = handleConditionalVariablesKernel(this,XTable);
            end 
        end

        function XTable = handleConditionalVariablesEnsemble(~,XTable)
        %HANDLECONDITIONALVARIABLESENSEMBLE  ensemble conditional variable
        %function
        %   XTABLE = HANDLECONDITIONALVARIABLESENSEMBLE(THIS,XTABLE)
        %   applies conditional constraints related to ensemble

            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'ensemble_Method', 'ensemble_LearnRate'})
                XTable.ensemble_LearnRate(XTable.ensemble_Method=='Bag') = NaN;
            end
            % Do not pass NumVariablesToSample if Method is not Bag
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'ensemble_Method', 'ensemble_NumVariablesToSample'})
                XTable.ensemble_NumVariablesToSample(XTable.ensemble_Method~='Bag') = NaN;
            end
            % Do not pass SplitCriterion if Method is LogitBoost or
            % GentleBoost, because they internally fit regression trees
            % and use 'mse'.
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'ensemble_Method', 'ensemble_SplitCriterion'})
                XTable.ensemble_SplitCriterion(XTable.ensemble_Method=='LogitBoost') = '<undefined>';
                XTable.ensemble_SplitCriterion(XTable.ensemble_Method=='GentleBoost') = '<undefined>';
            end
        end

        function XTable = handleConditionalVariablesNaiveBayes(~,XTable)
        %HANDLECONDITIONALVARIABLESNAIVEBAYES  naive Bayes conditional variable
        %function
        %   XTABLE = HANDLECONDITIONALVARIABLESNAIVEBAYES(THIS,XTABLE)
        %   applies conditional constraints related to naive Bayes

            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'nb_DistributionNames', 'nb_Kernel'})
                XTable.nb_Kernel(XTable.nb_DistributionNames ~= 'kernel') = '<undefined>';
            end
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'nb_DistributionNames', 'nb_Width'})
                XTable.nb_Width(XTable.nb_DistributionNames ~= 'kernel') = NaN;
            end
        end

        function XTable = handleConditionalVariablesDiscr(~,XTable)
        %HANDLECONDITIONALVARIABLESDISCR  discriminant conditional variable
        %function
        %   XTABLE = HANDLECONDITIONALVARIABLESDISCR(THIS,XTABLE)
        %   applies conditional constraints related to a discriminant
            % Do not pass Delta if discrim type is a quadratic

            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'discr_Delta', 'discr_DiscrimType'})
                XTable.discr_Delta(ismember(XTable.discr_DiscrimType, {'quadratic', ...
                    'diagQuadratic', 'pseudoQuadratic'})) = NaN;
            end
            % Gamma must be 0 if discrim type is a quadratic
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'discr_Gamma', 'discr_DiscrimType'})
                XTable.discr_Gamma(ismember(XTable.discr_DiscrimType, {'quadratic', ...
                    'diagQuadratic', 'pseudoQuadratic'})) = 0;
            end
        end

        function XTable = handleConditionalVariablesKNN(~,XTable)
        %HANDLECONDITIONALVARIABLESKNN  discriminant conditional variable
        %function
        %   XTABLE = HANDLECONDITIONALVARIABLESKNN(THIS,XTABLE)
        %   applies conditional constraints related to a KNN

            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'knn_Exponent', 'knn_Distance'})
                XTable.knn_Exponent(XTable.knn_Distance ~= 'minkowski') = NaN;
            end
        end

        function XTable = handleConditionalVariablesSVM(this,XTable)
        %HANDLECONDITIONALVARIABLESSVM  SVM conditional variable function
        %   XTABLE = HANDLECONDITIONALVARIABLESSVM(THIS,XTABLE)
        %   applies conditional constraints related to a SVM

            % KernelScale is irrelevant if KernelFunction~='rbf' or 'gaussian'
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'svm_KernelScale', 'svm_KernelFunction'})
                XTable.svm_KernelScale(~ismember(XTable.svm_KernelFunction, {'rbf','gaussian'})) = NaN;
            end
            % PolynomialOrder is irrelevant if KernelFunction~='polynomial'
            if bayesoptim.MultiModelOptim.hasVariables(XTable, {'svm_PolynomialOrder', 'svm_KernelFunction'})
                XTable.svm_PolynomialOrder(XTable.svm_KernelFunction ~= 'polynomial') = NaN;
            end
            
            % BoxConstraint must be 1 if NumClasses==1 in a fold
            if this.NumClasses == 1 && classreg.learning.paramoptim.BayesoptInfo.hasVariables(XTable, {'svm_BoxConstraint'})
                XTable.svm_BoxConstraint(:) = 1;
            end
        end

        function tf = XConstraintFunction(this,XTable)
        %XCONSTRAINTFUNCTION deterministic constraints function for the entire
        %optimization
        %   XTABLE = XCONSTRAINTFUNCTION(THIS,XTABLE)
        %   invokes individual deterministic constraint functions
            tf = handleXConstraintFunctionKNN(this,XTable);
        end

        function tf = handleXConstraintFunctionKNN(~,XTable)
        %HANDLEXCONSTRAINTFUNCTIONKNN  KNN deterministic function
        %optimization
        %   XTABLE = HANDLEXCONSTRAINTFUNCTIONKNN(THIS,XTABLE)
        %   applies deterministic constraints specific to KNN

            % When Standardize=true, prohibit seuclidean and mahalanobis, because the
            % result is the same when Standardize=false.
            if classreg.learning.paramoptim.BayesoptInfo.hasVariables(XTable, {'knn_Standardize', 'knn_Distance'})
                tf = ~(XTable.knn_Standardize=='true' & ismember(XTable.knn_Distance, {'seuclidean', 'mahalanobis'}));
            else
                tf = true(height(XTable),1);
            end
        end
        
        function checkNumClasses(this)
            if this.NumClasses <= 1
                error(message('stats:bayesoptim:multimodeloptimoptions:OneClassLearning'));
            end
        end
        
        function checkObservationsInRows(this)
            if ~this.DataSummary.ObservationsInRows
                error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:ObsInColsNotAllowed','fitcauto'));
            end
        end
        
        function checkForSparseData(this)
            if issparse(this.PredictorsMatrix)
                error(message('stats:bayesoptim:multimodeloptimoptions:SparseDataNotSupported'));
            end
        end
        
        function checkAndUpdateValidationMethod(this)
            if(strcmp(this.OptimizationOptions.ValidationMethod,'KFold') || ...
                    strcmp(this.OptimizationOptions.ValidationMethod,'Holdout'))
                if this.NumClasses == 0
                    % Do not stratify for regression models
                    this.OptimizationOptions.ValidationVal = ...
                        cvpartition(this.ResponseVector,this.OptimizationOptions.ValidationMethod,...
                        this.OptimizationOptions.ValidationVal,'Stratify',false);
                else
                    this.OptimizationOptions.ValidationVal = ...
                        cvpartition(this.ResponseVector,this.OptimizationOptions.ValidationMethod,...
                        this.OptimizationOptions.ValidationVal);
                end
                this.OptimizationOptions.ValidationMethod = 'CVPartition';
            end
        end
        
        function checkIfLearnersSupportMixedData(this)
        %CHECKIFLEARNERSSUPPORTMIXEDDATA  check if learners support mixed data
        %   CHECKIFLEARNERSSUPPORTMIXEDDATA(THIS) checks if
        %   the input data contains a combination of categorical and
        %   continuous predictors and skips the learners KNN and DISCR if
        %   necessary
            if any(ismember({'knn'},this.Learners))
                [~,nDims] = size(this.PredictorsMatrix);
                if ~ (isempty(this.DataSummary.CategoricalPredictors) || ...
                    (length(this.DataSummary.CategoricalPredictors) == nDims ...
                    && all(this.DataSummary.CategoricalPredictors==(1: nDims))))
                    if this.DataSummary.TableInput
                        warning(message('stats:bayesoptim:multimodeloptimoptions:KNNBadCategoricalTable')); 
                    else
                        warning(message('stats:bayesoptim:multimodeloptimoptions:KNNBadCategoricalPre')); 
                    end
                    this.Learners = setdiff(this.Learners,{'knn'});
                end   
            end

            if any(ismember({'discr'},this.Learners))
                 % Categorical predictors not allowed
                if ~isempty(this.DataSummary.CategoricalPredictors)
                    warning(message('stats:bayesoptim:multimodeloptimoptions:DiscrCategPred'));
                    this.Learners = setdiff(this.Learners,{'discr'});
                end

                if ~isfloat(this.PredictorsMatrix)
                    warning(message('stats:bayesoptim:multimodeloptimoptions:BadXType'));
                    this.Learners = setdiff(this.Learners,{'discr','linear','kernel',...
                        'nb','svm'});
                end
            end
        end
        function removeSVMOrKernelModel(this)
        % REMOVESVMORKERNEL  check to make sure that svm with gaussian
        % kernel and fitckernel with svm learner are not evaluated at the
        % same time since they solve the same optimization problem. On of
        % the two models are choosen depending on the dataset
        
            for i = 1:size(this.HyperparameterSearchRange,2)
                splits = split(this.HyperparameterSearchRange(i).Name,'_');
                if ~(numel(splits) == 2)
                    continue;
                end

                learner = splits{1};
                hyperparameter = splits{2};
                
                
                if strcmp(learner,'svm') && strcmp(hyperparameter,'KernelFunction') && ...
                        size(this.PredictorsMatrix,1)> 11000
                            this.HyperparameterSearchRange(i).Range = setdiff(this.HyperparameterSearchRange(i).Range,{'gaussian'});                            
                end
                
                if strcmp(learner,'kernel') && strcmp(hyperparameter,'Learner') && ...
                        size(this.PredictorsMatrix,1)<= 11000
                    this.HyperparameterSearchRange(i).Optimize = false;  
                    if this.OptimizationOptions.IsClassregRegressionFunction
                        addNonOptimizableNVPs(this,learner,hyperparameter,'leastsquares',0,...
                        []);
                    else
                        addNonOptimizableNVPs(this,learner,hyperparameter,'logistic',0,...
                        [])
                    end
                end       

            end
        end
        
        function removeLinearSvmModel(this)
        % REMOVESVMORKERNEL  check to make sure that svm with gaussian
        % kernel and fitckernel with svm learner are not evaluated at the
        % same time since they solve the same optimization problem. On of
        % the two models are choosen depending on the dataset
        
            for i = 1:size(this.HyperparameterSearchRange,2)
                splits = split(this.HyperparameterSearchRange(i).Name,'_');
                if ~(numel(splits) == 2)
                    continue;
                end

                learner = splits{1};
                hyperparameter = splits{2};
                
                
                if strcmp(learner,'svm') && strcmp(hyperparameter,'KernelFunction') 
                   RangeLength=numel(setdiff(this.HyperparameterSearchRange(i).Range,{'linear'}));
                   if RangeLength~=1
                    this.HyperparameterSearchRange(i).Range = setdiff(this.HyperparameterSearchRange(i).Range,{'linear'}); 
                   else
                    this.HyperparameterSearchRange(i).Optimize = false;
                    addNonOptimizableNVPs(this,learner,hyperparameter,'polynomial',0,...
                        []);
                   end
                end       
            end
        end
        
        
        function addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
               NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
               weakLearnerNoNOptimizableHyperparameters)

            if ~isfield(this.AutoSearchRange,learnerName)
                learnerSearchRange = bayesoptim.MultiModelOptim.defaultLearnerSearchRange();
            else
                learnerSearchRange = this.AutoSearchRange.(learnerName);
            end
                learnerSearchRange.OptimizableHyperparameters = bayesoptim.MultiModelOptim.nvpToStruct(learnerSearchRange.OptimizableHyperparameters,OptimizableHyperparameters);
                learnerSearchRange.NoNOptimizableHyperparameters = bayesoptim.MultiModelOptim.nvpToStruct(learnerSearchRange.NoNOptimizableHyperparameters,NoNOptimizableHyperparameters);

                if ~isempty(weakLearnerName)
                    learnerSearchRange.weakLearnerName = weakLearnerName;
                    learnerSearchRange.weakLearnerOptimizableHyperparameters = ...
                        bayesoptim.MultiModelOptim.nvpToStruct(learnerSearchRange.weakLearnerOptimizableHyperparameters,weakLearnerOptimizableHyperparameters);
                    learnerSearchRange.weakLearnerNoNOptimizableHyperparameters = ...
                        bayesoptim.MultiModelOptim.nvpToStruct(learnerSearchRange.weakLearnerNoNOptimizableHyperparameters,weakLearnerNoNOptimizableHyperparameters);
                end
                this.AutoSearchRange.(learnerName) = learnerSearchRange;
        end
        
        function updateSearchRangeForAll(this)
            % below function needs to be called when learners is set to all
            removeSVMOrKernelModel(this);
            removeLinearSvmModel(this);
        end
        
        function updateSearchRangeForAuto(this)
           for i = 1:size(this.HyperparameterSearchRange,2)
                splits = split(this.HyperparameterSearchRange(i).Name,'_');
                if ~(numel(splits) == 2)
                    continue;
                end

                learner = splits{1};
                hyperparameter = splits{2};

                if ~isfield(this.AutoSearchRange,learner)
                    continue;
                end
                
                if strcmp(learner,'ensemble') && strcmp(hyperparameter,'LearnRate') && ...
                        isfield(this.AutoSearchRange,'ensemble') && ...
                        isfield(this.AutoSearchRange.ensemble.OptimizableHyperparameters,'Method') && ...
                        (numel(this.AutoSearchRange.ensemble.OptimizableHyperparameters.Method) == 1) && ...
                        strcmp(this.AutoSearchRange.ensemble.OptimizableHyperparameters.Method,'Bag')
                            this.HyperparameterSearchRange(i).Optimize = false;
                end

                if strcmp(learner,'svm') && strcmp(hyperparameter,'PolynomialOrder')
                    if isfield(this.AutoSearchRange,'svm') && isfield(this.AutoSearchRange.svm.OptimizableHyperparameters,'KernelFunction')
                        if ~(numel(this.AutoSearchRange.svm.OptimizableHyperparameters.KernelFunction)> 1)
                            this.HyperparameterSearchRange(i).Optimize = false;
                        end
                    end
                end

                if strcmp(learner,'knn') && strcmp(hyperparameter,'Exponent')
                    if isfield(this.AutoSearchRange,'knn') && isfield(this.AutoSearchRange.knn.OptimizableHyperparameters,'Distance')
                        if ~any(strcmp(this.AutoSearchRange.knn.OptimizableHyperparameters.Distance,'minkowski'))
                            this.HyperparameterSearchRange(i).Optimize = false;
                        end
                    end
                end

                if ~strcmp(learner,'ensemble') && ~isfield(this.AutoSearchRange.(learner).OptimizableHyperparameters,hyperparameter)
                    continue;
                elseif strcmp(learner,'ensemble') && ~isfield(this.AutoSearchRange.(learner).OptimizableHyperparameters,hyperparameter) && ~isfield(this.AutoSearchRange.(learner).weakLearnerOptimizableHyperparameters,hyperparameter) 
                    continue;
                end

                isWeakLearner = false;

                if strcmp(learner,'ensemble') && isfield(this.AutoSearchRange.(learner).weakLearnerOptimizableHyperparameters,hyperparameter)
                    hyperparameterRange = this.AutoSearchRange.(learner).weakLearnerOptimizableHyperparameters.(hyperparameter);
                    isWeakLearner = true;
                else
                    hyperparameterRange = this.AutoSearchRange.(learner).OptimizableHyperparameters.(hyperparameter);
                end

                if numel(hyperparameterRange) > 1
                    this.HyperparameterSearchRange(i).Range = hyperparameterRange;
                    this.HyperparameterSearchRange(i).Optimize = true;
                else
                    this.HyperparameterSearchRange(i).Optimize = false;
                    addNonOptimizableNVPs(this,learner,hyperparameter,hyperparameterRange,isWeakLearner,this.AutoSearchRange.(learner).weakLearnerName);
                end
           end

           learners = fieldnames(this.AutoSearchRange);

           for i = 1:numel(learners)
                nonOptimizableNVPs = this.AutoSearchRange.(learners{i}).NoNOptimizableHyperparameters;
                hyperparameters = fieldnames(nonOptimizableNVPs);
                 for j = 1:numel(hyperparameters)
                     addNonOptimizableNVPs(this,learners{i},hyperparameters{j},nonOptimizableNVPs.(hyperparameters{j}),false,[]);
                 end

                 if ~isempty(this.AutoSearchRange.(learners{i}).weakLearnerName)
                     weakLearnerNonOptimizableNVPs = this.AutoSearchRange.(learners{i}).weakLearnerNoNOptimizableHyperparameters;
                     weakLearnerHyperparameters = fieldnames(weakLearnerNonOptimizableNVPs);
                     for j = 1:numel(weakLearnerHyperparameters)
                        addNonOptimizableNVPs(this,learners{i},weakLearnerHyperparameters{j},weakLearnerNonOptimizableNVPs.(weakLearnerHyperparameters{j}),true,this.AutoSearchRange.(learners{i}).weakLearnerName);
                     end
                 end
           end
        end

        function updateSearchRangeForAllLinear(this)
            for i = 1:size(this.HyperparameterSearchRange,2)
            splits = split(this.HyperparameterSearchRange(i).Name,'_');
            
            if ~(numel(splits) == 2)
                continue;
            end

            learner = splits{1};
            hyperparameter = splits{2};


           if strcmp(learner,'svm') && strcmp(this.OptimizeHyperparameters,'all')
                if strcmp(hyperparameter,'KernelFunction')
                    this.HyperparameterSearchRange(i).Optimize = false;
                end
                
                if strcmp(hyperparameter,'PolynomialOrder')
                    this.HyperparameterSearchRange(i).Optimize = false;
                end
           end
           
           if strcmp(learner,'discr')
                if strcmp(hyperparameter,'DiscrimType')
                    this.HyperparameterSearchRange(i).Range = {'linear','diagLinear','pseudoLinear'};
                    this.HyperparameterSearchRange(i).Optimize = true;
                end
           end
                     
           end
        end
        
        function updateSearchRangeForAllNoNLinear(this)
           % below function needs to be called when learners is set to
           % all-nonlinear
           removeSVMOrKernelModel(this);
           for i = 1:size(this.HyperparameterSearchRange,2)
            splits = split(this.HyperparameterSearchRange(i).Name,'_');
            if ~(numel(splits) == 2)
                continue;
            end

            learner = splits{1};
            hyperparameter = splits{2};


           if strcmp(learner,'svm')
                if strcmp(hyperparameter,'KernelFunction')
                    this.HyperparameterSearchRange(i).Optimize = true;
                    this.HyperparameterSearchRange(i).Range = {'gaussian','polynomial'};
                end
           end
           
           if strcmp(learner,'discr')
                if strcmp(hyperparameter,'DiscrimType')
                    this.HyperparameterSearchRange(i).Range = {'quadratic','diagQuadratic','pseudoQuadratic'};
                    this.HyperparameterSearchRange(i).Optimize = true;
                end
           end
                     
           end
        end
        
        
        function addNonOptimizableNVPs(this,learner,hyperparameter,hyperparameterValue,isWeakLearner,weakLearnerName)
           if iscell(hyperparameterValue)
               hyperparameterValue = hyperparameterValue{:};
           end
           if ~isWeakLearner
               if ~isfield(this.NonOptimizableNVPs,learner)
                    this.NonOptimizableNVPs.(learner) = {hyperparameter,hyperparameterValue};
               else
                this.NonOptimizableNVPs.(learner) = {this.NonOptimizableNVPs.(learner){:},...
                    hyperparameter,hyperparameterValue};
               end
           else
              if ~isfield(this.WeakLearnerNonOptimizableNVPs,weakLearnerName)
                    this.WeakLearnerNonOptimizableNVPs.(weakLearnerName) = {hyperparameter,hyperparameterValue};
              else
                this.WeakLearnerNonOptimizableNVPs.(weakLearnerName) = {this.WeakLearnerNonOptimizableNVPs.(weakLearnerName){:},...
                    hyperparameter,hyperparameterValue};
              end
           end
        end
        
        function addBaggingModels(this,isBig,surrogate)
            if isBig
                %Light RandomForest
                numLearningCycles = 50;
                predictorSelection = 'allsplits';
                addBaggingEnsembleModel(this,numLearningCycles,surrogate,predictorSelection)
            else
                %Magic RandomForest
                numLearningCycles = 200;
                predictorSelection = 'interaction-curvature';
                addBaggingEnsembleModel(this,numLearningCycles,surrogate,predictorSelection)
            end
        end
        
        function LearnerArgPairs =  AddSolverForLinear(this,optimizableHyperparams)
            if strcmp(this.OptimizeHyperparameters,'all')
                if strcmp(optimizableHyperparams{:,'Regularization'},'lasso')
                    LearnerArgPairs{1} = 'Solver';
                    LearnerArgPairs{2} = 'sparsa';
                else
                    LearnerArgPairs{1} = 'Solver';
                    LearnerArgPairs{2} = 'lbfgs';
                end
            elseif strcmp(this.OptimizeHyperparameters,'auto')  
                LearnerArgPairs{1} = 'Solver';
                LearnerArgPairs{2} = 'lbfgs'; % lbfgs is choosen since ridge is used by default
            else
                LearnerArgPairs = {};
            end
        end
    end
    
    methods(Static, Access = protected)
        % Utility functions
        function tf = hasVariables(Tbl, VarNames)
            % Return true if table Tbl has all variables VarNames.
            tf = all(ismember(VarNames, Tbl.Properties.VariableNames));
        end
        
        function argPairs = updateHyperparameterNVPs(argPairs,hyperparameterName,...
                relevantHyperparameter,optimizedHyperparameters)
            
            %extracts the hyperparameter name value pairs to be passed to the fit
            %functions
            value = optimizedHyperparameters.(relevantHyperparameter);
            
            if ismissing(value)
                return
            end
            
            if iscategorical(value)
                value = char(value);
            end
            
            if isempty(argPairs)
                argPairs{1} = hyperparameterName;
            else
                argPairs{end+1} = hyperparameterName;
            end
            
            argPairs{end+1} = value;
        end
        
        function LearnerArgPairs =  updateStandardize(LearnerArgPairs)
            standardizeIndex = [];
            for i = 1:size(LearnerArgPairs,2)
                if strcmp(LearnerArgPairs{i},'Standardize')
                    standardizeIndex = i+1;
                    break;
                end
            end
            
            if isempty(standardizeIndex)
                return;
            end
            
            if strcmp(LearnerArgPairs{standardizeIndex},'true')
                LearnerArgPairs{standardizeIndex} = true;
            else
                LearnerArgPairs{standardizeIndex} = false;
            end
        end
        
        function LearnerArgPairs =  removeLearnrateForBaggedEnsemble(parameterCombination,LearnerArgPairs)
            
            if ~any(ismember(parameterCombination.Properties.VariableNames,'ensemble_Method'))
                return;
            end
            
            if ~strcmp(char(parameterCombination.ensemble_Method),'Bag')
                return
            end
            
            LearnRateIndex = [];
            for i = 1:size(LearnerArgPairs,2)
                if strcmp(LearnerArgPairs{i},'LearnRate')
                    LearnRateIndex = i;
                    break;
                end
            end
            
            if isempty(LearnRateIndex)
                return;
            else
                LearnerArgPairs(LearnRateIndex:LearnRateIndex+1) = [];
            end
        end
        
        function [LearnerArgPairs,weakLearnerNonOptimizableArgPairs] =  AddNumBinsForBoostedEnsemble(hyperparamsTable,LearnerArgPairs,weakLearnerNonOptimizableArgPairs)
            
            if any(ismember(hyperparamsTable.Properties.VariableNames,'ensemble_Method'))
                method = hyperparamsTable.ensemble_Method;
            else
                MethodIndex = [];
                for i = 1:size(LearnerArgPairs,2)
                    if strcmp(LearnerArgPairs{i},'Method')
                        MethodIndex = i;
                        break;
                    end
                end
                
                if isempty(MethodIndex)
                    return;
                end
                method = LearnerArgPairs{MethodIndex+1};
            end
            if ismember(method,{'GentleBoost', 'LogitBoost', 'AdaBoostM1','AdaBoostM2','RUSBoost'})
                LearnerArgPairs{end+1} = 'NumBins';
                LearnerArgPairs{end+1} = 50;
                
                weakLearnerNonOptimizableArgPairs = bayesoptim.MultiModelOptim.removeInteractionCurvature(weakLearnerNonOptimizableArgPairs);
            end
        end     

        
        function weakLearnerArgPairs = removeInteractionCurvature(weakLearnerArgPairs)
            predictorSelectionIndex = [];
            for i = 1:size(weakLearnerArgPairs,2)
                if strcmp(weakLearnerArgPairs{i},'PredictorSelection')
                    predictorSelectionIndex = i;
                    break;
                end
            end
            
            if isempty(predictorSelectionIndex)
                return;
            end
            
            predictorSelection = weakLearnerArgPairs{predictorSelectionIndex+1};
            
            if strcmp(predictorSelection,'interaction-curvature')
                weakLearnerArgPairs(predictorSelectionIndex:predictorSelectionIndex+1) = [];
            end
        end
        
        function nvpStruct = nvpToStruct(nvpStruct,nvp)
            if isempty(fieldnames(nvpStruct))
                for i = 1:2:size(nvp,2)
                    if isnumeric(nvp{i+1}) || islogical(nvp{i+1})
                        nvpStruct.(nvp{i}) =  nvp{i+1};
                    else
                        nvpStruct.(nvp{i}) =  {nvp{i+1}};
                    end
                end
            else
                for i = 1:2:size(nvp,2)
                    if ~isfield(nvpStruct,nvp{i}) || isempty(nvpStruct.(nvp{i}))
                        if isnumeric(nvp{i+1}) || islogical(nvp{i+1})
                            nvpStruct.(nvp{i}) =  nvp{i+1};
                        else
                            nvpStruct.(nvp{i}) =  {nvp{i+1}};
                        end
                    else
                        if islogical(nvpStruct.(nvp{i}))
                            nvpStruct.(nvp{i}) =  unique([nvpStruct.(nvp{i}),nvp{i+1}]);
                        elseif isnumeric(nvpStruct.(nvp{i}))
                            range = unique([nvpStruct.(nvp{i}),nvp{i+1}]);
                            if numel(range) > 1
                                nvpStruct.(nvp{i}) = [min(range),max(range)];
                            else
                                nvpStruct.(nvp{i}) = range;
                            end
                        else
                            nvpStruct.(nvp{i}) =   unique({nvpStruct.(nvp{i}){:},nvp{i+1}});
                        end
                    end
                end
            end
        end
        
        function learnerSearchRange = defaultLearnerSearchRange()
            learnerSearchRange = struct();
            learnerSearchRange.OptimizableHyperparameters = struct();
            learnerSearchRange.NoNOptimizableHyperparameters = struct();
            learnerSearchRange.weakLearnerName = [];
            learnerSearchRange.weakLearnerOptimizableHyperparameters = struct();
            learnerSearchRange.weakLearnerNoNOptimizableHyperparameters = struct();
        end
        
        function template =  createTemplate(classifier,nonOptimizableLearnerArgPairs)
            switch classifier
                case 'svm'
                    template = templateSVM(nonOptimizableLearnerArgPairs{:});
                case 'linear'
                    template = templateLinear(nonOptimizableLearnerArgPairs{:});
                case 'kernel'
                    template = templateKernel(nonOptimizableLearnerArgPairs{:});
            end
        end
    end
end
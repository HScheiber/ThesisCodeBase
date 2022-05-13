classdef ClassificationMultiModelOptim < bayesoptim.MultiModelOptim
    %ClassificationMultiModelOptim This class performs bayesian
    %optimization across classification models and their hyperparameters

%   Copyright 2019-2020 The MathWorks, Inc.

    properties(GetAccess=public,Constant)
        % Cell array that stores all possible classification learner names
        AllLearners = {'svm','knn','nb','tree','discr','linear','ensemble',...
            'kernel'};
        
        % Cell array that stores all possible fitcfunction names
        FitFunctions = {'fitcsvm','fitcknn','fitcnb','fitctree','fitcdiscr',...
            'fitclinear','fitcensemble','fitckernel'};
    end
    
    methods
        % Constructor
        function this = ClassificationMultiModelOptim(X,Y,multiModelOptimOptions)          
            this.X = X;
            this.Y = Y;
            this.Learners = multiModelOptimOptions.Learners;
            this.OptimizeHyperparameters = multiModelOptimOptions.OptimizeHyperparameters;
            this.OptimizationOptions = multiModelOptimOptions.OptimizationOptions;
            this.FitFunctionArgs = multiModelOptimOptions.RemainingArgs;
            this.UserPassedMaxObjectiveEvaluations = multiModelOptimOptions.UserPassedMaxObjectiveEvaluations;
            this.FitFunctionMap = containers.Map(this.AllLearners,this.FitFunctions);
            this.IrrelevantHyperParamsMap = [];
            
            [this.PredictorsMatrix,this.ResponseVector,~,dataSummary,classSummary,~] = ...
                classreg.learning.classif.FullClassificationModel.prepareData(X, Y, this.FitFunctionArgs{:});
            
            this.NumClasses = numel(classSummary.ClassNames);
            this.DataSummary = dataSummary;
            updateLearners(this);
            checkNumClasses(this);
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
                classifier = char(parameterCombination.learner);
            else
                classifier = this.Learners{:};
            end

            % get the hyperparameters for the current algorithm
            optimizableHyperparams = getArgumentsTable(this,this.HyperparameterNames.(classifier), parameterCombination);
            nonOptimizableLearnerArgPairs = {};
            

            if isfield(this.NonOptimizableNVPs,classifier)
                nonOptimizableLearnerArgPairs = this.NonOptimizableNVPs.(classifier);
                if any(strcmp(classifier,{'svm','knn'}))
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                elseif strcmp(classifier,'ensemble')
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.removeLearnrateForBaggedEnsemble(parameterCombination,nonOptimizableLearnerArgPairs);
                end
            end
            
            if this.NumClasses > 2 && ismember(classifier, {'kernel', 'svm', 'linear'})
                codingIdx = find(strcmp(nonOptimizableLearnerArgPairs,'Coding'));
                codingNVP = {};
                if ~isempty(codingIdx)
                    codingNVP = {nonOptimizableLearnerArgPairs{codingIdx},nonOptimizableLearnerArgPairs{codingIdx+1}};
                    nonOptimizableLearnerArgPairs(codingIdx:codingIdx+1) = [];
                end
                template = bayesoptim.MultiModelOptim.createTemplate(classifier,nonOptimizableLearnerArgPairs);
                fitFunArgs = [{'Learners'}, {template},this.FitFunctionArgs(:)',codingNVP(:)'];
            elseif strcmp(classifier,'ensemble')
                weakLearnerNonOptimizableNVPs = {};
                if isfield(this.WeakLearnerNonOptimizableNVPs,'tree')
                    weakLearnerNonOptimizableNVPs = this.WeakLearnerNonOptimizableNVPs.('tree');
                    [nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs] = ...
                        bayesoptim.MultiModelOptim.AddNumBinsForBoostedEnsemble(parameterCombination,nonOptimizableLearnerArgPairs,weakLearnerNonOptimizableNVPs);
                end
                fitFunArgs = [{'Learners'},{templateTree(weakLearnerNonOptimizableNVPs{:})},this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            elseif strcmp(classifier,'linear')
                nonOptimizableLearnerArgPairs=AddSolverForLinear(this,optimizableHyperparams);                    
                fitFunArgs = [this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            else
                fitFunArgs = [this.FitFunctionArgs(:)',nonOptimizableLearnerArgPairs(:)'];
            end

            % get objective function of current algorithm
            objFcn = getObjectiveFcn(this,this.BayesOptInfo.(classifier), fitFunArgs);
            loss = objFcn(optimizableHyperparams);
        end
    end
    
    methods(Access = {?tClassificationMultiModelOptim, ?bayesoptim.MultiModelOptim})
        
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

                % Extract weakLearner specific arguments if the classifier is
                % ensemble
                if ismember(relevantHyperparameters{idx}, {'ensemble_MinLeafSize', ...
                        'ensemble_MaxNumSplits', 'ensemble_SplitCriterion', 'ensemble_NumVariablesToSample'})
                    weakLearnerArgPairs = bayesoptim.MultiModelOptim.updateHyperparameterNVPs(...
                        weakLearnerArgPairs,hyperparameterName,...
                        relevantHyperparameters{idx},optimizedHyperparameters);
                end

                % Extract arguments specific to the classifier
                if ~strcmp(hyperparameterName, 'Coding') ...
                    && ~ismember(relevantHyperparameters{idx}, {'ensemble_MinLeafSize', ...
                        'ensemble_MaxNumSplits', 'ensemble_SplitCriterion', 'ensemble_NumVariablesToSample'})
                    learnerArgPairs = bayesoptim.MultiModelOptim.updateHyperparameterNVPs(learnerArgPairs,...
                        hyperparameterName,relevantHyperparameters{idx},...
                        optimizedHyperparameters);
                end
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
           
            if this.NumClasses > 2 && ismember(learner, {'kernel', 'svm', 'linear'})
                model = trainEcocModels(this,learner,learnerArgPairs,nonOptimizableLearnerArgPairs,optimizedHyperparameters);
            else
                model = trainNonEcocModels(this,learner,learnerArgPairs,nonOptimizableLearnerArgPairs,...
                    weakLearnerArgPairs,weakLearnerNonOptimizableNVPs);
            end
        end
        
        function XTable = handleConditionalVariablesKernel(~,XTable)
            % Do nothing here. ClassificationKernel does not requite
            % conditional constraints but RegressionKernel does
        end
        
        function updateLearners(this)
        %UPDATELEARNERS  update the learner types for optimization
        %   UPDATELEARNERS(THIS) checks for the learner options
        %   'auto','all','all-linear' and 'all-nonlinear' and updates the
        %   learner types accordingly. If user passes specific learner
        %   types they are used directly without any change.
        
            learners = this.Learners{:};
            switch learners
                case 'all-linear'
                    this.LearnersIsAuto = false;
                    this.LearnersIsLinear = true;
                    this.Learners = {'linear','discr'}; % removed svm since linear got it covered
                case 'all-nonlinear'
                    this.LearnersIsAuto = false;
                    this.LearnersIsNoNLinear = true;
                    this.Learners = {'svm','knn','nb','tree','discr','ensemble','kernel'};
                case 'all'
                    this.LearnersIsAuto = false;
                    this.LearnersIsAll = true;   
                    this.Learners = this.AllLearners;
                case 'auto'
                     dataProp = dataProperties(this);
                     updateLearnersForAuto(this,dataProp);
                otherwise
                    this.LearnersIsAuto = false;
            end
            
            checkIfLearnersSupportMixedData(this);             
        end
        
        function model = trainNonEcocModels(this,learner,LearnerArgPairs,nonOptimizableLearnerArgPairs,...
                weakLearnerArgPairs,weakLearnerNonOptimizableNVPs)
        %TRAINNONECOCMODELS  trains a non ecoc model
        %   TRAINNONECOCMODELS(THIS,LEARNER,LEARNERARGPAIRS,WEAKLEARNERARGPAIRS) 
        %   trains and returns a nonecoc model. The LEARNERARGPAIRS and
        %   THIS.FITFUNCTIONARGS are passed to the fit function. Only
        %   ensemble models use the additional WEAKLEARNERARGPAIRS

            switch learner
                case 'linear'
                    model = fitclinear(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'svm'
                    LearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(LearnerArgPairs);
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                    model = fitcsvm(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'knn'
                    LearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(LearnerArgPairs);
                    nonOptimizableLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(nonOptimizableLearnerArgPairs);
                    model = fitcknn(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'kernel'
                    model= fitckernel(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'discr'
                    model = fitcdiscr(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'nb'
                    model = trainNaiveBayes(this,LearnerArgPairs,nonOptimizableLearnerArgPairs);
                case 'tree'
                    model = fitctree(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                case 'ensemble'
                    model = fitcensemble(this.X, this.Y, LearnerArgPairs{:}, nonOptimizableLearnerArgPairs{:},'Learners', ...
                        templateTree(weakLearnerArgPairs{:},weakLearnerNonOptimizableNVPs{:}),this.FitFunctionArgs{:});
            end
        end
        
        function model = trainEcocModels(this,learner,weakLearnerArgPairs,nonOptimizableLearnerArgPairs,optimizedHyperparameters)
        %TRAINECOCMODELS  trains a ecoc model
        %   TRAINECOCMODELS(THIS,LEARNER,WEAKLEARNERARGPAIRS,OPTIMIZEDHYPERPARAMETERS) trains
        %   and returns a ecoc model with a linear, SVM or Kernel templates. The WEAKLEARNERARGPAIRS and
        %   THIS.FITFUNCTIONARGS are passed to the template function

            codingIdx = find(strcmp(nonOptimizableLearnerArgPairs,'Coding'));
            codingNVP = {};
            if ~isempty(codingIdx)
                codingNVP = {nonOptimizableLearnerArgPairs{codingIdx},nonOptimizableLearnerArgPairs{codingIdx+1}};
                nonOptimizableLearnerArgPairs([codingIdx:codingIdx+1]) = [];
            else
                switch learner
                    case 'linear'
                        codingNVP = {'Coding', char(optimizedHyperparameters.linear_Coding)};
                    case 'svm'
                        codingNVP = {'Coding',char(optimizedHyperparameters.svm_Coding)};
                    case 'kernel'
                        codingNVP = {'Coding',char(optimizedHyperparameters.kernel_Coding)};
                end
            end

            switch learner
                case 'linear'
                    model = fitcecoc(this.X, this.Y, codingNVP{:}, 'Learner',...
                        templateLinear(weakLearnerArgPairs{:},nonOptimizableLearnerArgPairs{:}),this.FitFunctionArgs{:});
                case 'svm'
                    weakLearnerArgPairs = bayesoptim.MultiModelOptim.updateStandardize(weakLearnerArgPairs);
                    model = fitcecoc(this.X, this.Y, codingNVP{:}, 'Learners', ...
                        templateSVM(weakLearnerArgPairs{:},nonOptimizableLearnerArgPairs{:}),this.FitFunctionArgs{:});
                case 'kernel'
                    model= fitcecoc(this.X, this.Y, codingNVP{:}, 'Learners', ...
                        templateKernel(weakLearnerArgPairs{:},nonOptimizableLearnerArgPairs{:}),this.FitFunctionArgs{:});
            end
        end

        function model = trainNaiveBayes(this,LearnerArgPairs,nonOptimizableLearnerArgPairs)
        %TRAINNAIVEBAYES  trains a naive Bayes model
        %   TRAINNAIVEBAYES(THIS,LEARNERARGPAIRS) trains
        %   and returns a naive Bayes model with categorical predictors
        %   using 'mvmn' distribution

            if isempty(this.DataSummary.CategoricalPredictors)
                model = fitcnb(this.X, this.Y, LearnerArgPairs{:},nonOptimizableLearnerArgPairs{:},this.FitFunctionArgs{:});
                return;
            end

            distNameIndex = [];         
            for i = 1:size(LearnerArgPairs,2)
                if strcmp(LearnerArgPairs{i},'DistributionNames')
                    distNameIndex = i+1;
                    break;
                end
            end

            % X could be a table that includes response. So, use
            % PredictorsMatrix
            distNames = repmat(LearnerArgPairs(distNameIndex),[1,size(this.PredictorsMatrix,2)]);
            distNames(this.DataSummary.CategoricalPredictors) = {'mvmn'};
            LearnerArgPairs{distNameIndex} = distNames;
            model = fitcnb(this.X, this.Y, LearnerArgPairs{:},this.FitFunctionArgs{:});
        end

        function dataProp =  dataProperties(this)
            dataProp = struct();
            [N,D] = size(this.PredictorsMatrix);
            dataProp.hasCategorical = ~isempty(this.DataSummary.CategoricalPredictors);
            isValueMissing = ismissing(this.PredictorsMatrix);
            fracObsWithMissingVals = sum(any(isValueMissing,2))/N;
            dataProp.hasManyMissing = fracObsWithMissingVals > 0.05;
            dataProp.isLabelOrdinal = isa(labels(this.ResponseVector),'categorical') && isordinal(labels(this.ResponseVector));
            A = tabulate(labels(this.ResponseVector));
            if iscell(A)
                classSize = cell2mat(A(:,2));
            else
                classSize = A(:,2);
            end
            dataProp.smallestClassSize = min(classSize);
            largestClassSize = max(classSize);
            classRatio = largestClassSize/dataProp.smallestClassSize;
            dataProp.isImbalanced = classRatio > 5;
            dataProp.isBinary = this.NumClasses==2;
            dataProp.isWide = D >= N;
            dataProp.isHighD = D > 100;
            dataProp.isBig = N > 50000;
        end
        
        function updateLearnersForAuto(this,dataProp)
            hasManyMissing =  dataProp.hasManyMissing;
            isImbalanced = dataProp.isImbalanced;
            smallestClassSize = dataProp.smallestClassSize;
            isHighD = dataProp.isHighD;
            isBinary = dataProp.isBinary;
            hasCategorical = dataProp.hasCategorical;
            isWide = dataProp.isWide;
            isLabelOrdinal = dataProp.isLabelOrdinal;
            isBig = dataProp.isBig;
            
            if hasManyMissing
               surrogate = 'on';
               addTreeModel(this,surrogate);           
               addNaiveBayesModel(this);
               addBoostingModels(this,isImbalanced,smallestClassSize,isHighD,isBinary,surrogate);
               addBaggingModels(this,isBig,surrogate);
            elseif hasCategorical
               surrogate = 'off';
               addTreeModel(this,surrogate);
               % Naive bayes models have no hyperparameters to optimize when
               % all the predictors are categorical.
               if size(this.DataSummary.CategoricalPredictors,2) ~= size(this.PredictorsMatrix,2)
                addNaiveBayesModel(this);
               end
               addBoostingModels(this,isImbalanced,smallestClassSize,isHighD,isBinary,surrogate);
               addBaggingModels(this,isBig,surrogate);
               addSVMModels(this,'polynomial',isBinary,isLabelOrdinal,isBig)
            elseif isWide
                addDiscriminantModel(this,'linear');                
                addLinearModels(this,isBinary,isLabelOrdinal)
                if ~isBig
                    addKNNModel(this,'cosine',false)
                end
            else
                surrogate = 'off';
                addNaiveBayesModel(this);
                addTreeModel(this,surrogate);
                addBoostingModels(this,isImbalanced,smallestClassSize,isHighD,isBinary,surrogate);
                addBaggingModels(this,isBig,surrogate);
                addSVMModels(this,'rbf',isBinary,isLabelOrdinal,isBig)
               if ~isBig
                    if isHighD
                        addKNNModel(this,'cosine',false)
                   else
                        addKNNModel(this,'euclidean',true)
                   end
                end
            end
            
        end
                       
        function addBoostingModels(this,isImbalanced,smallestClassSize,isHighD,isBinary,surrogate)
            if  isImbalanced && smallestClassSize>=100                
                %Rough RUSBoost
                numLearningCycles = 500;
                maxNumSplits = 10;
                addBoostingEnsembleModel(this,'RUSBoost',numLearningCycles,maxNumSplits,surrogate);
                
                %Fine RUSBoost
                maxNumSplits = 5*this.NumClasses*(this.NumClasses-1);
                addBoostingEnsembleModel(this,'RUSBoost',numLearningCycles,maxNumSplits,surrogate);
           elseif ~isHighD
                if isBinary
                    %LogitBoost
                    numLearningCycles = 300;
                    maxNumSplits = 10;
                    addBoostingEnsembleModel(this,'LogitBoost',numLearningCycles,maxNumSplits,surrogate);
                else
                    %Rough AdaBoostM2
                    numLearningCycles = 300;
                    maxNumSplits = 10;
                    addBoostingEnsembleModel(this,'AdaBoostM2',numLearningCycles,maxNumSplits,surrogate);

                    %Fine AdaBoostM2
                    maxNumSplits = 5*this.NumClasses*(this.NumClasses-1);
                    addBoostingEnsembleModel(this,'AdaBoostM2',numLearningCycles,maxNumSplits,surrogate);
                end
            end
        end
        
        function addSVMModels(this,kernelFunction,isBinary,isLabelOrdinal,isBig)
            if isBinary
                addSVMModel(this,kernelFunction,[])
            else
                if strcmp(kernelFunction,'linear') || strcmp(kernelFunction,'polynomial') || ...
                        (strcmp(kernelFunction,'rbf') && ~isBig)
                    if isLabelOrdinal
                        addSVMModel(this,kernelFunction,'ordinal')
                    else
                        addSVMModel(this,kernelFunction,'onevsone')
                        addSVMModel(this,kernelFunction,'onevsall')
                    end
                end
            end
        end
        
        function addLinearModels(this,isBinary,isLabelOrdinal)
            if isBinary
                addLinearModel(this,[],'svm')
                addLinearModel(this,[],'logistic')                
            else
                if isLabelOrdinal
                    addLinearModel(this,'ordinal','svm')
                    addLinearModel(this,'ordinal','logistic')
                    
                else
                    addLinearModel(this,'onevsall','svm')
                    addLinearModel(this,'onevsone','logistic')                    
                end
            end
        end
        
        function addBoostingEnsembleModel(this,method,numLearningCycles,maxNumSplits,surrogate)
            learnerName = 'ensemble';
            NoNOptimizableHyperparameters = [];
            weakLearnerName = 'tree';
            
            if strcmp(this.OptimizeHyperparameters,'auto')
                OptimizableHyperparameters = {'Method',method,...
                    'NumLearningCycles',numLearningCycles...
                    'LearnRate',0.1};
                weakLearnerOptimizableHyperparameters = {'MaxNumSplits',maxNumSplits};
            else
                OptimizableHyperparameters = {'Method',method};
                if strcmp(method,'LogitBoost')
                    weakLearnerOptimizableHyperparameters = [];
                else
                    weakLearnerOptimizableHyperparameters = {'MaxNumSplits',maxNumSplits};
                end
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
        
        function addNaiveBayesModel(this)
            learnerName = 'nb';
            OptimizableHyperparameters = [];
            NoNOptimizableHyperparameters = [];
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
        end
        
        function addSVMModel(this,kernelFunction,coding)
            learnerName = 'svm';
            if isempty(coding)
                OptimizableHyperparameters = {'KernelFunction',kernelFunction,'Standardize',true};
            else
                OptimizableHyperparameters = {'KernelFunction',kernelFunction,'Standardize',true,'Coding',coding};
            end
            
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
        
        function addDiscriminantModel(this,discrimType)
            learnerName = 'discr';
            OptimizableHyperparameters = {'DiscrimType',discrimType};
            NoNOptimizableHyperparameters = {'FillCoeffs','off','SaveMemory','on'};
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters)
        end
        
        function addLinearModel(this,coding,learner)
            learnerName = 'linear';
            if isempty(coding)
                OptimizableHyperparameters = {'Learner',learner};
            else
                OptimizableHyperparameters = {'Learner',learner,'Coding',coding};
            end            
            % removing the solver set to lbfgs 
            NoNOptimizableHyperparameters = {};
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters)
        end
        
        function addKNNModel(this,distance,standardize)
            learnerName = 'knn';
            OptimizableHyperparameters = {'Distance',distance,'Standardize',standardize};
            NoNOptimizableHyperparameters = [];
            weakLearnerName = [];
            weakLearnerOptimizableHyperparameters = [];
            weakLearnerNoNOptimizableHyperparameters = [];
            addLearnerSearchRange(this,learnerName,OptimizableHyperparameters,...
                NoNOptimizableHyperparameters,weakLearnerName,weakLearnerOptimizableHyperparameters,...
                weakLearnerNoNOptimizableHyperparameters);
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
    end
end
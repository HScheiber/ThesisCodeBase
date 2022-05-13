classdef MultiModelOptimOptions < handle
    % MultiModelOptimOptions This class specifies the default options
    % for optimization across models and their hyperparameters

%   Copyright 2019-2020 The MathWorks, Inc.

    properties
        % Either 'auto' or 'all', represents the set of hyperparameters to
        % be optimized
        OptimizeHyperparameters = 'auto';   % Either 'auto' or 'all'
        
        % An object or table describing the results of hyperparameter optimization. 
        OptimizationOptions = [];
        
        % Either 'auto','all','all-linear' or 'all-nonlinear' or a cell array
        % of eligible learner names
        Learners = 'auto';
        
        % A cell array of additional arguments passed by the user
        RemainingArgs = {};
        
        % A boolean to indicate if user passed MaxObjectiveEvaluations
        UserPassedMaxObjectiveEvaluations = false;
    end
    
    methods
        %constructor
        function this = MultiModelOptimOptions(isClassification,inputArgs)
            remArgs = validatedAndUpdateLearners(this,isClassification,inputArgs);
            remArgs = validateAndUpdateHyperparameterOptions(this,isClassification,remArgs);
            this.RemainingArgs = remArgs;
        end
    end
    
    methods (Access = private)
        
        function remArgs = validatedAndUpdateLearners(this,isClassification,Args)
        %VALIDATEANDUPDATELEARNERS  validate and update the learners
        %   VALIDATEANDUPDATELEARNERS(THIS,ARGS) parses the input args for
        %   learners option, validates the learner names passed by the user
        %   and updates THIS.Learner based on the input. Returns the
        %   remaining arguments after parsing learners
            inputArgNames = {'learners'};
            defaultValues =  {'auto'};
            
            [learners,~,remArgs] = internal.stats.parseArgs(...
                inputArgNames,defaultValues, Args{:});
            learners = convertContainedStringsToChars(learners);
            learners = lower(learners);
            
            allLearners = {'svm','tree','linear','ensemble','kernel'};
            if isClassification
                allLearners = [allLearners{:},{'knn','nb','discr'}];
            else
                allLearners = [allLearners{:},{'gp'}];
            end
            
            % check for invalid learner datatypes
            if ~ischar(learners) && ~isstring(learners) && ~iscell(learners)
                error(message('stats:bayesoptim:multimodeloptimoptions:Learners'));
            end

            % convert to chars
            if isstring(learners) 
                learners = convertStringsToChars(learners);
                if ~iscell(learners)
                    learners = {learners};
                end
            elseif iscell(learners)
                [learners{:}] = convertStringsToChars(learners{:});
            elseif ischar(learners)
                learners = {learners};
            end

            %check if all, auto, all-linear and all-nonlinear are not passed
            %along with other learners
            
            if numel(learners) > 1 
                if any(ismember({'auto','all','all-linear','all-nonlinear'},learners))
                    error(message('stats:bayesoptim:multimodeloptimoptions:Learners'));
                end
            end
            
            if any(~ismember(learners,allLearners)) && any(~ismember(learners,{'auto','all','all-linear','all-nonlinear'}))
                    error(message('stats:bayesoptim:multimodeloptimoptions:Learners'));                
            end
            
            this.Learners = learners; 
        end
        
        function remArgs = validateAndUpdateHyperparameterOptions(this,isClassification,Args)
        %VALIDATEANDUPDATEHYPERPARAMETEROPTIONS  validate and update the 
        %HyperparameterOptimizationOptions
        %   VALIDATEANDUPDATEHYPERPARAMETEROPTIONS(THIS,ARGS) parses the input 
        %   args for HyperparameterOptimizationOptions, validates these options
        %   and updates the default values. Returns the remaining arguments after 
        %   parsing
        
            checkForInvalidOptimizationOptions(this,Args);
            
            [this.OptimizeHyperparameters, this.OptimizationOptions, remArgs] = ...
                classreg.learning.paramoptim.parseFitoptimizingArgs(Args, false);
            
            this.OptimizationOptions.Optimizer = 'bayesopt';
            this.OptimizationOptions.AcquisitionFunctionName = 'expected-improvement';
            this.OptimizationOptions.ModelType = 'MultiGaussianProcess';
            this.OptimizationOptions.IsClassregRegressionFunction = ~isClassification;

            updateValidationArgs(this);
            checkOptimizeHyperparameters(this);
        end

        function checkForInvalidOptimizationOptions(this,Args)
        %CHECKFORINVALIDOPTIMIZATIONOPTIONS  check if the user passed any
        %invalid hyperparameter optimization options
        %   CHECKFORINVALIDOPTIMIZATIONOPTIONS(THIS,ARGS) parses the input 
        %   args for any invalid HyperparameterOptimizationOptions. Some of
        %   the options like 'AcquisitionFunction', 'Optimizer', 
        %   'NumGridDivisions' are valid in regular fit functions but are
        %   not allowed in 'fitcauto'.
            [hyperOptimStruct,~,~] = internal.stats.parseArgs(...
                    'HyperparameterOptimizationOptions',[], Args{:});
                

            if ~isempty(hyperOptimStruct)
                if ~isstruct(hyperOptimStruct)
                    error(message('stats:bayesoptim:multimodeloptimoptions:OptimOptionsNotStruct'));
                end
                f = fieldnames(hyperOptimStruct);
                ArgList = cell(1,2*numel(f));
                ArgList(1:2:end) = f;
                ArgList(2:2:end) = struct2cell(hyperOptimStruct);
                Defaults = cell(1,numel(f));
                FieldNames = {'MaxObjectiveEvaluations','MaxTime', 'ShowPlots', 'SaveIntermediateResults',...
                    'Verbose','UseParallel','Repartition','CVPartition', 'Holdout', 'Kfold'};
                [values{1:numel(FieldNames)}, ~, extra] = internal.stats.parseArgs(...
                    FieldNames, Defaults, ArgList{:});
                if ~isempty(extra)
                    error(message('stats:bayesoptim:multimodeloptimoptions:BadStructField',extra{1}));
                end
                
                if ~isempty(values{1})
                    this.UserPassedMaxObjectiveEvaluations = true;
                end
                
                if ~isempty(values{5})
                    checkVerbose(this,values{5});
                end
            end
        end
        
        function updateValidationArgs(this)
        %UPDATEVALIDATIONARGS  update the validation args
        %   UPDATEVALIDATIONARGS(THIS) updates ValidationMethod and ValidationVal
        %   to be used in the bayesopt function.
            if ~isempty(this.OptimizationOptions.KFold)
                this.OptimizationOptions.ValidationMethod = 'KFold';
                this.OptimizationOptions.ValidationVal = this.OptimizationOptions.KFold;
            elseif ~isempty(this.OptimizationOptions.Holdout)
                this.OptimizationOptions.ValidationMethod = 'Holdout';
                this.OptimizationOptions.ValidationVal = this.OptimizationOptions.Holdout;
            elseif ~isempty(this.OptimizationOptions.CVPartition)
                this.OptimizationOptions.ValidationMethod = 'CVPartition';
                this.OptimizationOptions.ValidationVal = this.OptimizationOptions.CVPartition;
            end
        end
                
        function checkOptimizeHyperparameters(this)
        %CHECKOPTIMIZEHYPERPARAMETERS  checks the OptimizeHyperparameters
        %option
        %   CHECKOPTIMIZEHYPERPARAMETERS(THIS) throws error if the user
        %   passes any value other than 'auto' or 'all' along with 
        %   OptimizeHyperparameters option
            if ~ismember(this.OptimizeHyperparameters, {'auto', 'all'})
                error(message('stats:bayesoptim:multimodeloptimoptions:OptimizeHyperparameters'));
            end
        end
        
        function checkVerbose(~,verbose)
            if ~(isnumeric(verbose) && ismember(verbose,[0,1,2]))
                error(message('stats:bayesoptim:multimodeloptimoptions:Verbose'));
            end
        end
    end
end

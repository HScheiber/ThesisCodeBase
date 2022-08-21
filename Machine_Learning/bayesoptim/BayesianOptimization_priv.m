classdef BayesianOptimization_priv
    %BayesianOptimization    The results of a Bayesian Optimization.
    %
    %   BayesianOptimization Properties:
    %
    %     Problem definition properties:
    %       ObjectiveFcn   	- The ObjectiveFcn argument that was passed to
    %                       bayesopt.
    %       VariableDescriptions - The VariableDescriptions argument that
    %                       was passed to bayesopt.
    %       Options         - A read-only struct containing options that
    %                       were used to perform the optimization. It has
    %                       the following fields:
    %
    %                             AcquisitionFunctionName
    %                             IsObjectiveDeterministic
    %                             ExplorationRatio
    %
    %                             NumSeedPoints
    %                             MaxObjectiveEvaluations
    %                             MaxTime
    %
    %                             XConstraintFcn
    %                             ConditionalVariableFcn
    %                             NumCoupledConstraints
    %                             CoupledConstraintTolerances
    %                             AreCoupledConstraintsDeterministic
    %
    %                             Verbose
    %                             OutputFcn
    %                             SaveVariableName
    %                             SaveFileName
    %                             PlotFcn
    %
    %                             InitialX
    %                             InitialObjective
    %                             InitialConstraintViolations
    %                             InitialErrorValues
    %                             InitialObjectiveEvaluationTimes
    %                             InitialIterationTimes
    %                             InitialUserData
    %
    %     Solution properties:
    %    	MinObjective	- The minimum observed feasible Objective
    %                       function value, with feasibility defined by the
    %                       final constraint models (including the Error
    %                       constraint model).
    %     	XAtMinObjective - A 1-by-D table representing the value of X at
    %                       the minimum observed Objective function
    %                       value with feasibility defined by the final
    %                       constraint models (including the Error
    %                       constraint model).
    %    	MinEstimatedObjective - The minimum estimated feasible
    %                       Objective function value, according to the
    %                       model of Objective, with feasibility defined by
    %                       the final constraint models (including the
    %                       Error constraint model).
    %    	XAtMinEstimatedObjective - A 1-by-D table representing the
    %                       value of X at the minimum estimated feasible
    %                       Objective function value, according to the
    %                       model of Objective, with feasibility defined by
    %                       the final constraint models (including the
    %                       Error constraint model).
    %       NumObjectiveEvaluations - The number of Objective function
    %                       evaluations performed.
    %       TotalElapsedTime - The total elapsed time of the optimization.
    %       NextPoint       - A 1-by-D table representing the next point to
    %                       be evaluated were the optimization to continue.
    %
    %     Trace properties:
    %    	XTrace          - A T-by-D table of the variable values on
    %                       which ObjectiveFcn was evaluated, where T =
    %                       NumObjectiveEvaluations and D is the number of
    %                       variables.
    %     	ObjectiveTrace  - A T-by-1 array of Objective values resulting
    %                       from the T evaluations of ObjectiveFcn on
    %                       XTrace.
    %     	ConstraintsTrace - A T-by-K array of constraint values for K
    %                       coupled constraints resulting from the T
    %                       evaluations of ObjectiveFcn on XTrace.
    %     	UserDataTrace   - A T-by-1 cell array of the UserData outputs
    %                       of the Objective function evaluations.
    %     	ObjectiveEvaluationTimeTrace - A T-by-1 array of the runtime of
    %                       each of	the T evaluations of ObjectiveFcn on
    %                       XTrace.
    %     	IterationTimeTrace - A T-by-1 array of the total time taken for
    %                       each iteration, including both function
    %                       evaluation time and overhead.
    %    	ErrorTrace      - A T-by-1 vector of error values in {-1,1}
    %                       indicating on which points in XTrace the
    %                       objective function returned a nonfinite value.
    %     	FeasibilityTrace - A T-by-1 logical vector indicating which
    %                       points in XTrace are feasible according to the
    %                       final constraint models (including the Error
    %                       constraint model). The tolerances defined
    %                       in the 'CoupledConstraintTolerances' argument
    %                       are used to decide feasibility for each
    %                       constraint.
    %     	FeasibilityProbabilityTrace - A T-by-1 vector of the
    %                       probabilities that each point in XTrace is
    %                       feasible according to the final constraint
    %                       models (including the Error constraint model).
    %     	IndexOfMinimumTrace - A T-by-1 array of integer indices
    %                       indicating which element of the trace had the
    %                       minimum feasible Objective up to that
    %                       iteration, with feasibility defined by the
    %                       Constraint models at that iteration (including
    %                       the Error constraint model).
    %     	ObjectiveMinimumTrace - A T-by-1 array of the minimum feasible
    %                       Objective value found up to that iteration,
    %                       with feasibility defined by the Constraint
    %                       models at that iteration (including the Error
    %                       constraint model).
    %     	EstimatedObjectiveMinimumTrace - A T-by-1 array of the
    %                       estimated minimum Objective value found up to
    %                       that iteration, with feasibility defined by the
    %                       Constraint models at that iteration (including
    %                       the Error constraint model).
    %
    %   BayesianOptimization methods:
    %       bestPoint       - Return the best point according to a
    %                       criterion.
    %       plot            - Create plots from a set of plot functions.
    %       predictObjective - Predict the Objective at a set of points.
    %       predictObjectiveEvaluationTime - Predict the Objective function
    %                       evaluation time at a set of points.
    %       predictConstraints - Predict Constraint values at a set of
    %                       points.
    %       predictError    - Predict the Error value at a set of points.
    %       resume          - Resume a Bayesian optimization with modified
    %                       settings and stopping criteria.
    %
    %   See also: BAYESOPT, OPTIMIZABLEVARIABLE
    
    %   Copyright 2016-2020 The MathWorks, Inc.
    
    %% User API
    properties(Dependent, SetAccess=protected)
        % Problem definition
        ObjectiveFcn;
        VariableDescriptions;
        Options;
        % Solution / State of optimization
        MinObjective;
        XAtMinObjective;
        MinEstimatedObjective;
        XAtMinEstimatedObjective;
        NumObjectiveEvaluations;
        TotalElapsedTime;
        NextPoint;
        % Traces
        XTrace;
        ObjectiveTrace;
        ConstraintsTrace;
        UserDataTrace;
        ObjectiveEvaluationTimeTrace;
        IterationTimeTrace;
        ErrorTrace;
        FeasibilityTrace;
        FeasibilityProbabilityTrace;
        IndexOfMinimumTrace;
        ObjectiveMinimumTrace;
        EstimatedObjectiveMinimumTrace;
    end
    
    methods
        function [X, CriterionValue, Iteration] = bestPoint(this, varargin)
            % BESTPOINT Return the best point found in a Bayesian
            % optimization according to a criterion.
            %
            %   (1) X = bestPoint(RESULTS) returns the 1-row table X
            %       representing the best feasible point according to the
            %       default criterion.
            %
            %   (2) [X, CriterionValue] = bestPoint(RESULTS) also returns
            %       the value of the specified criterion at X.
            %
            %   (3) [X, CriterionValue, Iteration] = bestPoint(RESULTS)
            %       when Criterion is 'min-observed', 'min-visited-mean',
            %       or 'min-visited-upper-confidence-interval', also
            %       returns the iteration number at which the best point
            %       occurred.
            %
            %   (4) [...] = bestPoint(RESULTS, Name, Value,...)
            %       specifies additional name/value parameters:
            %
            %       Criterion   - Specifies the criterion used to choose
            %                   the best point. Possible values are the
            %                   following (case insensitive, hyphens are
            %                   optional, and unique prefixes are
            %                   recognized):
            %           'min-observed'  - X is the feasible point at which
            %                           the minimum observed Objective was
            %                           found.
            %           'min-mean'      - X is the feasible point at which
            %                           the mean of the Objective function
            %                           model is minimized.
            %           'min-upper-confidence-interval'
            %                           - X is the feasible point that
            %                           minimizes an upper confidence
            %                           interval according to the Objective
            %                           function model. See the ALPHA
            %                           parameter.
            %           'min-visited-mean'
            %                           - X is the feasible point at which
            %                           the mean of the Objective function
            %                           model is minimized, among those
            %                           points at which the Objective
            %                           function was evaluated.
            %           'min-visited-upper-confidence-interval'
            %                           - X is the feasible point that
            %                           minimizes an upper confidence
            %                           interval according to the Objective
            %                           function model, among those points
            %                           at which the Objective function was
            %                           evaluated. See the ALPHA parameter.
            %           Default: 'min-visited-upper-confidence-interval'
            %
            %       Alpha   - A probability value defining the upper
            %               confidence interval when CRITERION is
            %               'min-upper-confidence-interval' or
            %               'min-visited-upper-confidence-interval'. The
            %               upper confidence interval at X is the value Y
            %               such that Prob(mean(ObjectiveFcn(X))>Y)=Alpha
            %               according to the Objective function model.
            %               Default: .01
            [varargin{:}] = convertStringsToChars(varargin{:});
            C = bayesoptim_priv.suppressWarnings();
            % Parse Criterion and Alpha
            [Criterion, Alpha] = internal.stats.parseArgs({'Criterion', 'Alpha'}, ...
                {'minvisitedupperconfidenceinterval', .01}, varargin{:});
            % Parse Criterion value
            Crit = bayesoptim_priv.parseArgValue(Criterion, {'minobserved', ...
                'minmean', 'minupperconfidenceinterval', 'minvisitedmean',...
                'minvisitedupperconfidenceinterval'});
            % Check Alpha
            if ~(isscalar(Alpha) && isnumeric(Alpha) && isreal(Alpha) && Alpha>0 && Alpha<1)
                bayesoptim_priv.err('BadAlpha');
            end
            % Apply criterion
            Iteration = NaN;
            switch Crit
                case {'minobserved'}
                    [X, CriterionValue, Iteration] = minObservedPoint(this);
                case {'minmean'}
                    [X, CriterionValue] = minMeanPoint(this);
                case {'minvisitedmean'}
                    [X, CriterionValue, Iteration] = minVisitedMeanPoint(this);
                case {'minupperconfidenceinterval'}
                    [X, CriterionValue] = minUCIPoint(this, Alpha);
                case {'minvisitedupperconfidenceinterval'}
                    [X, CriterionValue, Iteration] = minVisitedUCIPoint(this, Alpha);
            end
        end
        
        function plot(this, varargin)
            %PLOT Plot Bayesian Optimization results
            %   (1) plot(RESULTS) or plot(RESULTS, 'all') calls all
            %       pre-defined plot functions on RESULTS.
            %
            %   (2) plot(RESULTS, plotFcn1, plotFcn2, ...) calls plot
            %       functions plotFcn1, plotFcn2, ... on RESULTS.
            %
            % Pre-defined plot functions:
            %     @plotObjectiveModel
            %     @plotAcquisitionFunction
            %     @plotObjectiveEvaluationTimeModel
            %     @plotConstraintModels
            %     @plotObjective
            %     @plotObjectiveEvaluationTime
            %     @plotMinObjective
            %     @plotElapsedTime
            [varargin{:}] = convertStringsToChars(varargin{:});
            C = bayesoptim_priv.suppressWarnings();
            if isempty(varargin) || isequal(varargin{1}, 'all')
                PlotFcn = {...
                    @plotObjectiveModel,...
                    @plotAcquisitionFunction,...
                    @plotObjectiveEvaluationTimeModel,...
                    @plotConstraintModels,...
                    @plotObjective,...
                    @plotObjectiveEvaluationTime,...
                    @plotMinObjective,...
                    @plotElapsedTime};
            else
                PlotFcn = varargin;
            end
            if ~all(cellfun(@(x)isa(x, 'function_handle'), PlotFcn))
                bayesoptim_priv.err('PlotArgs');
            end
            for i = 1:numel(PlotFcn)
                PlotFcn{i}(this, 'standalone');
            end
            drawnow;
        end
        
        function [ConstraintViolations, SDs] = predictConstraints(this, XTable)
            %PREDICTCONSTRAINTS Predict Constraint Violations at a
            %set of points.
            %   [ConstraintViolations, SDs] =
            %   predictConstraints(Results, XTable) returns the
            %   posterior means and standard deviations of all Constraints
            %   at the points in XTable, according to a Gaussian Process
            %   model of each constraint. ConstraintViolations and SDs are
            %   both N-by-K matrices, given N rows of XTable and K coupled
            %   constraints.
            C = bayesoptim_priv.suppressWarnings();
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictConstraints');
            if isempty(XTable) || this.PrivOptions.NumCoupledConstraints == 0
                ConstraintViolations = [];
                SDs = [];
            else
                for k = this.PrivOptions.NumCoupledConstraints : -1 : 1
                    if isempty(this.ConstraintModels{k})
                        ConstraintViolations(:,k) = NaN(height(XTable),1);
                        SDs(:,k) = NaN(height(XTable),1);
                    else
                        [ConstraintViolations(:,k), SDs(:,k)] = predict(...
                            this.ConstraintModels{k}, transformPoints(this, XTable));
                    end
                end
            end
        end
        
        function [Error, SD] = predictError(this, XTable)
            %PREDICTERROR Predict Error value at a set of points.
            %   [Error, SD] = predictError(Results, XTable) returns
            %   the posterior mean and standard deviation of the Error
            %   value at the points in XTable, according to a Gaussian
            %   Process model. Error and SD are both N-by-1 vectors, given
            %   N rows of XTable.
            C = bayesoptim_priv.suppressWarnings();
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictError');
            if isempty(XTable)
                Error = [];
                SD = [];
            elseif all(this.ErrorTrace < 0)
                Error = -ones(height(XTable),1);
                SD = zeros(height(XTable),1);
            elseif isempty(this.ErrorModel)
                Error = NaN(height(XTable),1);
                SD = NaN(height(XTable),1);
            else
                [Error, SD] = predict(this.ErrorModel, transformPoints(this, XTable));
            end
        end
        
        function [Objective, SD] = predictObjective(this, XTable)
            %PREDICTOBJECTIVE Predict the Objective function at a
            %set of points.
            %   [Objective, SD] = predictObjective(Results, XTable)
            %   returns the posterior mean and standard deviation of the
            %   objective function at the points in XTable, according to a
            %   Gaussian Process model. Objective and SD are both N-by-1
            %   vectors, given N rows of XTable.
            C = bayesoptim_priv.suppressWarnings();
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictObjective');
            if isempty(XTable)
                Objective = [];
                SD = [];
            elseif isempty(this.ObjectiveFcnModel)
                Objective = NaN(height(XTable),1);
                SD = NaN(height(XTable),1);
            else
                [Objective, SD] = predict(this.ObjectiveFcnModel, transformPoints(this, XTable));
            end
        end
        
        function ObjectiveEvaluationTime = predictObjectiveEvaluationTime(this, XTable)
            %PREDICTOBJECTIVEEVALUATIONTIME Predict the objective
            %function evaluation runtime at a set of points.
            %   ObjectiveEvaluationTime =
            %   predictObjectiveEvaluationTime(Results, XTable)
            %   returns the posterior mean of ObjectiveEvaluationTime at
            %   the points in XTable, according to a Gaussian Process
            %   model. ObjectiveEvaluationTime is an N-by-1 vector, given N
            %   rows of XTable.
            C = bayesoptim_priv.suppressWarnings();
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictObjectiveEvaluationTime');
            if isempty(XTable)
                ObjectiveEvaluationTime = [];
            elseif isempty(this.ObjectiveEvaluationTimeModel)
                ObjectiveEvaluationTime = NaN(height(XTable),1);
            else
                logObjectiveEvaluationTime = predict(this.ObjectiveEvaluationTimeModel, transformPoints(this, XTable));
                ObjectiveEvaluationTime = exp(logObjectiveEvaluationTime);
            end
        end
        
        function NewResults = resume(Results, varargin)
            %RESUME Resume a Bayesian Optimization.
            %   RESULTS = resume(RESULTS, 'PARAM1', val1, 'PARAM2'
            %   ,val2,...) continues running an optimization with modified
            %   options until new stopping criteria are met. All name/value
            %   parameters accepted by bayesopt are accepted, except
            %   'UseParallel', 'NumSeedPoints' and the 'Initial___'
            %   parameters. Important parameters for resuming are:
            %       'MaxObjectiveEvaluations'
            %                   - Specifies the additional number of
            %                     function evaluations to run. Default: 30
            %       'MaxTime'   - Specifies the additional number of
            %                     seconds to run. Default: Inf
            %       'VariableDescriptions'
            %                   - Specifies new variable descriptions to
            %                     use. The variable names to be optimized
            %                     must remain the same. Numeric variables
            %                     may have their Range, Type, and Transform
            %                     changed, but must remain numeric.
            %                     Categorical variables may not be changed.
            [varargin{:}] = convertStringsToChars(varargin{:});
            C = bayesoptim_priv.suppressWarnings();
            if iAnyIllegalResumeArgs(varargin)
                bayesoptim_priv.err('BadResumeArgs');
            else
                % Parse out new variable descriptions
                [PassedVariableDescrips, RemainingArgs] = iParseVariableDescriptionsArg(varargin);
                if isempty(PassedVariableDescrips)
                    VariableDescrips = Results.PrivOptions.VariableDescriptions;
                else
                    checkVariableDescriptions(PassedVariableDescrips, Results.PrivOptions.VariableDescriptions);
                    VariableDescrips = PassedVariableDescrips;
                end
                % Create NewOptions, install new variable descriptions
                NewOptions = Results.PrivOptions;
                NewOptions.VariableDescriptions = VariableDescrips;
                % Overwrite stopping criteria with defaults
                NewOptions.MaxObjectiveEvaluations = bayesoptim_priv.BayesoptOptions.DefaultMaxObjectiveEvaluations();
                NewOptions.MaxTime                 = bayesoptim_priv.BayesoptOptions.DefaultMaxTime();
                % Overwrite any options with contents of passed NVPs
                NewOptions = fillFromNVPs(NewOptions, RemainingArgs);
                % Go back and correct stopping criteria. Interpret them
                % as increments to the current state.
                NewOptions.MaxObjectiveEvaluations = NewOptions.MaxObjectiveEvaluations + Results.NumObjectiveEvaluations;
                NewOptions.MaxTime                 = NewOptions.MaxTime + Results.TotalElapsedTime;
                % Check updated options as a whole
                NewOptions = validateAndFillDefaults(NewOptions);
                % Save current XTrace
                CurXTrace = Results.XTrace;
                % Install new options (this would change XTrace)
                Results.PrivOptions = NewOptions;
                % Update XTrain because variableDescriptions may have changed
                Results.XTrain = XTraceToXTrain(Results, CurXTrace);
                % Run the optimization using the new options
                Results = initializeParallel(Results);
                Results = initializeTimers(Results, true);
                Resuming = true;
                NewResults = run(Results,Resuming);
            end
        end
    end
    
    %% Hidden API
    methods(Hidden)
        function this = BayesianOptimization_priv(Options)
            this.PrivOptions = Options;
            this = initializeParallel(this);
            this = initializeTimers(this, false);
            this = run(this);
        end
        
        % This is public so BayesoptParallel can call it
        function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                ObjectiveEvaluationTime, ObjectiveNargout] = callObjFcn(this, X)
            if ~isempty(this.PrivOptions.ObjectiveNargout)
                % nargout is known. Call ObjectiveFcn normally.
                ObjectiveNargout = this.PrivOptions.ObjectiveNargout;
                [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                    = callObjNormally(this, X);
            else
                % nargout is unknown. Try nargout=3:
                StartTic = tic;
                ErrorConstraintViolation = -1;
                try
                    [Objective, ConstraintViolations, UserData] = this.ObjectiveFcn(conditionalizeX(this, X));
                catch msg
                    % nargout may be < 3.
                    % Set nargout based on NumCoupledConstraints:
                    ObjectiveNargout = nargoutFromNumConstraints(this.PrivOptions.NumCoupledConstraints);
                    this.PrivOptions.ObjectiveNargout = ObjectiveNargout;
                    % nargout is now known. Call ObjectiveFcn normally.
                    [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                        = callObjNormally(this, X);
                    return;
                end
                % Calling with nargout=3 succeeded. Set nargout and finish up.
                ObjectiveNargout = 3;
                [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                    = finishObjEval(this, Objective, ConstraintViolations, -1, StartTic);
            end
        end
        
        % Output functions
        function Stop = assignInBase(bo, State)
            %assignInBase     OutputFcn for assigning a Bayesian
            %Optimization to a variable in the base workspace.
            %   Stop = assignInBase(bo, State) assigns the
            %   BayesianOptimization instance 'bo' to a variable in the
            %   base workspace. The variable name is obtained from
            %   bo.Options.SaveVariableName.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            narginchk(2,2);
            Stop = false;
            % Remove things that will cause a parallel pool to be opened
            % upon load
            bo.FantasyModels = [];
            bo.Parallel = [];
            % Do the assign
            try
                assignin('base', bo.PrivOptions.SaveVariableName, bo);
            catch me
                bayesoptim_priv.warn('AssignInBaseFailed');
                disp(me);
            end
        end
        
        function Stop = saveToFile(bo, State)
            %saveToFile     OutputFcn for saving a Bayesian Optimization to
            %a file.
            %   Stop = saveToFile(bo, State) saves the BayesianOptimization
            %   instance 'bo' to a file. The file name is obtained from
            %   bo.Options.SaveFileName.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            narginchk(2,2);
            % Remove things that will cause a parallel pool to be opened
            % upon load
            bo.FantasyModels = [];
            bo.Parallel = [];
            % Do the save
            Stop = false;
            SaveFilename = bo.PrivOptions.SaveFileName;
            prevFilename = sprintf('%s.PREV', SaveFilename);
            try
                if isfile(SaveFilename)
                    if isfile(prevFilename)
                        delete(prevFilename);
                    end
                    copyfile(SaveFilename, prevFilename);
                end
                BayesoptResults = bo;
                save(SaveFilename, 'BayesoptResults', '-v7.3');
            catch me
                if isfile(SaveFilename)
                    copyfile(SaveFilename,prevFilename);
                end
                save(SaveFilename, 'BayesoptResults', '-v7.3');
                bayesoptim_priv.warn('SaveToFileFailed', SaveFilename);
                disp(me);
                disp(['Current directory: ' pwd])
            end
        end
        
        % Predefined plot functions
        function Stop = plotAcquisitionFunction(bo, State)
            %plotAcquisitionFunction     PlotFcn for plotting the
            %Acquisition Function surface.
            %   Stop = plotAcquisitionFunction(bo, State) plots the
            %   Acquisition function surface. The value is set to NaN for
            %   points that violate the xConstraintFcn.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            if bo.NumVars > 2 || acquisitionFunctionIsGridOrRand(bo)
                return;
            end
            % Plot AF using fantasy models if UseParallel is true
            if bo.PrivOptions.UseParallel && ~isempty(bo.FantasyModels)
                bo = bo.FantasyModels;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updatePlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupPlot(bo);
                    updatePlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = W/3 + 1;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height],'Tag','bayesopt.AcqFcn');
                Axes = axes(f);
                title(Axes, bayesoptim_priv.infoString('AcqFcn', Results.PrivOptions.AcquisitionFunctionName));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updatePlot(Results, Axes, Curves)
                import bayesoptim_priv.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                FModel = Results.ObjectiveFcnModel;
                if acquisitionFunctionIsGridOrRand(Results)
                    return;
                end
                if isempty(FModel)
                    return;
                end
                ObjectiveEvaluationTimeModel = Results.ObjectiveEvaluationTimeModel;
                if isempty(ObjectiveEvaluationTimeModel) && ismember(Options.AcquisitionFunctionName, ...
                        {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'})
                    return;
                end
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = table2array(untransformPoints(Results, EvalGrid, false));
                        % Compute AF at gridpoints and at next point
                        switch Options.AcquisitionFunctionName
                            case 'probabilityofimprovement'
                                AFGrid = probabilityOfImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnModel.Sigma);
                                AFNext = probabilityOfImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnModel.Sigma);
                            case 'lowerconfidencebound'
                                kappa = Options.ExplorationRatio;
                                AFGrid = lowerConfidenceBound(LegalizedEvalGrid, ...
                                    FModel, kappa);
                                AFNext = lowerConfidenceBound(legalizePoints(Results, Results.XNext), ...
                                    FModel, kappa);
                            case {'expectedimprovement', 'expectedimprovementplus'}
                                AFGrid = expectedImprovement(LegalizedEvalGrid, ...
                                    FModel, ...
                                    Results.IncumbentF);
                                AFNext = expectedImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, ...
                                    Results.IncumbentF);
                            case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                                AFGrid = expectedImprovementPerSecond(LegalizedEvalGrid, ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                                AFNext = expectedImprovementPerSecond(legalizePoints(Results, Results.XNext), ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                        end
                        % NaN-out points that violate xConstraint
                        mask = double(satisfiesXConstraint(Results, EvalGrid));
                        mask(mask==0) = NaN;
                        AFGrid = AFGrid .* mask;
                        % Apply constraint weighting
                        AFGrid = AFGrid .* ProbAllConstraintsSatisfied(Results, EvalGrid);
                        % Plot surface
                        Curves{end+1} = plot(Axes, PlotGrid, AFGrid, 'r');
                        % Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, AFNext, MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, Options.AcquisitionFunctionName);
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = table2array(untransformPoints(Results, EvalGrid, false));
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % Compute AF at gridpoints and at next point
                        switch Options.AcquisitionFunctionName
                            case 'probabilityofimprovement'
                                AFGrid = probabilityOfImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnModel.Sigma);
                                AFNext = probabilityOfImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnModel.Sigma);
                            case 'lowerconfidencebound'
                                kappa = Options.ExplorationRatio;
                                AFGrid = lowerConfidenceBound(LegalizedEvalGrid, ...
                                    FModel, kappa);
                                AFNext = lowerConfidenceBound(legalizePoints(Results, Results.XNext), ...
                                    FModel, kappa);
                            case {'expectedimprovement', 'expectedimprovementplus'}
                                AFGrid = expectedImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF);
                                AFNext = expectedImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF);
                            case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                                AFGrid = expectedImprovementPerSecond(LegalizedEvalGrid, ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                                AFNext = expectedImprovementPerSecond(legalizePoints(Results, Results.XNext), ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                        end
                        % NaN-out points that violate xConstraint
                        mask = double(satisfiesXConstraint(Results, LegalizedEvalGrid));
                        mask(mask==0) = NaN;
                        AFGrid = AFGrid .* mask;
                        % Apply constraint weighting
                        AFGrid = AFGrid .* ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        
                        % Plot surface
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(AFGrid, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'b', 'FaceAlpha', .6);
                        % Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), AFNext, ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, Options.AcquisitionFunctionName);
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                end
            end
        end
        
        function Stop = plotConstraintModels(bo, State)
            %plotConstraintModels     PlotFcn for plotting Constraint-model
            %surfaces.
            %   Stop = plotConstraintModels(bo, State) plots each
            %   Constraint-model surface, including the built-in Error
            %   constraint if there were any errors. It also plots a
            %   Prob(Feasible) surface.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves NoErrorsYet;
            narginchk(2,2);
            Stop = false;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            NumConstraints = bo.PrivOptions.NumCoupledConstraints;
            % Order of plots is: {C1, ..., Cn, ProbFeas, Error}. Note that
            % the error plot will not exist until there is an error.
            switch State
                case 'initial'
                    % User constraints
                    for i = 1:NumConstraints
                        Axes{i} = setupPlot(bo, bayesoptim_priv.infoString('ConstraintNum', i), ...
                            bayesoptim_priv.infoString('DegreeOfViolation'), ...
                            sprintf('bayesopt.Constraint%d', i));
                        Curves{i} = {};
                    end
                    % Prob feasible
                    Axes{NumConstraints+1} = setupPlot(bo, bayesoptim_priv.infoString('ProbFeas'), ...
                        bayesoptim_priv.infoString('ProbFeas'), 'bayesopt.ProbFeas');
                    Curves{NumConstraints+1} = {};
                    % Error Constraint
                    NoErrorsYet = true;
                    if ~isempty(bo.ErrorModel)
                        NoErrorsYet = false;
                        Axes{NumConstraints+2} = setupPlot(bo, bayesoptim_priv.infoString('ErrorConstraint'), ...
                            bayesoptim_priv.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                        Curves{NumConstraints+2} = {};
                    end
                case {'iteration', 'done'}
                    % User constraints
                    for i = 1:NumConstraints
                        if isvalid(Axes{i})
                            [Axes{i}, Curves{i}] = updateConstraintPlot(bo, i, Axes{i}, Curves{i});
                        end
                    end
                    % Prob feasible
                    if isvalid(Axes{NumConstraints+1})
                        [Axes{NumConstraints+1}, Curves{NumConstraints+1}] = updateProbFeasPlot(bo, Axes{NumConstraints+1}, Curves{NumConstraints+1});
                    end
                    % Error Constraint
                    if ~isempty(bo.ErrorModel)
                        if NoErrorsYet
                            NoErrorsYet = false;
                            Axes{NumConstraints+2} = setupPlot(bo, bayesoptim_priv.infoString('ErrorConstraint'), ...
                                bayesoptim_priv.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                            Curves{NumConstraints+2} = {};
                        end
                        if isvalid(Axes{NumConstraints+2})
                            [Axes{NumConstraints+2}, Curves{NumConstraints+2}] = updateErrorPlot(...
                                bo, Axes{NumConstraints+2}, Curves{NumConstraints+2});
                        end
                    end
                case 'standalone'
                    % User constraints
                    for i = 1:NumConstraints
                        StandaloneAxes{i} = setupPlot(bo, bayesoptim_priv.infoString('ConstraintNum', i), ...
                            bayesoptim_priv.infoString('DegreeOfViolation'), ...
                            sprintf('bayesopt.Constraint%d', i));
                        updateConstraintPlot(bo, i, StandaloneAxes{i}, {});
                    end
                    % Prob feasible
                    StandaloneAxes{NumConstraints+1} = setupPlot(bo, bayesoptim_priv.infoString('ProbFeas'), ...
                        bayesoptim_priv.infoString('ProbFeas'), 'bayesopt.ProbFeas');
                    updateProbFeasPlot(bo, StandaloneAxes{NumConstraints+1}, {});
                    % Error Constraint
                    if ~isempty(bo.ErrorModel)
                        StandaloneAxes{NumConstraints+2} = setupPlot(bo, bayesoptim_priv.infoString('ErrorConstraint'), ...
                            bayesoptim_priv.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                        updateErrorPlot(bo, StandaloneAxes{NumConstraints+2}, {});
                    end
            end
            
            function Axes = setupPlot(Results, Title, ZLabel, Tag)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height],'Tag',Tag);
                Axes = axes(f);
                title(Axes, Title);
                zlabel(Axes, ZLabel)
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updateConstraintPlot(Results, ConstraintNum, Axes, Curves)
                import bayesoptim_priv.*
                if isempty(Results.ConstraintModels{ConstraintNum})
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ConstraintModels{ConstraintNum};
                XTrace = Results.XTrace;
                ValTrace = Results.ConstraintsTrace(:,ConstraintNum);
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points on all axes
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        Curves{end+1} = plot(Axes, [PlotGrid(1); PlotGrid(end)], [0;0], 'g');     % Plot line at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        ConstraintVals = predictConstraints(Results, Results.NextPoint);
                        Curves{end+1} = plot(Axes, XNextToPlot, ConstraintVals(ConstraintNum), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim_priv.infoString('ConstraintNum', ConstraintNum));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves{end+1} = surf(Axes, X1Plot, X2Plot, zeros(reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'g', 'FaceAlpha', .5, 'LineStyle','-');     % Plot surface at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        ConstraintVals = predictConstraints(Results, Results.NextPoint);
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ConstraintVals(ConstraintNum),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim_priv.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
            function [Axes, Curves] = updateErrorPlot(Results, Axes, Curves)
                import bayesoptim_priv.*
                if isempty(Results.ErrorModel)
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ErrorModel;
                XTrace = Results.XTrace;
                ValTrace = Results.ErrorTrace;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        Curves{end+1} = plot(Axes, [PlotGrid(1); PlotGrid(end)], [0;0], 'g');     % Plot line at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictError(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim_priv.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves{end+1} = surf(Axes, X1Plot, X2Plot, zeros(reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'g', 'FaceAlpha', .5, 'LineStyle','-');     % Plot surface at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            predictError(Results, Results.NextPoint),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim_priv.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
            function [Axes, Curves] = updateProbFeasPlot(Results, Axes, Curves)
                import bayesoptim_priv.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                XTrace = Results.XTrace;
                ValTrace = Results.FeasibilityProbabilityTrace;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        % (2) Plot training points on all axes
                        Curves{end+1} = plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, ...
                            ProbAllConstraintsSatisfied(Results, Results.XNext), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim_priv.infoString('ProbFeas'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            ProbAllConstraintsSatisfied(Results, Results.XNext),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim_priv.infoString('ProbFeas'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
        end
        
        function Stop = plotElapsedTime(bo, State)
            %plotElapsedTime     PlotFcn for plotting total elapsed time
            %versus the number of function evaluations.
            %   Stop = plotElapsedTime(bo, State) plots total elapsed time
            %   versus the number of function evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Title = bayesoptim_priv.infoString('ElapsedTimeTitle');
            XLabel = bayesoptim_priv.infoString('FEvals');
            YLabel = bayesoptim_priv.infoString('TotElapsed');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.ElapsedTime');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curves = {};
                case {'iteration', 'done'}
                    cellfun(@delete, Curves);
                    if isvalid(Axes)
                        obj = nancumsum(bo.ObjectiveEvaluationTimeTrace);
                        tot = bo.TotalElapsedTimeTrace;
                        Curves{end+1} = plotCurve(tot, Axes, 'b');
                        Curves{end+1} = plotCurve(obj, Axes, 'g');
                        if bo.PrivOptions.UseParallel
                            % Do not include modeling time, and relabel
                            % legend
                            legend(Axes, {bayesoptim_priv.infoString('RealTime'),...
                                bayesoptim_priv.infoString('FEvalTimeParallel')},...
                                'Location', 'best');
                        else
                            % Include modeling time
                            Curves{end+1} = plotCurve(tot-obj, Axes, 'r');
                            legend(Axes, {bayesoptim_priv.infoString('TotTime'),...
                                bayesoptim_priv.infoString('FEvalTime'),...
                                bayesoptim_priv.infoString('ModelingTime')},...
                                'Location', 'best');
                        end
                    end
                case 'standalone'
                    obj = nancumsum(bo.ObjectiveEvaluationTimeTrace);
                    tot = bo.TotalElapsedTimeTrace;
                    f = figure('Tag','bayesopt.ElapsedTime');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(tot, StandaloneAxes, 'b');
                    plotCurve(obj, StandaloneAxes, 'g');
                    if bo.PrivOptions.UseParallel
                        legend(StandaloneAxes, {bayesoptim_priv.infoString('RealTime'),...
                            bayesoptim_priv.infoString('FEvalTimeParallel')}, ...
                            'Location', 'best');
                    else
                        plotCurve(tot-obj, StandaloneAxes, 'r');
                        legend(StandaloneAxes, {bayesoptim_priv.infoString('TotTime'),...
                            bayesoptim_priv.infoString('FEvalTime'),...
                            bayesoptim_priv.infoString('ModelingTime')}, ...
                            'Location', 'best');
                    end
            end
            
            function Curve = plotCurve(Trace, Axes, color)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(Trace), Trace, ['-' color 'o'], 'MarkerSize', 2);
            end
        end
        
        function Stop = plotMinObjective(bo, State)
            %plotMinObjective     PlotFcn for plotting the minimum observed
            %Objective versus the number of function evaluations.
            %   Stop = plotMinObjective(bo, State) plots the minimum
            %   observed Objective versus the number of function
            %   evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes1 Curves1;
            narginchk(2,2);
            if bo.PrivOptions.FitcautoMultipleLearners || bo.PrivOptions.FitcAutoSingleLearner
                Title1 = bayesoptim_priv.infoString('MinObjTitleFitcauto');
                XLabel1 = bayesoptim_priv.infoString('FEvalsFitcautoPlot');
                if bo.PrivOptions.IsClassregRegressionFunction
                    YLabel = bayesoptim_priv.infoString('MinObjFitrauto');
                else
                    YLabel = bayesoptim_priv.infoString('MinObjFitcauto');
                end
            else
                Title1 = bayesoptim_priv.infoString('MinObjTitle');
                XLabel1 = bayesoptim_priv.infoString('FEvals');
                YLabel = bayesoptim_priv.infoString('MinObj');
            end
            
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.MinObjective');
                    Axes1 = axes(f);
                    title(Axes1, Title1);
                    xlabel(Axes1, XLabel1);
                    ylabel(Axes1, YLabel);
                    Axes1.YAxisLocation='right';
                    Curves1 = {};
                case {'iteration', 'done'}
                    cellfun(@delete, Curves1);
                    if isvalid(Axes1)
                        if bo.PrivOptions.FitcautoMultipleLearners || bo.PrivOptions.FitcAutoSingleLearner
                            Curves1{end+1} = plotCurveFitcauto(1:numel(bo.ObjectiveMinimumTrace), ...
                                bo.ObjectiveMinimumTrace, Axes1, '#0072bd',[0.00,0.45,0.74],[0.00,0.45,0.74]);
                        else
                            Curves1{end+1} = plotCurve(1:numel(bo.ObjectiveMinimumTrace), ...
                                bo.ObjectiveMinimumTrace, Axes1, 'b');
                        end
                        if all(isnan(bo.EstimatedObjectiveMinimumTrace))
                            if bo.PrivOptions.FitcautoMultipleLearners || bo.PrivOptions.FitcAutoSingleLearner
                                legend(Axes1, {bayesoptim_priv.infoString('MinObsObjFitcauto')}, 'Location', 'best');
                            else
                                legend(Axes1, {bayesoptim_priv.infoString('MinObsObj')}, 'Location', 'best');
                            end
                        else
                            if bo.PrivOptions.FitcautoMultipleLearners || bo.PrivOptions.FitcAutoSingleLearner
                                Curves1{end+1} = plotCurveFitcauto(1:numel(bo.ObjectiveMinimumTrace), ...
                                    bo.EstimatedObjectiveMinimumTrace, Axes1, '#a6dfff',[0.65,0.87,1.00],[0.65,0.87,1.00]);
                            else
                                Curves1{end+1} = plotCurve(1:numel(bo.ObjectiveMinimumTrace), ...
                                    bo.EstimatedObjectiveMinimumTrace, Axes1, 'g');
                            end
                            if bo.PrivOptions.FitcautoMultipleLearners || bo.PrivOptions.FitcAutoSingleLearner
                                legend(Axes1, {bayesoptim_priv.infoString('MinObsObjFitcauto'), ...
                                    bayesoptim_priv.infoString('EstMinObjFitcauto')}, 'Location', 'best');
                            else
                                 legend(Axes1, {bayesoptim_priv.infoString('MinObsObj'), ...
                                    bayesoptim_priv.infoString('EstMinObj')}, 'Location', 'best');
                            end
                                                            
                        end
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.MinObjective');
                    SAxes1 = axes(f);
                    title(SAxes1, Title1);
                    xlabel(SAxes1, XLabel1);
                    ylabel(SAxes1, YLabel);
                    SAxes1.YAxisLocation='right';
                    plotCurve(1:numel(bo.ObjectiveMinimumTrace), bo.ObjectiveMinimumTrace, SAxes1, 'b');
                    if all(isnan(bo.EstimatedObjectiveMinimumTrace))
                        legend(SAxes1, {bayesoptim_priv.infoString('MinObsObj')}, 'Location', 'best');
                    else
                        plotCurve(1:numel(bo.ObjectiveMinimumTrace), bo.EstimatedObjectiveMinimumTrace, SAxes1, 'g');
                        legend(SAxes1, {bayesoptim_priv.infoString('MinObsObj'), bayesoptim_priv.infoString('EstMinObj')},...
                            'Location', 'best');
                    end
            end
            
            function Curve = plotCurve(X, TraceToPlot, Axes, color)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, X, TraceToPlot, ['-' color 'o'], 'MarkerSize', 2);
            end
            
            function Curve = plotCurveFitcauto(X, TraceToPlot, Axes, color,markerFaceColor,markerEdgeColor)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, X, TraceToPlot, ['-' 'o'], 'MarkerSize', 2,...
                    'color',color,'MarkerFaceColor',markerFaceColor,'MarkerEdgeColor',markerEdgeColor);
            end
        end
        
        function Stop = plotObjective(bo, State)
            %plotObjective     PlotFcn for plotting each observed Objective
            %function value versus the number of function evaluations.
            %   Stop = plotObjective(bo, State) plots each observed
            %   Objective function value versus the number of function
            %   evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curve;
            narginchk(2,2);
            Title = bayesoptim_priv.infoString('ObjFunTitle');
            XLabel = bayesoptim_priv.infoString('FEvals');
            YLabel = bayesoptim_priv.infoString('ObjFun');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.Objective');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curve = [];
                case {'iteration', 'done'}
                    delete(Curve);
                    if isvalid(Axes)
                        Curve = plotCurve(bo.ObjectiveTrace, Axes);
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.Objective');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(bo.ObjectiveTrace, StandaloneAxes);
            end
            
            function Curve = plotCurve(ObjectiveTrace, Axes)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(ObjectiveTrace), ObjectiveTrace, '-bo', 'MarkerSize', 2);
            end
            
        end
        
        function Stop = plotObjectiveEvaluationTime(bo, State)
            %plotObjectiveEvaluationTime     PlotFcn for plotting each
            %observed Objective function evaluation time versus the number
            %of function evaluations.
            %   Stop = plotObjectiveEvaluationTime(bo, State) plots each
            %   observed Objective function evaluation time versus the
            %   number of function evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curve;
            narginchk(2,2);
            Title = bayesoptim_priv.infoString('FEvalTimeTitle');
            XLabel = bayesoptim_priv.infoString('FEvals');
            YLabel = bayesoptim_priv.infoString('FEvalTime');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.ObjectiveEvaluationTime');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curve = [];
                case {'iteration', 'done'}
                    delete(Curve);
                    if isvalid(Axes)
                        Curve = plotCurve(bo.ObjectiveEvaluationTimeTrace, Axes);
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.ObjectiveEvaluationTime');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(bo.ObjectiveEvaluationTimeTrace, StandaloneAxes);
            end
            
            function Curve = plotCurve(Trace, Axes)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(Trace), Trace, '-bo', 'MarkerSize', 2);
            end
        end
        
        function Stop = plotObjectiveEvaluationTimeModel(bo, State)
            %plotObjectiveEvaluationTimeModel     PlotFcn for plotting the
            %ObjectiveEvaluationTime-model surface.
            %   Stop = plotObjectiveEvaluationTimeModel(bo, State) plots
            %   the ObjectiveEvaluationTime-model surface.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updatePlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupPlot(bo);
                    updatePlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height], 'Tag','bayesopt.ObjectiveEvaluationTimeModel');
                Axes = axes(f);
                title(Axes, bayesoptim_priv.infoString('FEvalTimeModel'));
                zlabel(Axes, bayesoptim_priv.infoString('EstFEvalTime'));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updatePlot(Results, Axes, Curves)
                import bayesoptim_priv.*
                if isempty(Results.ObjectiveEvaluationTimeModel)
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ObjectiveEvaluationTimeModel;
                XTrace = Results.XTrace;
                ValTrace = Results.ObjectiveEvaluationTimeTrace;
                ValTrace(isnan(Results.ObjectiveTrace)) = NaN;      % Only plot valid function evaluations
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = exp(predict(Model, LegalizedEvalGrid));
                        % (2) Plot training points on all axes
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model's mean
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictObjectiveEvaluationTime(...
                            Results, Results.NextPoint), MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim_priv.infoString('FEvalTime'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = exp(predict(Model, LegalizedEvalGrid));
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot FModel's mean and minimum
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            predictObjectiveEvaluationTime(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
        end
        
        function Stop = plotObjectiveModel(bo, State)
            %plotObjectiveModel     PlotFcn for plotting the Objective
            %function model surface, along with the estimated location of
            %the minimum.
            %   Stop = plotObjectiveModel(bo, State) plots the Objective
            %   function model surface, along with the estimated location
            %   of the minimum.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            VarSpec = bo.PrivOptions.VarSpec;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupFPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updateFPlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupFPlot(bo);
                    updateFPlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupFPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                % (1) Plot True Objective on grid
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height], 'Tag','bayesopt.ObjectiveModel');
                Axes = axes(f);
                title(Axes, bayesoptim_priv.infoString('ObjModel'));
                zlabel(Axes, bayesoptim_priv.infoString('EstObj'));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updateFPlot(Results, Axes, Curves)
                import bayesoptim_priv.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                FModel = Results.ObjectiveFcnModel;
                if isempty(FModel)
                    return;
                end
                XTrace = Results.XTrace;
                ObjectiveTrace = Results.ObjectiveTrace;
                ModelMinXTable = untransformPoints(Results, Results.IncumbentX, true);
                ModelMinF = Results.IncumbentF;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        % Objective
                        [FPred, YSD, ~] = predict(FModel, LegalizedEvalGrid);
                        FSD = funStd(YSD, FModel.Sigma);
                        % (2) Plot training points on all axes
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        Curves{end+1} = plot(Axes, XPlottable, ObjectiveTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (2.5) Plot pending "fantasy" points if any
                        if ~isempty(Results.Parallel) && anyPendingX(Results.Parallel)
                            PendingXTable   = conditionalizeX(Results, Results.Parallel.PendingX);
                            PendingXToPlot  = makeXTablePlottable(Results, PendingXTable);
                            PendingObj      = Results.PendingObj;
                            Curves{end+1}   = plot(Axes, PendingXToPlot, PendingObj, 'go', 'MarkerFaceColor', 'green');
                        end
                        % (3) Plot FModel's mean, envelope
                        Curves{end+1} = plot(Axes, PlotGrid, FPred, 'r');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred + FSD, 'b');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred - FSD, 'b');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred + FModel.Sigma, 'c');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred - FModel.Sigma, 'c');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictObjective(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (6.1) Plot minimum if it exists
                        if ~isempty(ModelMinXTable)
                            XMinPlottable = makeXTablePlottable(Results, ModelMinXTable);
                            Curves{end+1} = plot(Axes, XMinPlottable, ModelMinF, 'r*');
                        end
                        % (6.2) Plot legend
                        if ~isempty(Results.Parallel) && anyPendingX(Results.Parallel)
                            LegendStrings = {...
                                bayesoptim_priv.infoString('Observed'),...
                                bayesoptim_priv.infoString('Evaluating'),...
                                bayesoptim_priv.infoString('ModelMean'),...
                                bayesoptim_priv.infoString('ModelErrorbars'),...
                                bayesoptim_priv.infoString('NoiseErrorbars'),...
                                bayesoptim_priv.infoString('NextPoint')};
                            CurveNums = [1,2,3,4,6,8];
                        else
                            LegendStrings = {...
                                bayesoptim_priv.infoString('Observed'),...
                                bayesoptim_priv.infoString('ModelMean'),...
                                bayesoptim_priv.infoString('ModelErrorbars'),...
                                bayesoptim_priv.infoString('NoiseErrorbars'),...
                                bayesoptim_priv.infoString('NextPoint')};
                            CurveNums = [1,2,3,5,7];
                        end
                        if ~isempty(ModelMinXTable)
                            LegendStrings{end+1} = bayesoptim_priv.infoString('ModelMin');
                            CurveNums(end+1) = CurveNums(end)+1;
                        end
                        legend(Axes, [Curves{CurveNums}], LegendStrings{:}, 'Location', 'best');
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim_priv.infoString('EstObj'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        % Objective
                        FPred = predict(FModel, LegalizedEvalGrid);
                        % (2) Plot training points
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ObjectiveTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (2.5) Plot pending "fantasy" points if any
                        if ~isempty(Results.Parallel) && anyPendingX(Results.Parallel)
                            PendingXTable   = conditionalizeX(Results, Results.Parallel.PendingX);
                            PendingXToPlot  = makeXTablePlottable(Results, PendingXTable);
                            PendingObj      = Results.PendingObj;
                            Curves{end+1} = plot3(Axes, PendingXToPlot(:,1), PendingXToPlot(:,2), PendingObj, ...
                                'go', 'MarkerFaceColor', 'green');
                        end
                        % (3) Plot FModel's mean
                        SurfAndContour = surfc(Axes, X1Plot, X2Plot, reshape(FPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves = [Curves, {SurfAndContour(1), SurfAndContour(2)}];
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), predictObjective(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (6.1) Plot model min if exists
                        if ~isempty(ModelMinXTable)
                            XMinPlottable = makeXTablePlottable(Results, ModelMinXTable);
                            Curves{end+1} = plot3(Axes, XMinPlottable(1), XMinPlottable(2), ModelMinF, 'r*');
                        end
                        
                        % (6.2) Plot legend
                        if ~isempty(Results.Parallel) &&  anyPendingX(Results.Parallel)
                            LegendStrings = {...
                                bayesoptim_priv.infoString('Observed'),...
                                bayesoptim_priv.infoString('Evaluating'),...
                                bayesoptim_priv.infoString('ModelMean'),...
                                bayesoptim_priv.infoString('NextPoint')};
                            CurveNums = [1,2,3,5];
                        else
                            LegendStrings = {...
                                bayesoptim_priv.infoString('Observed'),...
                                bayesoptim_priv.infoString('ModelMean'),...
                                bayesoptim_priv.infoString('NextPoint')};
                            CurveNums = [1,2,4];
                        end
                        if ~isempty(ModelMinXTable)
                            LegendStrings{end+1} = bayesoptim_priv.infoString('ModelMin');
                            CurveNums(end+1) = CurveNums(end)+1;
                        end
                        legend(Axes, [Curves{CurveNums}], LegendStrings{:}, 'Location', 'best');
                        
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                end
            end
        end
    end
    
    %% Internal
    properties(Hidden, SetAccess=protected)
        % Model training data
        XTrain;
        FTrain;
        ObjectiveEvaluationTimeTrain;
        ConstraintTrain;
        ErrorTrain;
        QueuedX;
        % Parallel
        Parallel;       % A bayesoptim.BayesoptParallel object.
        FantasyModels;
        PendingObj;     % Imputed values for pending X.
        % Models
        ModelStrategy;  % keep different models (GP and RF)
        ObjectiveFcnModel;
        ObjectiveEvaluationTimeModel;
        ConstraintModels;
        ErrorModel;
        % Old model properties that existed before 2020a for backward
        % loading compatibility:
        ObjectiveFcnGP
        ObjectiveEvaluationTimeGP
        ConstraintGPs
        ErrorGP
        % Current iteration info
        IncumbentF;
        IncumbentX;
        XNext;
        PlotFcnStop;
        OutputFcnStop
        % timing
        IterationStartTime;
        ObjectiveEvaluationTime;
        PreviousTime;
        CurrentRunStartTime;
        TotalElapsedTimeTrace;
        % other
        SampledGridIndices;
        PrivOptions;
        PrivIterationTimeTrace;
        PrivIndexOfMinimumTrace;
        PrivFMinTrace;
        PrivFMinEstimatedTrace;
        PrivUserDataTrace;
    end
    
    methods     % Dependent property getters
        
        function ObjectiveFcn = get.ObjectiveFcn(this)
            ObjectiveFcn = this.PrivOptions.ObjectiveFcn;
        end
        
        function VariableDescriptions = get.VariableDescriptions(this)
            VariableDescriptions = this.PrivOptions.VariableDescriptions;
        end
        
        function this = set.VariableDescriptions(this, VariableDescriptions)
            this.PrivOptions.VariableDescriptions = VariableDescriptions;
        end
        
        function s = get.Options(this)
            % Return all the visible public options properties in a struct.
            props = properties(this.PrivOptions);
            for i = 1:numel(props)
                s.(props{i}) = this.PrivOptions.(props{i});
            end
        end
        
        function XAtMinObjective = get.XAtMinObjective(this)
            XAtMinObjective = bestPoint(this, 'Criterion', 'min-observed');
        end
        
        function MinObjective = get.MinObjective(this)
            [~, MinObjective] = bestPoint(this, 'Criterion', 'min-observed');
        end
        
        function XAtMinEstimatedObjective = get.XAtMinEstimatedObjective(this)
            if  this.PrivOptions.FitcautoMultipleLearners || ...
                                     this.PrivOptions.FitcAutoSingleLearner
                XAtMinEstimatedObjective = bestPoint(this,'Criterion','min-visited-mean');
            else
                XAtMinEstimatedObjective = bestPoint(this);
            end
        end
        
        function MinEstimatedObjective = get.MinEstimatedObjective(this)
            X = this.XAtMinEstimatedObjective;
            if isempty(X)
                MinEstimatedObjective = NaN;
            else
                MinEstimatedObjective = predictObjective(this, X);
            end
        end
        
        function NumObjectiveEvaluations = get.NumObjectiveEvaluations(this)
            NumObjectiveEvaluations = size(this.XTrain,1);
        end
        
        function TotalElapsedTime = get.TotalElapsedTime(this)
            TotalElapsedTime = this.TotalElapsedTimeTrace(end);
        end
        
        function NextPoint = get.NextPoint(this)
            if isempty(this.XNext)
                NextPoint = table;
            else
                NextPoint = untransformPoints(this, this.XNext, true);
                NextPoint = applyConditionalVariableFcn(this, NextPoint);
                NextPoint = canonicalizePoints(this, NextPoint);
            end
        end
        
        function XTrace = get.XTrace(this)
            XTrace = conditionalizeX(this, this.XTrain);
        end
        
        function ObjectiveTrace = get.ObjectiveTrace(this)
            ObjectiveTrace = this.FTrain(:);
        end
        
        function ObjectiveEvaluationTimeTrace = get.ObjectiveEvaluationTimeTrace(this)
            ObjectiveEvaluationTimeTrace = this.ObjectiveEvaluationTimeTrain(:);
        end
        
        function ConstraintsTrace = get.ConstraintsTrace(this)
            ConstraintsTrace = this.ConstraintTrain;
        end
        
        function UserDataTrace = get.UserDataTrace(this)
            UserDataTrace = this.PrivUserDataTrace;
        end
        
        function ErrorTrace = get.ErrorTrace(this)
            ErrorTrace = this.ErrorTrain;
        end
        
        function FeasibilityTrace = get.FeasibilityTrace(this)
            FeasibilityTrace = allConstraintsSatisfied(this, this.XTrain, ...
                this.ConstraintModels, this.ErrorModel);
        end
        
        function FeasibilityProbabilityTrace = get.FeasibilityProbabilityTrace(this)
            FeasibilityProbabilityTrace(:,1) = ProbAllConstraintsSatisfied(this, this.XTrain);
        end
        
        function IndexOfMinimumTrace = get.IndexOfMinimumTrace(this)
            IndexOfMinimumTrace = this.PrivIndexOfMinimumTrace(:);
        end
        
        function ObjectiveMinimumTrace = get.ObjectiveMinimumTrace(this)
            ObjectiveMinimumTrace = this.PrivFMinTrace(:);
        end
        
        function EstimatedObjectiveMinimumTrace = get.EstimatedObjectiveMinimumTrace(this)
            EstimatedObjectiveMinimumTrace = this.PrivFMinEstimatedTrace(:);
        end
        
        function IterationTimeTrace = get.IterationTimeTrace(this)
            IterationTimeTrace = this.PrivIterationTimeTrace(:);
        end
    end
    
    methods(Access=protected)
        %% Main algorithm
        function this = run(this,varargin)
            if nargin > 1
                resuming = varargin{1};
            else
                resuming = false;
            end
            if this.PrivOptions.UseParallel
                this = runParallel(this);
            else
                this = runSerial(this,resuming);
            end
        end
        
        function this = initializeParallel(this)
            if this.PrivOptions.UseParallel
                this.Parallel = bayesoptim_priv.BayesoptParallel(this.ObjectiveFcn, this.PrivOptions.NumCoupledConstraints, this.PrivOptions.Verbose);
                this = checkAndSetMinWorkerUtilization(this);
                if ~isempty(this.Parallel.ObjectiveNargout)
                    this.PrivOptions.ObjectiveNargout = this.Parallel.ObjectiveNargout;
                end
            end
        end
        
        function this = checkAndSetMinWorkerUtilization(this)
            if isempty(this.PrivOptions.MinWorkerUtilization)
                this.PrivOptions.MinWorkerUtilization = floor(this.Parallel.NumWorkers * this.PrivOptions.DefaultMinWorkerUtilizationFraction);
            elseif this.PrivOptions.MinWorkerUtilization > this.Parallel.NumWorkers
                bayesoptim_priv.warn('MinWorkerUtilizationTooLarge', num2str(this.PrivOptions.MinWorkerUtilization), num2str(this.Parallel.NumWorkers));
                this.PrivOptions.MinWorkerUtilization = this.Parallel.NumWorkers;
            end
        end
        
        function this = initializeTimers(this, isResuming)
            this.CurrentRunStartTime = tic;
            if isResuming
                this.PreviousTime = this.TotalElapsedTime;
            else
                this.PreviousTime = 0;
            end
        end
        
        function this = runSerial(this,resuming)
            checkForOptimizableVariables(this);
            checkXConstraintFcnSatisfiability(this);
            this = processInitializationData(this);
            this.IterationStartTime = tic;
            this = callPlotFcn(this, 'initial');
            this = callOutputFcn(this, 'initial');
            this = fitModels(this);
            [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
            if resuming
                this.XNext = this.NextPoint{:,:};
            else
                this = chooseNextPoint(this);
            end
            iteration = this.NumObjectiveEvaluations + 1;          
            if iteration == 1 && this.PrivOptions.Verbose >= 1 && ...
                    (this.PrivOptions.FitcAutoSingleLearner || this.PrivOptions.FitcautoMultipleLearners)
                printFitcautoVerboseSummary(this);
            end
            while ~optimizationFinished(this, iteration)
                % Maybe print verbose line
                verboseXNextDisplay(this);
                % Do a function evaluation and record results
                [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                    ObjectiveFcnObjectiveEvaluationTime, ObjectiveNargout] = callObjFcn(this, this.XNext);
                this = recordFEvalResults(this, this.XNext, Objective, ConstraintViolations, ...
                    UserData, ErrorConstraintViolation, ObjectiveFcnObjectiveEvaluationTime, ObjectiveNargout);
                % Refit models
                FitModelsT = tic;
                this = fitModels(this);
                ModelFitTime = toc(FitModelsT);
                % Update traces
                this = updateTraces(this, iteration);
                % Find incumbent
                PointSelectionT = tic;
                [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
                % Choose next point
                this = chooseNextPoint(this);
                if this.PrivOptions.Verbose >= 2
                    bayesoptim_priv.printInfo('PointSelectionTime', num2str(toc(PointSelectionT)));
                end
                % Update timing, plots and output
                if this.PrivOptions.Verbose >= 2
                    bayesoptim_priv.printInfo('ModelFitTime', num2str(ModelFitTime));
                end
                this = updateTimingPlotsAndOutput(this, iteration);
                iteration = iteration + 1;
            end
            this = pruneGridIndices(this);
            this = callPlotFcn(this, 'done');
            this = callOutputFcn(this, 'done');
            showStoppingReason(this);
            showBestPoints(this);
            showVerboseDocHelp(this);
        end
        
        function Parallel = launchMultipleJobs(this, Xs)
            Parallel = this.Parallel;
            for Idx = 1:size(Xs,1)
                Parallel = launchJob(Parallel, this.PrivOptions.ObjectiveNargout, this.PrivOptions.NumCoupledConstraints, ...
                    Xs(Idx,:), conditionalizeX(this, Xs(Idx,:)));
            end
        end
        
        function this = runParallel(this)
            CleanupFutures = onCleanup(@()cancelAllJobs(this.Parallel));
            checkForOptimizableVariables(this);
            checkXConstraintFcnSatisfiability(this);
            this = processInitializationData(this);
            this.IterationStartTime = tic;
            this = callPlotFcn(this, 'initial');
            this = callOutputFcn(this, 'initial');
            this = fitModels(this);
            [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
            if this.NumObjectiveEvaluations == 0
                % We're not resuming, so add more points to the queue if
                % necessary to fill it...
                this = padQueuedX(this, this.Parallel.NumWorkers);
                % ... and launch the first NumWorkers points in the queue; Then remove them from queue.
                this.Parallel = launchMultipleJobs(this, this.QueuedX(1:this.Parallel.NumWorkers,:));
                this.QueuedX(1:this.Parallel.NumWorkers,:) = [];
            end
            this = chooseNextPointParallel(this);
            iteration = this.NumObjectiveEvaluations + 1;
            % Main loop
            if iteration == 1 && this.PrivOptions.Verbose >= 1 && ...
                    (this.PrivOptions.FitcAutoSingleLearner || this.PrivOptions.FitcautoMultipleLearners)
                printFitcautoVerboseSummary(this);
            end
            while ~optimizationFinished(this, iteration)
                % If there's a free worker and XNext, launch XNext.
                if anyFreeWorkers(this.Parallel) && ~isempty(this.XNext)
                    this.Parallel = launchJob(this.Parallel, this.PrivOptions.ObjectiveNargout, ...
                        this.PrivOptions.NumCoupledConstraints, this.XNext, conditionalizeX(this, this.XNext));
                end
                % If workers are under-utilized, launch jobs to fill all workers
                if numPendingX(this.Parallel) < this.PrivOptions.MinWorkerUtilization
                    NumFreeWorkers = this.Parallel.NumWorkers - numPendingX(this.Parallel);
                    [Xs, this] = randomOrGridXFeasiblePoints(this, NumFreeWorkers);
                    if this.PrivOptions.Verbose >= 2 && ~isempty(Xs)
                        bayesoptim_priv.printInfo('MultipleJobLaunch', num2str(numPendingX(this.Parallel)), num2str(size(Xs,1)));
                    end
                    this.Parallel = launchMultipleJobs(this, Xs);
                end
                % If there are no finished jobs and no free workers, wait for a job to finish.
                if ~anyFinishedJobs(this.Parallel) && ~anyFreeWorkers(this.Parallel)
                    waitForAnyJobToFinish(this.Parallel)
                end
                % If any jobs finished, record them all.
                NumFinishedJobs = 0;
                while anyFinishedJobs(this.Parallel)
                    NumFinishedJobs = NumFinishedJobs + 1;
                    [FinishedX, Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                        TotalJobTime, LaunchTime, ObjectiveNargout, this.Parallel] = getJobResults(this.Parallel);
                    this = recordFEvalResults(this, FinishedX, Objective, ConstraintViolations, ...
                        UserData, ErrorConstraintViolation, TotalJobTime, ObjectiveNargout);
                end
                FinishedIterationNums = iteration : iteration+NumFinishedJobs-1; % The iteration numbers of the finished fevals.
                % If there were any finished jobs, refit the models and
                % find the new incumbent
                PointSelectionT = [];
                if NumFinishedJobs > 0
                    FitModelsT = tic;
                    this = fitModels(this);
                    if this.PrivOptions.Verbose >= 2
                        bayesoptim_priv.printInfo('ModelFitTime', num2str(toc(FitModelsT)));
                    end
                    this = updateTraces(this, FinishedIterationNums);
                    PointSelectionT = tic;
                    [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
                end
                if isempty(PointSelectionT)
                    PointSelectionT = tic;
                end
                % Choose a new point to evaluate
                this = chooseNextPointParallel(this);
                if this.PrivOptions.Verbose >= 2
                    bayesoptim_priv.printInfo('PointSelectionTime', num2str(toc(PointSelectionT)));
                end
                % Update timing data for all newly-finished iterations, if any.
                this = updateTimingPlotsAndOutput(this, FinishedIterationNums);
                % Update iteration number
                if ~isempty(FinishedIterationNums)
                    iteration = FinishedIterationNums(end)+1;
                end
            end
            this = pruneGridIndices(this);
            this = callPlotFcn(this, 'done');
            this = callOutputFcn(this, 'done');
            showStoppingReason(this);
            showBestPoints(this);
            showVerboseDocHelp(this);
            % Remove things that will cause a parallel pool to be opened
            % upon load
            this.FantasyModels = [];
            this.Parallel = [];
        end
        
        function this = updateTimingPlotsAndOutput(this, FinishedIterationNums)
            % If running in parallel, there may be no new finished
            % iterations, but we still want to plot because there may be
            % new pending points.
            if ~isempty(FinishedIterationNums)
                this.PrivIterationTimeTrace(FinishedIterationNums) = toc(this.IterationStartTime);
                this.IterationStartTime = tic;
                arrayfun(@(i)printProgress(this, i, numel(FinishedIterationNums)), FinishedIterationNums);
                % Update total elapsed times
                this.TotalElapsedTimeTrace = [this.TotalElapsedTimeTrace;
                    repmat(this.PreviousTime + toc(this.CurrentRunStartTime), numel(FinishedIterationNums), 1)];
                this = callOutputFcn(this, 'iteration');
            end
            %  Update plots
            this = callPlotFcn(this, 'iteration');
        end
        
        function this = recordFEvalResults(this, X, Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                TotalJobTime, ObjectiveNargout)
            % Store the results of an feval in the training sets
            this.PrivOptions.ObjectiveNargout        = ObjectiveNargout;
            this.XTrain(end+1, :)                    = legalizePoints(this, X);
            this.FTrain(end+1)                       = Objective;
            this.PrivUserDataTrace{end+1,1}          = UserData;
            this.ErrorTrain(end+1,1)                 = ErrorConstraintViolation;
            this.ObjectiveEvaluationTimeTrain(end+1) = TotalJobTime;
            if this.PrivOptions.NumCoupledConstraints > 0
                this.ConstraintTrain(end+1,:)        = ConstraintViolations;
            end
        end
        
        function this = updateTraces(this, IterationNums)
            % Update Objective minimum traces and model-based traces
            [~, MinObsObjective, MinObsLoc]             = bestPoint(this, 'Criterion', 'min-observed');
            this.PrivIndexOfMinimumTrace(IterationNums) = MinObsLoc;
            this.PrivFMinTrace(IterationNums)           = MinObsObjective;
            this.PrivFMinEstimatedTrace(IterationNums)  = this.MinEstimatedObjective;
        end
        
        function this = padQueuedX(this, NumWorkers)
            % Add random points to this.QueuedX until its height is NumWorkers, or
            % until there are no more grid points.
            n = size(this.QueuedX,1);
            if n < NumWorkers
                if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                    % Add grid X-feasible points
                    for i = n+1 : NumWorkers
                        [X, this] = gridXFeasiblePoint(this);
                        if isempty(X)
                            % No more free grid points.
                            break;
                        else
                            this.QueuedX(i,:) = X;
                        end
                    end
                else
                    % Add random X-feasible points
                    this.QueuedX = [this.QueuedX; randomXFeasiblePoints(this, NumWorkers-n)];
                end
            end
        end
        
        function [Xs, this] = randomOrGridXFeasiblePoints(this, N)
            % Find up to N random [grid] points, or until there are no more
            % grid points. 'this' is updated to reflect sampling without
            % replacement for the grid.
            if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                % Add grid X-feasible points
                for i = 1:N
                    [X, this] = gridXFeasiblePoint(this);
                    if isempty(X)
                        % No more free grid points.
                        break;
                    else
                        Xs(i,:) = X;
                    end
                end
            else
                % Add random X-feasible points
                Xs = randomXFeasiblePoints(this, N);
            end
        end
        
        function this = processInitializationData(this)
            % Set this.QueuedX to be a set of points to evaluate first.
            if this.NumObjectiveEvaluations > 0
                % Only process initialization data if not resuming.
                return;
            end
            Opts = this.PrivOptions;
            if Opts.FitcautoMultipleLearners && any(strcmp(Opts.ModelType, {'MultiTreeBagger','GaussianProcess','MultiGaussianProcess'})) && Opts.NumSeedPoints > 0
                learnerIndex = 1;
                XPoints     = randomXFeasiblePointsForEachLearner(this, Opts.NumSeedPoints,learnerIndex);
                XTable      = untransformPoints(this, XPoints, true);
                idx = randperm(size(XTable,1));
                XTable = XTable(idx,:);
                this.QueuedX = transformPoints(this, XTable);
            elseif isempty(Opts.InitialX)
                % No initial points passed. Choose initial points.
                if isequal(Opts.AcquisitionFunctionName, 'grid')
                    % Random grid points
                    [X, this] = gridXFeasiblePoint(this);
                    XPoints = X;
                    while ~isempty(X) && size(XPoints,1) < Opts.NumSeedPoints
                        [X, this] = gridXFeasiblePoint(this);
                        if ~isempty(X)
                            XPoints(end+1, :) = X;
                        end
                    end
                else
                    % Random initial points
                    XPoints = initialXFeasiblePoints(this, this.PrivOptions.NumSeedPoints);
                end
                this.QueuedX = XPoints;
            elseif isempty(Opts.InitialObjective)
                % Initial X passed only. Recall that it is a table.
                this.QueuedX = transformPoints(this, Opts.InitialX);
            else
                % Initial evaluated points passed. use them directly.
                XTbl = canonicalizePoints(this, Opts.InitialX);
                this.XTrain = transformPoints(this, XTbl);
                this.FTrain = Opts.InitialObjective(:)';
                this.ObjectiveEvaluationTimeTrain = Opts.InitialObjectiveEvaluationTimes(:)';
                this.PrivIndexOfMinimumTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivFMinTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivFMinEstimatedTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivIterationTimeTrace = Opts.InitialIterationTimes;
                this = fillTotalElapsedTimeTrace(this);
                this.PrivUserDataTrace = Opts.InitialUserData;
                this.ConstraintTrain = Opts.InitialConstraintViolations;
                this.ErrorTrain = Opts.InitialErrorValues;
            end
        end
        
        function this = fillTotalElapsedTimeTrace(this)
            if ~isempty(this.PrivIterationTimeTrace)
                this.TotalElapsedTimeTrace = nancumsum(this.PrivIterationTimeTrace);
            end
        end
        
        function XBest = initialXFeasiblePoints(this, N)
            % Use random search to find N X-feasible points within bounds
            % for which the smallest Euclidean distance between points is
            % maximized.
            LB = this.PrivOptions.VarSpec.LBTrans;
            UB = this.PrivOptions.VarSpec.UBTrans;
            % Random search
            XBest = randomXFeasiblePoints(this, N);
            bsf = minDist(XBest);
            reps = 50;
            for r = 1:reps
                X = randomXFeasiblePoints(this, N);
                if minDist(X) > bsf
                    bsf = minDist(X);
                    XBest = X;
                end
            end
            function d = minDist(X)
                X = (X-LB)./(UB-LB);
                PointDists = pdist(X);
                d = min(PointDists(:));
            end
        end
        
        function this = chooseNextPoint(this)
            % Set this.Next to the next point to evaluate
            if this.NumObjectiveEvaluations < size(this.QueuedX,1)
                this.XNext = this.QueuedX(this.NumObjectiveEvaluations + 1, :);
            else
                [X, ChoseRandom, this] = chooseNextPointInternal(this);
                % Perform overexploitation loop. Overexploitation is impossible
                % if we have no IncumbentF, because we haven't found any
                % feasible points yet.
                if ~isnan(this.IncumbentF) && ~ChoseRandom && ...
                        ismember(this.PrivOptions.AcquisitionFunctionName,{'expectedimprovementplus','expectedimprovementpersecondplus'}) && ...
                        exploitingTooMuch(this, X)
                    % Save ObjectiveFcnGP
                    ObjFcnGP = this.ObjectiveFcnModel;
                    iteration = 1;
                    BullAdjustmentN = this.NumObjectiveEvaluations;
                    while exploitingTooMuch(this, X) && iteration <= this.PrivOptions.MaxExploitIterations
                        if this.PrivOptions.Verbose >= 2
                            fprintf('%s\n', overexploitDisplay(this, X));
                        end
                        % Refit GPs
                        this = fitBayesOptModels(this, true, BullAdjustmentN, false);
                        
                        % Choose the next point (maximize Acquisition Function)
                        [X, ~, this] = chooseNextPointInternal(this);
                        % Update loop
                        iteration = iteration + 1;
                        BullAdjustmentN = BullAdjustmentN*10;
                    end
                    % Restore ObjectiveFcnModel
                    this.ObjectiveFcnModel = ObjFcnGP;
                end
                this.XNext = X;
            end
        end
        
        function this = chooseNextPointParallel(this)
            % Choose the next point to evaluate, making use of predictions
            % of pending points.
            if this.NumObjectiveEvaluations < size(this.QueuedX,1)
                % There are queued initial points
                this.XNext = this.QueuedX(this.NumObjectiveEvaluations + 1, :);
                this.FantasyModels = this;
                this.FantasyModels.FantasyModels = [];
                % Impute pending objective for plotting purposes only.
                if anyPendingX(this.Parallel)
                    PendingXTable = conditionalizeX(this, this.Parallel.PendingX);
                    this.PendingObj = imputeObjective(this, PendingXTable, 'modelprediction');
                end
            elseif acquisitionFunctionIsGridOrRand(this)
                % Either grid or random. Choose point normally.
                this = chooseNextPoint(this);
                this.FantasyModels = this;
                this.FantasyModels.FantasyModels = [];
                % Impute pending objective for plotting purposes only.
                if anyPendingX(this.Parallel)
                    PendingXTable = conditionalizeX(this, this.Parallel.PendingX);
                    this.PendingObj = imputeObjective(this, PendingXTable, 'modelprediction');
                end
            elseif ~anyPendingX(this.Parallel)
                % No pending X. Choose point normally.
                this = chooseNextPoint(this);
                this.FantasyModels = this;
                this.FantasyModels.FantasyModels = [];
            else
                % We have pending points.
                % Predict objective, runtime, and constraints of pending points
                PendingXTable           = conditionalizeX(this, this.Parallel.PendingX);
                PendingObjective        = imputeObjective(this, PendingXTable, this.PrivOptions.ParallelMethod);
                PendingMeanRuntime      = predictObjectiveEvaluationTime(this, PendingXTable);
                PendingMeanConstraints  = predictConstraints(this, PendingXTable);
                PendingMeanError        = predictError(this, PendingXTable);
                % Make a copy of the BO object to do fantasy prediction
                BO2 = this;
                % Append predictions to training sets.
                BO2.XTrain = [BO2.XTrain; legalizePoints(BO2, BO2.Parallel.PendingX)];
                BO2.FTrain = [BO2.FTrain, PendingObjective'];
                BO2.ErrorTrain = [BO2.ErrorTrain; PendingMeanError];
                if BO2.PrivOptions.NumCoupledConstraints > 0
                    BO2.ConstraintTrain = [BO2.ConstraintTrain; PendingMeanConstraints];
                end
                BO2.ObjectiveEvaluationTimeTrain = [BO2.ObjectiveEvaluationTimeTrain, PendingMeanRuntime'];
                % Refit models to all points including fantasy points
                BO2 = fitModels(BO2, true);
                % Choose next point using the fantasy models
                BO2         = chooseNextPoint(BO2);
                this.XNext  = BO2.XNext;
                % Save fantasy models for Acquisition Function plotting.
                % Save PendingObj for Objective model plotting.
                this.FantasyModels = BO2;
                this.FantasyModels.FantasyModels = [];
                this.PendingObj    = PendingObjective;
            end
        end
        
        function ImputedObj = imputeObjective(this, PendingXTable, ParallelMethod)
            % Use the selected method to impute Objective values for the pending X points.
            switch ParallelMethod
                case 'clippedmodelprediction'
                    % Worst case of: Model mean at X, Model feasible min, Observed feasible min.
                    ImputedObj = max(max(predictObjective(this, PendingXTable), this.IncumbentF), this.MinObjective);
                case 'modelprediction'
                    % model predicted mean at X
                    ImputedObj = predictObjective(this, PendingXTable);
                case 'minobserved'
                    % min feasible observed
                    ImputedObj = repmat(this.MinObjective, height(PendingXTable), 1);
                case 'maxobserved'
                    % max feasible observed
                    FeasMask = double(this.FeasibilityTrace);
                    FeasMask(FeasMask==0) = NaN;
                    if isempty(this.FTrain)
                        ImputedObj = NaN(height(PendingXTable), 1);
                    else
                        ImputedObj = repmat(nanmax(this.FTrain(:).*FeasMask), height(PendingXTable), 1);
                    end
                otherwise
                    assert(ismember(ParallelMethod, {'clippedmodelprediction', ...
                        'modelprediction', 'minobserved', 'maxobserved'}));
            end
        end
        
        function [tf, ReasonKey] = shouldChooseRandomPoint(this)
            tf = false;
            ReasonKey = [];
            if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                return;
            elseif isequal(this.PrivOptions.AcquisitionFunctionName, 'random')
                % AF is 'random'
                tf = true;
                ReasonKey = 'RandomPointRandAcq';
            elseif sum(this.ErrorTrain ~= 1) < this.PrivOptions.NumSeedPoints
                % Not enough non-error points
                tf = true;
                ReasonKey = 'RandomPointInsuffSeedPoints';
            elseif isempty(this.ObjectiveFcnModel)
                % No Objective model
                tf = true;
                ReasonKey = 'RandomPointNoObjModel';
            else
                if shouldChooseRandomPoint(this.ObjectiveFcnModel, this.PrivOptions.SigmaFTol)
                    tf = true;
                    ReasonKey = 'RandomPointFlatObjModel';
                end
            end
        end
        
        function [XBest, ChooseRandom, this] = chooseNextPointInternal(this)
            [ChooseRandom, ReasonKey] = shouldChooseRandomPoint(this);
            if ChooseRandom
                if this.PrivOptions.Verbose >= 2
                    bayesoptim_priv.printInfo('ChoosingRandomPoint');
                    bayesoptim_priv.printInfo(ReasonKey);
                end
                XBest = randomXFeasiblePoints(this, 1);
            elseif isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                [XBest, this] = gridXFeasiblePoint(this);
            else
                % Choose the next point by maximizing the
                % constraint-weighted acquisition function.
                switch this.PrivOptions.AcquisitionFunctionName
                    case 'probabilityofimprovement'
                        AFcn = @(X)bayesoptim_priv.probabilityOfImprovement(X, this.ObjectiveFcnModel, ...
                            this.IncumbentF, this.ObjectiveFcnModel.Sigma);
                    case {'expectedimprovement', 'expectedimprovementplus'}
                        AFcn = @(X)bayesoptim_priv.expectedImprovement(X, this.ObjectiveFcnModel, ...
                            this.IncumbentF);
                    case 'lowerconfidencebound'
                        Kappa = this.PrivOptions.ExplorationRatio;
                        AFcn = @(X)bayesoptim_priv.lowerConfidenceBound(X, this.ObjectiveFcnModel, Kappa);
                    case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                        AFcn = @(X)bayesoptim_priv.expectedImprovementPerSecond(X, this.ObjectiveFcnModel, ...
                            this.ObjectiveEvaluationTimeModel, this.IncumbentF);
                    otherwise
                        bayesoptim_priv.err('AFUnknown', this.PrivOptions.AcquisitionFunctionName);
                end
                VarSpec = this.PrivOptions.VarSpec;
                XBest = iFminbndGlobal(@constraintWeightedNegAF, VarSpec.LBTrans, ...
                    VarSpec.UBTrans, this.PrivOptions.NumRestartCandidates, ...
                    this.PrivOptions.NumRestarts, this.PrivOptions.VerboseRestarts, ...
                    this.PrivOptions.MaxIterPerRestart, this.PrivOptions.RelTol);
            end
            if ~isempty(XBest)
                XBest = legalizePoints(this, XBest);
            end
            
            function NegAF = constraintWeightedNegAF(X)
                % Return the constraint-weighted negative AF value for all
                % feasible points in X (to be minimized). NegAF(i) = Inf if
                % not feasible
                NegAF = inf(size(X,1),1);
                [Xcanon, InputFeasible] = legalizePoints(this, X);
                if any(InputFeasible)
                    Xfeasible = Xcanon(InputFeasible,:);
                    if isnan(this.IncumbentF)
                        % We have constraints, but no incumbent, meaning that
                        % no feasible points were found. Maximize only the
                        % constraint probability.
                        NegAF(InputFeasible) = -ProbAllConstraintsSatisfied(this, Xfeasible);
                    else
                        % Constraints and an incumbent. Maximize constraint-weighted AF.
                        NegAF(InputFeasible) = -AFcn(Xfeasible) .*  ProbAllConstraintsSatisfied(this, Xfeasible);
                    end
                end
            end
        end
        
        function [XBest, this] = sampleGridWithoutReplacement(this)
            % Return XBest and update this.SampledGridIndices. Return [] if
            % the grid has already been fully sampled.
            Divs = this.PrivOptions.NumGridDivisions;
            if size(this.SampledGridIndices, 1) == prod(Divs)
                if this.PrivOptions.Verbose >= 2
                    bayesoptim_priv.printInfo('GridFullySampled');
                end
                this.SampledGridIndices = [];
            end
            Indices = iFindUnusedGridIndices(Divs, this.SampledGridIndices);
            XBest = XFromGridIndices(this, Indices);
            this.SampledGridIndices(end+1, :) = Indices;
        end
        
        function this = pruneGridIndices(this)
            if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid') && ~isempty(this.XNext)
                this.SampledGridIndices(end, :) = [];
            end
        end
        
        function X = XFromGridIndices(this, Indices)
            % Grid points include both endpoints. Indices are origin 1.
            GridRes = this.PrivOptions.NumGridDivisions;
            LB = this.PrivOptions.VarSpec.LBTrans;
            UB = this.PrivOptions.VarSpec.UBTrans;
            X = (Indices-1)./(GridRes-1) .* (UB-LB) + LB;
        end
        
        function TF = exploitingTooMuch(this, X)
            % X is a design matrix. TF is a logical vector. TF(i)=true if
            % the function std is less than alpha times the noise std at
            % the point X(i,:).
            %
            % "You can't overexploit an infeasible point": If X is being
            % evaluated for overexploitation, then it has been chosen by
            % the AF as the most desirable point. If this point is not
            % known to be feasible, then there is something to be gained by
            % evaluating it. We need to know whether it is feasible or not.
            if isempty(this.ObjectiveFcnModel)
                TF = false(size(X,1),1);
            else
                % Check ObjectiveFcnModel condition
                [~, YSD] = predict(this.ObjectiveFcnModel, X);
                FSD = funStd(YSD, this.ObjectiveFcnModel.Sigma);
                
                if strcmp(this.PrivOptions.ModelType,'GaussianProcess') || strcmp(this.PrivOptions.ModelType,'MultiGaussianProcess')
                    TF = FSD < this.PrivOptions.ExplorationRatio*this.ObjectiveFcnModel.Sigma;
                    % Check feasibility:
                    TF = TF && allConstraintsSatisfied(this, X, this.ConstraintModels, this.ErrorModel);
                else
                    TF = false(size(X,1),1);
                end
                 
            end
        end
        
        function tf = optimizationFinished(this, iteration)
            if iteration == 1
                tf = false;
            elseif isempty(this.XNext)
                tf = true;
            elseif this.PlotFcnStop || this.OutputFcnStop || ...
                    this.NumObjectiveEvaluations >= this.PrivOptions.MaxObjectiveEvaluations || ...
                    this.TotalElapsedTime >= this.PrivOptions.MaxTime
                tf = true;
            else
                tf = false;
            end
        end
        
        %% Function evaluation
        function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                = callObjNormally(this, X)
            % Call the objective fcn with the correct nargout, catching
            % errors and timing runtime
            StartTic = tic;
            ErrorConstraintViolation = -1;
            try
                switch this.PrivOptions.ObjectiveNargout
                    case 1
                        Objective = this.ObjectiveFcn(conditionalizeX(this, X));
                        ConstraintViolations = [];
                        UserData = [];
                    case 2
                        [Objective, ConstraintViolations] = this.ObjectiveFcn(conditionalizeX(this, X));
                        UserData = [];
                    otherwise
                        [Objective, ConstraintViolations, UserData] = this.ObjectiveFcn(conditionalizeX(this, X));
                end
            catch msg
                if isequal(msg.identifier, 'MATLAB:unassignedOutputs') && this.PrivOptions.NumCoupledConstraints > 0
                    bayesoptim_priv.err('ObjectiveNargoutWrong', this.PrivOptions.NumCoupledConstraints);
                else
                    rethrow(msg);
                end
            end
            [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                = finishObjEval(this, Objective, ConstraintViolations, ErrorConstraintViolation, StartTic);
        end
        
        function [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                = finishObjEval(this, Objective, ConstraintViolations, ErrorConstraintViolation, StartTic)
            % Check constraint violation dimensions
            if numel(ConstraintViolations) ~= this.PrivOptions.NumCoupledConstraints
                bayesoptim_priv.err('ObjectiveNumConstraintsWrong',...
                    numel(ConstraintViolations), this.PrivOptions.NumCoupledConstraints);
            end
            % Set illegal Objective and ConstraintViolation to NaN
            Objective = iNanIfBad(Objective);
            ConstraintViolations = arrayfun(@iNanIfBad, ConstraintViolations);
            if isnan(Objective)
                ErrorConstraintViolation = 1;
            end
            % Record runtime
            ObjectiveEvaluationTime = toc(StartTic);
        end
        
        %% Modeling
        function this = fitModels(this, useFitMethodNone)
            if nargin<2
                useFitMethodNone = false;
            end
            switch this.PrivOptions.ModelType
                case 'GaussianProcess'
                    if isempty(this.ModelStrategy)
                        this.ModelStrategy = bayesoptim_priv.GPStrategy();
                    end
                    this = fitBayesOptModels(this, false, 1, useFitMethodNone);
                case 'MultiGaussianProcess'
                    if isempty(this.ModelStrategy)
                        this.ModelStrategy = bayesoptim_priv.MultiGPStrategy();
                    end
                    this = fitBayesOptModels(this, false, 1, useFitMethodNone);
                case 'SingleTreeBagger'
                    if isempty(this.ModelStrategy)
                        this.ModelStrategy = bayesoptim_priv.SingleTreeBaggerStrategy();
                    end
                    this = fitBayesOptModels(this, useFitMethodNone);
                case 'MultiTreeBagger'
                    if isempty(this.ModelStrategy)
                        this.ModelStrategy = bayesoptim_priv.MultiTreeBaggerStrategy();
                    end
                    this = fitBayesOptModels(this, useFitMethodNone);
                otherwise
                    assert(false);
            end
        end
        
        function this = fitBayesOptModels(this, varargin)
            % Fit models for objective functions, objective function
            % evaluation time, constraints and errors using either Random
            % Forests or Gaussian Processes depending on the strategy
            % object
            if this.PrivOptions.FitModels
                this.ObjectiveFcnModel = fitObjectiveFcnModel(this.ModelStrategy, ...
                    this.XTrain, this.FTrain, this.PrivOptions, this.ObjectiveFcnModel, varargin);
                this.ObjectiveEvaluationTimeModel = fitObjectiveEvaluationTimeModel(this.ModelStrategy, ...
                    this.XTrain, this.FTrain, this.PrivOptions, ...
                    this.ObjectiveEvaluationTimeTrain, this.ObjectiveEvaluationTimeModel, varargin);
                this.ConstraintModels = fitConstraintModels(this.ModelStrategy, ...
                    this.XTrain, this.ConstraintTrain, this.ConstraintModels, ...
                    this.PrivOptions, varargin);
                this.ErrorModel = fitErrorModel(this.ModelStrategy, this.XTrain, ...
                    this.FTrain, this.PrivOptions, this.ErrorTrain, this.ErrorModel, varargin);
            else
                this.ObjectiveFcnModel = [];
                this.ObjectiveEvaluationTimeModel = [];
                this.ConstraintModels = cell(1, this.PrivOptions.NumCoupledConstraints);
                this.ErrorModel = [];
            end
        end
        
        %% Proxy optimization
        function [IncumbentF, IncumbentX] = findIncumbent(this, BOOptions)
            % Find the point that minimizes the mean of the Objective
            % function model. Optionally accept modified Options.
            if nargin < 2
                BOOptions = this.PrivOptions;
            end
            if isempty(this.ObjectiveFcnModel) || acquisitionFunctionIsGridOrRand(this) || (this.PrivOptions.NumCoupledConstraints>0 && any(cellfun(@isempty, this.ConstraintModels)))
                IncumbentF = NaN;
                IncumbentX = [];
            else
                VarSpec = this.PrivOptions.VarSpec;
                IncumbentX = iFminbndGlobal(@(X)GPFMeanOnPoints(this,X), VarSpec.LBTrans, VarSpec.UBTrans, ...
                    BOOptions.NumRestartCandidates, BOOptions.NumRestarts, BOOptions.VerboseRestarts, ...
                    BOOptions.MaxIterPerRestart, BOOptions.RelTol);
                % XAtMin is not discretized (although it was discretized during optimization)
                IncumbentX = legalizePoints(this, IncumbentX);
                IncumbentF = GPFMeanOnPoints(this, IncumbentX);
                % If no finite incumbent was found, set it to NaN
                if ~isfinite(IncumbentF)
                    IncumbentF = NaN;
                    IncumbentX = [];
                end
            end
        end
        
        function Objective = GPFMeanOnPoints(this, X)
            % Return the mean ObjectiveFcnGP value for all feasible points
            % in X. Return Inf for infeasible points. Note that we must
            % transform to native space in order to apply the user-supplied
            % conditionalVariableFcn and constraint functions, and after
            % that must canonicalize X and transform back before we can
            % apply the GP models.
            [Xcanon, InputFeasible] = legalizePoints(this, X);
            allSat = allConstraintsSatisfied(this, Xcanon, this.ConstraintModels, this.ErrorModel);
            feasible = InputFeasible & allSat;
            Objective = inf(size(X,1),1);
            Objective(feasible) = predict(this.ObjectiveFcnModel, Xcanon(feasible,:));
        end
        
        function TF = isInputFeasible(this, X)
            % TF(i) is true iff X(i,:) satisfies the bounds and XConstraintFcn.
            TF = all(bsxfun(@ge, X, this.PrivOptions.VarSpec.LBTrans), 2) & ...
                all(bsxfun(@le, X, this.PrivOptions.VarSpec.UBTrans), 2);
            TF = TF & satisfiesXConstraint(this, X);
        end
        
        function InputFeasible = applyXConstraint(this, XTable)
            if isempty(this.PrivOptions.XConstraintFcn)
                InputFeasible = true(height(XTable),1);
            else
                try
                    InputFeasible = this.PrivOptions.XConstraintFcn(XTable);
                catch me
                    bayesoptim_priv.warn('XConstraintFcnError');
                    rethrow(me);
                end
                if ~(islogical(InputFeasible) && isequal(size(InputFeasible), [height(XTable),1]))
                    bayesoptim_priv.err('XConstraintFcnOutput', num2str(height(XTable)), class(InputFeasible), num2str(size(InputFeasible)));
                end
            end
        end
        
        function checkForOptimizableVariables(this)
            if this.PrivOptions.VarSpec.NumVars == 0
                bayesoptim_priv.err('NoOptimizableVariables');
            end
        end
        
        function checkXConstraintFcnSatisfiability(this)
            % Check satisfiability of the XConstraint
            N = this.PrivOptions.XConstraintSatisfiabilitySamples;
            X = iUniformPointsWithinBounds(N, this.PrivOptions.VarSpec.LBTrans, ...
                this.PrivOptions.VarSpec.UBTrans);
            NumFeas = sum(isInputFeasible(this, X));
            if NumFeas == 0
                bayesoptim_priv.err('NoFeasiblePoints', N);
            elseif NumFeas < .01*N
                bayesoptim_priv.warn('FewFeasiblePoints', N, NumFeas);
            end
        end
        
        function TF = allConstraintsSatisfied(this, Xcanon, ConstraintModels, ErrorModel)
            % Return true for every row of Xcanon that satisfies all of the
            % constraints to their required tolerances. Empty constraints
            % are never satisfied. XCanon is canonicalized: It is
            % well-formed and has NaNs filled with LBs.
            if numel(ConstraintModels) == 0
                Sat = true;
            else
                for i = numel(ConstraintModels):-1:1
                    if isempty(ConstraintModels{i})
                        Sat(1:size(Xcanon,1), i) = false;
                    else
                        [means, ysds] = predict(ConstraintModels{i}, Xcanon);
                        fsds = funStd(ysds, ConstraintModels{i}.Sigma);
                        Sat(:,i) = normcdf(0, means, fsds) >= 1-this.PrivOptions.CoupledConstraintTolerances(i);
                    end
                end
            end
            if isempty(ErrorModel)
                % Assume no error if no error GP
                ErrorConstraintSat(1:size(Xcanon,1), 1) = true;
            else
                [means, ysds] = predict(ErrorModel, Xcanon);
                fsds = funStd(ysds, ErrorModel.Sigma);
                ErrorConstraintSat = normcdf(0, means, fsds) >= 1-this.PrivOptions.ErrorConstraintTol;
            end
            TF = all(Sat, 2) & ErrorConstraintSat;
        end
        
        %% Output and plotting
        % Plotting helper methods
        function printProgress(this, iteration, NumFinishedIterations)
            if this.PrivOptions.Verbose >= 1
                % Compute Action
                [~,~,MinObsLoc] = bestPoint(this, 'Criterion', 'min-observed');
                if isnan(this.FTrain(iteration))
                    Action = 'Error ';
                elseif ~isempty(this.ConstraintsTrace) && any(this.ConstraintsTrace(iteration,:) > 0)
                    Action = 'Infeas';
                elseif MinObsLoc == iteration
                    Action = 'Best  ';
                else
                    Action = 'Accept';
                end
                XTable = this.XTrace(iteration,:);
                nvars = this.NumVars;
                Names = this.PrivOptions.VarSpec.Names;
                VarNameFieldFitCAuto = 0;
                tblNames = this.XTrace.Properties.VariableNames;
                if this.PrivOptions.FitcautoMultipleLearners
                    learnerName = char(XTable{end,1});
                    ind = find(startsWith(tblNames(2:end),learnerName))+1;
                elseif this.PrivOptions.FitcAutoSingleLearner
                    splits = split(XTable.Properties.VariableNames{1},'_');
                    if size(splits,1) > 1
                        learnerName = splits{1};
                        ind = find(startsWith(tblNames(1:end),learnerName));
                    end
                end
                
                isfitcauto = this.PrivOptions.FitcautoMultipleLearners || ...
                    this.PrivOptions.FitcAutoSingleLearner;

                if isfitcauto
                    % For AutoML, create a small table of this learner's
                    % parameter names and values
                    fitcnames = cell(numel(ind),1);
                    fitcvals  = cell(numel(ind),1);
                    for idx = 1:numel(ind)
                        i = ind(idx);
                        c = XTable{end,i};
                        cname = char(erase(tblNames(i),[learnerName,'_']));
                        fitcnames{idx} = cname;
                        fitcvals{idx} = c;
                    end
                    smalltable = table(fitcnames,fitcvals);
                    XTable = XTable(:,1);
                    XTable.Parameters = {smalltable};
                    XTable.Properties.VariableNames = {'Learner','Hyperparameter:       Value'};
                    if this.PrivOptions.FitcAutoSingleLearner
                        XTable.Learner = categorical({learnerName});
                    end
                    Names = XTable.Properties.VariableNames;
                    nvars = 2;
                    VarNameFieldFitCAuto = 15;
                end
                
                % Display header.
                IterEvalWidth = 7+9;
                ParallelWidth = 10;
                VarNameFieldWidth = 12;

                IncludeEstimCol = ~acquisitionFunctionIsGridOrRand(this);
                MandatoryFields = 3+IncludeEstimCol;
                if rem(iteration, this.PrivOptions.DisplayHeaderInterval) == 1 || ...
                        (isfitcauto && this.PrivOptions.Verbose >= 2) || ...
                        (isfitcauto && (rem(iteration,10) == 1) ) 
                    if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        ObjectiveFieldLine1 = ' Validation ';
                        ObjectiveFieldLine2 = ' loss       ';
                    elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        ObjectiveFieldLine1 = ' log(1 + valLoss) ';
                        ObjectiveFieldLine2 = '                  ';
                    elseif this.PrivOptions.IsClassregRegressionFunction
                        ObjectiveFieldLine1 = ' Objective:  ';
                        ObjectiveFieldLine2 = ' log(1+loss) ';
                    else
                        ObjectiveFieldLine1 = ' Objective   ';
                        ObjectiveFieldLine2 = '             ';
                    end
                    % Header top bar
                    if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        NumEquals = 129;
                    elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        NumEquals = 137;
                    else
                        NumEquals = IterEvalWidth + MandatoryFields*(14) - 1 ...
                            + (this.PrivOptions.NumCoupledConstraints + nvars)*(VarNameFieldWidth+3) ...
                            + VarNameFieldFitCAuto;
                    end
                    if this.PrivOptions.UseParallel
                        NumEquals = NumEquals + ParallelWidth;
                    end
                    Bar = [sprintf('|'), sprintf(repmat('=', 1, NumEquals)), sprintf('|\n')];
                    fprintf(Bar);
                    % Header line 1
                    if this.PrivOptions.UseParallel
                        if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                            fprintf('| Iter | Active  | Eval   |%s| Time for training | Observed min    |', ObjectiveFieldLine1);
                        elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                            fprintf('| Iter | Active  | Eval   |%s| Time for training | Observed min     |', ObjectiveFieldLine1);
                        else
                            fprintf('| Iter | Active  | Eval   |%s| Objective   | BestSoFar   |', ObjectiveFieldLine1);
                        end
                    else
                        if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                            fprintf('| Iter | Eval   |%s| Time for training | Observed min    |', ObjectiveFieldLine1);
                        elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                            fprintf('| Iter | Eval   |%s| Time for training | Observed min     |', ObjectiveFieldLine1);
                       else
                            fprintf('| Iter | Eval   |%s| Objective   | BestSoFar   |', ObjectiveFieldLine1);
                        end
                    end
                    if IncludeEstimCol && isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' Estimated min   |');
                    elseif IncludeEstimCol && isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' Estimated min    |');
                    elseif IncludeEstimCol && ~isfitcauto
                        fprintf(' BestSoFar   |');
                    end
                    for c = 1:this.PrivOptions.NumCoupledConstraints
                        fprintf(' Constraint%d  |', c);
                    end
                    % Variable names in header line 1
                    if isfitcauto
                        fprintf(' Learner      | Hyperparameter:       Value |');
                    else 
                        for v = 1:nvars
                            Name = Names{v};
                            width = VarNameFieldWidth;
                            if numel(Name) == VarNameFieldWidth+1
                                ControlString = sprintf(' %%%ds|', VarNameFieldWidth+1);
                                width = VarNameFieldWidth+1;
                            elseif numel(Name) > VarNameFieldWidth
                                ControlString = sprintf(' %%%ds-|', VarNameFieldWidth);
                            else
                                ControlString = sprintf(' %%%ds |', VarNameFieldWidth);
                            end
                            if isfitcauto && (v == nvars)
                                ControlString = ' %21s |';
                                fprintf(ControlString, Name);
                            else
                                fprintf(ControlString, Name(1:min(numel(Name), width)));
                            end
                        end
                    end
                    fprintf('\n');
                    % Header line 2
                    if this.PrivOptions.UseParallel
                        if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                            fprintf('|      | workers | result |%s| & validation (sec)| validation loss |', ObjectiveFieldLine2);
                        elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                            fprintf('|      | workers | result |%s| & validation (sec)| log(1 + valLoss) |', ObjectiveFieldLine2);
                        else
                            fprintf('|      | workers | result |%s| runtime     | (observed)  |', ObjectiveFieldLine2);
                        end
                    else
                        if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                            fprintf('|      | result |%s| & validation (sec)| validation loss |', ObjectiveFieldLine2);
                        elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                            fprintf('|      | result |%s| & validation (sec)| log(1 + valLoss) |', ObjectiveFieldLine2);
                        else
                            fprintf('|      | result |%s| runtime     | (observed)  |', ObjectiveFieldLine2);
                        end
                    end
                    if IncludeEstimCol && isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' validation loss |');
                    elseif IncludeEstimCol && isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' log(1 + valLoss) |');
                    elseif IncludeEstimCol && ~isfitcauto
                        fprintf(' (estim.)    |');
                    end
                    for c = 1:this.PrivOptions.NumCoupledConstraints
                        fprintf('              |');
                    end
                    % Variable names in header line 2
                    for v = 1:nvars
                        Name = Names{v};
                        if numel(Name) > VarNameFieldWidth+1
                            str = Name(VarNameFieldWidth+1 : min(numel(Name), 2*VarNameFieldWidth));
                        else
                            str = '';
                        end
                        width = VarNameFieldWidth;
                        ControlString = sprintf(' %%-%ds |', width);
                        if ~(isfitcauto && (v == nvars))
                            fprintf(ControlString, str);
                        else
                            fprintf('%30s','|');
                        end
                    end
                    fprintf('\n');
                    % Header bottom bar
                    fprintf(Bar);
                end
                
                % Display a line of data
                if this.PrivOptions.UseParallel
                    % Number of jobs that will be running after we launch a new one:
                    assert(NumFinishedIterations > 0);  % Because printProgress should only be called when jobs finish.
                    if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        fprintf('| %4d | %7d | %s | %10.5g | %17.5g | %15.5g |', ...
                            iteration, 1+numPendingX(this.Parallel), Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        fprintf('| %4d | %7d | %s | %16.5g | %17.5g | %16.5g |', ...
                            iteration, 1+numPendingX(this.Parallel), Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    else
                        fprintf('| %4d | %7d | %s | %11.5g | %11.5g | %11.5g |', ...
                            iteration, 1+numPendingX(this.Parallel), Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    end
                else
                    if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        fprintf('| %4d | %s | %10.5g | %17.5g | %15.5g |', ...
                            iteration, Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        fprintf('| %4d | %s | %16.5g | %17.5g | %16.5g |', ...
                            iteration, Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    else
                        fprintf('| %4d | %s | %11.5g | %11.5g | %11.5g |', ...
                            iteration, Action, this.FTrain(iteration), this.ObjectiveEvaluationTimeTrain(iteration), ...
                            this.MinObjective);
                    end
                end
                if IncludeEstimCol
                    if isfitcauto && ~this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' %15.5g |', this.MinEstimatedObjective);
                    elseif isfitcauto && this.PrivOptions.IsClassregRegressionFunction
                        fprintf(' %16.5g |', this.MinEstimatedObjective);
                    else
                        fprintf(' %11.5g |', this.MinEstimatedObjective);
                    end
                end
                for c = 1:this.PrivOptions.NumCoupledConstraints
                    fprintf(' %12.3g |', this.ConstraintTrain(iteration,c));
                end
                fprintf('%s', commandLineXDisplay(this, XTable));
                fprintf('\n');
            end
        end
        
        function N = numCharsBeforeX(this)
            % The number of characters in the command-line display before
            % the X section.
            N = 1+7+9+14+14+14;  % vertical bar, iter, eval result, objective. obj runtime, best so far (observed).
            if this.PrivOptions.UseParallel
                N = N+10;        % parallel jobs
            end
            IncludeEstimCol = ~acquisitionFunctionIsGridOrRand(this);
            if IncludeEstimCol
                N = N+14;       % best so far (estim)
            end
            N = N + this.PrivOptions.NumCoupledConstraints*15;  % Constraint fields
        end
        
        function Line = commandLineXDisplay(this, XTable)
            nvars = size(XTable,2);
            Line = '';
            for v = 1:nvars
                val = XTable.(v);
                if iscell(val) && isscalar(val) && istable(val{1})
                    % For AutoML, read this learner's parameter names and
                    % values out of a table. Display them all in the final
                    % cell of the displayed table, across multiple lines.
                    tbl = val{1};
                    nparams = size(tbl,1);
                    block = repmat(' ',nparams,27);
                    for j=1:nparams
                        % Left-justify parameter name
                        pname = [tbl{j,1}{1},':'];
                        plength = min(size(block,2),numel(pname));
                        block(j,1:plength) = pname(1:plength);
                        
                        % Right-justify parameter value
                        pval = tbl{j,2}{1};
                        if isnumeric(pval)
                            pval = strip(sprintf('%11.5g',pval));
                        else
                            pval = strip(char(pval));
                        end
                        vlength = min(size(block,2),numel(pval));
                        block(j,end-vlength+1:end) = pval(1:vlength);
                    end
                    n = numel(Line);
                    Line = [Line, sprintf(' %12s |', block(1,:))];
                    if size(block,1)>1
                        % Indent all lines after the first
                        ncb = numCharsBeforeX(this);
                        
                        if this.PrivOptions.UseParallel && ~this.PrivOptions.IsClassregRegressionFunction
                            fitcautoPadding = '|      |         |        |            |                   |                 |                 |              ';
                        elseif this.PrivOptions.UseParallel && this.PrivOptions.IsClassregRegressionFunction
                            fitcautoPadding = '|      |         |        |                  |                   |                  |                  |              ';
                        elseif ~this.PrivOptions.IsClassregRegressionFunction
                            fitcautoPadding = '|      |        |            |                   |                 |                 |              ';
                        else
                            fitcautoPadding = '|      |        |                  |                   |                  |                  |              ';
                        end

                        for idx = 2:size(block,1)
                            Line = [Line, newline, fitcautoPadding, sprintf('| %12s |', block(idx,:))];
                        end
                    end
                elseif iscategorical(val)
                    if ismissing(val)
                        Name = '-';
                    else
                        Name = char(val);
                    end
                    Line = [Line, sprintf(' %12s |', Name(1:min(numel(Name),12)))];
                else
                    if ismissing(val)
                        Name = '-';
                        Line = [Line, sprintf(' %12s |', Name(1:min(numel(Name),12)))];
                    else
                        Line = [Line, sprintf(' %12.5g |', val)];
                    end
                end
            end
        end
        
        function Line = verboseXNextDisplay(this)
            % Iterative display
            if this.PrivOptions.Verbose >= 2 && ~isempty(this.XNext)
                XTable = conditionalizeX(this, this.XNext);
                if ~(this.PrivOptions.FitcautoMultipleLearners || ...
                        this.PrivOptions.FitcAutoSingleLearner)
                    str = bayesoptim_priv.infoString('IterationMessage');
                    spaces = repmat(' ', 1, numCharsBeforeX(this) - numel(str));
                    Line = [str, spaces, commandLineXDisplay(this, XTable), '\n'];
                    fprintf(Line);
                else
                    [learnerName,tbl] = simplifyAutoResults(this,XTable);
                    displayFitcautoResults(tbl,learnerName,'IterationMessageFitcauto');
                end
            end
        end
        
        function printFitcautoVerboseSummary(this)
            fprintf([bayesoptim_priv.infoString('LearnerTypesFitcauto'),' ']);
            if this.PrivOptions.FitcAutoSingleLearner
                learnerNameSplit = split(this.VariableDescriptions(1).Name,'_');
                fprintf('%s ',learnerNameSplit{1});
            elseif this.PrivOptions.FitcautoMultipleLearners
                for i = 1:length(this.Options.VariableDescriptions(1).Range)
                    fprintf('%s',this.Options.VariableDescriptions(1).Range{i});
                    if i ~= length(this.Options.VariableDescriptions(1).Range)
                        fprintf(', ');
                    end
                end
            end
            fprintf('\n');
            fprintf([bayesoptim_priv.infoString('TotalIterationsFitcauto'),' (MaxObjectiveEvaluations): %d\n'],this.Options.MaxObjectiveEvaluations);
            fprintf([bayesoptim_priv.infoString('TotalTimeFitcauto'),' (MaxTime): %d\n'],this.Options.MaxTime);
            fprintf('\n');
        end
        
        function Line = overexploitDisplay(this, X)
            % Iterative display
            XTable = conditionalizeX(this, X);
            str = bayesoptim_priv.infoString('Overexploit');
            spaces = repmat(' ', 1, numCharsBeforeX(this) - numel(str));
            Line = [spaces, str, commandLineXDisplay(this, XTable)];
        end
        
        function this = callPlotFcn(this, State)
            PlotFcn = this.PrivOptions.PlotFcn;
            this.PlotFcnStop = false;
            if ~isempty(PlotFcn)
                if ~iscell(PlotFcn)
                    PlotFcn = {PlotFcn};
                end
                for i = 1:numel(PlotFcn)
                    try
                        stop = PlotFcn{i}(this, State);
                    catch me
                        bayesoptim_priv.warn('PlotFcnError', i);
                        bayesoptim_priv.printInfo('ErrorDisplay');
                        disp(me.message);
                        stop = false;
                    end
                    if ~(isscalar(stop) && islogical(stop))
                        bayesoptim_priv.warn('PlotFcnOutput', i);
                    end
                    this.PlotFcnStop = this.PlotFcnStop || stop;
                end
                drawnow;
            end
        end
        
        function this = callOutputFcn(this, State)
            OutputFcn = this.PrivOptions.OutputFcn;
            this.OutputFcnStop = false;
            if ~isempty(OutputFcn)
                if ~iscell(OutputFcn)
                    OutputFcn = {OutputFcn};
                end
                for i = 1:numel(OutputFcn)
                    try
                        stop = OutputFcn{i}(this, State);
                    catch me
                        bayesoptim_priv.warn('OutputFcnError', i);
                        bayesoptim_priv.printInfo('ErrorDisplay');
                        disp(me.message);
                        stop = false;
                    end
                    if ~(isscalar(stop) && islogical(stop))
                        bayesoptim_priv.warn('OutputFcnOutput', i);
                    end
                    this.OutputFcnStop = this.OutputFcnStop || stop;
                end
            end
        end
        
        function showStoppingReason(this)
            if this.PrivOptions.Verbose >= 1
                fprintf('\n__________________________________________________________\n');
                bayesoptim_priv.printInfo('OptimizationCompleted');
                if this.PrivOptions.FitcAutoSingleLearner || this.PrivOptions.FitcautoMultipleLearners
                    bayesoptim_priv.printInfo('TotalEvaluationsFitcauto', this.NumObjectiveEvaluations);
                    bayesoptim_priv.printInfo('TotalTime', num2str(this.TotalElapsedTime));
                    bayesoptim_priv.printInfo('TotalEvalTimeFitcauto', num2str(nansum(this.ObjectiveEvaluationTimeTrain)));
                else
                    if this.OutputFcnStop
                        bayesoptim_priv.printInfo('OutputFcnStop');
                    end
                    if this.PlotFcnStop
                        bayesoptim_priv.printInfo('PlotFcnStop');
                    end
                    if this.NumObjectiveEvaluations >= this.PrivOptions.MaxObjectiveEvaluations
                        bayesoptim_priv.printInfo('MaxEvalsReached', this.PrivOptions.MaxObjectiveEvaluations);
                    end
                    if this.TotalElapsedTime >= this.PrivOptions.MaxTime
                        bayesoptim_priv.printInfo('MaxTimeReached', num2str(this.PrivOptions.MaxTime));
                    end
                    if isempty(this.XNext) && isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                        bayesoptim_priv.printInfo('GridSearched');
                    end
                    bayesoptim_priv.printInfo('TotalEvaluations', this.NumObjectiveEvaluations);
                    bayesoptim_priv.printInfo('TotalTime', num2str(this.TotalElapsedTime));
                    bayesoptim_priv.printInfo('TotalEvalTime', num2str(nansum(this.ObjectiveEvaluationTimeTrain)));
                end
            end
        end
        
        function showBestPoints(this)
            if this.PrivOptions.Verbose >= 1
                [MinObsXTable, MinObsObjective, MinObsLoc] = bestPoint(this, 'Criterion', 'min-observed');
                if ~isfinite(MinObsLoc)
                    bayesoptim_priv.printInfo('NoFeasibleResult');
                else
                    % Show min observed feasible point
                    fprintf('\n');
                    if ~(this.PrivOptions.FitcautoMultipleLearners || ...
                        this.PrivOptions.FitcAutoSingleLearner)
                        bayesoptim_priv.printInfo('BestObserved');
                        disp(MinObsXTable);
                    else
                        [learnerName,fitcautoResults] = simplifyAutoResults(this,MinObsXTable);
                        displayFitcautoResults(fitcautoResults,learnerName,'BestObservedFitcauto')
                    end
                    
                    if (this.PrivOptions.FitcautoMultipleLearners || ...
                        this.PrivOptions.FitcAutoSingleLearner)
                        if this.PrivOptions.IsClassregRegressionFunction
                            bayesoptim_priv.printInfo('ObservedObjectiveValueFitrauto', num2str(MinObsObjective));
                        else
                            bayesoptim_priv.printInfo('ObservedObjectiveValueFitcauto', num2str(MinObsObjective));
                        end
                        bayesoptim_priv.printInfo('ObservedEvaluationTimeFitcauto', num2str(this.ObjectiveEvaluationTimeTrace(MinObsLoc)));
                    else
                        bayesoptim_priv.printInfo('ObservedObjectiveValue', num2str(MinObsObjective));
                        if ~acquisitionFunctionIsGridOrRand(this) && ~isempty(MinObsXTable)
                            bayesoptim_priv.printInfo('EstimatedObjectiveValue', num2str(predictObjective(this, MinObsXTable)));
                        end
                        bayesoptim_priv.printInfo('ObservedEvaluationTime', num2str(this.ObjectiveEvaluationTimeTrace(MinObsLoc)));
                    end
                    

                    if ~isempty(this.ConstraintsTrace)
                        bayesoptim_priv.printInfoNoReturn('ObservedConstraintViolations');
                        fprintf('[ ');
                        fprintf('%f ', this.ConstraintsTrace(MinObsLoc,:));
                        fprintf(']\n');
                    end
                    % Unless using grid or random, show best point
                    % according to the default criterion
                    if ~acquisitionFunctionIsGridOrRand(this)
                        if  this.PrivOptions.FitcautoMultipleLearners || ...
                                     this.PrivOptions.FitcAutoSingleLearner
                                BestPoint = bestPoint(this,'Criterion','min-visited-mean');
                        else
                                BestPoint = bestPoint(this);
                        end
                        if ~isempty(BestPoint)
                            fprintf('\n');
                            if  ~(this.PrivOptions.FitcautoMultipleLearners || ...
                                     this.PrivOptions.FitcAutoSingleLearner)
                                bayesoptim_priv.printInfo('BestEstimated');
                                disp(BestPoint);
                            else
                                [learnerName,fitcautoResults] = simplifyAutoResults(this,BestPoint);
                                displayFitcautoResults(fitcautoResults,learnerName,'BestEstimatedFitcauto')
                            end
                            
                            if (this.PrivOptions.FitcautoMultipleLearners || ...
                                this.PrivOptions.FitcAutoSingleLearner)
                                if this.PrivOptions.IsClassregRegressionFunction
                                    bayesoptim_priv.printInfo('EstimatedObjectiveValueFitrauto', num2str(predictObjective(this, BestPoint)));
                                else
                                    bayesoptim_priv.printInfo('EstimatedObjectiveValueFitcauto', num2str(predictObjective(this, BestPoint)));
                                end
                                bayesoptim_priv.printInfo('EstimatedEvaluationTimeFitcauto', num2str(predictObjectiveEvaluationTime(this, BestPoint)));
                            else
                                bayesoptim_priv.printInfo('EstimatedObjectiveValue', num2str(predictObjective(this, BestPoint)));
                                bayesoptim_priv.printInfo('EstimatedEvaluationTime', num2str(predictObjectiveEvaluationTime(this, BestPoint)));
                            end
                           
                            if this.PrivOptions.NumCoupledConstraints > 0
                                ConstraintPred = cellfun(@(GP)predict(GP, transformPoints(this, BestPoint)), ...
                                    this.ConstraintModels);
                                bayesoptim_priv.printInfoNoReturn('EstimatedConstraintViolations');
                                fprintf('[ ');
                                fprintf('%f ', ConstraintPred);
                                fprintf(']\n');
                            end
                        end
                    end
                    fprintf('\n');
                end
            end
        end
        
        function showVerboseDocHelp(this)
            if (this.PrivOptions.FitcAutoSingleLearner || this.PrivOptions.FitcautoMultipleLearners) && ...
                    this.PrivOptions.Verbose >= 1  
                if this.PrivOptions.IsClassregRegressionFunction
                    veboseLinkText = bayesoptim_priv.infoString('veboseLinkTextFitrauto');
                    fprintf(['<a href="matlab:helpview(fullfile(docroot,''stats'',''stats.map''), ''fitrautoVerboseDoc'')">',veboseLinkText,'</a>\n']);
                else
                    veboseLinkText = bayesoptim_priv.infoString('veboseLinkTextFitcauto');
                    fprintf(['<a href="matlab:helpview(fullfile(docroot,''stats'',''stats.map''), ''fitcautoVerboseDoc'')">',veboseLinkText,'</a>\n']);
                end
            end
        end
        
        function N = NumVars(this)
            N = this.PrivOptions.VarSpec.NumVars;
        end
        
        function [Xcanon, InputFeasible] = legalizePoints(this, X)
            % Given points X in real, transformed space, return the corresponding
            % conditional-mapped, canonicalized points in the same space.
            XTable = untransformPoints(this, X, true);
            XTable = applyConditionalVariableFcn(this, XTable);
            InputFeasible = applyXConstraint(this, XTable);
            XTable = canonicalizePoints(this, XTable);
            Xcanon = transformPoints(this, XTable);
        end
        
        function XTrain = XTraceToXTrain(this, XTable)
            % Called in resume()
            XTable = applyConditionalVariableFcn(this, XTable);
            XTable = canonicalizePoints(this, XTable);
            XTrain = transformPoints(this, XTable);
        end
        
        function X = makeXTablePlottable(this, XTable)
            % Returns a matrix X that is a plottable version of XTable.
            % NaNs are canonicalized, Categoricals are int-coded.
            XTable = canonicalizePoints(this, XTable);
            for v = 1:width(XTable)
                if iscategorical(XTable.(v))
                    XTable.(v) = intCodeCategorical(XTable.(v));
                end
            end
            X = table2array(XTable);
        end
        
        function Prob = ProbAllConstraintsSatisfied(this, X)
            % The probability that ALL the constraints are satisfied at
            % each row of X. Constraint without models are never satisfied.
            N = size(X,1);
            % Coupled constraints
            if this.PrivOptions.NumCoupledConstraints == 0
                % There are no constraints at all.
                ProbSat = 1;
            else
                for i = numel(this.ConstraintModels):-1:1
                    if isempty(this.ConstraintModels{i})
                        % There is a constraint but it has no model yet.
                        ProbSat(1:N,i) = 0;
                    else
                        [means, ysds] = predict(this.ConstraintModels{i}, X);
                        fsds = funStd(ysds, this.ConstraintModels{i}.Sigma);
                        ProbSat(:,i) = normcdf(0, means, fsds);
                    end
                end
            end
            % Error constraint
            if isempty(this.ErrorModel)
                % There is no model of error yet.
                ProbErrorSat(1:N,1) = 1;
            else
                [means, ysds] = predict(this.ErrorModel, X);
                fsds = funStd(ysds, this.ErrorModel.Sigma);
                ProbErrorSat(:,1) = normcdf(0, means, fsds);
            end
            % XConstraint
            IsXConstraintSat = satisfiesXConstraint(this, X);
            
            % Calculate prob
            Prob = prod(ProbSat, 2) .* ProbErrorSat .* IsXConstraintSat;
        end
        
        function TF = satisfiesXConstraint(this, X)
            % Returns true for each row for which xConstraint is
            % satisfied. An empty Xconstraint is always satisfied.
            if isempty(this.PrivOptions.XConstraintFcn)
                TF = true(size(X,1),1);
            else
                XTable = untransformPoints(this, X, true);
                XTable = applyConditionalVariableFcn(this, XTable);
                TF = applyXConstraint(this, XTable);
            end
        end
        
        %% Choosing the best point (used by the bestPoint method)
        function [BestXTable, MinObserved, Iteration] = minObservedPoint(this)
            % Returns the feasible point BestXTable with the minimum
            % observed Objective value. Feasibility is judged by the
            % current constraint and error models.
            FeasMask = double(this.FeasibilityTrace);
            FeasMask(FeasMask==0) = NaN;
            [MinObserved, Iteration] = nanmin(this.ObjectiveTrace .* FeasMask);
            if isfinite(MinObserved)
                BestXTable = this.XTrace(Iteration,:);
            else
                BestXTable = [];
                MinObserved = NaN;
                Iteration = NaN;
            end
        end
        
        function [BestXTable, MinMean] = minMeanPoint(this)
            % Attempts to find the feasible point BestXTable that minimizes
            % the mean of the ObjectiveFcn model.
            if isempty(this.ObjectiveFcnModel) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintModels)))
                BestXTable = [];
                MinMean = NaN;
            else
                % Set tighter optimization stopping criteria
                BOOptions = this.PrivOptions;
                BOOptions.MaxIterPerRestart = 40;
                BOOptions.RelTol = 1e-7;
                [MinMean, BestX] = findIncumbent(this, BOOptions);
                BestXTable = untransformPoints(this, BestX, true);
            end
            % If none found, revert to the minVisitedMeanPoint
            if isempty(BestXTable)
                [BestXTable, MinMean] = minVisitedMeanPoint(this);
            end
        end
        
        function [BestXTable, MinUCI] = minUCIPoint(this, Alpha)
            % Attempts to find the feasible point BestXTable that minimizes
            % the 100(1-Alpha)% upper confidence interval (UCI) of the
            % ObjectiveFcn model.
            if isempty(this.ObjectiveFcnModel) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintModels)))
                MinUCI = NaN;
                BestXTable = [];
            else
                BOOptions = this.PrivOptions;
                BOOptions.MaxIterPerRestart = 40;
                BOOptions.RelTol = 1e-7;
                VarSpec = BOOptions.VarSpec;
                BestX = iFminbndGlobal(@(X)GPF_UCI_OnPoints(this, X, Alpha), VarSpec.LBTrans, VarSpec.UBTrans, ...
                    BOOptions.NumRestartCandidates, BOOptions.NumRestarts, BOOptions.VerboseRestarts, ...
                    BOOptions.MaxIterPerRestart, BOOptions.RelTol);
                % XAtMin is not discretized (although it was discretized during optimization)
                BestX = legalizePoints(this, BestX);
                MinUCI = GPF_UCI_OnPoints(this, BestX, Alpha);
                BestXTable = untransformPoints(this, BestX, true);
            end
            % If none found, revert to the minVisitedUCIPoint
            if isempty(BestXTable)
                [BestXTable, MinUCI] = minVisitedUCIPoint(this, Alpha);
            end
        end
        
        function [BestXTable, MinMean, Iteration] = minVisitedMeanPoint(this)
            % Attempts to find the feasible point BestXTable that minimizes
            % the mean of the ObjectiveFcn model, from among those points
            % visited (evaluated).
            if isempty(this.ObjectiveFcnModel) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintModels)))
                BestXTable = [];
                MinMean = NaN;
                Iteration = NaN;
            else
                ObjMeans = GPFMeanOnPoints(this, this.XTrain);
                [MinMean, Iteration] = nanmin(ObjMeans);
                BestXTable = this.XTrace(Iteration,:);
            end
        end
        
        function [BestXTable, MinUCI, Iteration] = minVisitedUCIPoint(this, Alpha)
            % Attempts to find the feasible point BestXTable that minimizes
            % the 100(1-Alpha)% upper confidence interval (UCI) of the
            % ObjectiveFcn model, from among those points visited
            % (evaluated).
            UCI = GPF_UCI_OnPoints(this, this.XTrain, Alpha);
            [MinUCI, Iteration] = nanmin(UCI);
            BestXTable = this.XTrace(Iteration,:);
        end
        
        function UCI = GPF_UCI_OnPoints(this, X, Alpha)
            % Return the Alpha% UCI of the ObjectiveFcnGP value for all
            % feasible rows of X. Return Inf for infeasible rows.
            UCI = inf(size(X,1),1);
            if isempty(this.ObjectiveFcnModel) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintModels)))
                return;
            else
                [Xcanon, InputFeasible] = legalizePoints(this, X);
                allSat = allConstraintsSatisfied(this, Xcanon, this.ConstraintModels, this.ErrorModel);
                feasible = InputFeasible & allSat;
                [FMean, YSD] = predict(this.ObjectiveFcnModel, Xcanon(feasible,:));
                FSD = funStd(YSD, this.ObjectiveFcnModel.Sigma);
                UCI(feasible) = norminv(1-Alpha, FMean, FSD);
            end
        end
        
        %% Misc
        function [X, this] = gridXFeasiblePoint(this)
            % A random grid point that is XFeasible, or [] if none exist.
            [X, this] = sampleGridWithoutReplacement(this);
            while ~isempty(X) && ~isInputFeasible(this, X)
                [X, this] = sampleGridWithoutReplacement(this);
            end
        end
        
        function X = randomXFeasiblePoints(this, N)
            % Generate N random initial points X that satisfy the
            % xConstraint. Return them in transformed space.
            X = iUniformPointsWithinBounds(this.PrivOptions.NumRestartCandidates, ...
                this.PrivOptions.VarSpec.LBTrans, this.PrivOptions.VarSpec.UBTrans);
            % Minimize the number of calls to isInputFeasible
            Idxs = [];
            for i = 1:size(X,1)
                if isInputFeasible(this, X(i,:))
                    Idxs(end+1) = i;
                    if numel(Idxs) == N
                        break;
                    end
                end
            end
            X = X(Idxs,:);
        end
        
          function X = randomXFeasiblePointsForEachLearner(this, N,learnerIdx)
            % Generate N random initial points X per learner that satisfy the
            % xConstraint. Return them in transformed space.
            X = [];
            learnerLB = this.PrivOptions.VarSpec.LBTrans(learnerIdx);
            learnerUB = this.PrivOptions.VarSpec.UBTrans(learnerIdx);
            learnerLBTrans = this.PrivOptions.VarSpec.LBTrans;
            learnerUBTrans = this.PrivOptions.VarSpec.UBTrans;
            for j = learnerLB:learnerUB
                learnerLBTrans(learnerIdx) = j;
                learnerUBTrans(learnerIdx) = j;
                XperLearner = iUniformPointsWithinBounds(this.PrivOptions.NumRestartCandidates, ...
                    learnerLBTrans, learnerUBTrans);
                % Minimize the number of calls to isInputFeasible
                Idxs = [];
                for i = 1:size(XperLearner,1)
                    if isInputFeasible(this, XperLearner(i,:))
                        Idxs(end+1) = i;
                        if numel(Idxs) == N
                            break;
                        end
                    end
                end
                X = [X;XperLearner(Idxs,:)];
            end
        end
                
        function XTable = conditionalizeX(this, X)
            XTable = untransformPoints(this, X, true);
            XTable = applyConditionalVariableFcn(this, XTable);
        end
        
        function NewXTable = applyConditionalVariableFcn(this, XTable)
            NewXTable = XTable;
            if ~isempty(this.PrivOptions.ConditionalVariableFcn)
                try
                    NewXTable = this.PrivOptions.ConditionalVariableFcn(XTable);
                catch me
                    bayesoptim_priv.warn('ConditionalVariableFcnError');
                    rethrow(me);
                end
                % Check output
                if ~(istable(NewXTable) && ...
                        isequal(size(NewXTable), size(XTable)) && ...
                        isequal(varfun(@class, NewXTable, 'OutputFormat','cell'), varfun(@class, XTable, 'OutputFormat','cell')) && ...
                        isequal(NewXTable.Properties.VariableNames, XTable.Properties.VariableNames))
                    bayesoptim_priv.err('ConditionalVariableFcnOutput');
                end
            end
        end
        
        function XTable = canonicalizePoints(this, XTable)
            % "Canonicalization" is how we handle parameters that are NaN
            % due to the action of the ConditionalVariableFcn. To give
            % these points a "canonical" location in X-space, we map
            % numeric missing values to the minimum of their range, and we
            % map categoricals to their first category.
            
            % Only canonicalize points when the model type is GP since RF
            % can handle NaN and undefined values
            if strcmp(this.PrivOptions.ModelType, 'GaussianProcess') ||...
                    strcmp(this.PrivOptions.ModelType, 'MultiGaussianProcess')
                VarSpec = this.PrivOptions.VarSpec;
                missing = ismissing(XTable);
                for v = 1:width(XTable)
                    if any(missing(:,v))
                        if VarSpec.isCat(v)
                            cats = VarSpec.Categories{v};
                            val = cats(1);
                        else
                            val = VarSpec.LBs(v);
                        end
                        XTable{missing(:,v), v} = val;
                    end
                end
            end
        end
        
        function X = transformPoints(this, XTable)
            % Convert a table in native space to a matrix of transformed points.
            VarSpec = this.PrivOptions.VarSpec;
            IsNumericLinear = ~VarSpec.isCat & ~VarSpec.isLog;
            IsNumericLog = ~VarSpec.isCat & VarSpec.isLog;
            if any(IsNumericLinear)
                X(:,IsNumericLinear) = XTable{:,IsNumericLinear};
            end
            if any(IsNumericLog)
                X(:,IsNumericLog) = log(XTable{:,IsNumericLog});
            end
            if any(VarSpec.isCat)
                X(:,VarSpec.isCat) = table2array(varfun(@intCodeCategorical, XTable(:,VarSpec.isCat)));
            end
        end
        
        function XTable = untransformPoints(this, X, HonorInts)
            % Convert a matrix of transformed points to a table in native
            % space. HonorInts=true means that the transformation will
            % round reals to ints (and then to cats if necessary). false
            % means the points are left as reals after transforming to
            % native space. 'false' is only used for plotting meshes.
            if isempty(X)
                XTable = table;
            else
                VarSpec = this.PrivOptions.VarSpec;
                Vars = struct;
                for v = 1:this.NumVars
                    Name = VarSpec.Names{v};
                    switch VarSpec.Types{v}
                        case 'real'
                            if isequal(VarSpec.Transforms{v}, 'log')
                                Vars.(Name) = exp(X(:,v));
                            else
                                Vars.(Name) = X(:,v);
                            end
                        case 'integer'
                            if isequal(VarSpec.Transforms{v}, 'log')
                                Vars.(Name) = exp(X(:,v));
                            else
                                Vars.(Name) = X(:,v);
                            end
                            if HonorInts
                                % Convert real to int
                                Vars.(Name) = round(Vars.(Name));
                            end
                        case 'categorical'
                            Vars.(Name) = X(:,v);
                            if HonorInts
                                % Convert real to int to cat
                                idx              = round(Vars.(Name));
                                nanIdx           = isnan(idx);
                                idx(nanIdx)      = 1;                       % Brief dummy value
                                catNames         = VarSpec.Categories{v}(idx);
                                catNames(nanIdx) = '<undefined>';
                                Vars.(Name)      = catNames;
                            end
                    end
                end
                XTable = struct2table(Vars);
            end
        end
        
        function T = checkAndPrepareTableForPrediction(this, XTable, caller)
            % Check input, then find the required variables by Name and put
            % them in the order matching VarSpec.
            if isempty(XTable)
                T = table;
                return;
            elseif ~istable(XTable)
                bayesoptim_priv.err('PredictArgNotTable', caller);
            else
                T = table;
                VarNames = this.PrivOptions.VarSpec.Names;
                for i = 1:numel(VarNames)
                    if ~ismember(VarNames{i}, XTable.Properties.VariableNames)
                        bayesoptim_priv.err('PredictVarMissing', VarNames{i});
                    end
                    T.(VarNames{i}) = XTable.(VarNames{i});
                end
                % Apply conditional variable function and canonicalize
                T = applyConditionalVariableFcn(this, T);
                T = canonicalizePoints(this, T);
            end
        end
        
        
        
        function tf = isPlottingObjFcnModel(this)
            tf = any(cellfun(@(x)isequal(x,@plotObjectiveModel), this.PrivOptions.PlotFcn));
        end
        
        
        % For fitcauto, remove irrelevant parameters
        function [learnerName,tbl] = simplifyAutoResults(this,tbl)
            tblNames = tbl.Properties.VariableNames;
            if this.PrivOptions.FitcautoMultipleLearners
                learnerName = char(tbl{end,1});
                tag = [learnerName,'_'];
                ind = find(startsWith(tblNames(2:end),tag))+1;
                tbl = tbl(:,[1,ind]);
            elseif this.PrivOptions.FitcAutoSingleLearner
                splits = split(tbl.Properties.VariableNames{1},'_');
                if size(splits,1) < 2
                    return;
                end
                learnerName = splits{1};
                ind = find(startsWith(tblNames(1:end),learnerName));
                tag = [learnerName,'_'];
            end
            
            if ~isempty(ind)
                tbl.Properties.VariableNames = erase(tbl.Properties.VariableNames,tag);
            end
        end
    end
    
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                obj = builtin('loadobj',s);
            else
                obj = s;
            end
            % TotalElapsedTimeTrace was added as a property in R2017b
            if isempty(obj.TotalElapsedTimeTrace)
                obj.TotalElapsedTimeTrace = nancumsum(obj.PrivIterationTimeTrace)';
            end
            % Changes in 2020a:
            %
            % The following properties were renamed, and the GP models were
            % packaged up inside objects of class
            % bayesoptim_priv.BayesOptGPModel:
            %
            %   ObjectiveFcnGP              --> ObjectiveFcnModel
            %   ObjectiveEvaluationTimeGP   --> ObjectiveEvaluationTimeModel
            %   ConstraintGPs               --> ConstraintModels
            %   ErrorGP                     --> ErrorModel
            %
            % Also, a new property was added: 'ModelStrategy', which should
            % be set to bayesoptim_priv.GPStrategy() for old models.
            if isempty(obj.ModelStrategy)
                obj.ModelStrategy = bayesoptim_priv.GPStrategy();
                if ~isempty(obj.ObjectiveFcnGP)
                    obj.ObjectiveFcnModel = bayesoptim_priv.BayesOptGPModel(...
                        obj.ObjectiveFcnGP.Sigma, ...
                        obj.ObjectiveFcnGP, ...
                        obj.ObjectiveFcnGP.KernelInformation, ...
                        obj.ObjectiveFcnGP.Beta);
                    obj.ObjectiveFcnGP = [];
                end
                if ~isempty(obj.ObjectiveEvaluationTimeGP)
                    obj.ObjectiveEvaluationTimeModel = bayesoptim_priv.BayesOptGPModel(...
                        obj.ObjectiveEvaluationTimeGP.Sigma, ...
                        obj.ObjectiveEvaluationTimeGP, ...
                        obj.ObjectiveEvaluationTimeGP.KernelInformation, ...
                        obj.ObjectiveEvaluationTimeGP.Beta);
                    obj.ObjectiveEvaluationTimeGP = [];
                end
                obj.ConstraintModels = cellfun(@(m)bayesoptim_priv.BayesOptGPModel(...
                    m.Sigma, m, m.KernelInformation, m.Beta), obj.ConstraintGPs, 'UniformOutput', false);
                obj.ConstraintGPs = {};
                if ~isempty(obj.ErrorGP)
                    obj.ErrorModel = bayesoptim_priv.BayesOptGPModel(...
                        obj.ErrorGP.Sigma, ...
                        obj.ErrorGP, ...
                        obj.ErrorGP.KernelInformation, ...
                        obj.ErrorGP.Beta);
                    obj.ErrorGP = [];
                end
            end
        end
    end
end

%% Local functions
function [VariableDescriptions, RemainingArgs] = iParseVariableDescriptionsArg(Args)
[VariableDescriptions,~,RemainingArgs] = internal.stats.parseArgs(...
    {'VariableDescriptions'}, {[]}, Args{:});
end

function val = iNanIfBad(val)
% val must be a finite real number. If not, set it to NaN.
if ~isscalar(val) || ~isnumeric(val) || ~isreal(val) || ~isfinite(val)
    val = NaN;
end
end

% function GP = iFitrgpRobust(X, Y, SigmaLowerBound, varargin)
% % Fit a GP. If it fails due to singularity, iteratively double
% % SigmaLowerBound until it succeeds (up to 10 times). If it never succeeds,
% % return [];
% if all(isnan(Y))
%     bayesoptim_priv.err('AllModelTargetsNaN');
% end
% SigmaLowerBound = max(SigmaLowerBound, 1e-6);
% GP = [];
% success = false;
% doublings = 0;
% while ~success && doublings <= 10
%     try
%         GPModel = compact(fitrgp(X, Y, varargin{:}, 'SigmaLowerBound', SigmaLowerBound));
%         GP = bayesoptim_priv.BayesOptGPModel(GPModel.Sigma, GPModel);
%         success = true;
%     catch me
%         SigmaLowerBound = 2*SigmaLowerBound;
%         doublings = doublings + 1;
%     end
% end
% end

% function RF = iFitrrfRobust(X, Y)
% % Fit a RF.
% if all(isnan(Y))
%     bayesoptim_priv.err('AllModelTargetsNaN');
% end
% treeBagger = TreeBagger(100, X, Y, 'Method', 'regression');
% RF = bayesoptim_priv.BayesOptRFModel(treeBagger);
% end

% function GP = refitGPWithFitMethodNone(GP, XTrain, FTrain)
% % Recompute alphas. This requires passing in the current parameters and
% % using 'FitMethod' 'none'.
% GP = iFitrgpRobust(XTrain, FTrain, 0,...
%     'KernelFunction', 'ardmatern52', ...
%     'KernelParameters', GP.KernelInformation.KernelParameters, ...
%     'Sigma', GP.Sigma, ...
%     'BasisFunction', 'constant', ...
%     'Beta', GP.Beta, ...
%     'FitMethod', 'none', ...
%     'PredictMethod', 'exact', ...
%     'Standardize', false);
% end

% % fixme: Do we need a special case for this
% function RF = refitRFWithFitMethodNone(RF, XTrain, FTrain)
% RF = iFitrrfRobust(XTrain, FTrain);
% end

% function s = iGetSigmaF(GPR)
% s = GPR.KernelInformation.KernelParameters(end);
% end

function [x, fval] = iFminbndGlobal(fcn, LB, UB, NumCandidates, BestK, ...
    Verbose, MaxIter, TolFun)
% Attempts to find the global minimum of fcn() within bounds given by the
% vectors LB and UB. Generates NumCandidates random points within bounds,
% runs fcn() on all of them and selects the BestK best. Then runs a bounded
% fminsearch on each of the K points and returns the best overall. The
% flags Verbose, MaxIter, and TolFun are options passed directly to
% fminsearch.

% Choose the best K points from N random candidates
[x0, fvals] = iMinKPointsFromN(fcn, LB, UB, NumCandidates, BestK);
% Do a local optimization on each of the K points
% if sum(isfinite(fvals)) < BestK
%     bayesoptim_priv.err('iMinKPointsFromNFailed');
% end
if Verbose
    VerboseOpts = {'Display', 'iter'};
else
    VerboseOpts = {'Display', 'off'};
end
for row = BestK:-1:1
    [xs(row,:), fvals(row)] = fminsearch(@boundedFcn, x0(row,:), ...
        optimset(...
        'MaxIter', MaxIter, ...
        'TolX', Inf, ...
        'TolFun', TolFun, ...
        VerboseOpts{:}));
end
[fval, loc] = min(fvals);
x = xs(loc,:);

    function y = boundedFcn(x)
        if any(x < LB | x > UB)
            y = Inf;
        else
            y = fcn(x);
        end
    end
end

function [X, fvals] = iMinKPointsFromN(fcn, LB, UB, N, K)
% Generate a random set of N points within bounds. Return the K points that
% best minimize fcn.
X = iUniformPointsWithinBounds(N, LB, UB);
fvals = fcn(X);
[~, rows] = sort(fvals, 'ascend');
X = X(rows(1:K), :);
fvals = fvals(rows(1:K));
end

function X = iUniformPointsWithinBounds(N, LB, UB)
Unifs = rand(N, size(LB,2));
X = Unifs .* repmat(UB - LB, N, 1);
X = X + repmat(LB, N, 1);
end

function tf = iAnyIllegalResumeArgs(NVPs)
Names = {'UseParallel', 'InitialX', 'InitialObjective', 'InitialConstraintViolations', ...
    'InitialErrorValues', 'InitialObjectiveEvaluationTimes', 'InitialIterationTimes', ...
    'InitialUserData', 'NumSeedPoints'};
Defaults = cell(1, numel(Names));
[~,~,~,~,~,~,~,~,~,setflag,~] = internal.stats.parseArgs(Names, Defaults, NVPs{:});
tf = bayesoptim_priv.anyFlagsSet(setflag);
end

function Indices = iFindUnusedGridIndices(NumGridDivisions, SampledGridIndices)
if isempty(SampledGridIndices)
    for p=numel(NumGridDivisions):-1:1
        Indices(p) = randi(NumGridDivisions(p));
    end
else
    Indices = SampledGridIndices(1,:);
    while ismember(Indices, SampledGridIndices, 'rows')
        for p=1:numel(NumGridDivisions)
            Indices(p) = randi(NumGridDivisions(p));
        end
    end
end
end

function checkVariableDescriptions(NewVariableDescrips, OldVariableDescrips)
% The variable names to be optimized must be the same. Numeric variables
% may have their Range, Type, and Transform changed, but must remain
% numeric. Categoricals cannot be changed.
Old = sortVDByName(OldVariableDescrips([OldVariableDescrips.Optimize]));
New = sortVDByName(NewVariableDescrips([NewVariableDescrips.Optimize]));
if ~isequal({Old.Name}, {New.Name})
    bayesoptim_priv.err('ResumeVarsChanged');
end
for i = 1:numel(Old)
    switch Old(i).Type
        case 'categorical'
            if ~isequal(New(i).Type, 'categorical') || ~isequal(Old(i).Range, New(i).Range)
                bayesoptim_priv.err('ResumeCatsChanged');
            end
        case {'integer','real'}
            if ~ismember(New(i).Type, {'integer','real'})
                bayesoptim_priv.err('ResumeNumericChanged');
            end
    end
end
end

function VD = sortVDByName(VD)
[~,I] = sort({VD.Name});
VD = VD(I);
end

function doubleVec = intCodeCategorical(catVec)
doubleVec = double(catVec);
end

function ObjectiveNargout = nargoutFromNumConstraints(NumCoupledConstraints)
% We know nargout < 3
if NumCoupledConstraints > 0
    ObjectiveNargout = 2;
else
    ObjectiveNargout = 1;
end
end

function S = nancumsum(X)
X(isnan(X)) = 0;
S = cumsum(X);
end

function FSD = funStd(YSD, Sigma)
% Return the standard deviation of the posterior over function values (F),
% given the std of the posterior over observed values (Y) and the noise
% standard deviation (Sigma). FSD, YSD, and Sigma are arrays.
FSD = sqrt(max(0, YSD.^2 - Sigma.^2)); % Use max to avoid complex numbers.
end

function Cleanup = holdOn(Axes)
% Turn hold on and return an oncleanup object to turn hold off.
hold(Axes, 'on');
Cleanup = onCleanup(@()holdOff(Axes));
    function holdOff(A)
        if isvalid(A)
            hold(A, 'off');
        end
    end
end

function tf = acquisitionFunctionIsGridOrRand(this)
tf = ismember(this.PrivOptions.AcquisitionFunctionName, {'grid','random'});
end

function eraseVerboseLine(VerboseLine)
fprintf(repmat('\b',1,numel(VerboseLine)))
end

function displayFitcautoResults(tbl,learnerName,displayMessage)

% remove learner from tbl to prevent duplicate printing
learnerIdx = find(strcmp('learner',tbl.Properties.VariableNames));

if ~isempty(learnerIdx)
    tbl(:,learnerIdx)= [];
end

if strcmp(learnerName,'nb')
    learnerName = 'naive Bayes';
end

codingIdx = find(strcmp('Coding',tbl.Properties.VariableNames));
if ~isempty(codingIdx)
    learnerName = ['multiclass ' learnerName];
end

if any(strcmp(learnerName,{'svm','ensemble'}))
    displayMessage = [displayMessage,'EnsembleSVM'];
end

if startsWith(displayMessage,'Iteration')
    fprintf('\n');
end

bayesoptim_priv.printInfo(displayMessage,learnerName);

for i=1:size(tbl,2)
    block = repmat(' ',1,27);
    
    % Left-justify parameter name
    pname = [tbl.Properties.VariableNames{i},':'];
    
    if strcmp(pname,'Coding:')
        pname = 'Coding (ECOC):';
    end
    
    plength = min(size(block,2),numel(pname));
    block(1:plength) = pname(1:plength);

    % Right-justify parameter value
    pval = tbl.(tbl.Properties.VariableNames{i});
    if isnumeric(pval)
        pval = strip(sprintf('%11.5g',pval));
    else
        pval = strip(char(pval));
    end
    
    vlength = min(size(block,2),numel(pval));
    block(end-vlength+1:end) = pval(1:vlength);
    
    if i == size(tbl,2) && startsWith(displayMessage,'Iteration')
        fprintf('\t%s',block);
    else
        fprintf('\t%s\n',block);
    end
    
end
end
function Results = bayesopt_priv(ObjectiveFcn, VariableDescriptions, varargin)
%BAYESOPT   Find the global minimum of a function using Bayesian
%Optimization.
%
%   (1) RESULTS = BAYESOPT(OBJECTIVEFCN, VARIABLEDESCRIPTIONS) attempts to
%   find values for the variables described in VARIABLEDESCRIPTIONS that
%   globally minimizes the function OBJECTIVEFCN. OBJECTIVEFCN may be
%   multivariate, expensive, non-deterministic, and lacking analytic
%   derivatives.
%
%     RESULTS   - An instance of class <a href="matlab:help BayesianOptimization">BayesianOptimization</a>.
%
%     OBJECTIVEFCN   - The objective function to minimize. This can 
%                    be a function handle, or if 'UseParallel' is true, it
%                    can also be a parallel.pool.Constant whose Value is a
%                    function handle. The function must have one of the
%                    signatures:
%
%                   Objective = objectiveFcn(X)
%
%                   [Objective, ConstraintViolation] = objectiveFcn(X)
%
%                   [Objective, ConstraintViolation, UserData] = objectiveFcn(X)
%
%           Input argument:
%               X           - A 1-by-D table of variable values, for D
%                           variables. The variable names of X match the
%                           names defined in the VariableDescriptions
%                           argument below. Variables declared 'real' or
%                           'integer' have MATLAB type 'double', and
%                           variables declared 'categorical' have MATLAB
%                           type 'categorical'.
%           Output arguments:
%               Objective	- The real, scalar function value evaluated at
%                           X. 
%               ConstraintViolation
%                           - A length-K vector of reals representing the
%                           degree of violation of a set of K "coupled"
%                           inequality constraints evaluated at X. A value
%                           less than or equal to 0 indicates that the
%                           corresponding constraint is satisfied.
%               UserData   	- Any additional object resulting from the
%                           function evaluation. A list of these values is
%                           stored in the Results object.
%
%     VARIABLEDESCRIPTIONS
%           - A vector of <a href="matlab:help optimizableVariable">optimizableVariable</a> objects describing the variables
%           over which to optimize. 
%
%   (2) RESULTS = BAYESOPT(__, 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%   optional name/value pairs:
%
%   Controlling the search:
%       AcquisitionFunctionName  - The name of the Acquisition Function to
%                                use in choosing the next point to
%                                evaluate. The following are accepted
%                                (hyphens are optional and unique prefixes
%                                are recognized):
%                   	'expected-improvement'
%                       'expected-improvement-plus'
%                       'expected-improvement-per-second'
%                       'expected-improvement-per-second-plus'
%                       'lower-confidence-bound'
%                       'probability-of-improvement'
%                     Default: 'expected-improvement-per-second-plus'
%
%       IsObjectiveDeterministic
%               	- Set to true if objectiveFcn is a deterministic
%               	function of X. If false, a noise level will be
%               	estimated during optimization. 
%                   Default: false
%
%       ExplorationRatio
%                   - When AcquisitionFunctionName is
%                   'expected-improvement-plus' or
%                   'expected-improvement-per-second-plus', this controls
%                   the exploration vs. exploitation tradeoff. Larger
%                   values result in more exploration, but limit the
%                   precision of the final solution. Specifically, if the
%                   posterior standard deviation of Objective at X is less
%                   than ExplorationRatio times the posterior standard
%                   deviation of the noise at X, then the algorithm
%                   concludes that it is "overexploiting" at X. It then
%                   makes an adjustment that encourages exploration.
%                   Default: 0.5
%
%   Parallel computing (Requires Parallel Computing Toolbox):
%       UseParallel - A logical scalar specifying whether to perform
%                   function evaluations on the current parallel pool.
%                   Default: false
%       ParallelMethod
%                   - A char array specifying the method used to impute
%                   Objective function values for points that are currently
%                   being evaluated. Choices are:
%                   'model-prediction', 'clipped-model-prediction', 
%                   'min-observed', 'max-observed'. 
%                   Default: 'clipped-model-prediction'
%       MinWorkerUtilization
%                   - A positive integer specifying the minimum acceptable
%                   number of busy parallel workers. When the number of
%                   busy workers falls below this number, all idle workers
%                   in the pool are given random points to evaluate.
%                   Default: The greatest integer not more than 80% of the
%                   number of workers in the current pool.
%
%   Model fitting:
%       GPActiveSetSize
%                   - A positive integer specifying the maximum number of
%                   points on which to fit the Gaussian Process models.
%                   This can be used to reduce the time needed to choose
%                   the next point to evaluate, at the cost of a less
%                   accurate model. 
%                   Default: 300
%
%   Starting and Stopping:
%       NumSeedPoints
%                   - The number of randomly-chosen points on which to
%                   evaluate the Objective function before using the
%                   Acquisition function. 
%                   Default: 4
%       MaxObjectiveEvaluations	
%                   - Optimization stops when this many function
%                   evaluations have been performed, provided all initial
%                   points have been evaluated. 
%                   Default: 30
%       MaxTime     - Optimization stops when the optimization has run for
%                   this many seconds of "wall clock time", provided all
%                   initial points have been evaluated. The actual runtime
%                   may be higher because function evaluations are not
%                   interrupted. 
%                   Default: Inf
%
%   Constrained optimization:
%       XConstraintFcn 
%                   - A function handle defining a deterministic constraint
%                   on the input, X, that solutions must satisfy. It must
%                   have signature:
%
%                           TF = xConstraint(X)
%
%                   where 'X' is a width-D table of arbitrary height. TF
%                   must be a logical column vector the same height as X,
%                   where TF(i) is true if the X constraint is satisfied
%                   for point X(i,:). This function should be fast to
%                   execute, because it will be called on each iteration of
%                   the optimization, on a table X whose height will be in
%                   the tens of thousands. Using vectorized operations is
%                   recommended. 
%                   Default: []
%       ConditionalVariableFcn
%                   - A function handle enforcing conditional variable
%                   dependencies. Two types of conditional variables are
%                   supported: (Type-1) A variable that is only defined
%                   when another variable has a particular value, and
%                   (Type-2) A variable that must take on a specific value,
%                   conditional on the values of other variable(s). The
%                   conditionalVariableFcn function must have signature:
%
%                           X = conditionalVariableFcn(X)
%
%                   where 'X' is a width-D table of arbitrary height. For
%                   each row of X, the function must set all type-1
%                   variables to NaN (or <undefined> if categorical) if not
%                   defined in that row, and set all type-2 variables to
%                   their correct values given the other variables in that
%                   row. 
%                   Default: []
%       NumCoupledConstraints
%                   - An integer specifying the number of coupled
%                   constraints. This argument is required if there are
%                   any coupled constraints, and must match the length of
%                   the second output of objectiveFcn. 
%                   Default: 0
%       CoupledConstraintTolerances
%                   - A vector of positive real numbers specifying a
%                   tolerance on the probability with which each coupled
%                   constraint must be satisfied. For a solution to be
%                   considered feasible, each coupled constraint C must be
%                   satisfied with probability at least
%                   1-CoupledConstraintTolerances(C). 
%                   Default: 0.01 for each coupled constraint
%       AreCoupledConstraintsDeterministic
%                   - A logical vector specifying whether each coupled
%                   constraint is deterministic (true) or nondeterministic
%                   (false). 
%                   Default: true for all coupled constraints
%
%   Reporting, Saving, and Interrupting execution:
%       Verbose     - A non-negative integer specifying the verbosity
%                   level of command-line display: 
%                       0: Do not display anything. 
%                       1: Each iteration, display the current point just 
%                          evaluated, the observed function value for that
%                          point, the function evaluation time, and the
%                          best observed function value so far.
%                       2: Same as (1), but also display the current X and
%                          current constraint values.
%                   Default: 1
%       OutputFcn  - A function handle or cell array of function handles
%                   that bayesopt calls after each iteration. This may be
%                   used to stop the optimization, produce command-line
%                   output, or perform arbitrary calculations. Choose from
%                   a predefined set or define your own. 
%                   Default: {}
%               Predefined output functions:
%                   * @assignInBase constructs a BayesianOptimization
%                   instance each iteration and assigns it to a variable in
%                   the base workspace.
%                   * @saveToFile constructs a BayesianOptimization
%                   instance each iteration and saves it to a file in the
%                   current directory.
%               User-defined output functions:
%                   Each output function must have the following signature:
%
%                       Stop = outputFcn(Results, State);
%
%                   Input Arguments:
%                       'Results'   - An object of class
%                                   BayesianOptimization describing the
%                                   current state of the optimization.
%                       'State'     - A character vector indicating the
%                                   state of the optimization. Possible
%                                   values are:
%                                   * 'initial': The algorithm is in the
%                                     initial state before the first
%                                     iteration.
%                                   * 'iteration': The algorithm is at the
%                                     end of an iteration.
%                                   * 'done': The algorithm is in the final
%                                     state after the last iteration.
%                   Output Argument:
%                       'Stop'      - A logical scalar indicating whether
%                                   the optimization should stop.
%       SaveVariableName
%                   - The name of the variable to use when using the
%                   @assignInBase output function. 
%                   Default: 'BayesoptResults'
%       SaveFileName
%                   - The name of the file to use when using the
%                   @saveToFile output function. 
%                   Default: 'BayesoptResults.mat'
%
%   Plotting:
%       PlotFcn    - A function handle, cell array of function handles, or
%                    'all', specifying a set of plotting functions to be
%                    called after each iteration. Choose from a predefined
%                    set or define your own. 'all' calls all predefined
%                    plot functions. To turn off plots, set to {}.
%                   Default: {@plotObjectiveModel, @plotMinObjective}
%             Predefined plot functions:
%               Model plots: These apply when D<=2.
%                   * @plotObjectiveModel plots the ObjectiveFcn model
%                     surface, the estimated location of the minimum, and
%                     the location of the next proposed point to evaluate.
%                     For 1-dimensional problems, envelopes are plotted 1
%                     credible interval above and below the mean function,
%                     as well as envelopes 1 noise SD above and below the
%                     mean.
%                   * @plotAcquisitionFunction plots the Acquisition
%                     Function surface.
%                   * @plotObjectiveEvaluationTimeModel plots the
%                     Function evaluation time model surface.
%                  	* @plotConstraintModels plots each Constraint model
%                  	  surface, including the built-in Error constraint. It
%                  	  also plots a Prob(Feasible) surface.
%               Trace plots: These apply for all D.
%                   * @plotObjective plots each observed Objective versus
%                     the number of function evaluations.
%                   * @plotObjectiveEvaluationTime plots each observed
%                     ObjectiveEvaluationTime versus the number of function
%                     evaluations.
%                   * @plotMinObjective plots the minimum observed
%                     Objective versus the number of function evaluations,
%                     and also versus elapsed time.
%                   * @plotElapsedTime plots the total elapsed time
%                     versus the number of function evaluations.
%           User-defined plot functions:
%                   Plot function syntax is the same as Output function
%                   syntax:
%
%                       Stop = plotFcn(Results, State);
%                   
%                   Typically, a plot function creates a plot window in the
%                   'initial' state and updates it in the 'iteration'
%                   state.
%
%   Initialization from prior function evaluations:
%       InitialX	- An N-by-D table of initial points to evaluate, for N
%                   points and D variables. Default: N = NumSeedPoints
%                   random initial points within bounds.
%       InitialObjective	- A length-N vector of initial Objective values
%                   corresponding to InitialX. Default: [].
%       InitialConstraintViolations
%                   - An N-by-K matrix of initial constraint values, for N
%                   points and K coupled constraints. Default: [].
%       InitialErrorValues	
%                   - A length-N vector of initial Error values {-1,1}
%                   corresponding to InitialX. Default: [].
%       InitialObjectiveEvaluationTimes
%                   - A length-N vector of initial ObjectiveEvaluationTime
%                   values corresponding to InitialX. Default: [].
%       InitialIterationTimes
%                   - A length-N vector of initial iteration times
%                   corresponding to InitialX. Default: [].
%       InitialUserData
%                   - A length-N cell vector of initial UserData
%                   corresponding to InitialX. Default: {}.
%       (Note: If only InitialX is provided, the objective function is
%        evaluated at InitialX. If any other initialization parameters
%        are provided, the objective function is not evaluated. Any
%        missing values are set to NaN.)
%
%   See also: BAYESIANOPTIMIZATION, OPTIMIZABLEVARIABLE

%   Copyright 2016-2017 The MathWorks, Inc.

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

C = bayesoptim.suppressWarnings();
Options = bayesoptim_priv.BayesoptOptions(ObjectiveFcn, VariableDescriptions, varargin);
Results = BayesianOptimization_priv(Options);
end

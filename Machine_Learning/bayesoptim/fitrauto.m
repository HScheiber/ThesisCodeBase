function[optimizedModel,optimizationResults] = fitrauto(X,Y,varargin)
% FITRAUTO performs Bayesian optimization across regression models and their 
% hyperparameters.
%	MODEL = FITRAUTO(TBL,Y) returns a compact regression model optimized for
%   the predictor data in the table TBL and the response variable Y.
%   Y can be any of the following:
%      1. A column vector of double or single numbers.
%      2. The name of a variable in TBL. This variable is used as the
%         response Y, and the remaining variables in TBL are used as
%         predictors.
%      3. A formula string such as 'y ~ x1 + x2 + x3' specifying that the
%         variable y is used as the response, and the other variables
%         in the formula are predictors. Any table variables not listed in
%         the formula are not used.
%
%       MODEL = FITRAUTO(X,Y) is an alternative syntax that accepts an
%       N-by-P matrix of predictors X with one row per observation and one
%       column per predictor. Y is the response vector.
%
%       [MODEL,OPTIMIZATIONRESULTS] = FITRAUTO() also returns an object
%       describing the results of the optimization.
%
%       MODEL = FITRAUTO (X,Y,Name,Value...) specifies optional parameter
%       name-value pair arguments:
%
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                          'CategoricalPredictors' as one of the following:
%                            * A numeric vector with indices between 1 and P,
%                              where P is the number of columns of X or
%                              variables in TBL.
%                            * A logical vector of length P, where a true
%                              entry means that the corresponding column of X
%                              or TBL is a categorical variable. 
%                            * 'all', meaning all predictors are categorical.
%                            * A string array or cell array of character
%                              vectors, where each element in the array is
%                              the name of a predictor variable. The names
%                              must match entries in the 'PredictorNames' 
%                              name-value pair argument.
%                          Default: For a matrix input X, no categorical
%                          predictors are selected; for a table TBL, predictors 
%                          are treated as categorical if they are strings, cell
%                          arrays of character vectors, logical, or
%                          categorical.
%       'OptimizeHyperparameters' 
%                        - Hyperparameters to optimize for the models specified 
%                          in the 'Learners' name-value pair argument, either 
%                          'auto' or 'all'. Use the 'HyperparameterOptimizationOptions' 
%                           name-value pair argument to control other aspects of the  
%                           optimization. The 'auto' value indicates to optimize 
%                           over a suitable subset of optimizable hyperparameters 
%                           for each model specified in 'Learners'. The 'all' value 
%                           indicates to optimize over all optimizable hyperparameters 
%                           for each model specified in 'Learners'. Default: 'auto'                          
%         'PredictorNames' - A cell array of names for the predictor
%                          variables.  For a matrix X, the order of the names 
%                          must match the order of the predictors in X.
%                          Matrix default: {'x1','x2',...}. For a table TBL, these
%                          names must be a subset of the variable names in
%                          TBL, and only the selected variables are used. Not
%                          allowed when Y is a formula. Table default: all
%                          variables other than Y.
%                          This name-value pair argument is not supported when Y is a formula.
%         'ResponseName' - Name of the response variable Y, a character vector or 
%                          string scalar. This name-value pair argument is not supported
%                          when Y is a formula. Default: 'Y'
%         'Weights'      - Vector of observation weights, with one weight per
%                          observation. fitrauto normalizes the weights to
%                          sum to 1. Default: ones(size(X,1),1).
%                          For an input table TBL, the 'Weights' value can be
%                          the name of a variable in TBL.
%         'Learners'     - Types of regression models to try during the 
%                          optimization process. Either 'auto','all','alllinear', 
%                          'allnonlinear' or a string or cell array of  
%                          eligible learner names. If you specify 'auto', a set of  
%                          learners suitable for the data set is selected and used  
%                          to find the optimized model. 'all' is equivalent to   
%                          {'svm','gp','tree','linear','ensemble','kernel'}.
%                          Default: 'auto'
%                          
%                          
%   Refer to the MATLAB documentation for information on parameters for
%       <a href="matlab:helpview(fullfile(docroot,'stats','stats.map'), 'fitrautoHyperparameterOptimizationOptions')">Hyperparameter Optimization</a>
%
%   Example: Create an optimized regression model for the carsmall data
%      load carsmall
%      X = [Weight,Displacement,Acceleration];
%      mdl = fitrauto(X,MPG);

%   Copyright 2020 The MathWorks, Inc.
if nargin > 1
    Y = convertStringsToChars(Y);
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

options = bayesoptim.MultiModelOptimOptions(false,varargin);
regressionMultiModelOptim = bayesoptim.RegressionMultiModelOptim(X,Y,options);
[optimizedModel,optimizationResults] = regressionMultiModelOptim.performOptimization;
end

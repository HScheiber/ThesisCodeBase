function[optimizedModel,optimizationResults] = fitcauto(X,Y,varargin)
% FITCAUTO performs Bayesian optimization across classification models and their 
% hyperparameters.
%	MODEL = FITCAUTO(TBL,Y) returns a compact classification model optimized for  
%   the predictor data in the table TBL and response variable Y.
%   Y can be any of the following:
%      1. An array of class labels. Y can be a categorical array, logical
%         vector, numeric vector, string array or cell array of character
%         vectors.
%      2. The name of a variable in TBL. This variable is used as the
%         response Y, and the remaining variables in TBL are used as
%         predictors.
%      3. A formula string such as 'y ~ x1 + x2 + x3' specifying that the
%         variable y is used as the response, and the other variables
%         in the formula are predictors. Any table variables not listed in
%         the formula are not used.
%
%       MODEL = FITCAUTO(X,Y) is an alternative syntax that accepts an 
%       N-by-P matrix of predictors X with one row per observation and one column 
%       per predictor. The response Y is an array of N class labels.
%
%       [MODEL,OPTIMIZATIONRESULTS] = FITCAUTO(...) also
%       returns an object describing the results of the optimization.
%
%       MODEL = FITCAUTO (X,Y,'PARAM1',val1,'PARAM2',val2,...) 
%       specifies optional parameter name-value pair arguments:
%
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                          'CategoricalPredictors' as one of the following:
%                            * A numeric vector with indices between 1 and P,
%                              where P is the number of columns of X or
%                              variables in TBL.
%                            * A logical vector of length P, where a true
%                              entry means that the corresponding column of X
%                              or T is a categorical variable. 
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
%       'ClassNames'   -   Array of class names. Use the data type that
%                          exists in Y. You can use this name-value pair argument 
%                          to order the classes or select a subset of classes for
%                          training. Default: All class names in Y
%       'Cost'         -   Square matrix, where COST(I,J) is the
%                          cost of classifying a point into class J if its
%                          true class is I. Alternatively, COST can be a
%                          structure S with two fields: S.ClassificationCosts
%                          containing the cost matrix C, and S.ClassNames
%                          containing the class names and defining the
%                          order of classes used for the rows and columns
%                          of the cost matrix. For S.ClassNames use the data
%                          type that exists in Y. Default: COST(I,J)=1 if
%                          I~=J, and COST(I,J)=0 if I=J  
%       'OptimizeHyperparameters' 
%                        - Hyperparameters to optimize for the models specified 
%                          in the 'Learners' name-value pair argument, either 
%                          'auto' or 'all'. Use the HyperparameterOptimizationOptions 
%                           name-value pair argument to control other aspects of the  
%                           optimization. The 'auto' value indicates to optimize 
%                           over a suitable subset of optimizable hyperparameters 
%                           for each model specified in 'Learners'. The 'all' value 
%                           indicates to optimize over all optimizable hyperparameters 
%                           for each model specified in 'Learners'. Default: 'auto'                          
%         'PredictorNames' - A cell array of names for the predictor
%                          variables.For a matrix X, the order of the names must 
%                          match the order of the predictors in X.
%                          Matrix default: {'x1','x2',...}. 
%                          For a table TBL, these names must be a subset of the 
%                          variable names in TBL, and only the selected variables  
%                          are used.Table default: all variables other than Y
%                          allowed when Y is a formula. Default: all
%                          variables other than Y
%                          This name-value pair argument is not supported when Y 
%                          is a formula.
%         'Prior'        - Prior probabilities for each class. Specify as one
%                          of the following: 
%                           * A character vector or string scalar:
%                             - 'empirical' determines class probabilities
%                               from class frequencies in Y
%                             - 'uniform' sets all class probabilities equal
%                           * A vector (one scalar value for each class)
%                           * A structure S with two fields: S.ClassProbs
%                             containing a vector of class probabilities, and
%                             S.ClassNames containing the class names and
%                             defining the order of classes used for the
%                             elements of this vector
%                           If you pass numeric values, fitcauto normalizes
%                           them to sum to one. Default: 'empirical'
%       'ScoreTransform' - Function handle for transforming scores, or
%                          a character vector or string scalar representing a 
%                          built-in transformation function. The available  
%                          functions: 'symmetric','invlogit', 'ismax', 'symmetricismax', 
%                          'none', 'logit', 'doublelogit', 'symmetriclogit',
%                          and 'sign'. Default: 'none'
%         'ResponseName' - Name of the response variable Y, a character vector or 
%                          string scalar. This name-value pair argument is not supported 
%                          when Y is a name or formula. 
%                          Default: 'Y'
%         'Weights'      - Vector of observation weights, with one weight per
%                          observation. fitcauto normalizes the weights to
%                          sum to the value of the prior probability in
%                          each respective class. Default: ones(size(X,1),1)
%                          For an input table TBL, the 'Weights' value can be
%                          the name of a variable in TBL.
%         'Learners'     - Types of classification models to try during the 
%                          optimization process, 'auto', 'all', 'all-linear', 
%                          'all-nonlinear' or a string or cell array of  
%                          eligible learner names. If you specify 'auto', a set of  
%                          learners suitable for the data set is selected and used  
%                          to find the optimized model. 'all' is equivalent to   
%                          {'svm','knn','nb','tree','discr','linear','ensemble',
%                           'kernel'}. Default: 'auto'
%                          
%                          
%   Refer to the MATLAB documentation for information on parameters for
%       <a href="matlab:helpview(fullfile(docroot,'stats','stats.map'), 'fitcautoHyperparameterOptimizationOptions')">Hyperparameter Optimization</a>
%
%   Example: Create an optimized classifier for the Fisher iris data, and compare
%            the actual and predicted species in a confusion matrix.
%      t = readtable('fisheriris.csv');
%      mdl = fitcauto(t,'Species');
%      confusionmat(t.Species,predict(mdl,t))

%   Copyright 2019 The MathWorks, Inc.
if nargin > 1
    Y = convertStringsToChars(Y);
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

options = bayesoptim.MultiModelOptimOptions(true,varargin);
classificationMultiModelOptim = bayesoptim.ClassificationMultiModelOptim(X,Y,options);
[optimizedModel,optimizationResults] = classificationMultiModelOptim.performOptimization;
end

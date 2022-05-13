function obj = fitrmultirf(X,y,numtrees,varargin)
%FITRMULTIRF Fit a multilevel random forest for regression
%   RF=FITRMULTIRF(X,Y,NUMTREES) returns model RF with several random
%   forest models, one per unique value in the first column of the input
%   matrix X, for regression response Y. FITRMULTIRF splits data into
%   subsets, each corresponding to a unique value in X(:,1) and fits a
%   random forest regression model with NUMTREES trees for each subset.
%
%   Pass X as a floating-point matrix with one row per observation and one
%   column per predictor. Pass Y as a floating-point vector with as many
%   rows as there are in X.
%   
%   RF=fitrmultirf(X,Y,NUMTREES,'PARAM1',val1,'PARAM2',val2,...) accepts
%   all optional parameters for TreeBagger.
%
%   See also TreeBagger.

%   Copyright 2019 The MathWorks, Inc.

D = size(X,2) - 1; % dimensionality minus the 1st column

% Get the unique groups
grp = X(:,1);
groups = unique(grp);
G = numel(groups);

% Initialize arrays for models and predictors
models = cell(G,1);
predictors = false(D,G);

% Catch categorical predictors
[catpreds,~,extra] = internal.stats.parseArgs(...
    {'categoricalpredictors'},{[]},varargin{:});

% Validate CategoricalPredictors
if isempty(catpreds)
    catpreds = false(1,D);
else
    if islogical(catpreds)
        catpreds(1) = [];
    elseif isnumeric(catpreds)
        catpreds = ismember(2:D+1,catpreds);
    else
        error(message('stats:bayesoptim:fitrmultirf:BadCatPred',D))
    end
end

% Fit a TreeBagger model to each unique group
for g=1:G
    isg = grp==groups(g);
    Xg = X(isg,2:end);
    notnan = ~all(isnan(Xg),1);
    predictors(:,g) = notnan';
    if ~all(isnan(y(isg)))
        models{g} = compact(TreeBagger(numtrees,Xg(:,notnan),y(isg),extra{:},...
            'CategoricalPredictors',catpreds(notnan),'Method','regression'));
    end
end

% Pass data and models to MultiTreeBagger
ymean = mean(y,'omitnan');
ystd = std(y,'omitnan');
obj = bayesoptim.MultiTreeBagger(groups,predictors,models,ymean,ystd);
end

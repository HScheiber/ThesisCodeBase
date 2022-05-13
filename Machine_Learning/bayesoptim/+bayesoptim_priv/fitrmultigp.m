function obj = fitrmultigp(X,y,SigmaLowerBound,varargin)
%FITRMULTIGP Fit a multilevel gaussian process regression model
%   GP=FITRMULTIGP(X,Y,SIGMALOWERBOUND) returns model GP with several GPR
%   models, one per unique value in the first column of the input matrix X,
%   for regression response Y. FITRMULTIGP splits data into subsets, each
%   corresponding to a unique value in X(:,1) and fits a gaussian process
%   regression model for each subset.
%
%   Pass X as a floating-point matrix with one row per observation and one
%   column per predictor. Pass Y as a floating-point vector with as many
%   rows as there are in X.
%   
%   GP=FITRMULTIGP(X,Y,SIGMALOWERBOUND,'PARAM1',val1,'PARAM2',val2,...)
%   accepts all optional parameters for GP.
%
%   See also GP.

%   Copyright 2020 The MathWorks, Inc.

D = size(X,2) - 1; % dimensionality minus the 1st column

% Get the unique groups
grp = X(:,1);
groups = unique(grp);
G = numel(groups);

% Initialize arrays for models and predictors
models = cell(G,1);
predictors = false(D,G);

% Fit a GP model to each unique group
for g=1:G
    isg = grp==groups(g);
    Xg = X(isg,2:end);
    notnan = ~all(isnan(Xg),1);
    predictors(:,g) = notnan';
    varargin = updateSizeForKernelOptions(varargin,notnan);
    if ~all(isnan(y(isg)))
        models{g} = getModel(Xg(:,notnan),y(isg),SigmaLowerBound,varargin{:});
    end
end

% Pass data and models to MultiGP
ymean = mean(y,'omitnan');
ystd = std(y,'omitnan');
obj = bayesoptim.MultiGP(groups,predictors,models,ymean,ystd);
end

function GP = getModel(X,Y,SigmaLowerBound,varargin)
% Fit a GP model. For reference, check the GPStrategy.iFitRobust method
SigmaLowerBound = max(SigmaLowerBound, 1e-6);
GP = [];
success = false;
doublings = 0;
while ~success && doublings <= 10
    try
        GP = compact(fitrgp(X, Y, varargin{:}, 'SigmaLowerBound', SigmaLowerBound));
        success = true;
    catch me
        SigmaLowerBound = 2*SigmaLowerBound;
        doublings = doublings + 1;
    end
end
end

function args = updateSizeForKernelOptions(args,notnan)
% From KernelParameters N-V pair, remove columns that 
idx = find(contains(args(1:2:end),'KernelParameters'));
idx = idx*2;
for ind = 1:numel(idx)
    param = args{idx(ind)};
    paramCopy = param(2:end-1);
    if numel(paramCopy) == numel(notnan)
        args{idx(ind)} = [paramCopy(notnan);param(end)];
    end
end
end
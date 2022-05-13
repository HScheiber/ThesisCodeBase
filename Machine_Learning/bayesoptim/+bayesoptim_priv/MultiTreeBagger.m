classdef MultiTreeBagger
%

%   Copyright 2019 The MathWorks, Inc.

    properties(GetAccess = 'public', SetAccess = 'private')
        % A NumGroups-by-1 vector that stores the unique categories in the
        % first column of X as passed to fitrmultirf
        Groups
        
        % A NumPredictors-by-NumGroups logical array with ones(true)
        % indicating predictors that do not contain all NaNs and
        % zeros(false) indicating predictors containing only NaNs for each
        % group
        Predictors
        
        % A NumGroups-by-1 cell array of compact TreeBagger models. One per
        % group
        Models
        
        % Scalar that stores the mean of the response y
        ResponseMean
        
        % Scalar that stores the standard deviation of the response y
        ResponseStandardDeviation
    end
    
    methods(Hidden)
        function obj = MultiTreeBagger(groups,predictors,models,ymean,ystd)
            obj.Groups = groups;
            obj.Predictors = predictors;
            obj.Models = models;
            obj.ResponseMean = ymean;
            obj.ResponseStandardDeviation = ystd;
        end
    end
       
    methods
        function obj = refit(obj,group,X,y,numtrees,varargin)
        %REFIT Refit random forest for a specific group
        %   RF=REFIT(RF,GROUP,X,Y,NUMTREES) refits random forest model RF for group
        %   index GROUP. Use the same syntax as for FITRMULTIRF for all other
        %   inputs.
            
            isg = X(:,1)==group;
            X = X(isg,2:end);
            
            if ~any(isg)
                error(message('stats:bayesoptim:fitrmultirf:MissingGroup',group))
            end
            
            if all(isnan(y(isg)))
                return
            else
                [found,loc] = ismember(group,obj.Groups);
                if ~found
                    obj.Groups(end+1) = group;
                    obj.Models(end+1) = {[]};
                    obj.Predictors(:,end+1) = false;
                    loc = numel(obj.Groups);
                end

                notnan = ~all(isnan(X),1);
                obj.Predictors(:,loc) = notnan';

                % Catch categorical predictors
                [catpreds,~,extra] = internal.stats.parseArgs(...
                    {'categoricalpredictors'},{[]},varargin{:});
                D = size(obj.Predictors,1);
                if isempty(catpreds)
                    catpreds = false(1,D);
                else
                    if islogical(catpreds)
                        catpreds(1) = [];
                    elseif isnumeric(catpreds)
                        catpreds = ismember(2:D+1,catpreds);
                    else
                        error(message('stats:bayesoptim:fitrmultirf:BadCatPred',D));
                    end
                end
            
                % Update mean and standard deviation
                obj.ResponseMean = mean(y,'omitnan');
                obj.ResponseStandardDeviation = std(y,'omitnan');

                % Fit
                obj.Models{loc} = compact(TreeBagger(numtrees,...
                    X(:,notnan),y(isg),extra{:},...
                    'CategoricalPredictors',catpreds(notnan),'Method','regression'));
            end
        end
        
        function [yfit,stdevs] = predict(obj,X)
        %PREDICT Prediction by regression random forest.
        %   [YFIT,STDEVS]=PREDICT(RF,X) returns predicted response YFIT and
        %   standard deviation of response STDEVS for random forest RF and
        %   floating-point matrix X. 
        %
        %   Pass RF as a random forest model trained by FITRMULTIRF. Pass X as a
        %   floating-point matrix with first column filled with group numbers used
        %   for training this model.
            
            grp = X(:,1);
            
            G = numel(obj.Groups); % number of groups
            N = size(X,1); % number of observations
            
            yfit = repmat(obj.ResponseMean,N,1);
            stdevs = repmat(obj.ResponseStandardDeviation,N,1);
            
            for g=1:G
                isg = grp==obj.Groups(g);
                preds = [false; obj.Predictors(:,g)];
                if any(isg) && ~isempty(obj.Models{g})
                    [yfit(isg),stdevs(isg)] = predict(obj.Models{g},X(isg,preds));
                elseif isempty(obj.Models{g})
                    yfit = NaN(N,1);
                    stdevs = NaN(N,1);
                end
            end
        end
    end
    
end

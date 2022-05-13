classdef MultiGP
%

%   Copyright 2020 The MathWorks, Inc.

    properties(GetAccess = 'public', SetAccess = 'private')
        % A NumGroups-by-1 vector that stores the unique categories in the
        % first column of X as passed to fitrmultigp
        Groups
        
        % A NumPredictors-by-NumGroups logical array with ones(true)
        % indicating predictors that do not contain all NaNs and
        % zeros(false) indicating predictors containing only NaNs for each
        % group
        Predictors
        
        % A NumGroups-by-1 cell array of compact GP models. One per
        % group
        Models
        
        % Scalar that stores the mean of the response y
        ResponseMean
        
        % Scalar that stores the standard deviation of the response y
        ResponseStandardDeviation
    end
    
    methods(Hidden)
        function obj = MultiGP(groups,predictors,models,ymean,ystd)
            obj.Groups = groups;
            obj.Predictors = predictors;
            obj.Models = models;
            obj.ResponseMean = ymean;
            obj.ResponseStandardDeviation = ystd;
        end
    end
       
    methods
        function obj = refit(obj,group,X,y,sigmaLB,varargin)
        %REFIT Refit gaussian process for a specific group
        %   GP=REFIT(GP,GROUP,X,Y) refits gaussian process model GP for
        %   group index GROUP. Use the same syntax as for FITRMULTIGP for
        %   all other inputs.
            
            % Set X to be only the observations for the current learner
            % ('group'), and remove the first predictor, which is the
            % learner:
            isg = X(:,1)==group;
            X = X(isg,2:end);
            % Remove the learner from the NVP arguments 'KernelParameters'
            % and 'ConstantKernelParameters':
            idx = find(cellfun(@(x)isequal(x,'KernelParameters'), varargin), 1, 'last');
            if ~isempty(idx)
                varargin{idx+1} = varargin{idx+1}(2:end);
            end
            idx = find(cellfun(@(x)isequal(x,'ConstantKernelParameters'), varargin), 1, 'last');
            if ~isempty(idx)
                varargin{idx+1} = varargin{idx+1}(2:end);
            end
            
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
            
                % Update mean and standard deviation
                obj.ResponseMean = mean(y,'omitnan');
                obj.ResponseStandardDeviation = std(y,'omitnan');

                % Fit
                obj.Models{loc} = iFitRobust(obj, X(:,notnan),y(isg),...
                    sigmaLB,varargin{:});
            end
        end
        
        function [yfit,stdevs,yint] = predict(obj,X)
        %PREDICT Prediction by regression gaussian process.
        %   [YFIT,STDEVS]=PREDICT(GP,X) returns predicted response YFIT and
        %   standard deviation of response STDEVS for gaussian process GP
        %   and floating-point matrix X.
        %
        %   Pass GP as a gaussian process model trained by FITRMULTIGP.
        %   Pass X as a floating-point matrix with first column filled with
        %   group numbers used for training this model.
            
            grp = X(:,1);
            
            G = numel(obj.Groups); % number of groups
            N = size(X,1); % number of observations
            
            yfit = repmat(obj.ResponseMean,N,1);
            stdevs = repmat(obj.ResponseStandardDeviation,N,1);
            yint = NaN(N,2);
            
            for g=1:G
                isg = grp==obj.Groups(g);
                preds = [false; obj.Predictors(:,g)];
                if any(isg) && ~isempty(obj.Models{g})
                    [yfit(isg),stdevs(isg),yint(isg,:)] = predict(obj.Models{g},X(isg,preds));
                elseif isempty(obj.Models{g})
                    yfit = NaN(N,1);
                    stdevs = NaN(N,1);
                end
            end
        end
    end
    
    methods(Access=protected)
        function GP = iFitRobust(~, X, Y, SigmaLowerBound, varargin)
            % This is equivalent to the iFitRobust method on GPStrategy.
            % There are occasions where the SigmaLowerBound estimate is too
            % small and we need to attempt to fit the GP with a larger
            % estimate.
            
            % Fit a GP. If it fails due to singularity, iteratively double
            % SigmaLowerBound until it succeeds (up to 10 times). If it never succeeds,
            % return [];
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
    end
end

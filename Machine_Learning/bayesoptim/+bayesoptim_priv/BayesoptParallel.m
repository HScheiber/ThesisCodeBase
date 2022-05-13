classdef BayesoptParallel
    %BayesoptParallel  A class to manage parallel bayesopt.
    
        %   Copyright 2016-2017 The MathWorks, Inc.
    
    properties(SetAccess=protected)
        NumWorkers
        PendingX
        ObjectiveNargout = [];        % Needed by BayesianOptimization class
    end
    
    properties(Access=protected)
        FuturesArray
        ObjFcnPPC       % A parallel.pool.Constant
        LaunchTime      % The number of seconds it took to launch each job.
        Cleanup
    end
    
    properties(Dependent, Access=protected)
        ParPool
    end
    
    methods
        % Getters
        function Pool = get.ParPool(~)
            Pool = gcp;
            if isempty(Pool)
                bayesoptim.err('NoParallelPool');
            end
        end
        
        function tf = poolIsOpen(~)
            tf = ~isempty(gcp('nocreate'));
        end
        
        % Constructor, initialization
        function this = BayesoptParallel(objFcn, NumCoupledConstraints, Verbose)
            % Start the parallel pool by accessing it.
            pool = this.ParPool;
            this.NumWorkers = pool.NumWorkers;
            % Make sure the objFcn and supporting code is on all the workers
            this = ensureObjFcnIsOnWorkers(this, objFcn, Verbose);
            % Get nargout if possible
            this = checkObjectiveFcnNargout(this, objFcn, NumCoupledConstraints);
            % Init FuturesArray to a 1x0 Future
            this.FuturesArray = parallel.FevalFuture;
            this.FuturesArray(1) = [];
            % Suppress serialization warning about futures
            this = suppressSerializationWarning(this);
        end
        
        function this = ensureObjFcnIsOnWorkers(this, objFcn, Verbose)
            % objFcn is either a function handle or a parallel.pool.Constant
            if isa(objFcn, 'function_handle')
                % Copy the objFcn and supporting code to all the workers.
                if Verbose>0
                    bayesoptim.printInfo('CopyObjToWorkers');
                end
                try
                    this.ObjFcnPPC = copyFunctionHandleToWorkers(objFcn);
                catch ME
                    if isequal(ME.identifier, 'parallel:cluster:FileDoesNotExist')
                        bayesoptim.err('ParallelclusterFileDoesNotExist', ME.message);
                    else
                        rethrow(ME);
                    end
                end
                if Verbose>0
                    bayesoptim.printInfo('DoneCopyObjToWorkers');
                end
            elseif isa(objFcn, 'parallel.pool.Constant')
                % It's already a Constant. Verify that we can access its Value property.
                parfor i=1:this.ParPool.NumWorkers
                    try
                        v = objFcn.Value;
                        failed(i) = false;
                    catch
                        failed(i) = true;
                    end
                end
                if any(failed)
                    bayesoptim.err('ObjFcnPPCHasNoValue');
                else
                    this.ObjFcnPPC = objFcn;
                end
            end
        end
        
        % Job management
        function this = launchJob(this, ObjectiveNargout, NumCoupledConstraints, X, conditionalizedX)
            t                          = tic;
            this.FuturesArray(end+1)   = parfeval(this.ParPool, @bayesoptim.callObjFcnParallel, 6, this.ObjFcnPPC, ...
                ObjectiveNargout, NumCoupledConstraints, conditionalizedX);
            this.PendingX(end+1,:)     = X;
            this.LaunchTime(end+1)     = toc(t);
        end
        
        function waitForAnyJobToFinish(this)
            assert(~isempty(this.FuturesArray));
            fetchNext(this.FuturesArray);   % This line blocks until any future finishes.
        end
        
        function [FinishedX, Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                TotalJobTime, LaunchTime, ObjectiveNargout, this] = getJobResults(this)
            % Fetch the results of the first finished job and remove it
            % from FuturesArray, PendingX, LaunchTime.
            Idx = find(arrayfun(@(f)isequal(f.State,'finished'), this.FuturesArray), 1, 'first');
            assert(~isempty(Idx));      % Assert this because getJobResults should only be called after anyFinishedJobs returns true.
            FinishedX = this.PendingX(Idx,:);
            [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime, ObjectiveNargout] ...
                = fetchOutputs(this.FuturesArray(Idx));
            TotalJobTime = ObjectiveEvaluationTime + this.LaunchTime(Idx);
            LaunchTime = this.LaunchTime(Idx);
            this.FuturesArray(Idx)   = [];
            this.PendingX(Idx,:)     = [];
            this.LaunchTime(Idx) = [];
        end
        
        function cancelAllJobs(~)
            p=gcp;
            arrayfun(@cancel, p.FevalQueue.RunningFutures);
        end
        
        % Queries
        function tf = anyPendingX(this)
            tf = ~isempty(this.PendingX);
        end
        
        function n = numPendingX(this)
            n = size(this.PendingX,1);
        end
        
        function tf = anyFreeWorkers(this)
            tf = numel(this.FuturesArray) < this.NumWorkers;
        end
        
        function tf = anyFinishedJobs(this)
            rethrowJobError(this);
            tf = any(arrayfun(@(f)isequal(f.State,'finished'), this.FuturesArray));
        end
    end
    
    methods(Access=protected)
        
        function this = suppressSerializationWarning(this)
            % Suppress serialization warning about futures
            s = warning;
            this.Cleanup = onCleanup(@()warning(s));
            warning('off','parallel:cluster:CannotSaveCorrectly');
        end
        
        function rethrowJobError(this)
            % Rethrow the first error found
            Idx = find(arrayfun(@(f)~isempty(f.Error), this.FuturesArray),1);
            if ~isempty(Idx)
                rethrow(this.FuturesArray(Idx).Error.remotecause{1});
            end
        end
        
        function this = checkObjectiveFcnNargout(this, objFcn, NumCoupledConstraints)
            if isa(objFcn, 'parallel.pool.Constant')
                C = objFcn;
                % Check that all workers have a function handle and get nargouts
                p = gcp;
                parfor i=1:p.NumWorkers
                    Nargout{i} = checkObjFcnTypeAndNargout(C.Value, NumCoupledConstraints);
                end
                Nargouts = cell2mat(Nargout(:));
                if ~all(Nargouts==Nargouts(1))
                    bayesoptim.err('ObjFcnsHaveDifferentNargouts');
                else
                    % Nargout may be negative, in which case do not set it.
                    if Nargouts(1) > 0
                        this.ObjectiveNargout = Nargouts(1);
                    end
                end
            end
        end
    end
end

function Nargout = checkObjFcnTypeAndNargout(ObjFcn, NumCoupledConstraints)
if ~isa(ObjFcn, 'function_handle')
    bayesoptim.err('ObjectiveFcnTypeInPPC');
else
    try
        Nargout = nargout(ObjFcn);
    catch me
        switch me.identifier
            case 'MATLAB:nargin:isScript'
                % this.ObjectiveFcn is a function defined
                % inside a script.
                return;
            case 'MATLAB:narginout:functionDoesnotExist'
                bayesoptim.err('ObjFcnDoesNotExist', char(ObjFcn));
            otherwise
                rethrow(me);
        end
    end
    if Nargout == 0
        bayesoptim.err('ObjectiveNargoutZero');
    elseif Nargout == 1 && NumCoupledConstraints > 0
        bayesoptim.err('NumCoupledConstraintsMismatch');
    end
end
end
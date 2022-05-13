function C = copyFunctionHandleToWorkers(fHandle)
%copyFunctionHandleToWorkers Copy objective function handle to parallel
%workers, sending any required code.
%
%   C = copyFunctionHandleToWorkers(fHandle) copies the function handle
%   fHandle to all workers in the current parallel pool, and sends any code
%   required to execute the function on the workers. C is an object of
%   class parallel.pool.Constant.

%   Copyright 2017-2019 The MathWorks, Inc.

% IMPORTANT: This function must attach the files BEFORE sending the function handle.
attachMissingFunctionHandleFilesToPool(fHandle, gcp);
% Use an anonymous function that will be evaluated on workers to set the
% function handle into the constant.
C = parallel.pool.Constant(@() fHandle);
end

function attachMissingFunctionHandleFilesToPool(fHandle, pool)
% The command
% parallel.internal.apishared.AttachedFiles.convertFunctionHandleForDependencyAnalysis(fHandle)
% is internal and may not be supported in the future.
NeededFunctionsAndFiles = parallel.internal.apishared.AttachedFiles.convertFunctionHandleForDependencyAnalysis(fHandle);
Paths = cellfun(@which, NeededFunctionsAndFiles, 'UniformOutput',false);
Builtin = cellfun(@isBuiltin, Paths);
ToAttach = setdiff(Paths(~Builtin), pool.AttachedFiles);
if ~isempty(ToAttach)
    pool.addAttachedFiles(ToAttach);
end
end

function tf = isBuiltin(path)
tf = strncmp(path, 'built-in', 8);
end

function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime, ObjectiveNargout] = ...
    callObjFcnParallel(ObjectiveFcnPPC, ObjectiveNargout, NumCoupledConstraints, conditionalizedX)
% ObjectiveFcnPPC.Value is either a function handle or a cell containing a
% function handle.

    %   Copyright 2017 The MathWorks, Inc.

ObjectiveFcn = ObjectiveFcnPPC.Value;
if ~isempty(ObjectiveNargout)
    % nargout is known. Call ObjectiveFcn normally.
    [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
        = callObjNormallyParallel(ObjectiveFcn, ObjectiveNargout, NumCoupledConstraints, conditionalizedX);
else
    % nargout is unknown. Try nargout=3:
    StartTic = tic;
    ErrorConstraintViolation = -1;
    try
        [Objective, ConstraintViolations, UserData] = ObjectiveFcn(conditionalizedX);
    catch msg
        % nargout may be < 3.
        % Set nargout based on NumCoupledConstraints:
        ObjectiveNargout0 = nargoutFromNumConstraints(NumCoupledConstraints);
        ObjectiveNargout = ObjectiveNargout0;
        % nargout is now known. Call ObjectiveFcn normally.
        [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
            = callObjNormallyParallel(ObjectiveFcn, ObjectiveNargout, NumCoupledConstraints, conditionalizedX);
        return;
    end
    % Calling with nargout=3 succeeded. Set nargout and finish up.
    ObjectiveNargout = 3;
    [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
        = finishObjEvalParallel(Objective, ConstraintViolations, -1, StartTic, NumCoupledConstraints);
end
end

function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
    = callObjNormallyParallel(ObjectiveFcn, ObjectiveNargout, NumCoupledConstraints, conditionalizedX)
% Call the objective fcn with the correct nargout, catching
% errors and timing runtime
StartTic = tic;
ErrorConstraintViolation = -1;
try
    switch ObjectiveNargout
        case 1
            Objective = ObjectiveFcn(conditionalizedX);
            ConstraintViolations = [];
            UserData = [];
        case 2
            [Objective, ConstraintViolations] = ObjectiveFcn(conditionalizedX);
            UserData = [];
        otherwise
            [Objective, ConstraintViolations, UserData] = ObjectiveFcn(conditionalizedX);
    end
catch msg
    if isequal(msg.identifier, 'MATLAB:unassignedOutputs') && NumCoupledConstraints > 0
        bayesoptim.err('ObjectiveNargoutWrong', NumCoupledConstraints);
    else
        rethrow(msg);
    end
end
[Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
    = finishObjEvalParallel(Objective, ConstraintViolations, ErrorConstraintViolation, StartTic, NumCoupledConstraints);
end

function [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
    = finishObjEvalParallel(Objective, ConstraintViolations, ErrorConstraintViolation, StartTic, NumCoupledConstraints)
% Check constraint violation dimensions
if numel(ConstraintViolations) ~= NumCoupledConstraints
    bayesoptim.err('ObjectiveNumConstraintsWrong',...
        numel(ConstraintViolations), NumCoupledConstraints);
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

function val = iNanIfBad(val)
% val must be a finite real number. If not, set it to NaN.
if ~isscalar(val) || ~isnumeric(val) || ~isreal(val) || ~isfinite(val)
    val = NaN;
end
end

function ObjectiveNargout = nargoutFromNumConstraints(NumCoupledConstraints)
% We know nargout < 3
if NumCoupledConstraints > 0
    ObjectiveNargout = 2;
else
    ObjectiveNargout = 1;
end
end


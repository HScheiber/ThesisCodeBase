function [stop,options,optchanged]  = outputFcn_patternsearch_meltingpoint(optimvalues,options,flag,Settings)
%PSOUTPUTFCNTEMPLATE Template to write custom OutputFcn for PATTERNSEARCH.
%   [STOP,OPTIONS,OPTCHANGED] = PSOUTPUTFCNTEMPLATE(OPTIMVALUES,OPTIONS,FLAG) 
%   where OPTIMVALUES is a structure containing information about the state
%   of the optimization:
% 
%       PATTERNSEARCH:
%                   x: current point X
%           iteration: iteration count
%                fval: function value
%            meshsize: current mesh size
%           funccount: number of function evaluations
%              method: method used in last iteration
%              TolFun: tolerance on function value in last iteration
%                TolX: tolerance on X value in last iteration
%          nonlinineq: nonlinear inequality constraints at X
%            nonlineq: nonlinear equality constraints at X
% 
%       PARETOSEARCH:
%                   x: current Pareto set in parameter space
%                fval: current Pareto set in functional space
%           iteration: iteration count
%           funccount: number of function evaluations
%          nonlinineq: nonlinear inequality constraints at X
%            nonlineq: nonlinear equality constraints at X
%              volume: total volume of polytope defined by Pareto set
%     averagedistance: weighted average distance between Pareto points
%              spread: measure of how far apart Pareto points are from outliers in the Pareto set
%
%   OPTIONS: Options object used by PATTERNSEARCH.
%
%   FLAG: Current state in which OutPutFcn is called. Possible values are:
%         'init': initialization state 
%         'iter': iteration state
%    'interrupt': subproblem for nonlinear constraints state
%         'done': final state
% 		
%   STOP: A boolean to stop the algorithm.
%
%   OPTCHANGED: A boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, PARETOSEARCH, GA, OPTIMOPTIONS, SEARCHFCNTEMPLATE

%   Copyright 2003-2018 The MathWorks, Inc.

stop = false;
optchanged = false;

switch flag
    case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end
  
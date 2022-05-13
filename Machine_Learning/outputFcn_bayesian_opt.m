function stop = outputFcn_bayesian_opt(x,optimValues,state)
% outputFcn_global()
%
% OutputFun for optimizers (fminunc, fmincon etc),  saving intermediate results 
% in global variable "intermediate_data" for later access. 
%
% It is not supposed for live updates during optimization but for 
% later inspection, which is much more efficient. 
%
% Usage 
%   options = optimoptions( ... , 'OutputFcn',@outputFcn_global ); 
%   [XOpt,fval,exitflag,output] = fminunc(@fun, X0, options); 
%   intermediate_data(k).x 
%
% See also the supplied example file. 
%
% Last Changes
%   Daniel Frisch, ISAS, 10.2020: created example
%   Daniel Frisch, ISAS, 11.2019: improved documentation 
% Created
%   Daniel Frisch, ISAS, 10.2019 
%
stop = false;
global intermediate_data
switch state
  case 'init'
    intermediate_data = struct();
    intermediate_data.x = x;
    intermediate_data.optimValues = optimValues;
    intermediate_data.timerVal = tic;
  case 'iter'
    ind = length(intermediate_data)+1;
    intermediate_data(ind).x = x;
    intermediate_data(ind).optimValues = optimValues;
    intermediate_data(ind).timerVal = toc(intermediate_data(1).timerVal);
  case 'done'
    %
  otherwise
    error('wrong switch')
end

% Save to file
filename = 'intermediate_bayesian_opt.mat';
save(filename,'intermediate_data')
end
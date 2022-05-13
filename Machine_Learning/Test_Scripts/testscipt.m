clc
N = 300;
% - let's create 300 points with gaussian shape
x = linspace(0,10,N)';
b = [1; 3; 2; 2]; % model parameters
mf = @modelfun;
y = mf(b,x);
% - let's create search parameters
num(1) = optimizableVariable('base',    [1,9],'Type','real');
num(2) = optimizableVariable('height',  [1,9],'Type','real');
num(3) = optimizableVariable('location',[1,9],'Type','real');
num(4) = optimizableVariable('width',   [1,9],'Type','real');
% - define loss function
bf=@(num)bayesfun(num,x,y);
%  - call bayesopt with initial parameters close to model values
results = bayesopt_priv(bf,num, 'MaxObjectiveEvaluations',2000, 'isobj',true,...
    'AcquisitionFunctionName','lower-confidence-bound','ExplorationRatio',0.5,...
    'plotfcn','all','KernelFunction','matern52',...
    'OutputFcn',@assignInBase,'UseParallel',true,'GPActiveSetSize',2000);
b0 = table2array(results.bestPoint);
results.bestPoint
original_objective = exp(bayesfun(results.bestPoint,x,y)) - 1;  % Undo log1p


% let's see plot of model data and fit data
%------------------------
plot(x,y,'.b')
hold on
plot(x,mf(b0,x),'.r')
hold off
legend('y','fit')

function out = bayesfun(t,x,y)
    y1 = modelfun(table2array(t),x);
    out = sum((y1 - y).^2);
    out = log1p(out);
end

function out = modelfun(b,x)
    out = b(1) + b(2)*exp(-((b(3)-x)/b(4)).^2);
end
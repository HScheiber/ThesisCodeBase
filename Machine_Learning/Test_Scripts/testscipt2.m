clc
N = 2000;

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
bf = @(num)bayesfun(num,x,y);

%  - call bayesopt with initial parameters close to model values


options = optimoptions('surrogateopt','Display','iter',...
    'MaxFunctionEvaluations',500,...
    'MinSurrogatePoints',50,...
    'PlotFcn','surrogateoptplot',...
    'UseParallel',false,...
    'MinSampleDistance',1e-1);
ub = [10 10 10 10];
lb = [0 0 0 0];
[b0,fval,exitflag,output,trials] = surrogateopt(bf,lb,ub,options);



% options = optimoptions('ga','Display','iter',...
%     'CreationFcn','gacreationuniform',...
%     'PlotFcn',{'gaplotbestf' 'gaplotselection'},...
%     'UseParallel',false,'FunctionTolerance',1e-8);
% ub = [100 100 100 100];
% lb = [1 1 1 1];
% [b0,fval,exitflag,output,population,scores] = ga(bf,4,[],[],[],[],lb,ub,[],options);




% options = optimoptions('particleswarm','Display','iter',...
%     'PlotFcn','pswplotbestf',...
%     'UseParallel',false,'FunctionTolerance',1e-8);
% ub = [100 100 100 100];
% lb = [1 1 1 1];
% [b0,fval,exitflag,output] = particleswarm(bf,4,lb,ub,options);


original_objective = exp(bayesfun(b0,x,y)) - 1;  % Undo log1p


% let's see plot of model data and fit data
%------------------------
% plot(x,y,'-b')
% hold on
% plot(x,mf(b0,x),'.r')
% hold off
% legend('y','fit')

function out = bayesfun(t,x,y)
    y1 = modelfun(t,x);
    out = sum((y1 - y).^2);
    out = log1p(out);
end

function out = modelfun(b,x)
    out = b(1) + b(2)*exp(-((b(3)-x)/b(4)).^2);
end
N = 1e8;

mu = 1000;
Delta_T_min = 1;
step_size = 10;
gamma = 2;
T_UB = nan(1,N);
T_LB = nan(1,N);
MP =  nan(1,N);
for idx = 1:N
    T0 = mu + (rand - 0.5)*1000;
    [T_LB(idx),T_UB(idx),MP(idx)] = sim(mu,T0,step_size,gamma,Delta_T_min);
end

figure
histogram(T_LB,(mu-0.6250):0.025:mu,'Normalization','probability')
hold on
histogram(T_UB,mu:0.025:(mu+0.6250),'Normalization','probability')
xlim([mu-Delta_T_min*1.5 mu+Delta_T_min*1.5])
set(gca,'fontsize',25)
legend({'$T_{LB}$' '$T_{UB}$'},'interpreter','latex','fontsize',25)
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
ylabel('$\rho(T)$','interpreter','latex')
xlabel('$T$','interpreter','latex')
xticks(gca,[mu-0.6250 mu mu+0.6250])
xticklabels(gca,{'$T_{m}-0.625$' '$T_{m}$' '$T_{m}+0.625$'})
figure
histogram(MP,(mu-0.3125):0.025:(mu+0.3125),'Normalization','probability')
xlim([mu-Delta_T_min*0.75 mu+Delta_T_min*0.75])
set(gca,'fontsize',25)
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
ylabel('$\rho(T)$','interpreter','latex')
xlabel('$T$','interpreter','latex')
xticks(gca,[mu-0.3125 mu mu+0.3125])
xticklabels(gca,{'$T_{m}-0.3125$' '$T_{m}$' '$T_{m}+0.3125$'})

% When gamma neq 1
% x is the number of steps taken to reach the point where lb and ub are defined
% y is the number of iterations to finish after that point
% initial_size = step_size*(gamma^x);
% step_size*(gamma^x)/(2^y) < Delta_T_min
% 
% log(step_size/Delta_T_min) <  ylog(2) - xlog(gamma)
% 

% log(step_size*(gamma^x)/Delta_T_min) <  ylog(2)

width = nan(1,20);
for x = 1:20
    y = ceil((log(step_size) + x*log(gamma) - log(Delta_T_min))/log(2));
    width(x) = (step_size*(gamma^x))/(2^y);
end
unique(width)

% % When gamma = 1
% x = ceil(log(step_size/Delta_T_min)/log(2));
% width = step_size/(2^x)

% exportgraphics(gca,'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Manuscript_4\SI_Figures\Simulated_MP_Dist.pdf','Resolution',300)
% exportgraphics(gca,'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Manuscript_4\SI_Figures\Simulated_TLB_TUB_Dist.pdf','Resolution',300)


function [lb,ub,mp] = sim(mu,T0,step_size,gamma,Delta_T_min)
    lb = 0;
    ub = inf;
    T = T0;
    
    %ceil(log(1 - (1-gamma)*abs(T0 - mu)/step_size)/log(gamma));
    steps = 0;
    while ub - lb > Delta_T_min

        if lb > 0 && ~isinf(ub)
            if T > mu
                ub = T;
            elseif T < mu
                lb = T;
            else
                error('oh no')
            end
            T = (ub + lb)/2;
        else
            if T >= mu
                ub = T;
                if lb > 0 && ~isinf(ub)
                    T = (ub + lb)/2;
                else
                    T = T - step_size;
                end
            elseif T < mu
                lb = T;
                if lb > 0 && ~isinf(ub)
                    T = (ub + lb)/2;
                else
                    T = T + step_size;
                end
            end
        end
        step_size = step_size*gamma;
        steps = steps + 1;
    end
    mp = (ub + lb)/2;
end

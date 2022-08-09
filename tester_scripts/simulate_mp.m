N=1e6;

mu = 1000;
Delta_T_min = 1;
step_size = 1;
gamma = 3;
T_UB = nan(1,N);
T_LB = nan(1,N);
MP =  nan(1,N);
for idx = 1:N
    T0 = mu + (rand - 0.5)*1000;
    [T_LB(idx),T_UB(idx),MP(idx)] = sim(mu,T0,step_size,gamma,Delta_T_min);
end

figure
histogram(T_LB,100)
hold on
histogram(T_UB,100)
xlim([mu-Delta_T_min*1.5 mu+Delta_T_min*1.5])
figure
histogram(MP,100)
xlim([mu-Delta_T_min*1.5 mu+Delta_T_min*1.5])


% When gamma neq 1
% y is the number of steps taken to reach the point where lb and ub are
% defined
% x is the number of iterations to finish after that point
% initial_size = step_size*(gamma^y);
% step_size*(gamma^y)/(2^x) < Delta_T_min
% 
% log(step_size/Delta_T_min) <  xlog(2) - ylog(gamma)
% 

% log(step_size*(gamma^y)/Delta_T_min) <  xlog(2)

width = nan(1,20);
for y = 1:20
    x = ceil((log(step_size) + y*log(gamma) - log(Delta_T_min))/log(2));
    width(y) = (step_size*(gamma^y))/(2^x);
end
unique(width)

% % When gamma = 1
% x = ceil(log(step_size/Delta_T_min)/log(2));
% width = step_size/(2^x)

function [lb,ub,mp] = sim(mu,T0,step_size,gamma,Delta_T_min)
    lb = 0;
    ub = inf;
    T = T0;

    while ub - lb > Delta_T_min

        if lb > 0 && ~isinf(ub)
            ub - lb;
            if T > mu
                ub = T;
            elseif T < mu
                lb = T;
            else
                error('oh no')
            end
            T = (ub + lb)/2;
        else
            if T > mu
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
    end
    ub - lb;
    mp = (ub + lb)/2;

end
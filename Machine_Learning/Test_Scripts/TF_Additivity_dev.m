
epsilon = 10; % Well depth
gamma = 48/7-6; % Steepness parameter
r0 = 1; % Well minima location
r = 0:0.01:5;


%% Form A
k = (2 - (gamma/6) - (gamma/8));
f = epsilon.*((2/k).*exp(gamma.*(1 - (r./r0))) - (gamma/(6*k)).*((r0./r).^6) - (gamma/(8*k)).*((r0./r).^8) );

%Df = -(2*epsilon/k)*exp(gamma*(1 - r./r0))*(gamma/r0) + (gamma*epsilon*(r0^6)/k).*(1./(r.^7)) + (gamma*epsilon*(r0^8)/k).*(1./(r.^9));
%Df = -2*exp(gamma*(1 - r./r0))*(gamma/r0) + (gamma*(r0^6)).*(1./(r.^7)) + (gamma*(r0^8)).*(1./(r.^9));
%Df = -2*gamma*exp(gamma)*exp( - gamma.*r./r0) + gamma*((r0./r).^7) + gamma*((r0./r).^9);
%Df = -2*exp(gamma)*exp( - gamma.*r./r0) + ((r0./r).^7) + ((r0./r).^9);
% Df = -2*exp(gamma)*exp( - gamma.*r./r0) + ((r0./r).^7) + ((r0./r).^9);
% 
% syms gamma r0 r
% gamma = 48/7 + 0.1; % Steepness parameter
% r0 = 1; % Well minima location
% eqn = -2*exp(gamma)*exp( - gamma.*r./r0) + ((r0./r).^7) + ((r0./r).^9) == 0;
% solve(eqn,r)



plot(r,f)
ylim([-15 15])
yline(0,':k')

%% Form B
alpha = gamma/r0                                        % TF Steepness parameter
B = 48*epsilon*exp(gamma)/(48 - 4*gamma - 3*gamma)      % Repulsive prefactor
C6 = 4*gamma*(r0^6)*epsilon/(48 - 4*gamma - 3*gamma)    % 1/r6 dispersion
C8 = (3/4)*(r0^2)*C6                                    % 1/r8 dispersion

% r = 0:0.01:5;
% alpha = 1;                                       % TF Steepness parameter
% B = 10;                                          % Repulsive prefactor
% C6 = 1;                                          % 1/r6 dispersion
% C8 = 1;                                          % 1/r8 dispersion


g = B.*exp(-alpha.*r) - C6./(r.^6) - C8./(r.^8);


%% Back to form A

epsilon = B*exp(-sqrt( (3/4)*(alpha^2)*(C6/C8) ))*( 1 - (7/96)*alpha*sqrt(3*C6/C8) );
gamma   = (1/2)*alpha*sqrt(3*C6/C8);
r0      = 2*sqrt((1/3)*(C8/C6));

k = (2 - (gamma/6) - (gamma/8));
f = epsilon.*((2/k).*exp(gamma.*(1 - (r./r0))) - (gamma/(6*k)).*((r0./r).^6) - (gamma/(8*k)).*((r0./r).^8) );




%% Back to form B
% alpha = gamma/r0;                                       % TF Steepness parameter
% B = 48*epsilon*exp(gamma)/(48 - 4*gamma - 3*gamma);     % Repulsive prefactor
% C6 = 4*gamma*(r0^6)*epsilon/(48 - 4*gamma - 3*gamma);   % 1/r6 dispersion
% C8 = (3/4)*(r0^2)*C6;                                   % 1/r8 dispersion



plot(r,g,'-k',r,f,':r','linewidth',3)

ylim([-15 15])



%% Further development 1


B  = 10;
a  = 1;
C6 = 10;
C8 = 10;

r = 0:0.01:20;
f = B*exp(-a.*r) - C6./(r.^6) - C8./(r.^8);
figure;
plot(r,f)
ylim([-0.2 0.2])
xlim([0 20])

[fmax,ridx] = max(f);
rmax = r(ridx);

rmin = 


dr = -a*B*exp(-a.*r) + 6*C6./(r.^7) - 8*C8./(r.^9);

syms r
rmin = vpasolve(-a*B*exp(-a.*r) + 6*C6./(r.^7) - 8*C8./(r.^9) == 0, r);
fmin = double(B*exp(-a.*rmin) - C6./(rmin.^6) - C8./(rmin.^8));

hold on
scatter(rmin,fmin)
ylim([min(-0.2,fmin) 0.2])
xlim([0 40])


%% Further development 2

syms r
a_vec = 1:0.1:2;
S = zeros(size(a_vec));
for idx = 1:length(a_vec)
    B  = 1;
    a  = a_vec(idx);
    C6 = 1;
    C8 = 1;
    S(idx) = vpasolve(B*exp(-a*r) - C6/(r^6) - C8/(r^8), r);
end
plot(a_vec,S)

idx = 5;

B  = 1;
a  = a_vec(idx);
C6 = 1;
C8 = 1;

rmin = S(idx);
fmin = B*exp(-a.*rmin) - C6./(rmin.^6) - C8./(rmin.^8);


r = -10:0.01:20;
f = B*exp(-a.*r) - C6./(r.^6) - C8./(r.^8);
plot(r,f)
hold on
scatter(rmin,fmin)
ylim([min(-1,fmin) 1])
xlim([-10 20])

%% Further development 3

epsilon = 10; % Well depth
gamma = 48/7+0.01; % Steepness parameter
r0 = 1; % Well minima location
r = 0:0.001:5;


% Loose form
k = (2 - (gamma/6) - (gamma/8));
f = epsilon.*((2/k).*exp(gamma.*(1 - (r./r0))) - (gamma/(6*k)).*((r0./r).^6) - (gamma/(8*k)).*((r0./r).^8) );


% Tight form
alpha = gamma/r0;                                        % TF Steepness parameter
B = 48*epsilon*exp(gamma)/(48 - 4*gamma - 3*gamma);      % Repulsive prefactor
C = 4*gamma*(r0^6)*epsilon/(48 - 4*gamma - 3*gamma);    % 1/r6 dispersion
D = (3/4)*(r0^2)*C;                                    % 1/r8 dispersion

% r = 0:0.01:5;
% alpha = 1;                                       % TF Steepness parameter
% B = 10;                                          % Repulsive prefactor
% C6 = 1;                                          % 1/r6 dispersion
% C8 = 1;                                          % 1/r8 dispersion


g = B.*exp(-alpha.*r) - C./(r.^6) - D./(r.^8);
plot(r,g,'-k',r,f,':r','linewidth',3)
ylim([-15 15])

% Back to Loose form
gamma = (2/3)*sqrt(3*D/C)*alpha; % Steepness parameter
r0 = (2/3)*sqrt(3*D/C); % Well minima location
%epsilon = B*exp(-gamma)*(1 - (7/48)*gamma); % Well depth
%epsilon = B*exp(-sqrt( (4/3)*(alpha^2)*(D/C) ))*( 1 - (7/72)*alpha*sqrt(3*D/C) );
epsilon = 12*C/((128/2187)*alpha*((3*D/C)^(7/2))) - 7*C/((256/729)*((3*D/C)^3));
%epsilon = C*(12/(gamma*(r0^6)) - 7/(4*(r0^6)));

k = (2 - (gamma/6) - (gamma/8));
f = epsilon.*((2/k).*exp(gamma.*(1 - (r./r0))) - (gamma/(6*k)).*((r0./r).^6) - (gamma/(8*k)).*((r0./r).^8) );


plot(r,g,'-k',r,f,':r','linewidth',3)
ylim([-15 15])


%% Further development 4

epsilon = 10;
gamma_vec = 0.1:0.1:48/7;
r0_vec = 0.1:0.1:10;
Qmin = 38.8722633361816;

rt = zeros(length(gamma_vec),length(r0_vec));


for kdx = 1:length(gamma_vec)
    gamma = gamma_vec(kdx);
    k = 1/(2 - (gamma/6) - (gamma/8));
    for idx = 1:length(r0_vec)
        r0 = r0_vec(idx);
        ffun = @(r)2*epsilon.*k.*exp(gamma.*(1 - (r./r0))) - (epsilon*gamma*k/6).*((r0./r).^6) - (epsilon*gamma*k/8).*((r0./r).^8);
        fpfun = @(r)-2*epsilon*gamma*k.*exp(gamma.*(1 - (r./r0)))./r0 + epsilon*gamma*k*((r0.^6)./(r.^7)) + epsilon*gamma*k*((r0.^8)./(r.^9));

        rt(kdx,idx) = fzero(fpfun,2*r0+20/gamma);
        jdx = 1;
        while abs(ffun(rt(kdx,idx)) - epsilon) < 1e-10
            rt(kdx,idx) = fzero(fpfun,Qmin*r0/gamma + 10*jdx);
            jdx = jdx + 1;
        end
    end
end

[X,Y] = meshgrid(r0_vec,gamma_vec);

surf(X,Y,rt)
xlabel('r0')
ylabel('gamma')


surf(X,Y,Qmin*X./Y - rt)
xlabel('r0')
ylabel('gamma')

Q = 30:60;
fun = @(Q) sum((Q*(X./Y) - rt).^2,'all');
fQ = arrayfun(fun,Q);

plot(Q,fQ)

Qmin = fminsearch(fun,39)


%% Further development 5

lambertw(k,x)



gamma = 7.2;
epsilon = 10;
r0 = 1;
k = 1/(gamma - 6);

f1 = @(r) 6*epsilon*k*exp(gamma*(1 - r./r0)) - epsilon*gamma*k*((r0./r).^6);



alpha = gamma/r0;
B = 6*epsilon*k*exp(gamma);
C = epsilon*gamma*k*(r0^6);

f2 = @(r) B*exp(-alpha.*r) - C./(r.^6);

r = 0.1:0.001:5;

plot(r,f1(r),r,f2(r))
ylim([-10 10])


exp(gamma)/(gamma^7)



Q = B/(6*C*(alpha^6));

syms gamma A
assume(gamma,{'positive','real'})
assume(A,{'positive','real'})

solve(exp(gamma)/(gamma^7) - A == 0,gamma,'ReturnConditions',true)


-7*lambertw(0, -1/(7*Q^(1/7)))



-7*lambertw( 0, (-1/7)*( 6*C*(alpha^6)/B  )^(1/7) )




% No solution to loose form of exp-C6-C8
syms r A B C
assume(r,{'positive','real'})
assume(A,{'positive','real'})
assume(B,{'positive','real'})
assume(C,{'positive','real'})

soln = solve((r^7)*exp(-A*r) - (r^2)*B - C == 0,r,'ReturnConditions',true)



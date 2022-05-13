
syms A a C D r positive real

A = 50;
a = 1;
C = 5e3;
D = 5e3;

eqn = (r^9)*a*A*exp(-a*r) - 6*(r^2)*C - 8*D == 0;

%6*C/(r^7) + 8*D/(r^9) - B*A*exp(-B*r) == 0;

sol = solve(eqn, r,'Real',true)

pretty(sol)



r = 0:0.01:20;

A = 50;
a = 1;
C = 5e3;
D = 5e3;

eqn =      A.*exp(-a.*r) -    C./(r.^6) -    D./(r.^8);
deqn = -a*A.*exp(-a.*r) + 6.*C./(r.^7) + (8*D)./(r.^9);

hold on
plot(r,eqn,'r',r,deqn,'b')
ylim([-1e-1,1e-1])
yline(0,':k')


% syms f(r) A a C D
% f(r) = A.*exp(-a.*r) -    C./(r.^6) -    D./(r.^8);
% Df = diff(f,r)




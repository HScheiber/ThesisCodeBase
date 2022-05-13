r = 0.01:0.01:4;
C6 = 1;
d = 20;
s = 1;
R_vdW = 0.6; % nm
gamma = 10;

% Fermi-type damping function: d and s are empirical parameters
% R_vdW is sum of van der Waals radii
f_r = 1./(1 + exp(-d.*( (r./(s*R_vdW)) - 1)) );

% Chai and Head-Gordon damping (D3(0) dispersion): gamma and s are empirical parameters
f_r = 1./(1 + 6.*(r./(s.*R_vdW)).^(-gamma));


% Becke_Johnson damping function (D3(BJ))
f_r = (r.^(6))./(r.^(6) + R_vdW.^(6));

E_disp = -C6.*f_r./(r.^6);

plot(r,E_disp)
ylim([-10 5])
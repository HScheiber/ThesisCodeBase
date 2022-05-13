%% NaCl calculation of enthalpy change from experimental heat capacities
% Data source: NIST_JANAF Thermochemical Tables 4th Edition, page 789
% Low temperature Cp data: https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1964.0090

%% Solid at NaCl experimental MP
Tm = 1073.8;
T  = [0.00      2.50	5.00	7.50	10.00	15.00	20.00	25.00	30.00	40.00	50.00	60.00	70.00	80.00	90.00	100.00	140.00	160.00	200.00	250.00	260.00	298.15	300.00	400.00	500.00	600.00	700.00	800.00	900.00	1000.00	1073.80]; % K
Cp = [0.0000    0.0018	0.0149	0.0515	0.1265	0.4874	1.3071	2.7029	4.6442	9.6441	15.1754	20.4137	25.0036	28.9198	32.3005	34.9320	41.8316	44.1747	46.8730	49.2289	49.5595	50.5090	50.5430	52.3500	53.9400	55.4760	57.2040	59.3120	61.8690	64.8650	67.3710]; % J / mol K

T_H =   [0 100    200    298.15 300    400    500    600    700     800     900     1000    1073.8]; % K
H_exp = [0 1.5690 5.8040 10.611 10.704 15.853 21.168 26.639 32.270  38.092  44.147  50.479  55.358]; % kJ /mol
S_exp = [0 23.723 52.605 72.115 72.428 87.226 99.080 109.05 117.727 125.498 132.627 139.296 144.002]; % J / mol K

T_intp = 0:0.001:Tm; % K
Cp_intp = interp1(T,Cp,T_intp,'spline'); % J / mol K

hold on
axh = gca;
plot(axh,T_intp,Cp_intp,'linewidth',3)
scatter(axh,T,Cp,50,'filled','linewidth',2)
scatter(axh,Tm,Cp(T==Tm),100,'rd','filled','linewidth',2)
xlabel(axh,'$T$ [K]','interpreter','latex','FontSize',24)
ylabel(axh,'$C_{p}$ [$\textrm{J mol}^{-1}\textrm{K}^{-1}$]','interpreter','latex','FontSize',24)
set(axh,'TickLabelInterpreter','latex','FontSize',24)
title(axh,'Experimental $C_{P}$ vs $T$ Curve for NaCl','interpreter','latex','FontSize',32)

% Calculating change in enthalpy w.r.t. T = 0
Delta_H = cumtrapz(T_intp,Cp_intp)./1000; % kJ / mol

yyaxis(axh,'right')
plot(axh,T_intp,Delta_H)
scatter(axh,T_H,H_exp)
ylabel(axh,'$\Delta H^{\circ} = H^{\circ}(T) - H^{\circ}(0)$ [$\textrm{kJ mol}^{-1}$]','interpreter','latex','FontSize',24)

% Calculating entropy w.r.t. T = 0
Cp_over_T = Cp_intp./T_intp; % J / mol K^2
Cp_over_T(1) = 0;
Delta_S = cumtrapz(T_intp,Cp_over_T); % J / mol K

hold on
plot(T_intp,Delta_S)
scatter(T_H,S_exp)


%% Enthalpy and entropy of solid at the melting temperature
Enthalpy_mp_solid = trapz(T_intp,Cp_intp)./1000; % kJ / mol
Entropy_mp_solid = trapz(T_intp,Cp_over_T); % J / mol K

%% Solid above MP at T = 1289 K: JC MP
% To calculate entropy of solid above the melting point, we will
% extrapolate the heat capacity of the solid into the liquid, but do not
% add the entropy of melting

% Once you have the enthalpy, estimate the entropy of the solid by
% integrating the extrapolated heat capacities.
Tm = 1289;
T  = [0.00 2.50 5.00 7.50 10.00 15.00 20.00 25.00 30.00 40.00 50.00 60.00 70.00 80.00 90.00 100.00 140.00 160.00 200.00 250.00 260.00 298.15 300.00 400.00 500.00 600.00 700.00 800.00 900.00 1000.00 1073.80 1100.00 1200.00 1300.00 1400.00 1500.00]; % K
Cp = [0.0000 0.0018 0.0149 0.0515 0.1265 0.4874 1.3071 2.7029 4.6442 9.6441 15.1754 20.4137 25.0036 28.9198 32.3005 34.9320 41.8316 44.1747 46.8730 49.2289 49.5595 50.5090 50.5430 52.3500 53.9400 55.4760 57.2040 59.3120 61.8690 64.8650 67.3710 68.3250 71.9650 75.1450 77.5300 79.2450]; % J / mol K

T_H =   [0 100    200    298.15 300    400    500    600    700     800     900     1000    1073.8  1100.00	1200.00	1300.00	1400.00	1500.00]; % K
H_exp = [0 1.5690 5.8040 10.611 10.704 15.853 21.168 26.639 32.270  38.092  44.147  50.479  55.358 57.135	64.153	71.516	79.155	86.999]; % kJ /mol
S_exp = [0 23.723 52.605 72.115 72.428 87.226 99.080 109.05 117.727 125.498 132.627 139.296 144.002 145.637	151.742	157.633	163.293	168.704]; % J / mol K

T_intp = 0:0.001:Tm; % K
Cp_intp = interp1(T,Cp,T_intp,'spline'); % J / mol K

hold on
plot(T_intp,Cp_intp)
scatter(T,Cp)

Delta_H = cumtrapz(T_intp,Cp_intp)./1000; % kJ / mol

hold on
plot(T_intp,Delta_H)
scatter(T_H,H_exp)

% Calculating entropy w.r.t. T = 0
Cp_over_T = Cp_intp./T_intp; % J / mol K^2
Cp_over_T(1) = 0;
Delta_S = cumtrapz(T_intp,Cp_over_T); % J / mol K

hold on
plot(T_intp,Delta_S)
scatter(T_H,S_exp)

%% Enthalpy and entropy of solid at T = 1289 K: JC MP
Enthalpy_mp_solid = trapz(T_intp,Cp_intp)./1000; % kJ / mol
Entropy_mp_solid = trapz(T_intp,Cp_over_T); % J / mol K


% Calculating entropy change
Delta_S_Fusion_tm = 26.2230582976346; % J / mol K

Delta_S_liquid_to_tm = trapz(T_intp(T_intp >= 1073.8 & T_intp <= Tm),Cp_over_T(T_intp >= 1073.8 & T_intp <= Tm))



%% Solid above MP at T = 1081 K: TF MP
% To calculate entropy of solid above the melting point, we will
% extrapolate the heat capacity of the solid into the liquid, but do not
% add the entropy of melting

% Once you have the enthalpy, estimate the entropy of the solid by
% integrating the extrapolated heat capacities.
Tm = 1081;
T  = [0.00 2.50 5.00 7.50 10.00 15.00 20.00 25.00 30.00 40.00 50.00 60.00 70.00 80.00 90.00 100.00 140.00 160.00 200.00 250.00 260.00 298.15 300.00 400.00 500.00 600.00 700.00 800.00 900.00 1000.00 1073.80 1100.00 1200.00 1300.00 1400.00 1500.00]; % K
Cp = [0.0000 0.0018 0.0149 0.0515 0.1265 0.4874 1.3071 2.7029 4.6442 9.6441 15.1754 20.4137 25.0036 28.9198 32.3005 34.9320 41.8316 44.1747 46.8730 49.2289 49.5595 50.5090 50.5430 52.3500 53.9400 55.4760 57.2040 59.3120 61.8690 64.8650 67.3710 68.3250 71.9650 75.1450 77.5300 79.2450]; % J / mol K

T_H =   [0 100    200    298.15 300    400    500    600    700     800     900     1000    1073.8  1100.00	1200.00	1300.00	1400.00	1500.00]; % K
H_exp = [0 1.5690 5.8040 10.611 10.704 15.853 21.168 26.639 32.270  38.092  44.147  50.479  55.358 57.135	64.153	71.516	79.155	86.999]; % kJ /mol
S_exp = [0 23.723 52.605 72.115 72.428 87.226 99.080 109.05 117.727 125.498 132.627 139.296 144.002 145.637	151.742	157.633	163.293	168.704]; % J / mol K

T_intp = 0:0.001:Tm; % K
Cp_intp = interp1(T,Cp,T_intp,'spline'); % J / mol K

hold on
plot(T_intp,Cp_intp)
scatter(T,Cp)

Delta_H = cumtrapz(T_intp,Cp_intp)./1000; % kJ / mol

hold on
plot(T_intp,Delta_H)
scatter(T_H,H_exp)

% Calculating entropy w.r.t. T = 0
Cp_over_T = Cp_intp./T_intp; % J / mol K^2
Cp_over_T(1) = 0;
Delta_S = cumtrapz(T_intp,Cp_over_T); % J / mol K

hold on
plot(T_intp,Delta_S)
scatter(T_H,S_exp)

%% Enthalpy and entropy of solid at T = 1289 K: JC MP
Enthalpy_mp_solid = trapz(T_intp,Cp_intp)./1000; % kJ / mol
Entropy_mp_solid = trapz(T_intp,Cp_over_T); % J / mol K

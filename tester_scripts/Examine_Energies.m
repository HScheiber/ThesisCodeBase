WorkDir = 'C:\Users\Hayden\AppData\Local\Temp\tp8527fc7b_41e0_4886_bfd6_d3c880cd5871\Test6';

data = import_xvg('energy.xvg');

% @ s0 legend "Potential"
% @ s1 legend "Pressure"
% @ s2 legend "Volume"
% @ s3 legend "Enthalpy"
% @    xaxis  label "Time (ps)"
% @    yaxis  label "(kJ/mol), (bar), (nm^3)"

% %[ps] time constant for coupling T. Should be 20*Nsttcouple*timestep

figure
plot(data(:,1),data(:,2)./nmol_solid) % potential (kj/mol ion pairs)
hold on
plot(data(:,1),data(:,3))              % Pressure (bar)
plot(data(:,1),data(:,4)./nmol_solid) % Volume (nm^3 / ion pair)
plot(data(:,1),data(:,5)./nmol_solid) % Enthalpy ( kJ/mol ion pairs)

ylim([-2000 100])

midpoint = ceil(length(data)/2);


P = data(midpoint:end,3);
V = data(midpoint:end,4).*(10^3)./nmol_solid;
H = data(midpoint:end,5)./nmol_solid;

% Statistics on Pressure
SEM = std(P)/sqrt(length(P));               % Standard Error
ts = tinv([0.025  0.975],length(P)-1);      % T-Score
disp(['Pressure = ' num2str(mean(P),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])

% Statistics on Volume
SEM = std(V)/sqrt(length(V));               % Standard Error
ts = tinv([0.025  0.975],length(V)-1);      % T-Score
disp(['Volume = ' num2str(mean(V),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' A^3 / Ion pair'])

% Statistics on Enthalpy
SEM = std(H)/sqrt(length(H));               % Standard Error
ts = tinv([0.025  0.975],length(H)-1);      % T-Score
disp(['Enthalpy = ' num2str(mean(H),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' kJ / mol'])


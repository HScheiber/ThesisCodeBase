Energy_file = 'Test_MP.edr';
system(['wsl source ~/.bashrc; echo 9 11 15 18 0 ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o energy.xvg'])
%system(['wsl source ~/.bashrc; echo "4 0" ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o energy.xvg'])
%     En_xvg_file = fullfile(Settings.WorkDir,'Prep_Liq.xvg');
%     Data = import_xvg(En_xvg_file);
%     plot(Data(:,1),Data(:,2)./nmol_liquid) % Potential
%     ylim([-1000 1000])

Data = import_xvg('energy.xvg');

% @ s0 legend "Potential"
% @ s1 legend "Pressure"
% @ s2 legend "Volume"
% @ s3 legend "Enthalpy"
% @    xaxis  label "Time (ps)"
% @    yaxis  label "(kJ/mol), (bar), (nm^3)"

% %[ps] time constant for coupling T. Should be 20*Nsttcouple*timestep

nmol_solid = 2000;

hold on %figure
plot(Data(:,1),Data(:,2)./nmol_solid,'Linewidth',4,'Color','r') % potential (kj/mol ion pairs)
plot(Data(:,1),Data(:,2),'Linewidth',4,'Color','r') % potential (kj/mol ion pairs)
% std(Data(:,3)./Data(:,2))
% range(Data(:,3)./Data(:,2))
% min(Data(:,3)./Data(:,2))
% max(Data(:,3)./Data(:,2))

dev_from_init = Data(:,3)./Data(:,2) - Data(1,3)./Data(1,2);
max(dev_from_init)
min(dev_from_init)

ylim([-0.2 0.1])
hold on
plot(Data(:,1),Data(:,2)./max(abs(Data(:,2)))-1,'Linewidth',4,'Color','r')              % Pressure (bar)
xticks([])
yticks([])



plot(Data(:,1),Data(:,4)./nmol_solid) % Volume (nm^3 / ion pair)
plot(Data(:,1),Data(:,5)./nmol_solid) % Enthalpy ( kJ/mol ion pairs)

ylim([-2000 100])

startpoint = 50; %ceil(length(Data)/2);
nmol_solid = 1372;

T = Data(startpoint:end,2);
P = Data(startpoint:end,3);
V = Data(startpoint:end,4).*(10^3)./nmol_solid;
H = Data(startpoint:end,5)./nmol_solid;





startpoint = round(numel(Energy.Time)*0.2);
col = 'g';
NF = Energy.NF;
t = Energy.Time(startpoint:end);
T = Energy.Temperature(startpoint:end);
P = Energy.Pressure(startpoint:end);
%
if isfield(Energy,'Total_Energy')
    E = Energy.Total_Energy(startpoint:end)./NF;
    U = Energy.Potential(startpoint:end)./NF;
    hold on
    plot(t,E,'Linewidth',4,'Color',col)
else
    V = Energy.Volume(startpoint:end).*(10^3)./NF;
    H = Energy.Enthalpy(startpoint:end)./NF;
    hold on
    plot(t,T,'Linewidth',4,'Color',col)
end




%plot(conv(T,ones(3000,1)./3000),'Linewidth',4,'Color','r')


% Statistics on Temperature
SEM = std(T(1:1:end))/sqrt(length(T(1:1:end)));               % Standard Error
ts = tinv([0.025  0.975],length(T(1:1:end))-1);      % T-Score
disp(['Temperature = ' num2str(mean(T(1:1:end)),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' K'])

% Statistics on Pressure
SEM = std(P)/sqrt(length(P));               % Standard Error
ts = tinv([0.025  0.975],length(P)-1);      % T-Score
disp(['Pressure = ' num2str(mean(P),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])

% Statistics on Total Energy
SEM = std(E)/sqrt(length(E));               % Standard Error
ts = tinv([0.025  0.975],length(E)-1);      % T-Score
disp(['Total Energy = ' num2str(mean(E),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])

% Statistics on Potential
SEM = std(U)/sqrt(length(U));               % Standard Error
ts = tinv([0.025  0.975],length(U)-1);      % T-Score
disp(['Potential = ' num2str(mean(U),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])



% Statistics on Volume
SEM = std(V)/sqrt(length(V));               % Standard Error
ts = tinv([0.025  0.975],length(V)-1);      % T-Score
disp(['Volume = ' num2str(mean(V),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' A^3 / Ion pair'])

% Statistics on Enthalpy
SEM = std(H)/sqrt(length(H));               % Standard Error
ts = tinv([0.025  0.975],length(H)-1);      % T-Score
disp(['Enthalpy = ' num2str(mean(H),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' kJ / mol'])


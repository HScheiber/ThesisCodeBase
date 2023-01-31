Energy_file = 'ChkLiqStable_R_BF_Model_Alexandria_pol_NPT.edr';
%system(['wsl source ~/.bashrc; echo 0 ^| gmx_d energy -f ' windows2unix(Energy_file)],'-echo')
system(['wsl source ~/.bashrc; echo 6 9 11 15 18 0 ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o energy.xvg'])
%system(['wsl source ~/.bashrc; echo "4 0" ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o energy.xvg'])
%     En_xvg_file = fullfile(Settings.WorkDir,'Prep_Liq.xvg');
%     Data = import_xvg(En_xvg_file);
%     plot(Data(:,1),Data(:,2)./nmol_liquid) % Potential
%     ylim([-1000 1000])

Data = import_xvg('energy.xvg');

% @ s0 legend "Temperature"
% @ s1 legend "Pressure"
% @ s2 legend "Volume"
% @ s3 legend "Enthalpy"
% %[ps] time constant for coupling T. Should be 20*Nsttcouple*timestep

nmol_solid = 11664;

hold on %figure
%plot(Data(:,1),Data(:,2)./nmol_solid,'Linewidth',4,'Color','r') % potential (kj/mol ion pairs)
plot(Data(:,3),Data(:,2)./nmol_solid,'Linewidth',4,'Color','r') % T vs potential (kj/mol ion pairs)
plot(Data(:,3),Data(:,5)./nmol_solid,'Linewidth',4,'Color','b') % T vs Volume (nm^3/ion pair)

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

startpoint = 1; %ceil(length(Data)/2);
nmol_solid = 1000;
step = 1;
Energy.NF           = nmol_solid;
Energy.Time         = Data(startpoint:step:end,1);
%Energy.Potential    = Data(startpoint:step:end,2);
%Energy.Total_Energy       = Data(startpoint:step:end,3);
Energy.Temperature  = Data(startpoint:step:end,2);
Energy.Pressure     = Data(startpoint:step:end,3);
Energy.Volume       = Data(startpoint:step:end,4);
Energy.Enthalpy     = Data(startpoint:step:end,5);

% @ s0 legend "Temperature"
% @ s1 legend "Pressure"
% @ s2 legend "Volume"
% @ s3 legend "Enthalpy"



startpoint = 1;%round(numel(Energy.Time)*0.1);
step = 1;
col = 'c';
NF = Energy.NF;
t = Energy.Time(startpoint:step:end);
T = Energy.Temperature(startpoint:step:end);
P = Energy.Pressure(startpoint:step:end);
%
if isfield(Energy,'Total_Energy')
    E = Energy.Total_Energy(startpoint:step:end)./NF;
    U = Energy.Potential(startpoint:step:end)./NF;
    hold on
    plot(t,E,'Linewidth',4,'Color',col)
else
    V = Energy.Volume(startpoint:step:end).*(10^3)./NF;
    H = Energy.Enthalpy(startpoint:step:end)./NF;
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

% % Statistics on Total Energy
% SEM = std(E)/sqrt(length(E));               % Standard Error
% ts = tinv([0.025  0.975],length(E)-1);      % T-Score
% disp(['Total Energy = ' num2str(mean(E),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])
% 
% % Statistics on Potential
% SEM = std(U)/sqrt(length(U));               % Standard Error
% ts = tinv([0.025  0.975],length(U)-1);      % T-Score
% disp(['Potential = ' num2str(mean(U),'%.3f') ' -+' num2str(ts(2)*SEM,'%.3f') ' bar'])



% Statistics on Volume
SEM = std(V)/sqrt(length(V));               % Standard Error
ts = tinv([0.025  0.975],length(V)-1);      % T-Score
disp(['Volume = ' num2str(mean(V),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' A^3 / Ion pair'])

% Statistics on Enthalpy
SEM = std(H)/sqrt(length(H));               % Standard Error
ts = tinv([0.025  0.975],length(H)-1);      % T-Score
disp(['Enthalpy = ' num2str(mean(H),'%.5f') ' -+' num2str(ts(2)*SEM,'%.5f') ' kJ / mol'])


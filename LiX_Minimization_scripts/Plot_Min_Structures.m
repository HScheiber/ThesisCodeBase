Settings = Initialize_MD_Settings;
%Data = load(fullfile(Settings.home,'data','MX_JCTF_Min_Data.mat'),'Data').Data;
Data = load(fullfile(Settings.home,'data','MX_Alexandria_Min_Data.mat'),'Data').Data;
%Data = load(fullfile(Settings.home,'data','MX_Alexandria_PointQ_Min_Data.mat'),'Data').Data;

fs = 24;

Theory = 'Mie'; % {'BF' 'BH' 'JC' 'Mie'}
%Salts = {'NaCl'};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Structures = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' 'BetaBeO' 'CsCl'};
Structures_legend = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' '$\beta$-BeO' 'CsCl'};
Prop_of_intr = 'E';

X = 1:length(Salts);
Y = nan(length(Salts),length(Structures));
for idx = 1:length(Salts)
    Salt = Salts{idx};
    for jdx = 1:length(Structures)
        Structure = Structures{jdx};
        if isfield(Data.(Salt),Theory)
            Y(idx,jdx) = Data.(Salt).(Theory).Rocksalt.(Prop_of_intr) - Data.(Salt).(Theory).(Structure).(Prop_of_intr);
        end
    end
end

Colours = cbrewer('qual','Set3',length(Structures));
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'position',[0.1 0.1 0.85 0.85]);
hold(axh,'on')

colormap(Colours)
p = bar(axh,X,Y,1,'FaceColor','flat','EdgeColor','k','LineWidth',1);

for jdx = 1:length(Structures)
    p(jdx).CData = repmat(Colours(jdx,:),length(Salts),1);
end


for idx = X(1:end-1)
    xline(axh,idx+0.5,'--k','Linewidth',1)
end
    
switch Theory
    case 'TF'
        Thertxt = 'Tosi-Fumi';
    case 'JC'
        Thertxt = 'Jeung-Cheatham (SPC/E)';
    case 'JC3P'
        Thertxt = 'Jeung-Cheatham (TIP3P)';
    case 'JC4P'
        Thertxt = 'Jeung-Cheatham (TIP4P$_{\mathrm{EW}}$)';
    case 'JCSD'
        Thertxt = 'Smith-Dang';
    case 'BF'
        Thertxt = 'Alexandria (WBK)';
    otherwise
        Thertxt = Theory;
end

%title(['Comparison of Alkali Halide Crystal Energies: ' Thertxt ' Model.'],...
%    'Interpreter','latex','FontSize',fs)
xlim(axh,[0.5 length(X)+0.5])
ylim(axh,'padded')
xticks(axh,X);
xticklabels(axh,Salts);
set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
axh.XAxis.TickLength = [0 0];
axh.YGrid = 'on';
axh.YMinorGrid = 'On';
ylabel(axh,'$E_{\textrm{RS}} - E_{\textrm{Struc}}$ [kJ mol$^{-1}$]','Interpreter','latex');
legend(p,Structures_legend,'FontSize',fs,'Box','On','Interpreter','latex',...
    'NumColumns',4)
ylim(axh,[-60 10])
%xticks(axh,[]);
%xticklabels(axh,[]);

%exportgraphics(axh ,['Min_Structures_' Theory '.png'],'ContentType','image','BackgroundColor','none','Resolution',600)
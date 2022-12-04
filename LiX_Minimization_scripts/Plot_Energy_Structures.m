Settings = Initialize_MD_Settings;
%Data = load(fullfile(Settings.home,'data','MX_JCTF_Min_Data.mat'),'Data').Data;
%Data = load(fullfile(Settings.home,'data','MX_Alexandria_PointQ_Min_Data.mat'),'Data').Data;
Data = load(fullfile(Settings.home,'data','MX_Alexandria_Min_Data.mat'),'Data').Data;
%Data = load(fullfile(Settings.home,'data','MX_Alexandria_Polarized_Min_Data.mat'),'Data').Data;

fs = 24;

Theory = 'JC'; % {'BF' 'BH' 'JC' 'Mie'}
%Salts = {'NaCl'};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
% Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
%          'NaF' 'NaCl' 'NaBr' 'NaI' ...
%          'KF' 'KCl' 'KBr' 'KI' ...
%          'RbF' 'RbCl' 'RbBr' 'RbI' ...
%          'CsF' 'CsCl' 'CsBr' 'CsI'};
Structures_legend = {'Wurtzite' 'NiAs' 'Sphalerite' 'FiveFive' 'AntiNiAs' '$\beta$-BeO' 'CsCl'};
Prop_of_intr = 'E';

% Load experimental energy
refstructure = 'Wurtzite';
Exp = Load_Experimental_Data;
DFT = Load_Best_DFT_Data;

if false
    X = 1:length(Salts);
    Y = nan(length(Salts),1);
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        if isfield(Data.(Salt),Theory)
            if strcmp(Prop_of_intr,'a')
                Y(idx) = DFT.(Salt).(refstructure).a_zero - Data.(Salt).(Theory).(refstructure).(Prop_of_intr);
            else
                Y(idx) = DFT.(Salt).(refstructure).(Prop_of_intr) - Data.(Salt).(Theory).(refstructure).(Prop_of_intr);
            end
        end
    end
else
X = 1:length(Salts);
    Y = nan(length(Salts),1);
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        if isfield(Data.(Salt),Theory)
            Y(idx) = Exp.(Salt).Rocksalt.(Prop_of_intr) - Data.(Salt).(Theory).Rocksalt.(Prop_of_intr);
        end
    end
end

Colours = cbrewer('qual','Set3',length(Salts));
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'position',[0.1 0.1 0.85 0.85]);
hold(axh,'on')

colormap(Colours)
p = bar(axh,X,Y,1,'FaceColor','flat','EdgeColor','k','LineWidth',1);
p.CData = Colours;


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

% title(['Comparison of Crystal Energies: ' Thertxt ' Model.'],...
%    'Interpreter','latex','FontSize',fs)
xlim(axh,[0.5 length(X)+0.5])
ylim(axh,'padded')
xticks(axh,X);
xticklabels(axh,Salts);
set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
axh.XAxis.TickLength = [0 0];
axh.YGrid = 'on';
axh.YMinorGrid = 'On';

if strcmp(Prop_of_intr,'E')
    ylabel(axh,['$' Prop_of_intr '_{\textrm{RS}}$(Exp) - $' Prop_of_intr '_{\textrm{RS}}$ (Theory) [kJ mol$^{-1}$]'],'Interpreter','latex');
else
    ylabel(axh,['$' Prop_of_intr '_{\textrm{RS}}$(Exp) - $' Prop_of_intr '_{\textrm{RS}}$ (Theory) [\AA]'],'Interpreter','latex');
end
ylim(axh,'padded')
%xticks(axh,[]);
%xticklabels(axh,[]);

%exportgraphics(axh ,['Min_Structures_' Theory '.png'],'ContentType','image','BackgroundColor','none','Resolution',600)
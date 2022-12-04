% Predictive ionic Radii
Radius.Li = 73; % pm
Radius.Na = 105;
Radius.K = 137;
Radius.Rb = 152;
Radius.Cs = 175;

Radius.F = 129;
Radius.Cl = 179;
Radius.Br = 194;
Radius.I = 216;

% % Effective ionic Radii
% Radius.Li = 76; % pm
% Radius.Na = 102;
% Radius.K = 138;
% Radius.Rb = 152;
% Radius.Cs = 167;
% 
% Radius.F = 133;
% Radius.Cl = 181;
% Radius.Br = 196;
% Radius.I = 220;

% % Crystal ionic Radii
% Radius.Li = 90; % pm
% Radius.Na = 116;
% Radius.K = 152;
% Radius.Rb = 166;
% Radius.Cs = 167;
% 
% Radius.F = 119;
% Radius.Cl = 181;
% Radius.Br = 182;
% Radius.I = 206;

fs = 28;
RR_Ideal_Lines = [0.225 0.414 0.732 1.3661];
RR_Crossover_Lines = [0.326 0.717 1.3947];

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
     
N_Salts = numel(Salts);
RR = nan(N_Salts,1);
for idx = 1:N_Salts
    Salt = Salts{idx};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    
    RR(idx) = Radius.(Metal)/Radius.(Halide);
end

X = 1:N_Salts;

% Setup plot
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'position',[0.1 0.1 0.85 0.85]);
hold(axh,'on')
aobj = area(axh,[0 0 0; N_Salts+0.5 N_Salts+0.5 N_Salts+0.5],fliplr([1.3947-0.717 0.717-0.326 0.326; 1.3947-0.717 0.717-0.326 0.326]));

Cols = cbrewer('qual','Pastel1',3);
for idx = 1:length(aobj)
    aobj(idx).FaceColor = Cols(idx,:);
    aobj(idx).FaceAlpha = 0.5;
    aobj(idx).EdgeColor='none';
end

% % Draw lines separating the salts grouped by metal
% for XL = 0.5+(4:4:(length(X)-1))
%     xline(XL,':k','Linewidth',2)
% end

% Plot RRs
plot(axh,X,RR,'MarkerSize',15,'Marker','o','MarkerFaceColor','r',...
                'MarkerEdgeColor','k','linewidth',2,'Color','k','LineStyle','-');

% Plot lines of RR stability
for YL = RR_Crossover_Lines
    yline(YL,'--k','Linewidth',2)
end
for YL = RR_Ideal_Lines
    yline(YL,'--r','Linewidth',2)
end

xlim(axh,[0.5 length(X)+0.5])
ylim(axh,[0.225 1.3947])
xticks(axh,X);
xticklabels(axh,Salts);
yticks(axh,[0 0.326 0.717 1 1.3947])
set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex','Box','On')
ylabel(axh,'$r_{\textrm{M}^{+}}/r_{\textrm{X}^{-}}$','Interpreter','latex');
text(axh,10,1.2,'Cubic or Octahedral','FontSize',fs,'Interpreter','latex','HorizontalAlignment','center')
text(axh,10,0.5304,'Octahedral','FontSize',fs,'Interpreter','latex','HorizontalAlignment','center')
text(axh,10,0.27,'Tetrahedral','FontSize',fs,'Interpreter','latex','HorizontalAlignment','center')

exportgraphics(axh ,'Radius_Ratios.png','Resolution',300)
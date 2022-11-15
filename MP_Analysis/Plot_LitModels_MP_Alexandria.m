% Script to plot all melting points: experiment vs models

% Settings
fs = 32;
mkrsz = 14;
linestyles = {'-',':','-.','--'};
mkrshapes = {'s' 'd' '^' 'o'};
%DataSetName = 'Alexandria_Melting_Point_Data.mat';
DataSetName = 'Alexandria_Polarized_Melting_Point_Data.mat';
Models = {'BF' 'BH' 'JC' 'Mie'};
Legend_Labels = {'Experiment' 'WBK' 'BK' 'Lennard Jones [12-6]' 'Lennard Jones [8-6]'};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
if strcmp(DataSetName,'Alexandria_Melting_Point_Data.mat')
    titleg = 'Alexandria Models - No Polarization';
elseif strcmp(DataSetName,'Alexandria_Polarized_Melting_Point_Data.mat')
    titleg = 'Alexandria Models - With Polarization';
end
DefStructure = 'Rocksalt';
Structures = {'Rocksalt'};

Exclude = {};

Settings = Initialize_MD_Settings;

Exp = Load_Experimental_Data;
Lit = Load_Literature_Model_MPs;

SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Data = load(SaveDataDir).Data;

% This is the X coordinate to map to the salts
X = 1:length(Salts);
Y = nan(length(Salts),length(Models));
Y_rs = nan(length(Salts),length(Models));
Y_wz = nan(length(Salts),length(Models));
Y_cs = nan(length(Salts),length(Models));
Y_ff = nan(length(Salts),length(Models));
Y_Er = 2.*ones(length(Salts),length(Models));
Y_Exp = nan(length(Salts),1);
Y_Exp_Er = nan(length(Salts),1);
Y_indet = nan(length(Salts),length(Models),length(Structures));
N_indet = nan(length(Salts),length(Models),length(Structures));

% Y_Lit = nan(length(Salts),length(Models),length(Literature_Sets));
% Y_Lit_Er = nan(length(Salts),length(Models),length(Literature_Sets));

for idx = X
    Salt = Salts{idx};
    for jdx = 1:length(Models)
        Model = Models{jdx};
        
        % Grab MY newly-calculated melting points
        if isfield(Data,Salt)
            if isfield(Data.(Salt),Model)
                
                cMP = 0;
                for sdx = 1:length(Structures)
                    Structure = Structures{sdx};
                    if isfield(Data.(Salt).(Model),Structure)
                        if cMP < mean(Data.(Salt).(Model).(Structure).dT)
                            cMP = mean(Data.(Salt).(Model).(Structure).dT);
                            Y(idx,jdx) = cMP;
                        end
                        
                        if any(cellfun(@(x) all(strcmpi(x,{Salt Model Structure})),Exclude))
                            continue
                        end
                        if Data.(Salt).(Model).(Structure).Alt_Structure
                            continue
                        end
                         
                        switch Structure
                            case 'Rocksalt'
                                Y_rs(idx,jdx) = mean(Data.(Salt).(Model).(Structure).dT);
                                disp([Salt ' ' Model ' ' Structure ' : ' num2str(Y_rs(idx,jdx),'%.0f') ])
                            case 'Wurtzite'
                                Y_wz(idx,jdx) = mean(Data.(Salt).(Model).(Structure).dT);
                                disp([Salt ' ' Model ' ' Structure ' : ' num2str(Y_wz(idx,jdx),'%.0f') ])
                            case 'CsCl'
                                Y_cs(idx,jdx) = mean(Data.(Salt).(Model).(Structure).dT);
                                disp([Salt ' ' Model ' ' Structure ' : ' num2str(Y_cs(idx,jdx),'%.0f') ])
                            case 'FiveFive'
                                Y_ff(idx,jdx) = mean(Data.(Salt).(Model).(Structure).dT);
                                disp([Salt ' ' Model ' ' Structure ' : ' num2str(Y_ff(idx,jdx),'%.0f') ])
                        end
                        
                        MPs_idx = ~Data.(Salt).(Model).(Structure).Freeze_Trace & ...
                            ~Data.(Salt).(Model).(Structure).Melt_Trace;
                        Tms = Data.(Salt).(Model).(Structure).T_Trace(MPs_idx);
                        Fz = Data.(Salt).(Model).(Structure).T_Trace(logical(Data.(Salt).(Model).(Structure).Freeze_Trace));
                        Mt = Data.(Salt).(Model).(Structure).T_Trace(logical(Data.(Salt).(Model).(Structure).Melt_Trace));
                        Y_indet(idx,jdx,sdx) = diff([max([Fz 0]) min(Mt)]);
                        N_indet(idx,jdx,sdx) = length(Tms);
                    end
                end
            end
        end
        
%         % Grab literature melting points
%         Structure = DefStructure;
%         for kdx = 1:length(Literature_Sets)
%             LitSet = Literature_Sets{kdx};
%             
%             if isfield(Lit.(LitSet),Model)
%                 if isfield(Lit.(LitSet).(Model),Salt)
%                     if isfield(Lit.(LitSet).(Model).(Salt),Structure)
%                         Y_Lit(idx,jdx,kdx) = Lit.(LitSet).(Model).(Salt).(Structure).mp;
%                         Y_Lit_Er(idx,jdx,kdx) = Lit.(LitSet).(Model).(Salt).(Structure).dmp;
%                     end
%                 end
%             end
%         end
    end
    % Grab Experimental melting points
    Y_Exp(idx) = Exp.(Salt).mp;
    Y_Exp_Er(idx) = Exp.(Salt).dmp;
end

% Setup plot
p = gobjects(1,length(Models));
colours = max(min(cbrewer('qual','Set1',length(Models),'spline'),1),0);
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'position',[0.1 0.1 0.85 0.85]);
hold(axh,'on')

% Draw lines separating the salts grouped by metal
for XL = 0.5+(4:4:(length(X)-1))
    xline(XL,':k','Linewidth',2)
end

% Plot experimental data
plot(axh,X,Y_Exp,'linewidth',3,'Color','k')
p(1) = plot(axh,X(1:18),Y_Exp(1:18),'MarkerSize', mkrsz,'Marker','s','MarkerFaceColor','k',...
                'MarkerEdgeColor','k','linewidth',2,'Color','k','LineStyle','none');
plot(axh,X(19:end),Y_Exp(19:end),'MarkerSize', mkrsz,'Marker','o','MarkerFaceColor','k',...
                'MarkerEdgeColor','k','linewidth',2,'Color','k','LineStyle','none');

% Plot literature data
% for kdx = 1:length(Literature_Sets)
%     linestyle = linestyles{kdx};
%     for jdx = 1:length(Models)
%         col = colours(jdx,:);
%         mkrshape = mkrshapes{jdx};
%         
%         %plot(axh,X,Y_Lit(:,jdx,kdx),'linewidth',2,'Color',col,'LineStyle',linestyle)
%         errorbar(axh,X,Y_Lit(:,jdx,kdx),Y_Lit_Er(:,jdx,kdx),Y_Lit_Er(:,jdx,kdx),...
%             'linewidth',2,'Color',col,'LineStyle','none');
%         scatter(axh,X,Y_Lit(:,jdx,kdx),50,'Marker',mkrshape,'linewidth',1,'MarkerFaceColor',col,...
%             'MarkerEdgeColor','k')
%         
%     end
% end

% Plot my data
for jdx = 1:length(Models)
    col = colours(jdx,:);
    
    plot(axh,X,Y(:,jdx),'linewidth',3,'Color',col,'LineStyle','--');

    
    for sdx = 1:length(Structures)
        Structure = Structures{sdx};
        mkrshape = mkrshapes{sdx};
        switch Structure
            case 'Rocksalt'
                p(jdx+1) = plot(axh,X,Y_rs(:,jdx),'MarkerSize', mkrsz,'Marker',mkrshape,'MarkerFaceColor',col,...
                    'MarkerEdgeColor','k','linewidth',2,'Color',col,'LineStyle','none');
            case 'Wurtzite'
                plot(axh,X,Y_wz(:,jdx),'MarkerSize', mkrsz,'Marker',mkrshape,'MarkerFaceColor',col,...
                    'MarkerEdgeColor','k','linewidth',2,'Color',col,'LineStyle','none');
            case 'CsCl'
                plot(axh,X,Y_cs(:,jdx),'MarkerSize', mkrsz,'Marker',mkrshape,'MarkerFaceColor',col,...
                    'MarkerEdgeColor','k','linewidth',2,'Color',col,'LineStyle','none');
            case 'FiveFive'
                plot(axh,X,Y_ff(:,jdx),'MarkerSize', mkrsz,'Marker',mkrshape,'MarkerFaceColor',col,...
                    'MarkerEdgeColor','k','linewidth',2,'Color',col,'LineStyle','none');
        end
    end

end


xlim(axh,[0.5 length(X)+0.5])
ylim(axh,[400 1800])
xticks(axh,X);
xticklabels(axh,Salts);
set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
ylabel(axh,'T [K]','Interpreter','latex');
legend(p,Legend_Labels,'FontSize',fs,'Box','On','Interpreter','latex','NumColumns',2,'location','southeast')
grid(axh,'on')

title(axh,titleg,'FontSize',fs,'Interpreter','latex')

% exportgraphics(axh ,'C:\Users\Hayden\Documents\Patey_Lab\Thesis_Projects\Manuscript_4\Figures\MP_Alkali_Halides.pdf',...
%     'ContentType','vector','BackgroundColor','none')

% Initialize stuff
figtitle = 'Vibrational Analysis of Lithium Halides';
app.home = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase';
app.axpos = [0.065625 0.111111111111111 0.916145833333333 0.785046728971962];
app.FS = 24; % Font size
app.LW = 2; % Line width
app.figh = figure('WindowState','maximized','NumberTitle','off',...
	'Name',figtitle,'Visible','on');
app.axh = axes(app.figh,'Position',app.axpos);

% Turn off some toolbar options
set(0,'Showhidden','on')
ch = findall(app.figh,'tag','FigureToolBar');
UT = get(ch,'children');
set(UT([1 3:5]),'Visible','Off','Separator','Off')
set(UT(1),'Separator','On')

%-------------------------------------------------------------------------%
% Add crosshair togglebutton to toolbar:
%-------------------------------------------------------------------------%
% Read an image
[cdata,map] = imread(fullfile(app.home,'analysis','AppData','crosshair.gif'));

% Convert white pixels into a transparent background
map(map(:,1)+map(:,2)+map(:,3)==3) = NaN;

% Convert into 3D RGB-space
cdata = ind2rgb(cdata,map);

% Set the button icon
tbh = findall(app.figh,'Type','uitoolbar');
app.CrosshairButton = uitoggletool(tbh,'CData',cdata,...
    'Separator','off','Enable','off',...
    'TooltipString','Data Crosshair','OffCallback',@CrosshairOff,...
    'BusyAction','Cancel','tag','Exploration.Crosshair');

% Update crosshair button
%set(CrosshairButton,'OnCallback',{@crosshairfcn,X_Axis,Y_Axis,plot_extra,fs});
%CrosshairButton.Enable = 'on';
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

app.SaltLabel = uicontrol(app.figh,'Style','text','String',...
    'Salt(s)','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0.094866071428571,0.917168674698795,0.043154761904762,0.052116297857204]);
app.SaltBox = uicontrol(app.figh,'Style','listbox','String',...
    {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.142782738095238,0.899273104880582,0.053571428571429,0.096573208722741],...
    'Min',1,'Max',5);

app.StructureLabel = uicontrol(app.figh,'Style','text','String',...
    'Structures','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0.233333333333333,0.917168674698795,0.071354166666667,0.052116297857204]);
app.StructureBox = uicontrol(app.figh,'Style','listbox','String',...
    {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.3046875,0.899273104880582,0.074479166666667,0.096573208722741],...
    'Min',1,'Max',7);

app.TheoryLabel = uicontrol(app.figh,'Style','text','String',...
    'Theory','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0.420386904761904,0.917168674698795,0.050967261904762,0.052116297857204]);
app.TheoryBox = uicontrol(app.figh,'Style','listbox','String',...
    {'PBE-D3' 'PW1PW-D3' 'B3LYP-D3'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.4703125,0.899273104880582,0.090104166666667,0.096573208722741],...
    'Min',1,'Max',3);

app.CellSizeLabel = uicontrol(app.figh,'Style','text','String',...
    'Supercell Size','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0.601636904761905,0.917168674698795,0.107142857142857,0.052116297857204]);
app.CellSizeBox = uicontrol(app.figh,'Style','listbox','String',...
    {'1x1x1' '2x2x2' '3x3x3' '4x4x4'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.702157738095238,0.899273104880582,0.059300595238095,0.096573208722741],...
    'Min',1,'Max',4);

app.CalcTypeLabel = uicontrol(app.figh,'Style','text','String',...
    'Plot Type','FontSize',app.FS-5,'Units','Normalized','Position',...
    [0.808854166666666,0.917168674698795,0.070238095238095,0.052116297857204]);
app.CalcTypeBox = uicontrol(app.figh,'Style','listbox','String',...
    {'Total DOS' 'Metal DOS' 'Halide DOS' 'ZPE' 'PV vs T' 'Vib E vs T' 'S vs T' 'T vs V (Helmholtz)' 'T vs a (Helmholtz)' 'A vs T' 'G vs T (Harmonic)'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.87455357142857,0.899273104880582,0.088988095238096,0.096573208722741],...
    'Min',1,'Max',1);

app.PlotButton = uicontrol(app.figh,'Style','pushbutton','String',...
    {'Plot'},'FontSize',app.FS-5,'Units','Normalized',...
    'Position',[0.0171875,0.924195223260644,0.062499999999999,0.053997923156802],...
    'Callback',{@LoadPlots app});

set(app.axh,'box','on','TickLabelInterpreter','latex');
set(app.axh,'XMinorTick','on','YMinorTick','on','FontSize',app.FS);
xticks(app.axh,[]);
yticks(app.axh,[]);

% Load Plots function
function LoadPlots(~,~,app)
    % Clear previous plot
    cla(app.axh);
    legend(app.axh, 'off');
    drawnow;
    clearvars X Y Leg_txt
    ylabel(app.axh,'');
    xlabel(app.axh,'');
    set(app.axh, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
    set(app.axh, 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
    
    % Load Inputs
    Salts = app.SaltBox.String(app.SaltBox.Value);
    Structures = app.StructureBox.String(app.StructureBox.Value);
    Basis = 'pob-TZVP';
    Theories = app.TheoryBox.String(app.TheoryBox.Value);
    CellSizes = app.CellSizeBox.Value;
    CalcType = app.CalcTypeBox.String{app.CalcTypeBox.Value};

    % Load Data
    NSalt = length(Salts);
    NStrct = length(Structures);
    NTheory = length(Theories);
    NSuperCell = length(CellSizes);

    Data_dir = fullfile(app.home,'data','CRYSTAL');
    
    % Constants
    R = 0.0083144598; % Gas constant in kJ/K mol
    Ha_kJ = 4.35974e-21; % KiloJoules per Hartree
    NA = 6.0221409e23; % Molecules per mole
    Ha_kJ_NA = Ha_kJ*NA; % KiloJoule per Hartree per mole

    %% Load quantum data
    for i=1:NStrct
        Structure = Structures{i};
        Data.(Structure) = struct();

        for j=1:NSalt
            Salt = Salts{j};
            X = load(fullfile(Data_dir,[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
            Data.(Structure) = catstruct(Data.(Structure),X.Data);
            clearvars X
        end
    end


    %% Plot Data
    x = 1; % Keep track of successful data
    X = {};
    for i=1:NSalt
        Salt = Salts{i};
        for j=1:NStrct
            Structure = Structures{j};
            Labend = Label_replace(Structure);
            for k=1:NTheory
                Theory = strrep(Theories{k},'-','_');
                Label = [Salt Labend '_' Theory];
                DOI = Data.(Structure).(Label);

                % Select basis set
                Vib_DOI = [];
                for m=1:size(DOI,1)
                    if strcmp(DOI{m,1},Basis)
                        Vib_DOI = DOI(m,:);
                        break
                    end
                end
                if isempty(Vib_DOI)
                    continue
                end

                % Select supercell size
                for m=1:NSuperCell
                    SuperCellSize = CellSizes(m);
                    Vib_SS_DOI = [];
                    for n=1:size(Vib_DOI{1,5},1)
                        if Vib_DOI{1,5}{n} == SuperCellSize
                            Vib_SS_DOI = Vib_DOI{1,5}(n,:);
                            break
                        end
                    end
                    QHA_SS_DOI = [];
                    for n=1:size(Vib_DOI{1,6},1)
                        if Vib_DOI{1,6}{n} == SuperCellSize
                            QHA_SS_DOI = Vib_DOI{1,6}(n,:);
                            break
                        end
                    end
                    
                    if isempty(Vib_SS_DOI) && isempty(QHA_SS_DOI)
                        continue
                    end

                    NS = num2str(SuperCellSize);
                    Leg_txt{x} = [Salt ' ' Structure ' ' Theories{k} ' (' NS '$\times$' NS '$\times$' NS ')'];
                    
                    % Calc types
                    P = Vib_SS_DOI{1,5};
                    P_ind = (abs(P) < 0.101326);
                    
                    if strcmp(CalcType,'Total DOS')
                        X{x} = Vib_SS_DOI{1,3}(:,1); % wavenumber (cm-1)
                        Y{x} = Vib_SS_DOI{1,3}(:,2); % PDOS
                    elseif strcmp(CalcType,'Metal DOS')
                        X{x} = Vib_SS_DOI{1,17}(:,1); % wavenumber (cm-1)
                        Y{x} = Vib_SS_DOI{1,17}(:,2); % PDOS
                    elseif strcmp(CalcType,'Halide DOS')
                        X{x} = Vib_SS_DOI{1,18}(:,1); % wavenumber (cm-1)
                        Y{x} = Vib_SS_DOI{1,18}(:,2); % PDOS
                    elseif strcmp(CalcType,'ZPE')
                        X{x} = Vib_SS_DOI{1,2}; % Zero point energy (kJ/mol)
                    elseif strcmp(CalcType,'PV vs T')
                        X{x} = Vib_SS_DOI{1,6}(P_ind); % Temperature (K)
                        Y{x} = Vib_SS_DOI{1,9}(P_ind); % PV term (kJ/mol)
                    elseif strcmp(CalcType,'Vib E vs T')
                        X{x} = Vib_SS_DOI{1,6}(P_ind); % Temperature (K)
                        Y{x} = Vib_SS_DOI{1,4}(P_ind); % Vibrational E (kJ/mol)
                    elseif strcmp(CalcType,'S vs T')
                        X{x} = Vib_SS_DOI{1,6}(P_ind); % Temperature (K)
                        Y{x} = Vib_SS_DOI{1,7}(P_ind); % Entropy (J/mol K)
                    elseif strcmp(CalcType,'T vs V (Helmholtz)')
                        if ~isempty(QHA_SS_DOI)
                            X{x} = QHA_SS_DOI{1,20}(:,1); % Temperature (K)
                            Y{x} = QHA_SS_DOI{1,20}(:,2)/(SuperCellSize^3); % Volume (A^3)
                        else
                            continue
                        end
                    elseif strcmp(CalcType,'T vs a (Helmholtz)')
                        if ~isempty(QHA_SS_DOI)
                            X{x} = QHA_SS_DOI{1,21}(:,1); % Temperature (K)
                            if strcmp(Structure,'Rocksalt')
                                Y{x} = sqrt(2)*QHA_SS_DOI{1,21}(:,2)/SuperCellSize; % a (A)
                            else
                                Y{x} = QHA_SS_DOI{1,21}(:,2)/SuperCellSize; % a (A)
                            end
                        else
                            continue
                        end
                    elseif strcmp(CalcType,'A vs T')
                        X{x} = Vib_SS_DOI{1,6}(P_ind); % Temperature (K)
                        
                        % Get Entropy
                        Entropy = Vib_SS_DOI{1,7}(P_ind)/1000; % Entropy (kJ/mol K)
                        
                        % Get Lattice Energy
                        ELattice = min(Vib_DOI{1,2}{8}).*Ha_kJ_NA; % kJ/mol;
                        
                        % Get ZPE
                        ZPE = Vib_SS_DOI{1,2};
                        
                        % Get Vib E
                        EVib = Vib_SS_DOI{1,4}(P_ind);
                        
                        % Helmholtz free energy
                        Y{x} = ELattice + ZPE + EVib - Entropy.*X{x};
                    elseif strcmp(CalcType,'G vs T (Harmonic)')
                        T = Vib_SS_DOI{1,6}(P_ind); % Temperature (K)
                        
                        % Get Entropy
                        S = Vib_SS_DOI{1,7}(P_ind)/1000; % Entropy (kJ/mol K)
                        
                        % Get Lattice Energy
                        U = Vib_SS_DOI{1,16}.*Ha_kJ_NA; % kJ/mol;
                        
                        % Get ZPE
                        EV_ZPE = Vib_SS_DOI{1,2};
                        
                        % Get Vib E
                        EV = Vib_SS_DOI{1,4}(P_ind);
                        
                        % Get PV
                        PV = Vib_SS_DOI{1,9}(P_ind);
                        
                        TS = S.*T; % TS in kJ / mol
                        H = (U + EV_ZPE + EV) + PV; % Enthalpy kJ / mol
                        G = H - TS; % Gibbs free energy
                        
                        X{x} = T;
                        
                        % Helmholtz free energy
                        Y{x} = G;
                    end
                    x=x+1;
                end
            end
        end
    end

    if isempty(X)
        warndlg('No Data Selected.')
        return
    elseif strcmp(CalcType,'ZPE')
        set(app.axh,'Position',[0.065625 0.168224299065421 0.916145833333333 0.727933541017652])
        
        % Determine plot colours
        NC = length(X);
        if NC < 3
            NCol = 3;
        else
            NCol = NC;
        end
        
        Colours = cbrewer('qual','Set3',NCol);
        
        Leg_txt_Bar = cell(1,NC);
        for i = 1:NC
            Leg_txt_Bar{i} = ['\begin{tabular}{c}' strrep(Leg_txt{i},' ','\\ ') '\end{tabular}'];
            disp([strrep(Leg_txt{i},'$\times$','x') ' ZPE: ' num2str(X{i},'%6.8f') ' kJ/mol']);
        end
        
        p = bar(app.axh,cell2mat(X),'FaceColor','flat','Visible','on','BarWidth',1,...
            'LineWidth',2);
        ylim(app.axh,[0 max([X{:}]*1.05)])
        
        % Apply colours
        for i=1:NC
            p.CData(i,:) = Colours(i,:);
        end

        ylabel(app.axh,'Zero-Point Vibrational Energy (kJ mol$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
        set(app.axh,'box','on','TickLabelInterpreter','latex');
        set(app.axh,'XMinorTick','off','YMinorTick','on','FontSize',app.FS);
        
        % Modify x tick labels
        
        xticks(app.axh,1:length(X))
        xticklabels(app.axh,Leg_txt_Bar)

    else
        set(app.axh,'Position',app.axpos)
        % Determine plot colours
        NC = length(X);
        if NC < 3
            NCol = 3;
        else
            NCol = NC;
        end
        Colours = cbrewer('qual','Set1',NCol,'pchip');
        
        hold on
        if strcmp(CalcType,'Total DOS')
            for i=1:NC
                plot(app.axh,X{i},Y{i},'Color',Colours(i,:),...
                    'LineWidth',app.LW,'LineStyle','-',...
                    'Visible','on');
            end
            xlabel(app.axh,'Wavenumber (cm$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Total DOS (Normalized)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'Metal DOS')
            for i=1:NC
                plot(app.axh,X{i},Y{i},'Color',Colours(i,:),...
                    'LineWidth',app.LW,'LineStyle','-',...
                    'Visible','on');
            end
            xlabel(app.axh,'Wavenumber (cm$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Metal Contribution to DOS (Normalized)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'Halide DOS')
            for i=1:NC
                plot(app.axh,X{i},Y{i},'Color',Colours(i,:),...
                    'LineWidth',app.LW,'LineStyle','-',...
                    'Visible','on');
            end
            xlabel(app.axh,'Wavenumber (cm$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Halide Contribution to DOS (Normalized)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'PV vs T')
        	for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
        	end

            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'PV (at 1 atm) (kJ mol$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'Vib E vs T')
            for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
                if ~isempty(X{i})
                    E_298 = interp1(X{i},Y{i},298.15,'spline');
                    disp([strrep(Leg_txt{i},...
                        '$\times$','x') ' Vib E at 298.15K: ' num2str(E_298,'%6.8f') ' kJ/mol']);
                end
            end
            
            % Plot classical limit
            T_class = 0:0.01:1000;
            EV_Class = 3*2*R.*T_class; % classical limit
            plot(T_class,EV_Class,'-k','Visible','on','LineWidth',2)
            Leg_txt{end+1} = 'Classical Limit';
            
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Vibrational Contribution to Total E (kJ mol$^{-1}$)',...
                'fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'S vs T')
            for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
            end
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Entropy (J mol$^{-1}$ K$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'T vs V (Helmholtz)')
            for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
            end
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Volume Per Ion Pair (\AA$^{3}$)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'T vs a (Helmholtz)')
            for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
            end
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'a (\AA)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'A vs T')
            for i=1:NC
                scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                	'Marker','o','Visible','on');
            end
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Helmholz Free Energy (kJ mol$^{-1}$ K$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        elseif strcmp(CalcType,'G vs T (Harmonic)')
            for i=1:NC
                if ~isempty(X{i})
                    scatter(app.axh,X{i},Y{i},40,Colours(i,:),...
                        'Marker','o','Visible','on');
                    G_298 = interp1(X{i},Y{i},298.15,'spline');
                    disp([strrep(Leg_txt{i},'$\times$','x') ' G at 298.15K: ' num2str(G_298,'%6.8f') ' kJ/mol']);
                end
            end
            xlabel(app.axh,'Temperature (K)','fontsize',app.FS,'Interpreter','latex');
            ylabel(app.axh,'Gibbs Free Energy (kJ mol$^{-1}$ K$^{-1}$)','fontsize',app.FS,'Interpreter','latex');
            ylim(app.axh,[-inf inf])
        end
        set(app.axh,'box','on','TickLabelInterpreter','latex');
        set(app.axh,'XMinorTick','on','YMinorTick','on','FontSize',app.FS);
        
        app.Legendh = clickableLegend(app.axh,Leg_txt,'Interpreter','latex');
    end
    
    
    if ~strcmp(CalcType,'ZPE')
        ylim(app.axh,[-inf inf])
    end
    xlim(app.axh,[-inf inf])
end



home = find_home;
DataDir = fullfile(home,'data','Melting_Point_Data.mat');
Salt = 'NaCl';
SolStructure = 'Rocksalt';
Model = 'JC';
fs = 24;
% Set1: varying XY with fixed Z=10 nm, tmax = 1 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set2: varying Z with fixed XY = 6 nm, tmax = 1 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set3: droplet with varying R, tmax = 1 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set4: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 1 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set5: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set6: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 2 ns with $\Delta x_{\mathrm{liquid}}$ = 0.15
% Set7: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with $\Delta x_{\mathrm{liquid}}$ = 0.10
% Set8: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set9: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set10: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 10 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set11: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 10 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25

% Set12: Repeated Measurements [6x6x7 nm]   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set13: Repeated Measurements [6x6x4.5 nm] to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set14: Repeated Measurements 4000 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set15: Repeated Measurements 2000 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set16: Repeated Measurements 1024 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25

%% Comparisons of note:
Custom_Legend = true;
% Custom legend parameters:
% 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P' 'R_C' 'DispCorr' 'MaxTime'
% 'MeltFreezeThreshold' 'XYZ' 'XY' 'Z' 'N' 'Ewald' 'Fourier_Spacing'
% 'Ewald_rtol' 'ScaleCompressibility' 'PreEquilibrate'

% Sets = {'Set1'};
% Legend_Params = {'Z' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set2'};
% Legend_Params = {'XY' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set3'};
% Legend_Params = {'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set4' 'Set5'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set5' 'Set6'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set7' 'Set5' 'Set8'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set8' 'Set9' 'Set10'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set10' 'Set11'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set13' 'Set12' 'Set9'};
% Legend_Params = {'XYZ' 'N' 'MaxTime' 'MeltFreezeThreshold'};

% Sets = {'Set16' 'Set15' 'Set14'};
% Legend_Params = {'XYZ' 'N' 'R_C' 'DispCorr'};
% 
% Sets = {'Set16' 'Set15' 'Set14' 'Set13' 'Set12' 'Set9'};
% Legend_Params = {'XYZ' 'N' 'R_C' 'DispCorr'};
% 
% Sets = {'Set9' 'Set17' 'Set18' 'Set19'};
% Legend_Params = {'XYZ' 'N' 'R_C' 'DispCorr'};

% Sets = {'Set9' 'Set20' 'Set21' 'Set22' 'Set23' 'Set24'};
% Legend_Params = {'N' 'R_C' 'DispCorr' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P' 'Ewald'};

% Sets = {'Set17' 'Set19'};
% Legend_Params = {'XYZ' 'N' 'R_C' 'DispCorr'};

% Sets = {'Set9' 'Set20'};
% Legend_Params = {'N' 'R_C' 'DispCorr' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P' 'Ewald'};

% Sets = {'Set20' 'Set21'};
% Legend_Params = {'N' 'R_C' 'Ewald' 'Ewald_rtol' 'Fourier_Spacing'};

% Sets = {'Set21' 'Set26'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Barostat' 'ScaleCompressibility' 'PreEquilibrate'};

% Sets = {'Set9' 'Set25'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Barostat' 'PreEquilibrate'};

% Sets = {'Set27' 'Set28' 'Set29' 'Set30' 'Set31'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P'};

% Sets = {'Set32' 'Set33' 'Set34' 'Set35' 'Set36' 'Set37' 'Set38'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P'};

% Sets = {'Set39' 'Set40' 'Set41' 'Set42' 'Set43' 'Set44' 'Set45'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P'};

% Sets = {'Set46' 'Set47' 'Set48' 'Set49' 'Set50'};
% Legend_Params = {'N' 'R_C' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P'};

% Sets = {'Set51' 'Set55' 'Set56' 'Set57' 'Set58'};
% Legend_Params = {'N' 'R_C' 'DispCorr' 'Ewald_rtol' 'Fourier_Spacing'};


% Sets = {'Set51' 'Set52' 'Set55' 'Set56' 'Set57' 'Set58' 'Set59'};
% Legend_Params = {'N' 'MaxTime' 'MeltFreezeThreshold'};


Sets = {'Set32' 'Set33' 'Set34' 'Set35' 'Set36' 'Set37' 'Set38'};
Legend_Params = {'N' 'R_C' 'Thermostat' 'Tau_T' 'Barostat' 'Tau_P'};


N_Sets = length(Sets);
Colours = cbrewer('qual','Dark2',max(N_Sets,3));
Legend_txt = cell(1,N_Sets);
erh = gobjects(1,N_Sets);

figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
axh = axes(figh,'FontSize',24,'TickLabelInterpreter','latex',...
    'XMinorGrid','on','YLimitMethod','padded',...
    'YMinorGrid','on');
hold(axh,'on')

Exp = Load_Experimental_Data;
Tm_Exp = 0; %Exp.(Salt).mp;

Data = load(DataDir).Data.(Salt).(Model).(SolStructure);
AllJobs = fieldnames(Data);
for sidx = 1:N_Sets
    Set = Sets{sidx};
    col = Colours(sidx,:);
    
    [~,Dependent_var,Legend_txt{sidx},x_label,xlimits,x_ticks,Parameters] =  MP_Set_Info(Set);
    if Custom_Legend
        Legend_txt{sidx} = Parameters.(Legend_Params{1});
        for lidx = 2:length(Legend_Params)
            Legend_txt{sidx} = [Legend_txt{sidx} ', ' Parameters.(Legend_Params{lidx})];
        end
    end
    
    JobsOfInterest = AllJobs(contains(AllJobs,Set));
    X = zeros(length(JobsOfInterest),1);
    Y = zeros(length(JobsOfInterest),1);
    Y_pos_err = zeros(length(JobsOfInterest),1);
    Y_neg_err = zeros(length(JobsOfInterest),1);

    Xscatter = [];
    Yscatter = [];

    for jdx = 1:length(JobsOfInterest)

         if Data.(JobsOfInterest{jdx}).MP_confirmed
            MPs_idx = ~Data.(JobsOfInterest{jdx}).Freeze_Trace & ~Data.(JobsOfInterest{jdx}).Melt_Trace;
            Tms = Data.(JobsOfInterest{jdx}).T_Trace(MPs_idx);

            mean_Tm = mean([min(Tms) max(Tms)]);

            dTms = Tms - Tm_Exp;
            switch Dependent_var
                case 'Rep'
                    repno = regexp(JobsOfInterest{jdx},'Rep_([0-9]+)','once','tokens');
                    XPoints = repmat(str2double(repno{1}),1,length(dTms));
                    X(jdx) = str2double(repno{1});
                otherwise
                    XPoints = repmat(norm(Data.(JobsOfInterest{jdx}).(Dependent_var)),1,length(dTms));
                    X(jdx) = norm(Data.(JobsOfInterest{jdx}).(Dependent_var));
            end

            Y(jdx) = mean([min(dTms) max(dTms)]);

            for kdx = 1:length(Tms)
                Xscatter(end+1) = XPoints(kdx);
                Yscatter(end+1) = dTms(kdx);
            end

            Tfr_idx = Data.(JobsOfInterest{jdx}).T_Freeze_Trace < min(Tms);
            Tml_idx = Data.(JobsOfInterest{jdx}).T_Melt_Trace > max(Tms);

            T_freeze_max = max( Data.(JobsOfInterest{jdx}).T_Freeze_Trace(Tfr_idx) );
            T_melt_min = min( Data.(JobsOfInterest{jdx}).T_Melt_Trace(Tml_idx) );

            Y_pos_err(jdx) = T_melt_min - mean_Tm;
            Y_neg_err(jdx) = mean_Tm - T_freeze_max;

         else
            switch Dependent_var
                case 'Rep'
                    repno = regexp(JobsOfInterest{jdx},'Rep_([0-9]+)','once','tokens');
                    X(jdx) = str2double(repno{1});
                otherwise
                    X(jdx) = norm(Data.(JobsOfInterest{jdx}).(Dependent_var));
            end

            Tm = mean(Data.(JobsOfInterest{jdx}).dT);
            Y(jdx) = Tm - Tm_Exp;

            Y_pos_err(jdx) = Data.(JobsOfInterest{jdx}).dT(2) - Tm;
            Y_neg_err(jdx) = Tm - Data.(JobsOfInterest{jdx}).dT(1);
         end
    end

    [X,sortidx] = sort(X);
    Y = Y(sortidx);
    Y_pos_err = Y_pos_err(sortidx);
    Y_neg_err = Y_neg_err(sortidx);
    
    hold(axh,'on')
    erh(sidx) = errorbar(axh,X,Y,Y_neg_err,Y_pos_err,'linewidth',2,'Color',col);
    scatter(axh,Xscatter,Yscatter,50,'Marker','o','linewidth',1,'MarkerFaceColor','r',...
        'MarkerEdgeColor',col)
    set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
    xlim(axh,xlimits)
    ylim(axh,'padded')
    xticks(axh,x_ticks);
    xlabel(axh,x_label,'Interpreter','latex');
end
    
ylabel(axh,'$\Delta T_{m}$ [K]','Interpreter','latex')
grid(axh,'minor')
legend(axh,erh,Legend_txt,'location','Southeast','Interpreter','latex');


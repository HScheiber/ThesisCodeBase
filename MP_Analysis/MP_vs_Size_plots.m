home = find_home;
DataDir = fullfile(home,'data','Melting_Point_Data.mat');
Salt = 'NaCl';
SolStructure = 'Rocksalt';
% Set1: varying XY with fixed Z=10 nm, tmax = 1 ns with MeltFreezeThreshold = 0.15
% Set2: varying Z with fixed XY = 6 nm, tmax = 1 ns with MeltFreezeThreshold = 0.15
% Set3: droplet with varying R, tmax = 1 ns with MeltFreezeThreshold = 0.15
% Set4: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 1 ns with MeltFreezeThreshold = 0.15
% Set5: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with MeltFreezeThreshold = 0.15
% Set6: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 2 ns with MeltFreezeThreshold = 0.15
% Set7: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with MeltFreezeThreshold = 0.10
% Set8: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 2 ns with MeltFreezeThreshold = 0.25

% Set9: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Set10: Repeated Measurements [6x6x10 nm] to check for consistency at tmax = 10 ns with MeltFreezeThreshold = 0.25
% Set11: Repeated Measurements [8x8x10 nm] to check for consistency at tmax = 10 ns with MeltFreezeThreshold = 0.25

% Set12: Repeated Measurements [6x6x7 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25
% Set13: Repeated Measurements [6x6x4.5 nm] to check for consistency at tmax = 5 ns with MeltFreezeThreshold = 0.25

% Set12: Repeated Measurements [6x6x7 nm]   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set13: Repeated Measurements [6x6x4.5 nm] to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set14: Repeated Measurements 4000 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set15: Repeated Measurements 2000 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25
% Set16: Repeated Measurements 1024 atoms   to check for consistency at tmax = 5 ns with $\Delta x_{\mathrm{liquid}}$ = 0.25

Set = 'Set1';
fs = 24;

figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
Exp = Load_Experimental_Data;
Tm_Exp = Exp.(Salt).mp;

[Models,Dependent_var,Legend_txt,x_label,xlimits,x_ticks,Parameters] =  MP_Set_Info(Set);

N = length(Models);
t = tiledlayout(figh,N,1,'TileSpacing','compact');
Colours = cbrewer('qual','Set1',max(3,N));

title(t, ['\bf{' Legend_txt '}'],'Fontsize',fs+4,'Interpreter','latex')
for idx = 1:N
    
    axh = nexttile(t);
    hold(axh,'on')
    
    Model = Models{idx};
    title(axh,[Salt ' ' Model],'Fontsize',fs,'Interpreter','latex')
    
    Data = load(DataDir).Data.(Salt).(Model).(SolStructure);

    AllJobs = fieldnames(Data);
    JobsOfInterest = AllJobs(contains(AllJobs,[Set '_']));
    
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
    
    %scatter(axh,X,Y,100,'MarkerEdgeColor','k','MarkerFaceColor',Colours(idx,:),...
    %    'linewidth',2);
%     erh = errorbar(axh,X,Y,Y_neg_err,Y_pos_err,'linewidth',2,'Marker','o',...
%         'MarkerSize',10,'MarkerFaceColor',Colours(idx,:),...
%         'MarkerEdgeColor','k','Color',Colours(idx,:));
    erh = errorbar(axh,X,Y,Y_neg_err,Y_pos_err,'linewidth',2,'Color','k');
    scatter(axh,Xscatter,Yscatter,50,Colours(idx,:),'Marker','o','linewidth',1,'MarkerFaceColor',Colours(idx,:),...
        'MarkerEdgeColor','k')
    set(axh,'FontSize',fs,'Box','On','TickLabelInterpreter','latex')
    xlim(axh,xlimits)
    ylim(axh,'padded')
    xticks(axh,x_ticks);
    if idx < N
        xticklabels(axh,[])
    else
        xlabel(axh,x_label,'Interpreter','latex');
    end
    ylabel(axh,'$\Delta T_{m}$ [K]','Interpreter','latex')
    grid(axh,'minor')
end

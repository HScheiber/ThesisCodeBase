%% Plot result to visualize
Settings = Initialize_LiX_BO_Settings;
Settings.Salt = 'LiF';
Settings.Theory = 'BH';
Settings.SigmaEpsilon = true;
Settings.Comb_rule = 'Hogervorst';
N_chk = 10000;

Parinfo = bayesopt_params(Settings);
ParNames = {Parinfo.Name};
N_Vars = numel(ParNames);
ParTransforms = {Parinfo.Transform};
ParRanges = cell(1,N_Vars);
for idx = 1:N_Vars
    ParRanges{idx} = Parinfo(idx).Range;
end

% Build random table within parameter range
Param = table('Size',[N_chk N_Vars],'VariableTypes',repmat({'double'},1,N_Vars),'VariableNames',ParNames);
for idx = 1:N_Vars    
    if strcmp(ParTransforms{idx},'log')
        a = log(ParRanges{idx}(1));
        b = log(ParRanges{idx}(2));
        r = a + (b-a).*rand(N_chk,1);
        Param{:,idx} = exp(r);
    else
        a = ParRanges{idx}(1);
        b = ParRanges{idx}(2);
        Param{:,idx} = a + (b-a).*rand(N_chk,1);
    end
end

tf = LiX_Constraint_Fcn(Settings,Param);

figh = figure('WindowState','Maximized');
T = tiledlayout(figh,N_Vars,N_Vars,'TileSpacing','tight','padding','tight');
for idx = 1:N_Vars
    Y_VarName = ParNames{idx};
    Y_Data = Param.(Y_VarName);
    Y_Transform = ParTransforms{idx};
    Y_Range = ParRanges{idx};
    
    for jdx = 1:N_Vars
        ax = nexttile(T);
        hold(ax,'on')
        if jdx > idx
            set(ax,'visible','off')
            continue
        end
        if idx == jdx % Histogram plot
            X_feas = Y_Data(tf);
            X_infeas = Y_Data(~tf);
            
            if strcmp(Y_Transform,'log')
                [h_infeas, edges] = histcounts(log(X_infeas),30);
                [h_feas, ~] = histcounts(log(X_feas),edges);
                centers = exp((edges(1:end-1)+edges(2:end))/2);
                area(ax,centers, h_feas./(h_infeas+h_feas), 'FaceColor','r','Edgecolor','k','LineWidth',1)
                set(ax, 'XScale', 'log')
            else
                [h_infeas, edges] = histcounts(X_infeas,30);
                [h_feas, ~] = histcounts(X_feas,edges);
                centers = (edges(1:end-1)+edges(2:end))/2;
                area(ax,centers, h_feas./(h_infeas+h_feas), 'FaceColor','r','Edgecolor','k','LineWidth',1)
            end
            xlim(ax,Y_Range);
            set(ax,'Yticklabel',[])
            if idx == 1
                ylabel(ax,par_name_map(Y_VarName),'interpreter','latex');
                set(ax,'Xticklabel',[])
            elseif idx ~= N_Vars
                set(ax,'Xticklabel',[])
            else
                xlabel(ax,par_name_map(Y_VarName),'interpreter','latex');
            end
        else % pair correlation plots
            X_VarName = ParNames{jdx};
            X_Data = Param.(X_VarName);
            X_Transform = ParTransforms{jdx};
            X_Range = ParRanges{jdx};
            scatter(ax,X_Data(~tf),Y_Data(~tf),20,'k','x')
            scatter(ax,X_Data(tf),Y_Data(tf),20,'r','filled','linewidth',1,'MarkerEdgeColor','k')
            if strcmp(X_Transform,'log')
                set(ax, 'XScale', 'log')
            end
            if strcmp(Y_Transform,'log')
                set(ax, 'YScale', 'log')
            end
            xlim(ax,X_Range);
            ylim(ax,Y_Range);
            
            if jdx == 1
                ylabel(ax,par_name_map(Y_VarName),'interpreter','latex');
            else
                set(ax,'Yticklabel',[])
            end
            if idx == N_Vars
                xlabel(ax,par_name_map(X_VarName),'interpreter','latex');
            else
                set(ax,'Xticklabel',[])
            end
        end
    end
end
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
Theory = 'BH';
ModelID = 'MG';
BestOnly = false;
SelectOnly = [];
Reps = [1:5];
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';
PlotType = 'lj';
lw = 2;
fs = 28;
xmax = 6; % Angstrom
ylims = [-50 50];

% Find reps of models
N_Salts = numel(Salts);
Models = cell(1,length(Salts));
for idx = 1:N_Salts
    files = {dir(fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' ModelID '*_data.mat'])).name};
    Models{idx} = unique(regexp(files,[ModelID '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models{idx});
    Models{idx} = Models{idx}(nonempty_idx);
end
Models = unique(horzcat(Models{:}));
N_Reps = numel(Reps);
Models_of_Interest = cell(1,N_Reps);
for idx = 1:N_Reps
    Models_of_Interest{idx} = [ModelID num2str(Reps(idx))];
end
Models = intersect(Models_of_Interest,Models);

% Make sure the models are sorted
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end

% Instantiate figure
figh = figure('WindowState','Maximized');
T = tiledlayout(figh,3,4,'TileSpacing','compact','padding','tight');

ints = {'MM' 'XX' 'MX'};
N_Models = numel(Models);
Cols = cbrewer('qual','Set1',N_Models,'spline');
Cols = max(min(Cols,1),0);
Total_loss = nan(N_Salts,N_Models);
U = struct;
legtext= struct;
for idx = 1:N_Salts
    Salt = Salts{idx};
    for iidx = 1:N_Models
        Model = Models{iidx};

        % Find the fully optimized model
        dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);

        if ~isfile(dat_file)
            disp(['Could not load results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end

        % Load data
        try
            data = load(dat_file).full_data;

            if isfield(data,'secondary_result')
                optimvals = nan(1,length(data.secondary_result));
                for jdx = 1:length(data.secondary_result)
                    optimvals(jdx) = [data.secondary_result(jdx).optimValues.fval];
                end
                Total_loss(idx,iidx) = min(optimvals);
                params = data.full_opt_point;
            else
                [Total_loss(idx,iidx),midx] = min(data.bayesopt_results.ObjectiveTrace);
                params = data.bayesopt_results.XTrace{midx,:};
            end
        catch
            disp(['Could not obtain objective function data for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end

        Settings = data.Settings;
        Ptypes = bayesopt_params(Settings);
        Param_names = {Ptypes.Name};

        Param = table();
        for jdx = 1:numel(Param_names)
            Param.(Param_names{jdx}) = params(jdx);
        end
        Settings = Get_Scaling_Params(Settings,Param);
        Settings.MDP.vdw_modifier = 'None';
        Settings.Table_StepSize = 0.0005; % nm
        Settings.Table_Length = 4; % nm
        Settings_WBK = Settings;
        Settings_WBK.S = Init_Scaling_Object;

        if Settings.Additivity
            Param.r0_MX = (Param.r0_MM + Param.r0_XX)/2;
            Param.epsilon_MX = sqrt(Param.epsilon_MM * Param.epsilon_XX);
            Param.gamma_MM = Param.gamma_MX;
            Param.gamma_XX = Param.gamma_MX;
            Param_names{end+1} = 'r0_MX';
            Param_names{end+1} = 'epsilon_MX';
            Param_names{end+1} = 'gamma_MM';
            Param_names{end+1} = 'gamma_XX';
        end

        for kdx = 1:numel(ints)
            int = ints{kdx};
            r0_BH = Param.(['r0_' int]);
            gamma_BH = Param.(['gamma_' int]);
            if gamma_BH < 6
                epsilon_BH = -Param.(['epsilon_' int]);
            else
                epsilon_BH = Param.(['epsilon_' int]);
            end
            k = 1/(gamma_BH - 6);

            f_BH = @(r) 6*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) - epsilon_BH*gamma_BH*k*((r0_BH./r).^6);
            df_BH = @(r) -6*(gamma_BH/r0_BH)*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) + 6*epsilon_BH*gamma_BH*k*(r0_BH^6)./(r.^7);
            df_BH_neg = @(r) 6*(gamma_BH/r0_BH)*epsilon_BH*k*exp(gamma_BH*(1 - r./r0_BH)) - 6*epsilon_BH*gamma_BH*k*(r0_BH^6)./(r.^7);

            % Find minima and inflection points
            if gamma_BH < 7 % find the well minimum
                sigma_WBK = fzero(df_BH,r0_BH+1);
                eps_WBK = abs(f_BH(sigma_WBK));
            else % In this domain, r0 and epsilon describe the well minima
                sigma_WBK = r0_BH;
                eps_WBK = abs(epsilon_BH);
            end

            % find the inflection points
            r_infl = fminsearchbnd(df_BH,0.99*sigma_WBK,0,sigma_WBK);
            r_infl_outer = fminsearchbnd(df_BH_neg,1.01*sigma_WBK,sigma_WBK,10);

            % Define the WBK function with unknown gamma
            f_WBK_gam = @(r,gamma_WBK) (2*eps_WBK./(1 - (3./(gamma_WBK+3)))).*((sigma_WBK^6)./((sigma_WBK^6) + (r.^6))).*((3./(gamma_WBK+3)).*exp(gamma_WBK.*(1 - (r./sigma_WBK))) - 1);

            % Domain to compare BK and WBK functions
            r_check = r_infl:0.01:r_infl_outer;

            % Define a new function that computes the squared difference between BK and WBK functions for a given input gamma
            diff_BH_WBK = @(gamma_WBK) sum((f_BH(r_check) - f_WBK_gam(r_check,gamma_WBK)).^2);

            % Optimize the value of WBK gamma
            gamma_WBK = fminsearchbnd(diff_BH_WBK,gamma_BH,0,100);

            disp(repmat('*',1,50))
            disp([Salt ', ' Theory ', Model ' Model  '. Interaction: ' int '.'])
            disp(repmat('*',1,50))
            disp(['r0_' int ':' sprintf('\t\t') num2str(r0_BH) ' nm']);
            disp(['epsilon_' int ':' sprintf('\t') num2str(epsilon_BH) ' kJ/mol']);
            disp(['gamma_' int ':' sprintf('\t') num2str(gamma_BH)]);
            disp(repmat('*',1,50))
            disp('Best fit WBK Parameters:')
            disp(repmat('*',1,50))
            disp(['sigma_' int ':' sprintf('\t') num2str(sigma_WBK) ' nm']);
            disp(['epsilon_' int ':' sprintf('\t') num2str(eps_WBK) ' kJ/mol']);
            disp(['gamma_' int ':' sprintf('\t') num2str(gamma_WBK)]);
            Settings_WBK.S.S.(int) = sigma_WBK;
            Settings_WBK.S.E.(int) = eps_WBK;
            Settings_WBK.S.G.(int) = gamma_WBK;
            legtext.(Salt).(Model).(int) = [Model ' $\gamma = ' num2str(gamma_WBK,'%.2f') '$, ' ...
                '$\epsilon =$ ' num2str(eps_WBK,'%.2E') ' kJ/mol, $\sigma = ' num2str(sigma_WBK,'%.2f') '$ nm'];
        end
        U_BH.(Salt).(Model) = BH_Potential_Generator(Settings);
        U_WBK.(Salt).(Model) = BF_Potential_Generator(Settings_WBK);
        
    end
end

% Plot everything
for idx = 1:numel(ints)
    int = ints{idx};
    for jdx = 1:N_Salts
        Salt = Salts{jdx};
        ax = nexttile(T); % New tile
        hold(ax,'on')
        if idx == 1
            title(ax,Salt,'Fontsize',fs,'interpreter','latex')
        end
        subtitle(ax,int,'Fontsize',fs-4,'interpreter','latex')
        leg_plots = gobjects(1,N_Models);
        leg_texts = cell(1,N_Models);
        for kdx = 1:N_Models
            Model = Models{kdx};
            
            V_BH  = U_BH.(Salt).(Model).(int);
            r_BH  = U_BH.(Salt).(Model).r.*10;  % Angstrom
            V_WBK = U_WBK.(Salt).(Model).(int);
            r_WBK = U_WBK.(Salt).(Model).r.*10; % Angstrom
            
            switch lower(PlotType)
                case 'full'
                    plot(ax,r_BH,V_BH.f + V_BH.g + V_BH.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.f + V_WBK.g + V_WBK.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'full-derivative'
                    plot(ax,r_BH,V_BH.df + V_BH.dg + V_BH.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.df + V_WBK.dg + V_WBK.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'lj'
                    plot(ax,r_BH,V_BH.g + V_BH.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.g + V_WBK.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'lj-derivative'
                    plot(ax,r_BH,V_BH.dg + V_BH.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.dg + V_WBK.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'dispersion'
                    plot(ax,r_BH,V_BH.g,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.g,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'dispersion-derivative'
                    plot(ax,r_BH,V_BH.dg,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.dg,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'repulsive'
                    plot(ax,r_BH,V_BH.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.h,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
                case 'repulsive-derivative'
                    plot(ax,r_BH,V_BH.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','-');
                    leg_plots(kdx) = plot(ax,r_WBK,V_WBK.dh,'Color',Cols(kdx,:),'LineWidth',lw,'LineStyle','--');
            end
            leg_texts{kdx} = legtext.(Salt).(Model).(int);
        end
        % Update plot
        yline(ax,0,':k',"LineWidth",1)            
        ylim(ax,ylims)
        xlim(ax,[0,xmax])
        legend(ax,leg_plots,leg_texts,'interpreter','latex');
    end
end
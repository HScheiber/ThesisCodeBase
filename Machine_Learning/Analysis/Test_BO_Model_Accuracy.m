% Test_BO_Model_Accuracy
clear;

% Analysis parameters
Salt = 'LiI';
Theory = 'JC';
Model = 'CG1';
N_samples = 100; % How many samples to check
axpos = [0.10 0.100 0.800 0.850];  % [left bottom width height]
fs = 24; % fontsize

% Load data
[~,project] = find_home;
A = load([project filesep 'BO_Models' filesep Salt '_' Theory '_Model_' Model '_bayesopt.mat']);
results = A.results;
varnames = {results.VariableDescriptions.Name};
best_fx = results.MinObjective;

% Get random samples of parameters within their domain
N_Pars = length(results.VariableDescriptions);
Rdm_Samples = nan(N_samples,N_Pars);
for idx = 1:N_samples
    for jdx = 1:N_Pars
        Par_range = results.VariableDescriptions(jdx).Range;
        Rdm_Samples(idx,jdx) = Par_range(1) + rand*(Par_range(2) - Par_range(1));
    end
end

Rdm_Samples_t = array2table(Rdm_Samples,'VariableNames',varnames);

[Z_pred,Sigma_pred] = results.predictObjective(Rdm_Samples_t);

% Setup parallel mode
Parcores = feature('numcores');
Parcores = min(Parcores,N_samples);
ParStartTimer = tic; % Start timer
if ~isempty(gcp('nocreate'))
    Cur_Pool = gcp;
    Cur_Workers = Cur_Pool.NumWorkers;
    
    % Start the parallel pool
    if Cur_Workers ~= Parcores
        delete(Cur_Pool);
        % Create a "local" cluster object
        local_cluster = parcluster('local');

        % Modify the JobStorageLocation to a temporary directory
        tmp = tempname;
        if ~isfolder(tmp)
            mkdir(tmp)
        end
        local_cluster.JobStorageLocation = tempname;

        ppool = parpool(local_cluster,Parcores);
    else
        ppool = Cur_Pool;
    end
else
    ppool = parpool(Parcores); 
end
Time_Used = seconds(toc(ParStartTimer));
disp(['Time to Open Parallel Pool: ' datestr(Time_Used,'HH:MM:SS')]);

% Sample points
f(1:N_samples) = parallel.FevalFuture;
fun = results.ObjectiveFcn;

for idx = 1:N_samples   
    f(idx) = parfeval(ppool,fun,3,Rdm_Samples(idx,:));
end

wait(f)

% Collect outputs
Z_actual = nan(N_samples,1);
Minimization_Data = cell(N_samples,1);
for idx = 1:N_samples   
    [Z_actual(idx),~,Minimization_Data{idx}] = fetchOutputs(f(idx));
end

% Calculate error in predictions 
Z_error = Z_pred - Z_actual;


%% plot error in predictions vs GP uncertainty
figh = figure('WindowState','Maximized');
colormap(figh,'turbo')
ax = axes(figh,'Position',axpos,'YGrid','On','XGrid','On',...
    'FontSize',fs,'TickLabelInterpreter','latex');
hold(ax,'on')
xlabel(ax,'GP Uncertainty $\sigma(\mathbf{x})$','Interpreter','latex','FontSize',fs);
ylabel(ax,'$\hat{f}_{GP}(\mathbf{x}) - f(\mathbf{x})$','Interpreter','latex','FontSize',fs);


scatter(Sigma_pred,Z_error,50,Sigma_pred,'filled',...
    'MarkerEdgeColor','k')


%% plot error in predictions vs GP value

figh2 = figure('WindowState','Maximized');
colormap(figh2,'turbo')
ax2 = axes(figh2,'Position',axpos,'YGrid','On','XGrid','On',...
    'FontSize',fs,'TickLabelInterpreter','latex','DataAspectRatio',[1 1 1]);
hold(ax2,'on')
xlabel(ax2,'$f(\mathbf{x})$','Interpreter','latex','FontSize',fs);
%ylabel(ax2,'$\hat{f}_{GP}(\mathbf{x}) - f(\mathbf{x})$','Interpreter','latex','FontSize',fs);
ylabel(ax2,'$\hat{f}_{GP}(\mathbf{x}) - f(\mathbf{x})$','Interpreter','latex','FontSize',fs);


scatter(ax2,Z_actual,Z_error,50,Sigma_pred,'filled',...
    'MarkerEdgeColor','k') % Z_pred

scatter(ax2,best_fx,0,100,'r','filled',...
    'MarkerEdgeColor','k') % Z_pred

drawnow;

% Plot x = y lines
x1 = 0:0.1:ax2.YLim(2);
y1 = x1;

y2 = ax2.YLim(1):0.1:0;
x2 = abs(y2);

plot(ax2,x1,y1,':r','Linewidth',3)
plot(ax2,x2,y2,':r','Linewidth',3)

% Add colorbar
hCB = colorbar(ax2,'Position',[0.94 0.15 0.025 0.75],... % [left bottom width height]
    'TickLabelInterpreter','latex','FontSize',fs);
set(hCB.Title,{'String','Interpreter','FontSize'},{'$\hat{\sigma}_{GP}(\mathbf{x})$','latex',fs})
caxis(ax2,[0 ceil(max(Sigma_pred,[],'all'))])

%set(ax2,'xscale','log')

%% calculate average errors
av_error = mean(Z_error);
av_abs_error = mean(abs(Z_error));

%% Save the model
save([Salt '_' Theory '_Model_' Model '_Errors_N' num2str(N_samples) '.mat'],...
    'Z_pred','Sigma_pred','Z_actual','Z_error','N_samples','Minimization_Data');

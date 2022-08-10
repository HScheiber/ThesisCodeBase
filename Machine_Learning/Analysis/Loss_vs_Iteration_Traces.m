% Script for plotting the bayesian optimization trace of multiple models

%% Analysis Parameters
fs = 24;
Salt = 'LiBr';
Theory = 'BH';
Basenum = 'F';
Midnum = 'B';
min_iter = 900; % Set the first iterations to show, set to 1 for all numbers
max_iter = Inf; % Set the max to show, set to Inf for all numbers
ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

% Find reps of models
files = {dir(fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Basenum Midnum '*'])).name};
Models = unique(regexp(files,[Basenum Midnum '([0-9])+'],'match','once'));
nonempty_idx = ~cellfun(@isempty,Models);
Models = Models(nonempty_idx);
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end
Model_Labels = Models;

filename = [Salt '_' Theory '_' Models{1} '-' Models{end} '_Traces.png'];

N_Models = length(Models);
Colours = cbrewer('qual','Set1',N_Models,'spline');

plot_min = true; % If true plots the minimum loss vs step
plot_iter = true; % If true, plots the current loss value vs step
plot_params = true; % If true, plots the values of each parameter vs step. Note, all models should have the same number of parameters

%% Other parameters
savefile = false; % switch to save the final plots to file
fullopt_dat = true;
% Find model in results
warning('off','MATLAB:cellfun:NotACell')

% Preallocate data arrays
p = gobjects(N_Models,1);

% Loop over chosen models
for idx = 1:N_Models
    Model = Models{idx};
    
    % Find the bayesian optimized model and full opt history
    data_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);

	if ~isfile(data_file)
        %error(['No results found for: ' Salt ', ' Theory ', Model ' Model '.']);
        continue
	end

    % Load bayesopt data
    data = load(data_file).full_data;
    data_bayes = data.bayesopt_results;
    
    % Load fullopt data
    if fullopt_dat
        data_full  = data.secondary_result;
    end
    
    % Check how many parameters exist
    if idx == 1
        if plot_params
            N_Pars = size(data_bayes.XTrace,2);
            ParNames = data_bayes.XTrace.Properties.VariableNames;
        else
            N_Pars = 0;
        end
        
        N_Plots = sum([plot_min plot_iter N_Pars]);
        ydat = cell(N_Models,N_Plots);
        xdat = cell(N_Models,1);
        titles = cell(1,N_Plots);
        yaxislabel = cell(1,N_Plots);
        yaxislog = false(1,N_Plots);
    end
    
    if fullopt_dat
        xdat{idx} = 1:(length(data_bayes.ObjectiveTrace) + length(data_full));
    else
        xdat{idx} = 1:length(data_bayes.ObjectiveTrace);
    end
    pdx = 1;
    
    if plot_min
        if fullopt_dat
            ydat{idx,pdx} = [data_bayes.ObjectiveMinimumTrace; ...
                extract_fullopt_hist(data_full,'optimValues','fval')];
        else
            ydat{idx,pdx} = data_bayes.ObjectiveMinimumTrace;
        end
        
        if idx == 1
            titles{pdx} = 'Lowest known $f(\mathbf{x})$';
            yaxislog(pdx) = true;
            yaxislabel{pdx} = '$f(\mathbf{x^{*}})$';
        end
        pdx = pdx + 1;
    end
    if plot_iter
        if fullopt_dat
            ydat{idx,pdx} = [data_bayes.ObjectiveTrace; ...
                extract_fullopt_hist(data_full,'optimValues','fval')];
        else
            ydat{idx,pdx} = data_bayes.ObjectiveTrace;
        end
        
        if idx == 1
            titles{pdx} = '$f(\mathbf{x})$';
            yaxislog(pdx) = true;
            yaxislabel{pdx} = '$f(\mathbf{x})$';
        end
        pdx = pdx + 1;
    end
    if plot_params
        if fullopt_dat
            fulldat = extract_fullopt_hist(data_full,'optimValues','x');
        end
        
        for qq = 1:N_Pars
            if fullopt_dat
                ydat{idx,pdx} = [data_bayes.XTrace{:,qq}; ...
                    fulldat(:,qq)];
            else
                ydat{idx,pdx} = data_bayes.XTrace{:,qq};
            end
            titles{pdx} = param_name_map(ParNames{qq},Salt);
            yaxislabel{pdx} = param_name_map(ParNames{qq},Salt);
            pdx = pdx + 1;
        end
    end
end

% Create figure and axis
figh = figure('WindowState','maximized','NumberTitle','off',...
    'Name','','Visible','On');
t = tiledlayout(figh,N_Plots,1,'TileSpacing','tight');

for idx = 1:N_Plots
    yy = ydat(:,idx); % Select the plots for each model
    axobj = nexttile(t);
    hold(axobj,'on')
    
    for jdx = 1:N_Models
        if max_iter < length(xdat{jdx})
            maxn = max_iter;
        else
            maxn = length(xdat{jdx});
        end
        p(jdx) = plot(axobj,xdat{jdx}(min_iter:maxn),ydat{jdx,idx}(min_iter:maxn),'-','linewidth',2,'color',Colours(jdx,:));
    end
    
    % Add plot title and y axis
    set(axobj,'box','on','TickLabelInterpreter','latex');
    set(axobj,'XMinorTick','on','YMinorTick','on','FontSize',fs);
    grid(axobj,'minor');
    ylim(axobj,'padded')
    xlim(axobj,[min_iter max_iter])
    if idx ~= N_Plots
        xticklabels(axobj,[])
    end
    %set(axobj,'XLim',[0 N_Salts+0.5])
    
    if yaxislog(idx)
        set(axobj, 'YScale', 'log')
    end
    
    title(axobj, titles{idx},'Fontsize',fs,'Interpreter','latex')
    ylabel(axobj, yaxislabel{idx},'Fontsize',fs-4,'Interpreter','latex')
end

xlabel(axobj, 'Optimization Iteration','Fontsize',fs,'Interpreter','latex')



legh = legend(axobj,p,Model_Labels,'Position',[0.9334 0.4222 0.0622 0.2041],'Orientation','Vertical',...
    'Interpreter','latex','Box','off','fontsize',fs,'NumColumns', 1);
title(legh,'\bf{\underline{Models}}','Fontsize',fs,'Interpreter','latex')
% title(legh,['\bf{\underline{' Legend_title '}}'],'Fontsize',fs,'Interpreter','latex')

if savefile
    set(figh,'renderer','opengl')
    exportgraphics(figh,filename,'Resolution',300)
end



function output = param_name_map(p_name,Salt)
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    switch p_name
        case {'SDMM' 'SD6MM'}
            output = ['$S_D$[' Metal '-' Metal ']'];
        case {'SDXX' 'SD6XX'}
            output = ['$S_D$[' Halide '-' Halide ']'];
        case {'SDMX' 'SD6MX'}
            output = ['$S_D$[' Metal '-' Halide ']'];
        case 'SRMM'
            output = ['$S_R$[' Metal '-' Metal ']'];
        case 'SRXX'
            output = ['$S_R$[' Halide '-' Halide ']'];
        case 'SRMX'
            output = ['$S_R$[' Metal '-' Halide ']'];
        case 'SNMM'
            output = ['$S_N$[' Metal '-' Metal ']'];
        case 'SNXX'
            output = ['$S_N$[' Halide '-' Halide ']'];
        case 'SNMX'
            output = ['$S_N$[' Metal '-' Halide ']'];
        case 'SQ'
            output = '$S_Q$';
            
        case 'r0_MM'
            output = ['$r_{0}$[' Metal '-' Metal '] [nm]'];
        case 'r0_XX'
            output = ['$r_{0}$[' Halide '-' Halide '] [nm]'];
        case 'r0_MX'
            output = ['$r_{0}$[' Metal '-' Halide '] [nm]'];
        case 'gamma_MM'
            output = ['$\gamma$[' Metal '-' Metal ']'];
        case 'gamma_XX'
            output = ['$\gamma$[' Halide '-' Halide ']'];
        case 'gamma_MX'
            output = ['$\gamma$[' Metal '-' Halide ']'];
            
        case 'epsilon_MM'
            output = ['$\epsilon$[' Metal '-' Metal '] [kJ mol$^{-1}$]'];
        case 'epsilon_XX'
            output = ['$\epsilon$[' Halide '-' Halide '] [kJ mol$^{-1}$]'];
        case 'epsilon_MX'
            output = ['$\epsilon$[' Metal '-' Halide '] [kJ mol$^{-1}$]'];
            
        otherwise
            output = p_name;
    end
end

function output = extract_fullopt_hist(data,field,subfield)
    if strcmpi(subfield,'x')
        output = nan(length(data),length(data(1).(subfield)));
        for idx = 1:length(data)
            output(idx,:) = data(idx).(subfield);
        end
    else
        output = nan(length(data),length(data(1).(field).(subfield)));
        for idx = 1:length(data)
            output(idx,:) = data(idx).(field).(subfield);
        end
    end
end
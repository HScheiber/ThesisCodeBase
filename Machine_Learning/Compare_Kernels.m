%% Data options
KernelFunctions = {'ardsquaredexponential' 'ardexponential' 'ardmatern32' 'ardmatern52' 'ardrationalquadratic'};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; %  'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Theory = 'BF';
ModelID = 'PS';
IsDeterministic = false;
Testfits = 4;

Log_likelihood = nan(length(Salts),length(KernelFunctions),Testfits);

for salt_idx = 1:length(Salts)
    Salt = Salts{salt_idx};

    % Find reps of models
    ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\Model_Building\Completed';
    files = {dir(fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' ModelID '*_data.mat'])).name};
    Models = unique(regexp(files,[ModelID '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models);
    Models = Models(nonempty_idx);
    N_Models = numel(Models);

    % Make sure the models are sorted
    if ~isempty(Models)
        [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
        Models = Models(idx);
    end

    % Find the targets
    for jdx = 1:N_Models
        dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Models{jdx} '_data.mat']);
        if ~isfile(dat_file)
            continue
        else
            % Load data
            try
                data = load(dat_file).full_data;
                Bayesopt_model = data.Settings;
                Bayesopt_Loss_Options = Bayesopt_model.Loss_Options;
                params = data.bayesopt_results.Options.VariableDescriptions;
                Settings = data.Settings;
                break
            catch
                continue
            end
        end
    end

    % Load all data
    [Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
        Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
    Settings.Initialize_From_Model = {ModelID};
    Prev_Data = LoadPrevModelData(Settings,params,Inf);

    % Remove errors from data
    nan_idx = isnan(Prev_Data.InitialObjective);
    Data_X = Prev_Data.InitialX(~nan_idx,:);
    Data_Y = Prev_Data.InitialObjective(~nan_idx);

    % Generate parameters for GPR and transform X Data
    LBTrans = zeros(length(params),1);
    UBTrans = zeros(length(params),1);
    Data_X_Trans = Data_X;
    for idx = 1:length(params)
        LBTrans(idx) = params(idx).Range(1);
        UBTrans(idx) = params(idx).Range(2);
        if strcmpi(params(idx).Transform,'log')
            LBTrans(idx) = log(LBTrans(idx));
            UBTrans(idx) = log(UBTrans(idx));
            Data_X_Trans{:,idx} = log(Data_X{:,idx});
        end
    end
    isCat = false(length(params),1);

    BasisFunction = 'constant';
    Standardize = false;
    FitMethod = 'sd'; % 'exact' or 'sd'
    PredictMethod = 'sd'; % 'exact' or 'sd'
    ActiveSetSize = 1000;
    ActiveSetMethod = 'random';

    if IsDeterministic
        SigmaLowerBound = max(1e-8, 1e-4*nanstd(Data_Y));
        Sigma0 = SigmaLowerBound;
        ConstantFSigma = true;
    else
        Sigma0Divisor = 5;
        SigmaLowerBound = max(1e-8, nanstd(Data_Y)*1e-2);
        Sigma0 = max(SigmaLowerBound, nanstd(Data_Y)/Sigma0Divisor);
        ConstantSigma = false; % only set to true for deterministic, and combine with small sigma
    end
    
    sigmaF = nanstd(Data_Y)/sqrt(2);
    if isnan(sigmaF) || sigmaF == 0
        sigmaF = 1;
    end
    
    for kernel_idx = 1:length(KernelFunctions)
        KernelFunction = KernelFunctions{kernel_idx};
        switch KernelFunction
            case {'exponential' 'squaredexponential' 'matern32' 'matern52'}
                KernelParams0 = [mean((UBTrans(:)-LBTrans(:))/2); sigmaF];    % Last one is the default kernel amplitude
                ConstantKernelParameters = [false; false];
            case 'rationalquadratic'
                KernelParams0 = [mean((UBTrans(:)-LBTrans(:))/2); 1; sigmaF];    % Last one is the default kernel amplitude
                ConstantKernelParameters = [false; false; false];
            case {'ardexponential' 'ardsquaredexponential' 'ardmatern32' 'ardmatern52'}
                KernelParams0 = [(UBTrans(:)-LBTrans(:))/2; sigmaF];    % Last one is the default kernel amplitude
                ConstantKernelParameters = [isCat(:); false];
            case {'ardrationalquadratic'}
                KernelParams0 = [(UBTrans(:)-LBTrans(:))/2; 1; sigmaF];    % Last one is the default kernel amplitude
                ConstantKernelParameters = [isCat(:); false; false];
        end
        % Fit Gaussian Process
        %HyperparameterOptimizationOptions = struct('AcquisitionFunctionName','expected-improvement-plus','Holdout',0.8);
        
        for tf_idx = 1:Testfits
            gprMdl = fitrgp(Data_X_Trans,Data_Y,'BasisFunction',BasisFunction,'ConstantKernelParameters',ConstantKernelParameters,...
                'ConstantSigma',ConstantSigma,'KernelFunction',KernelFunction,'KernelParameters',KernelParams0,...
                'Sigma',Sigma0,'Standardize',Standardize,'FitMethod',FitMethod,'PredictMethod',PredictMethod,...
                'ActiveSetSize',ActiveSetSize,'ActiveSetMethod',ActiveSetMethod);
                %'OptimizeHyperparameters',{'KernelFunction'},'HyperparameterOptimizationOptions',HyperparameterOptimizationOptions);
            Log_likelihood(salt_idx,kernel_idx,tf_idx) = gprMdl.LogLikelihood;
        end
        disp(['Test Complete for ' Salt ' - ' KernelFunction])
        disp(['Measured Log Likelihood over ' num2str(Testfits) ' tests with active set size = ' num2str(ActiveSetSize)])
        disp(['Average Log Likelihood: ' num2str(mean(Log_likelihood(salt_idx,kernel_idx,:)),'%.4f')])
        disp('*****************************************************************')
    end
end

Average_Log_Likelihood = sum(mean(Log_likelihood,3),1);
% -8905.95258331056 -8944.86017690041 -8888.00960831035 -8871.70868254116 -8861.26987213706




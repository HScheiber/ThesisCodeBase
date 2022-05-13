classdef GPStrategy < bayesoptim.ModelStrategy
    %BAYESIANOPTIMIZATIONGP This class contains methods in which Gaussian
    %Processes are used as estimators and cases specific to Gaussian
    %Processes
    
    %   Copyright 2020 The MathWorks, Inc.

    methods (Access = private)
        
        function Args = GPFitMethodArgs(this, XTrain, YTrain, PrivOptions)
            NumFiniteY = sum(isfinite(YTrain));
            if PrivOptions.GPActiveSetSize <= NumFiniteY
                Args = {'FitMethod', 'sd', 'ActiveSetSize', PrivOptions.GPActiveSetSize};
            else
                Args = {'FitMethod', 'exact'};
            end
        end
        
        function [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, TrainingData, PrivOptions)
            VarSpec = PrivOptions.VarSpec;
            sigmaF = nanstd(TrainingData)/sqrt(2);
            if isnan(sigmaF) || sigmaF == 0
                sigmaF = 1;
            end
            
            switch PrivOptions.KernelFunction
                case {'exponential' 'squaredexponential' 'matern32' 'matern52'}
                    KernelParams0 = [mean((VarSpec.UBTrans(:)-VarSpec.LBTrans(:))/2); sigmaF];    % Last one is the default kernel amplitude
                    ConstantKernelParameters = [false; false];
                case 'rationalquadratic'
                    KernelParams0 = [mean((VarSpec.UBTrans(:)-VarSpec.LBTrans(:))/2); 1; sigmaF];    % Last one is the default kernel amplitude
                    ConstantKernelParameters = [false; false; false];
                case {'ardexponential' 'ardsquaredexponential' 'ardmatern32' 'ardmatern52'}
                    KernelParams0 = [(VarSpec.UBTrans(:)-VarSpec.LBTrans(:))/2; sigmaF];    % Last one is the default kernel amplitude
                    KernelParams0(VarSpec.isCat) = PrivOptions.CatKernelScale;
                    ConstantKernelParameters = [VarSpec.isCat(:); false];
                case {'ardrationalquadratic'}
                    KernelParams0 = [(VarSpec.UBTrans(:)-VarSpec.LBTrans(:))/2; 1; sigmaF];    % Last one is the default kernel amplitude
                    KernelParams0(VarSpec.isCat) = PrivOptions.CatKernelScale;
                    ConstantKernelParameters = [VarSpec.isCat(:); false; false];
            end
        end
        
        function ConstraintGP = fitConstraintModel(this, GP, ConstraintTrain, isDeterministic, useFitMethodNone, PrivOptions, XTrain)
            if all(isnan(ConstraintTrain))
                ConstraintGP = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitConstraint');
                end
            elseif useFitMethodNone && ~isempty(GP)
                ConstraintGP = refitWithFitMethodNone(this, GP, XTrain, ConstraintTrain, PrivOptions);
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, ConstraintTrain, PrivOptions);
                if isDeterministic
                    std = nanstd(ConstraintTrain);
                    if std == 0
                        SigmaLowerBound = 1e-4;
                    else
                        SigmaLowerBound = 1e-4*std;
                    end
                    sigma0 = SigmaLowerBound;
                    ConstantSigma = true;
                    FitMethodArgs = GPFitMethodArgs(this, XTrain, ConstraintTrain, PrivOptions);
                    ConstraintGP = iFitRobust(this, GP, XTrain, ConstraintTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'ConstantSigma', ConstantSigma,...
                        'KernelFunction', PrivOptions.KernelFunction, ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        FitMethodArgs{:});
                else
                    std = nanstd(ConstraintTrain);
                    if std == 0
                        SigmaLowerBound = 1e-2;
                    else
                        SigmaLowerBound = 1e-2*std;
                    end
                    FitMethodArgs = GPFitMethodArgs(this, XTrain, ConstraintTrain, PrivOptions);
                    ConstraintGP = iFitRobust(this, GP, XTrain, ConstraintTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'KernelFunction', PrivOptions.KernelFunction, ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Standardize', false,...
                        FitMethodArgs{:});
                end
            end
        end
    end
    
    methods
        
        function ObjectiveFcnGP = fitObjectiveFcnModel(this, XTrain, FTrain, PrivOptions, Model, varargin)
            BullAdjustment = varargin{1}{1};
            BullAdjustmentN = varargin{1}{2};
            useFitMethodNone = varargin{1}{3};
            if all(isnan(FTrain))
                ObjectiveFcnGP = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitObjective');
                end
            elseif useFitMethodNone && ~isempty(this)
                ObjectiveFcnGP = refitWithFitMethodNone(this, Model, XTrain, FTrain, PrivOptions);
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, FTrain, PrivOptions);
                if PrivOptions.IsObjectiveDeterministic
                    % Deterministic objective function
                    SigmaLowerBound = max(1e-8, 1e-4*nanstd(FTrain));
                    sigma0 = SigmaLowerBound;
                    ConstantFSigma = true;
                    FitMethodArgs = GPFitMethodArgs(this, XTrain, FTrain, PrivOptions);
                    ObjectiveFcnGP = iFitRobust(this, Model, XTrain, FTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'ConstantSigma', ConstantFSigma,...
                        'KernelFunction', PrivOptions.KernelFunction, ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        FitMethodArgs{:});
                else
                    % Nondeterministic objective function
                    SigmaLowerBound = max(1e-8, nanstd(FTrain)*1e-2);
                    sigma0 = max(SigmaLowerBound, nanstd(FTrain)/PrivOptions.Sigma0Divisor);
                    FitMethodArgs = GPFitMethodArgs(this, XTrain, FTrain, PrivOptions);
                    ObjectiveFcnGP = iFitRobust(this, Model, XTrain, FTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'KernelFunction', PrivOptions.KernelFunction, ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        FitMethodArgs{:});
                end
                if BullAdjustment
                    % Apply Bull's adjustment to the global kernel amplitude, and
                    % recompute alphas. This requires passing in all the current values
                    % of other parameters and using 'FitMethod' 'none'.
                    KernelParameters =  ObjectiveFcnGP.KernelInformation.KernelParameters;
                    KernelParameters(end) = KernelParameters(end)*BullAdjustmentN;
                    % recompute alphas only
                    ObjectiveFcnGP = iFitRobust(this, Model, XTrain, FTrain, SigmaLowerBound,...
                        'KernelFunction', PrivOptions.KernelFunction, ...
                        'KernelParameters', KernelParameters, ...
                        'Sigma', ObjectiveFcnGP.Sigma, ...
                        'BasisFunction', 'constant', ...
                        'Beta', ObjectiveFcnGP.Beta, ...
                        'FitMethod', 'none', ...
                        'PredictMethod', 'exact', ...
                        'Standardize', false);
                end
            end
        end
        
        function ObjectiveEvaluationTimeGP = fitObjectiveEvaluationTimeModel(this, XTrain, FTrain, PrivOptions, ObjectiveEvaluationTimeTrain, Model, varargin)
            useFitMethodNone = varargin{1}{1};
            % Only fit model to successful function evaluations
            FitIndices = ~isnan(ObjectiveEvaluationTimeTrain) & ~isnan(FTrain);
            if ~any(FitIndices)
                ObjectiveEvaluationTimeGP = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitObjectiveEvaluationTime');
                end
            elseif useFitMethodNone && ~isempty(this)
                ObjectiveEvaluationTimeGP = refitWithFitMethodNone(this, Model, ...
                                                                     XTrain(FitIndices,:), ...
                                                                     log(ObjectiveEvaluationTimeTrain(FitIndices)),...
                                                                     PrivOptions);
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, ObjectiveEvaluationTimeTrain, PrivOptions);
                SigmaLowerBound = std(log(ObjectiveEvaluationTimeTrain(FitIndices)))*1e-2;
                FitMethodArgs = GPFitMethodArgs(this, XTrain(FitIndices,:), log(ObjectiveEvaluationTimeTrain(FitIndices)), PrivOptions);
                ObjectiveEvaluationTimeGP = iFitRobust(this, Model, XTrain(FitIndices,:), ...
                    log(ObjectiveEvaluationTimeTrain(FitIndices)), ...
                    SigmaLowerBound, ...
                    'BasisFunction', 'constant', ...
                    'ConstantKernelParameters', ConstantKernelParameters, ...
                    'KernelFunction', PrivOptions.KernelFunction, ...
                    'KernelParameters', KernelParams0, ...
                    'PredictMethod', 'exact', ...
                    'Standardize', false,...
                    FitMethodArgs{:});
            end
        end
        
        function ErrorGP = fitErrorModel(this,  XTrain, FTrain, PrivOptions, ErrorTrain, Model, varargin)
            useFitMethodNone = varargin{1}{1};
            if all(isnan(ErrorTrain))
                ErrorGP = [];
                if PrivOptions.Verbose >= 3
                    bayesoptim.printInfo('CantFitError');
                end
            elseif all(ErrorTrain < 0)
                ErrorGP = [];
            elseif useFitMethodNone && ~isempty(this)
                ErrorGP = refitWithFitMethodNone(this, Model, XTrain, ErrorTrain, PrivOptions);
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, ErrorTrain, PrivOptions);
                SigmaLowerBound = nanstd(ErrorTrain)*1e-4;
                FitMethodArgs = GPFitMethodArgs(this, XTrain, ErrorTrain, PrivOptions);
                ErrorGP = iFitRobust(this, Model, XTrain, ErrorTrain, SigmaLowerBound, ...
                    'BasisFunction', 'constant', ...
                    'ConstantKernelParameters', ConstantKernelParameters, ...
                    'KernelFunction', PrivOptions.KernelFunction, ...
                    'KernelParameters', KernelParams0, ...
                    'PredictMethod', 'exact', ...
                    'Standardize', false,...
                    FitMethodArgs{:});
            end
        end
        
        function ConstrGPs = fitConstraintModels(this, XTrain, ConstraintTrain, ConstraintGPs, PrivOptions, varargin)
            useFitMethodNone = varargin{1}{1};
            ConstrGPs = cell(1,PrivOptions.NumCoupledConstraints);
            if ~isempty(ConstraintTrain)
                for i = PrivOptions.NumCoupledConstraints:-1:1
                    if useFitMethodNone && ~isempty(ConstraintGPs{i})
                        ConstrGPs{i} = fitConstraintModel(this, ConstraintGPs{i},...
                            ConstraintTrain(:,i), PrivOptions.AreCoupledConstraintsDeterministic(i), ...
                            useFitMethodNone, PrivOptions, XTrain);
                    else
                        ConstrGPs{i} = fitConstraintModel(this, [],...
                            ConstraintTrain(:,i), PrivOptions.AreCoupledConstraintsDeterministic(i), ...
                            useFitMethodNone, PrivOptions, XTrain);
                    end
                end
            end
        end
        
        function GP = iFitRobust(this, Model, X, Y, SigmaLowerBound, varargin)
            % Fit a GP. If it fails due to singularity, iteratively double
            % SigmaLowerBound until it succeeds (up to 10 times). If it never succeeds,
            % return [];
            if all(isnan(Y))
                bayesoptim.err('AllModelTargetsNaN');
            end
            SigmaLowerBound = max(SigmaLowerBound, 1e-6);
            GP = [];
            success = false;
            doublings = 0;
            while ~success && doublings <= 10
                try
                    GPModel = compact(fitrgp(X, Y, varargin{:}, 'SigmaLowerBound', SigmaLowerBound));
                    GP = bayesoptim.BayesOptGPModel(GPModel.Sigma, GPModel, ...
                        GPModel.KernelInformation, GPModel.Beta);
                    success = true;
                catch me
                    SigmaLowerBound = 2*SigmaLowerBound;
                    doublings = doublings + 1;
                end
            end
        end
        
        function GP = refitWithFitMethodNone(this, GP, XTrain, FTrain, PrivOptions)
            % Recompute alphas. This requires passing in the current parameters and
            % using 'FitMethod' 'none'.
            GP = iFitRobust(this, GP, XTrain, FTrain, 0,...
                'KernelFunction', PrivOptions.KernelFunction, ...
                'KernelParameters', GP.KernelInformation.KernelParameters, ...
                'Sigma', GP.Sigma, ...
                'BasisFunction', 'constant', ...
                'Beta', GP.Beta, ...
                'FitMethod', 'none', ...
                'PredictMethod', 'exact', ...
                'Standardize', false);
        end
    end
end


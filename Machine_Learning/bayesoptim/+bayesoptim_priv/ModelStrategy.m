classdef ModelStrategy
    %MODELSTRATEGY Abstract class for Model Strategies
    %   This class defines the methods that need to be implemented by child
    %   classes. The methods contain fitting models for different objects.

%   Copyright 2020 The MathWorks, Inc.
   
    methods (Abstract)
        fitObjectiveFcnModel(this, XTrain, FTrain, PrivOptions, ...
            Model, varargin);
        
        fitObjectiveEvaluationTimeModel(this, XTrain, FTrain, PrivOptions, ...
            ObjectiveEvaluationTimeTrain, Model, varargin);
        
        fitErrorModel(this,  XTrain, FTrain, PrivOptions, ErrorTrain, ...
            Model, varargin);
        
        fitConstraintModels(this, XTrain, ConstraintTrain, ConstraintGPs, ...
            PrivOptions, varargin);
        
        iFitRobust(this, X, Y, SigmaLowerBound, varargin);
        
        refitWithFitMethodNone(this, GP, XTrain, FTrain);
    end
    
end


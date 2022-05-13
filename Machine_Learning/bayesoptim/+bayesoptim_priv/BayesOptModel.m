classdef BayesOptModel
    %BAYESOPTMODEL This is the parent class for different models used in
    %esitmation differnt aspects in Bayesian Optimization and contains
    %common properties and methods of child classes.

%   Copyright 2020 The MathWorks, Inc.
    
    properties
        Sigma
        Model
    end
    
    methods
        function obj = BayesOptModel(Sigma, Model)
            obj.Sigma = Sigma;
            obj.Model = Model;
        end
    end
end


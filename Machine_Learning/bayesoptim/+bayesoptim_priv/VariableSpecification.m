classdef VariableSpecification
    % One instance represents all the variables.
    
%   Copyright 2016 The MathWorks, Inc.

    properties
        Names           % cellstr
        LBs             % vector
        UBs             % vector
        Categories      % cell array of categorical arrays
    end
    
    properties(SetAccess = private)
        Types           % cellstr
        Transforms      % cellstr
        isReal          % A logical vector used for efficiency
        isCat           % A logical vector used for efficiency
        isLog           % A logical vector used for efficiency
    end
    
    properties(Dependent)
        NumVars
        % Transformed lower and upper bounds
        LBTrans     	% vector
        UBTrans      	% vector
    end
    
    methods
        function this = VariableSpecification(Names, Types, Transforms)
            this.Names = Names;
            % Set types
            this.Types = Types;
            this.isReal = ismember(this.Types, 'real');
            this.isCat = ismember(this.Types, 'categorical');
            % Set transforms
            this.Transforms = Transforms;
            this.isLog = ismember(this.Transforms, 'log');
        end
            
        function N = get.NumVars(this)
            N = length(this.Names);
        end
        
        function LBT = get.LBTrans(this)
            LBT = this.LBs;
            LBT(this.isLog) = log(LBT(this.isLog));
        end
        
        function UBT = get.UBTrans(this)
            UBT = this.UBs;
            UBT(this.isLog) = log(UBT(this.isLog));
        end
    end
end

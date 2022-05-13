classdef optimizableVariable
    % optimizableVariable A variable that can be optimized using bayesopt
    % or other optimizers.
    %
    %   Properties:
    %       Name        - A character vector naming the variable.
    %     	Range       - Either a length-2 vector of the form
    %                     [LowerBound, UpperBound] defining the
    %                     range of a numeric variable, or a cell array
    %                     of character vectors defining the possible values
    %                     of a categorical variable.
    %       Type        - A character vector, one of 'real', 'integer', or
    %                     'categorical' indicating the type of the
    %                     variable. If Type is 'integer', then the Range
    %                     argument must be whole numbers. Default: 'real'
    %                     if the Range argument is a length-2 numeric
    %                     vector, or 'categorical' if Range is a cell array
    %                     of character vectors.
    %     	Transform   - A character vector: either 'log' or 'none'. 'log'
    %                     indicates that the variable should be searched
    %                     and modeled on a log scale.
    %       Optimize    - A logical scalar indicating whether the
    %                     variable should be optimized.
    %   
    %   Construction:
    %       (1) Variable = optimizableVariable(Name, Range) creates a
    %       variable with the given Name and Range of values.
    %     
    %       (2) Variable = optimizableVariable(..., Param, value,...)
    %       specifies additional variable properties via optional name
    %       value pairs:
    %           Type        - A character vector, one of 'real', 'integer',
    %                         or 'categorical'. If Type is 'integer', then
    %                         the Range argument must be whole numbers.
    %                         Default: 'real' if Range is a length-2
    %                         numeric vector; 'categorical' if Range is a
    %                         cell array of character vectors.
    %           Transform   - A character vector, either 'log' or 'none'.
    %                         Default: 'none'
    %           Optimize    - A logical scalar indicating whether the
    %                         variable should be optimized.
    %                         Default: true
    %
    %   Examples:
    %       Variable1 = optimizableVariable('Var1', [1 20])
    %     
    %       Variable2 = optimizableVariable('Var2', [1 20], ...
    %                                    'Type', 'integer', ...
    %                                    'Transform', 'log')
    %     
    %       Variable3 = optimizableVariable('Var3', {'cat1', 'cat2'},...
    %                                    'Optimize', false)
    %     
    %       AllVariables = [Variable1, Variable2, Variable3]
    %
    %   See also: BAYESOPT, HYPERPARAMETERS

    %   Copyright 2016-2017 The MathWorks, Inc.
        
    properties(Dependent)
        Name
        Range
        Type
        Transform
        Optimize
    end
    
    properties(Access=protected)        
        % Whether 'Type' was explicitly set by the user.
        TypeSet = false;
        
        % Dependent property values 
        PrivName = '';
        PrivRange = [0 1];
        PrivType = 'real';
        PrivTransform = 'none';
        PrivOptimize = false;
    end
    
    methods
        function this = optimizableVariable(Name, Range, varargin)
            if nargin > 0
                % Convert String inputs to char/cellstr
                Name = convertStringsToChars(Name);
                Range = convertStringsToChars(Range);
                [varargin{:}] = convertStringsToChars(varargin{:});
                
                this.PrivName = Name;
                this.PrivRange = Range;
                ArgNames = {'Type', 'Transform', 'Optimize'};
                Defaults = {'real', 'none', true};
                [Type, Transform, Optimize, SetFlags] = internal.stats.parseArgs(...
                    ArgNames, Defaults, varargin{:});
                this.TypeSet = SetFlags.Type;
                this.PrivType = Type;
                this.PrivTransform = Transform;
                this.PrivOptimize = Optimize;
                this = checkAndFill(this);
            end
        end
        
        % Dependent property getters
        function Name = get.Name(this)
            Name = this.PrivName;
        end
        
        function Range = get.Range(this)
            Range = this.PrivRange;
        end
        
        function Type = get.Type(this)
            Type = this.PrivType;
        end
        
        function Transform = get.Transform(this)
            Transform = this.PrivTransform;
        end
        
        function Optimize = get.Optimize(this)
            Optimize = this.PrivOptimize;
        end

        % Dependent setters
        function this = set.Name(this, Name)
            Name = convertStringsToChars(Name);
            this.PrivName = Name;
            this = checkAndFill(this);
        end
        
        function this = set.Range(this, Range)
            Range = convertStringsToChars(Range);
            this.PrivRange = Range;
            this = checkAndFill(this);
        end
        
        function this = set.Type(this, Type)
            Type = convertStringsToChars(Type);
            this.PrivType = Type;
            this.TypeSet = true;
            this = checkAndFill(this);
        end
        
        function this = set.Transform(this, Transform)
            Transform = convertStringsToChars(Transform);
            this.PrivTransform = Transform;
            this = checkAndFill(this);
        end
        
        function this = set.Optimize(this, Optimize)
            this.PrivOptimize = Optimize;
            this = checkAndFill(this);
        end
    end
    
    methods(Access=protected)
        function this = checkAndFill(this)
            checkName(this);
            checkRange(this);
            this = checkType(this);
            checkTransform(this);
            checkOptimize(this);
        end
        
        function checkName(this)
            if ~isa(this.PrivName, 'char')
                bayesoptim.err('VarName');
            end
        end
        
        function checkRange(this)
            if numel(this.PrivRange) < 2
                bayesoptim.err('RangeType');
            elseif iscellstr(this.PrivRange)
                return;
            elseif ~(isnumeric(this.PrivRange) && isreal(this.PrivRange) && ...
                    isvector(this.PrivRange) && numel(this.PrivRange)==2)
                bayesoptim.err('RangeType');
            elseif ~(this.PrivRange(1) < this.PrivRange(2))
                bayesoptim.err('RangeOrder');
            end
        end
        
        function this = checkType(this)
            if this.TypeSet
                if isempty(this.PrivType)
                    bayesoptim.err('TypeType')
                elseif ~isa(this.PrivType, 'char') || ~ismember(this.PrivType, {'real', 'integer', 'categorical'})
                    bayesoptim.err('TypeType');
                elseif isequal(this.PrivType, 'real') && iscellstr(this.PrivRange)
                    bayesoptim.err('TypeCatIfRangeCellstr');
                elseif isequal(this.PrivType, 'integer') && iscellstr(this.PrivRange)
                    bayesoptim.err('TypeCatIfRangeCellstr');
                elseif isequal(this.PrivType, 'integer') && iAnyNonIntegral(this.PrivRange)
                    bayesoptim.err('TypeIntRangeReal');
                elseif isequal(this.PrivType, 'categorical') && ~iscellstr(this.PrivRange)
                    bayesoptim.err('TypeCatRangeBad');
                end
            elseif iscellstr(this.PrivRange)
                this.PrivType = 'categorical';
            end
        end
        
        function checkTransform(this)
            if isempty(this.PrivTransform)
                bayesoptim.err('TransformType');
            elseif ~isa(this.PrivTransform, 'char') || ...
                    ~ismember(this.PrivTransform, {'none', 'log'})
                bayesoptim.err('TransformType');
            elseif isequal(this.PrivTransform, 'log') && iscellstr(this.PrivRange)
                bayesoptim.err('TransformLogRangeCellstr');
            elseif isequal(this.PrivTransform, 'log') && any(this.PrivRange <= 0)
                bayesoptim.err('TransformLogRangeNonpositive');
            end
        end
        
        function checkOptimize(this)
            if ~(isscalar(this.PrivOptimize) && islogical(this.PrivOptimize))
                bayesoptim.err('Optimize');
            end
        end
    end
end

function tf = iAnyNonIntegral(Range)
tf = any(floor(Range) ~= Range);
end

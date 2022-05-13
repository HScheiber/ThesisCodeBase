function c = elcm(varargin)
%elcm extends the functionality of lcm to multiple inputs.
% Any size, shape, or dimension input can be used.
% This function maintains 100% compatibility with the function vlcm,
% updated Nov-4-2014. Any code that uses vlcm can use elcm.
% 
% Note that the variable-precision function of vlcm can be retained using the following:
% c = vpa(elcm(inputs),digits)
inputs = [];
% Input handling
    if nargin == 0
        % No inputs were found.
            fprintf('No inputs were found for elcm.\n');
            return
    elseif nargin == 1
        % Only one input, potentially multiple elements.
            inputs = varargin{1};
            
        % Shortcut if there is only one input.
            if length(inputs) == 1
                c = inputs;
            end
    elseif nargin > 1
        % Multiple inputs, potentially with multiple elements.
        for i = 1:nargin
            inputvec = reshape(varargin{i},[1 numel(varargin{i})]);
            inputs = [inputs inputvec];
        end
    end
    
    
% The LCM must be a multiple of all numbers. Using that, it is simple to step through all
% multiples of the last number, stopping when the current multiple is divisible by all the
% other numbers. For prime number inputs, the lcm is the product of all the inputs, so this
% represents an upper bound on the LCM.
    for i = inputs(end):inputs(end):prod(inputs)
        
        % Sum the remainders when dividing the current multiple by all the other inputs.
            sum = 0;
            for j = 1:length(inputs)
                sum = sum + mod(i,inputs(j));
            end
        
        % Stop when the first common multiple is found.
            if sum == 0
                c = i;
                return
            end
    end
    
end
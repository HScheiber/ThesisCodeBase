function[output] = elements(varargin)

% USAGE:
% output = elements(input_specifier, input [, output_specifier])
%
% WHAT DOES IT DO:
% This function searches for the element specified, and
% returns specific information about that element.  Only one element can be
% looked up at a time.  The element is specified with it's name, symbol, or
% atomic number.  The user must specify, using the input_specifier
% variable, which of these she plans to use. 
% The user can also specify the information about the element they would
% like returned.  They can specify that they want the name, the symbol, the
% atomic mass, the atomic number, or all of these in a struct (See
% 'all_struct' notes below).  If the user does not specify a third
% paramater, the return value will be the struct of full information about
% the element.
%
% The input_specifier and output_specifier fields are always chars.  They
% can be any of the following.  
%
%   specifier             example input or output
%   -----------------------------------------------------------------------
%   name                        oxygen, iron
%   Symbol, Sym, symbol, sym    O, Fe  (first letter in caps)
%   SYMBOL, SYM                 O, FE  (all caps)
%   atomic_number               8, 26
%   atomic_mass(output only)    55.845, 15.9994   
%
%
% There is one additional output_specifer: 'all_struct' 'all_struct'
% returns a struct with keys: name, Symbol, SYMBOL, atomic_number, and
% atomic_mass. If you need more than once piece of information about an
% element, it is faster to use the 'all_struct' option
% 
% EXAMPLES:
% Example uses:
% output = elements(input_specifier, input, output_specifier)
% name = elements('Symbol', 'Fe', 'name');
% mass = elements('atomic_number', 26, 'atomic_mass');
% standardSymbol = elements('SYMBOL', 'FE', 'Symbol');
% 
% Example results:
% >> mass = elements('atomic_number', 26, 'atomic_mass')
% mass =  55.8450
%
%>>  name = elements('Symbol', 'Fe', 'name')
% name = iron
%
%
% NOTES:
% The associated file, elements.txt is required for this method.  You may
% modify this file if you have additional elements (!!).  
% 
% atomic_mass is not unique across all elements and hence can not be used
% as an input specifier
%
% all data taken from
% http://www.webelements.com/
%
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 

% the following variable is a global variable so that the file need only be
% read once in your matlab session
% The first time this method is called, it reads the file into an array and
% stores it in the following global variable.  The next time this script is
% called, it checks the existance of this variable.  If it exists, it
% allows will not read the file in, but instead access the data already in
% memory
global chem_elements_full_data_array;

% open file and read in data only if not already done
if isnumeric(chem_elements_full_data_array)
    chem_elements_full_data_array = open_file();
end

%%%%%%%%%%% specify arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  check that user sent correct number of arguments in
%  assign arguments a variable name to make code easier to read
if nargin==3
    output_specifier = varargin{3};
elseif (nargin == 2)
    output_specifier = 'all_struct';
else
    error('USAGE: output = elements(input_specifier, input [, output_specifier])');
end
input_specifier = varargin{1};
input = varargin{2};

% In determining the input_sepcifier and output_specifier, switches with
% list of possible symbols are used. The switch function is quite fast in
% Matlab.  This is faster than using string comparison methods.  

%%%%%%%%%%% get the information based on the input %%%%%%%%%%%%%%%%%%%%%%%%

% based on the input specifier, find the element
% and put information about it in the element struct
switch input_specifier
    case {'name','Name'}
        element_struct = fromName(lower(input));
    case {'Symbol','Sym','symbol','sym','SYMBOL','SYM'}
        element_struct = fromSymbol(toStandardSymbol(input));
    case {'atomic_number',...
            'atomic_Number',...
            'Atomic_number',...
            'Atomic_Number',...
            'atomicnumber',...
            'Atomicnumber',...
            'atomicNumber',...
            'AtomicNumber',...
            'z','Z'}
        element_struct = fromAtomicNumber(input);
    otherwise
        error('invalid input_specifier %s', input_specifier);
end

%%%%%%%% get the correct output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch output_specifier
    case {'name','Name'}
        output = element_struct.name;
    case {'Symbol','Sym','symbol','sym'}
        output = element_struct.Symbol;
    case {'SYMBOL','SYM'}
        output = upper(element_struct.Symbol) ;
    case {'atomic_number',...
            'atomic_Number',...
            'Atomic_number',...
            'Atomic_Number',...
            'atomicnumber',...
            'Atomicnumber',...
            'atomicNumber',...
            'AtomicNumber',...
            'z','Z'}
        output = element_struct.atomic_number;
    case {'atomic_mass',...
            'atomic_Mass',...
            'Atomic_mass',...
            'Atomic_Mass',...
            'atomicmass',...
            'Atomicmass',...
            'atomicMass',...
            'AtomicMass',...
            'm','M'}
        output = element_struct.atomic_mass;
    case{'all_struct'}
        % fix the struct slightly, remove the found field
        % which is only used programmtically
        element_struct = rmfield(element_struct, 'found');
        % go ahead and add in the SYMBOL field
        element_struct.SYMBOL = upper(element_struct.Symbol);
        output = element_struct;
    otherwise
        error('invalid input_specifier %s', output_specifier);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions below this point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function element_struct = fromAtomicNumber(atomicNumber)
    atmNumber_cell = chem_elements_full_data_array{3};
    element_struct.found = false;
    for i=1:length(atmNumber_cell)
        if  atmNumber_cell(i) == atomicNumber
            element_struct.name = chem_elements_full_data_array{1}{i};
            element_struct.Symbol = chem_elements_full_data_array{2}{i};
            element_struct.atomic_number = atomicNumber;
            element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
            element_struct.found = true;
            i=length(atmNumber_cell);
        end
    end
    if element_struct.found == false
        error('Unable to find element with atomic Number %d', atomicNumber);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function element_struct = fromName(name)    
    names_cell = chem_elements_full_data_array{1};
    element_struct.found = false;
    for i=1:length(names_cell)        
        if strcmp(names_cell{i}, name)            
            element_struct.name = name;
            element_struct.Symbol = chem_elements_full_data_array{2}{i};
            element_struct.atomic_number = chem_elements_full_data_array{3}(i);
            element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
            element_struct.found = true;
            i=length(names_cell);
        end
    end
    if element_struct.found == false
        error('Unable to find element with name %s', name);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function element_struct = fromSymbol(symbol)
    symbol_cell = chem_elements_full_data_array{2};
    element_struct.found = false;
    for i=1:length(symbol_cell)
        if strcmp(symbol_cell{i}, symbol)
            element_struct.name = chem_elements_full_data_array{1}{i};
            element_struct.Symbol = symbol;
            element_struct.atomic_number = chem_elements_full_data_array{3}(i);
            element_struct.atomic_mass = chem_elements_full_data_array{4}(i);
            element_struct.found = true;
            i=length(symbol_cell);
        end
    end
    if element_struct.found == false
        error('Unable to find element with symbol %s', symbol);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% NOTE: THIS FUNCTION WILL RETURN STRANGE RESULTS
% AtomicMass is not a unique identifier
% since several elements contain identical atomicMass
% if you insist on modifying this code to use this function,
% make sure the tolerance (variable "tol") is reasonable for your needs
function element_struct = fromAtomicMass(atomicMass)        
    mass_cell = chem_elements_full_data_array{4};
    element_struct.found = false;    
    
    tol = .01;
    for i=1:length(mass_cell)
              if  ( mass_cell(i) - tol <= atomicMass ) ...
              & ( mass_cell(i)  + tol >= atomicMass )
            
            element_struct.name = chem_elements_full_data_array{1}{i};
            element_struct.Symbol = chem_elements_full_data_array{2}{i};
            element_struct.atomic_number = chem_elements_full_data_array{3}(i);
            element_struct.atomic_mass = atomicMass;
            element_struct.found = true;
            i=length(mass_cell); 
        end
    end
     if element_struct.found == false
        error('Unable to find element with atomic mass %f', atomicMass);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function element_cell_array = open_file()
    fid = fopen('elements.txt', 'r');
    element_cell_array = textscan(fid, '%s %s %d %f', -1);

    if (fclose(fid) == -1)
        error('Failed to close the input data file');
    end
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function standardSymbol = toStandardSymbol(symbol)    
  
   % uppercase the first letter and lowercase everything else
   standardSymbol = upperCaseFirstLetter(symbol);
      
   % check the element name is valid (from data file)
   valid_element_symbols = chem_elements_full_data_array{2};
   found = false; 
   for j = 1:length(valid_element_symbols)        
       if strcmp(valid_element_symbols{j}, standardSymbol)
           found = true;
           j = length(valid_element_symbols);           
       end
   end

   if (found == false) 
      error('standardSymbol %s is not a valid element', standardSymbol);
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Word = upperCaseFirstLetter(word)
    % uppercase the first letter and lowercase everything else
    letters = regexpi(word, '(\w)(.*)', 'tokens');
    Word = strcat(upper(letters{1}{1}),lower(letters{1}{2}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% end elements function
end

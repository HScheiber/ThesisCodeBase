%% IMPORT DATA
FF_fields = fieldnames(FF_Data);
F_fields = fieldnames(F_Data);

FF_Lattice = FF_Data;
%% First Loop: Loop through FF Data
% Loop through levels of theory
for i = 1:numel(FF_fields)
    
    % Current calculation type
    FF_Calculation = regexprep(FF_fields{i},'^(.+?)_','');
    
    N = size(FF_Data.(FF_fields{i}),1); % Number of different basis sets in FF data

    for j = 1:N % Loop through FF Basis sets
        
        % Current FF CPU time and Basis set
        FF_Basis_set = FF_Data.(FF_fields{i}){j,1};
        
        % Number of different FF sub-methods
        M = size(FF_Data.(FF_fields{i}){j,2},1) - 1;
        
        for k = 1:M % Loop through FF theory types
            
            % Current Theory Level
            FF_Theory = FF_Data.(FF_fields{i}){j,2}{k,1};
            
            % Wipe FF lattice energy clear
            FF_Lattice.(FF_fields{i}){j,2}{k,3} = [];
            
            %% Looping through F data
            for x = 1:numel(F_fields)
                
                % Current F calculation type
                F_Calculation = regexprep(F_fields{x},'^(.+?)_','');
                
                % Match calculation types
                if ~strcmp(FF_Calculation,F_Calculation)
                    continue
                end
                
                % Number of different basis sets in F data
                NF = size(F_Data.(F_fields{x}),1);
            
                for y = 1:NF % Loop through F Basis sets

                    % Current F Basis set
                    F_Basis_set = F_Data.(F_fields{x}){y,1};
                    
                    % Match Basis sets
                    if ~strcmp(FF_Basis_set,F_Basis_set)
                        continue
                    end
                    
                    % Number of different F sub-methods
                    MF = size(F_Data.(F_fields{x}){y,2},1) - 1;
                
                    for z = 1:MF % Loop through F sub-methods
            
                        % Current sub-methods
                        F_Theory = F_Data.(F_fields{x}){y,2}{z,1};
                        
                        % Match sub-methods
                        if strcmp(FF_Theory,F_Theory)
                            
                            FF_Lattice.(FF_fields{i}){j,2}{k,3} = ...
                                FF_Data.(FF_fields{i}){j,2}{k,3} - ...
                                2*F_Data.(F_fields{x}){y,2}{z,2};
                        end
                    end
                end
            end
        end
    end
end

% Remove empty data points
for i = 1:numel(FF_fields)
    N = size(FF_Data.(FF_fields{i}),1); % Number of different basis sets in FF data
    for j = 1:N % Loop through FF Basis sets
        % Number of different FF sub-methods
        M = size(FF_Data.(FF_fields{i}){j,2},1) - 1;
        for k = M:-1:1 % Loop through FF theory types
            % Wipe FF lattice energy clear
            if isempty( FF_Lattice.(FF_fields{i}){j,2}{k,3} )
                FF_Lattice.(FF_fields{i}){j,2}(k,:) = [];
            end
        end
    end
end



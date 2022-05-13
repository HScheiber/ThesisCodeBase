%% IMPORT DATA
LiLi_fields = fieldnames(LiLi_Data);
Li_fields = fieldnames(Li_Data);

LiLi_Lattice = LiLi_Data;
%% First Loop: Loop through LiLi Data
% Loop through levels of theory
for i = 1:numel(LiLi_fields)
    
    % Current calculation type
    LiLi_Calculation = regexprep(LiLi_fields{i},'^(.+?)_','');
    
    N = size(LiLi_Data.(LiLi_fields{i}),1); % Number of different basis sets in LiLi data

    for j = 1:N % Loop through LiLi Basis sets
        
        % Current LiLi CPU time and Basis set
        LiLi_Basis_set = LiLi_Data.(LiLi_fields{i}){j,1};
        
        % Number of different LiLi sub-methods
        M = size(LiLi_Data.(LiLi_fields{i}){j,2},1) - 1;
        
        for k = 1:M % Loop through LiLi theory types
            
            % Current Theory Level
            LiLi_Theory = LiLi_Data.(LiLi_fields{i}){j,2}{k,1};
            
            % Wipe LiLi lattice energy clear
            LiLi_Lattice.(LiLi_fields{i}){j,2}{k,3} = [];
            
            %% Looping through F data
            for x = 1:numel(Li_fields)
                
                % Current F calculation type
                F_Calculation = regexprep(Li_fields{x},'^(.+?)_','');
                
                % Match calculation types
                if ~strcmp(LiLi_Calculation,F_Calculation)
                    continue
                end
                
                % Number of different basis sets in F data
                NF = size(Li_Data.(Li_fields{x}),1);
            
                for y = 1:NF % Loop through F Basis sets

                    % Current F Basis set
                    F_Basis_set = Li_Data.(Li_fields{x}){y,1};
                    
                    % Match Basis sets
                    if ~strcmp(LiLi_Basis_set,F_Basis_set)
                        continue
                    end
                    
                    % Number of different F sub-methods
                    MF = size(Li_Data.(Li_fields{x}){y,2},1) - 1;
                
                    for z = 1:MF % Loop through F sub-methods
            
                        % Current sub-methods
                        F_Theory = Li_Data.(Li_fields{x}){y,2}{z,1};
                        
                        % Match sub-methods
                        if strcmp(LiLi_Theory,F_Theory)
                            
                            LiLi_Lattice.(LiLi_fields{i}){j,2}{k,3} = ...
                                LiLi_Data.(LiLi_fields{i}){j,2}{k,3} - ...
                                2*Li_Data.(Li_fields{x}){y,2}{z,2};
                        end
                    end
                end
            end
        end
    end
end

% Remove empty data points
for i = 1:numel(LiLi_fields)
    N = size(LiLi_Data.(LiLi_fields{i}),1); % Number of different basis sets in LiLi data
    for j = 1:N % Loop through LiLi Basis sets
        % Number of different LiLi sub-methods
        M = size(LiLi_Data.(LiLi_fields{i}){j,2},1) - 1;
        for k = M:-1:1 % Loop through LiLi theory types
            % Wipe LiLi lattice energy clear
            if isempty( LiLi_Lattice.(LiLi_fields{i}){j,2}{k,3} )
                LiLi_Lattice.(LiLi_fields{i}){j,2}(k,:) = [];
            end
        end
    end
end



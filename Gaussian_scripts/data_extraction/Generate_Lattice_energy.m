%% Parameters
function Data = Generate_Lattice_energy(Salt,crystal_type,datadir)

% Preprocessing
[match,~] = regexp(Salt,{'Li' 'F' 'Cl' 'Br' 'I'},'match','tokens');
matches = match(~cellfun('isempty',match));
metal = matches{1}{1};
halide = matches{2}{1};

%% Load Data
try
if strcmp(crystal_type,'Wurtzite')
    Label = [Salt '_W'];
    X = load([datadir filesep Salt '_Wurtzite_Total_Energies.mat'],'Data');
    Salt_Data = X.Data;
    Y = load([datadir filesep metal '_Ion_Total_Energies_DFT.mat'],'Data');
    Metal_Data = Y.Data;
    Z = load([datadir filesep halide '_Ion_Total_Energies_DFT.mat'],'Data');
    Halide_Data = Z.Data;
    clearvars X Y Z
elseif strcmp(crystal_type,'Rocksalt')
    Label = [Salt '_R'];
    X = load([datadir filesep Salt '_Rocksalt_Total_Energies.mat'],'Data');
    Salt_Data = X.Data;
    Y = load([datadir filesep metal '_Ion_Total_Energies_DFT.mat'],'Data');
    Metal_Data = Y.Data;
    Z = load([datadir filesep halide '_Ion_Total_Energies_DFT.mat'],'Data');
    Halide_Data = Z.Data;
    clearvars X Y Z
end
catch
    error(['Unable to load selected data, check ' datadir])
end

Salt_fields = fieldnames(Salt_Data);
Halide_fields = fieldnames(Halide_Data);
Metal_fields = fieldnames(Metal_Data);

Salt_Lattice = Salt_Data;
%% First Loop: Loop through LiF Data
% Loop through levels of theory
for i = 1:numel(Salt_fields)
    
    % Current calculation type
    Salt_Calculation = strrep(Salt_fields{i},[Label '_'],'');
    
    N = size(Salt_Data.(Salt_fields{i}),1); % Number of different basis sets in LiF data

    for j = 1:N % Loop through LiF Basis sets
        
        % Current LiF CPU time and Basis set
        Salt_Basis_set = Salt_Data.(Salt_fields{i}){j,1};
        
        % Number of different LiF sub-methods
        M = size(Salt_Data.(Salt_fields{i}){j,2},1) - 1;
        
        for k = 1:M % Loop through LiF theory types
            
            % Current Theory Level
            Salt_Theory = Salt_Data.(Salt_fields{i}){j,2}{k,1};
            
            % Wipe LiF lattice energy clear
            Salt_Lattice.(Salt_fields{i}){j,2}{k,3} = [];
            
            %% Looping through F data
            for x = 1:numel(Halide_fields)
                
                % Current F calculation type
                Halide_Calculation = regexprep(Halide_fields{x},'^(.+?)_','');
                
                % Match calculation types
                if ~strcmp(Salt_Calculation,Halide_Calculation)
                    continue
                end
                
                % Number of different basis sets in F data
                NF = size(Halide_Data.(Halide_fields{x}),1);
            
                for y = 1:NF % Loop through F Basis sets

                    % Current F Basis set
                    Halide_Basis_set = Halide_Data.(Halide_fields{x}){y,1};
                    
                    % Match Basis sets
                    if ~strcmp(Salt_Basis_set,Halide_Basis_set)
                        if strcmp(Salt_Basis_set,'pob-TZVP') && strcmp(Halide_Basis_set,'Def2TZVP')
                        else
                            continue
                        end
                    elseif strcmp(Salt_Basis_set,'pob-TZVP') && strcmp(Halide_Basis_set,'pob-TZVP')
                        continue
                    end
                    
                    % Number of different F sub-methods
                    MF = size(Halide_Data.(Halide_fields{x}){y,2},1) - 1;
                
                    for z = 1:MF % Loop through F sub-methods
            
                        % Current sub-methods
                        Halide_Theory = Halide_Data.(Halide_fields{x}){y,2}{z,1};
                        
                        % Match sub-methods
                        if strcmp(Salt_Theory,Halide_Theory)
                            
                            if strcmp(crystal_type,'Pair')
                                Salt_Lattice.(Salt_fields{i}){j,2}{k,3} = ...
                                    Salt_Data.(Salt_fields{i}){j,2}{k,3} - ...
                                    Halide_Data.(Halide_fields{x}){y,2}{z,2} - ...
                                    Metal_Data.(Metal_fields{x}){y,2}{z,2};
                            elseif strcmp(crystal_type,'Wurtzite')
                                Salt_Lattice.(Salt_fields{i}){j,2}{k,3} = ...
                                    (Salt_Data.(Salt_fields{i}){j,2}{k,3} - ...
                                    2*Halide_Data.(Halide_fields{x}){y,2}{z,2} - ...
                                    2*Metal_Data.(Metal_fields{x}){y,2}{z,2})./2;
                            elseif strcmp(crystal_type,'Rocksalt')
                                Salt_Lattice.(Salt_fields{i}){j,2}{k,3} = ...
                                    (Salt_Data.(Salt_fields{i}){j,2}{k,3} - ...
                                    4*Halide_Data.(Halide_fields{x}){y,2}{z,2} - ...
                                    4*Metal_Data.(Metal_fields{x}){y,2}{z,2})./4;
                            else
                                error('Unknown Crystal Type.')
                            end
                        end
                    end
                end
            end
        end
    end
end

% Remove empty data points
for i = 1:numel(Salt_fields)
    N = size(Salt_Data.(Salt_fields{i}),1); % Number of different basis sets in LiF data
    for j = 1:N % Loop through LiF Basis sets
        % Number of different LiF sub-methods
        M = size(Salt_Data.(Salt_fields{i}){j,2},1) - 1;
        for k = M:-1:1 % Loop through LiF theory types
            % Wipe LiF lattice energy clear
            if isempty( Salt_Lattice.(Salt_fields{i}){j,2}{k,3} )
                Salt_Lattice.(Salt_fields{i}){j,2}(k,:) = [];
            end
        end
    end
end

%% Save
Data = Salt_Lattice;

save([datadir filesep Salt '_' crystal_type '_Lattice_Energies.mat'],'Data')

end

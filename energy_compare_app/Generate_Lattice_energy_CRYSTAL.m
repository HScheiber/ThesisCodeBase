%% Parameters
function Salt_Lattice = Generate_Lattice_energy_CRYSTAL(Salt,Label,Structure,datadir,Calc_type)
% Default reference ion basis set
General_Reference_Basis = 'def2-TZVPD';

% Reference ion basis when gCP is active
gCP_Reference_Basis = 'def2-TZVPD';

General_Crystal_Basis = {'pob-TZVP' 'STO-3G' '7-311G' 'aug-cc-pV5Z' 'cc-pV5Z' 'def2-TZVPD',...
    'Mpob-TZVP_Hpob-TZVPD' 'Mpob-TZVPP_Hpob-TZVPD' 'Mpob-TZVP-rev2_Hpob-TZVP' ...
    'Mpob-TZVPP_Hpob-TZVP'};

% Define axuiliary basis set for I atoms involving the ECP28MDF type basis set
I_Atom_Aux_Basis = 'ECP28MDF-AVTZ';
LiI_Crystal_Aux_Basis = {'Mpob-TZVP_HECP28MDF-VTZ' 'Mpob-TZVP_HECP28MDF-AVTZ' ...
    'Mpob-TZVP_HECP28MDF-V5Z' 'Mpob-TZVP_HECP28MDF-AV5Z' ...
    'Mpob-TZVPP_HECP28MDF-VTZ' 'Mpob-TZVPP_HECP28MDF-AVTZ' ...
    'Mpob-TZVPP_HECP28MDF-V5Z' 'Mpob-TZVPP_HECP28MDF-AV5Z'};

% Define tertiary basis set for atoms involving the def2-QZVPD type basis set
Tert_Basis = 'def2-QZVPD';
Crystal_Tert_Basis = {'Mpob-TZVP_Hdef2-QZVP'};

% Preprocessing
[Metal,Halide] = Separate_Metal_Halide(Salt);

if isempty(Metal) && strcmp(Structure,'Pair')
    SkipMetal = true;
    SkipHalide = false;
    QN = 2;
elseif isempty(Halide) && strcmp(Structure,'Pair')
    SkipMetal = false;
    SkipHalide = true;
    QN = 2;
else
    SkipMetal = false;
    SkipHalide = false;
    QN = 1;
end

%% Load Data
try
    X = load([datadir filesep Salt '_' Structure '_Total_Energies.mat'],'Data');
    Salt_Data = X.Data;
    
    if strcmp(Calc_type,'Lattice')
        if ~isempty(Metal)
            Y = load([datadir filesep Metal '_Ion_Total_Energies.mat'],'Data');
        else
            Y.Data = struct();
        end
        if ~isempty(Halide)
            Z = load([datadir filesep Halide '_Ion_Total_Energies.mat'],'Data');
        else
            Z.Data = struct();
        end
    elseif strcmp(Calc_type,'Cohesive')
        if ~isempty(Metal)
            Y = load([datadir filesep Metal '_Atom_Total_Energies.mat'],'Data');
        else
            Y.Data = struct();
        end
        if ~isempty(Halide)
            Z = load([datadir filesep Halide '_Atom_Total_Energies.mat'],'Data');
        else
            Z.Data = struct();
        end
    else
        error('Error. Unable to load Single Atom/Ion Energies.');
    end
    
    Metal_Data = Y.Data;
    Halide_Data = Z.Data;
    clearvars X Y Z
catch
    disp(['Unable to load ' Salt ' ' Structure ' data for ' Calc_type ...
        ' energy calculation. Skipping.'])
    Salt_Lattice = [];
    return
end

Salt_fields = fieldnames(Salt_Data);
Halide_fields = fieldnames(Halide_Data);
Metal_fields = fieldnames(Metal_Data);

Salt_Lattice = Salt_Data;
%% First Loop: Loop through levels of theory
K = numel(Salt_fields);
for i = 1:K
    
    % Current calculation type
    Salt_Calculation = strrep(Salt_fields{i},[Label '_'],'');
    
    % Note if gCP corrected calculation
    gCP_Corrected = ~isempty(regexp(Salt_Calculation,'gCP','ONCE'));
    
    Salt_Calculation = regexprep(Salt_Calculation,'_D[2-4].*','');
    Salt_Calculation = regexprep(Salt_Calculation,'_gCP','');
    
    N = size(Salt_Data.(Salt_fields{i}),1); % Number of different basis sets in LiF data

    for j = N:-1:1 % Loop through Basis sets
        
        % Current Salt Basis set
        Salt_Basis_set = Salt_Data.(Salt_fields{i}){j,1};
        Salt_M_Basis_set = regexprep(regexprep(Salt_Basis_set,'^M',''),'_H.+$','');
        Salt_H_Basis_set = regexprep(Salt_Basis_set,'^M.+?_H','');

        matched_Halide = SkipHalide;
        matched_Metal = SkipMetal;

        %% Looping through Halide data
        if ~SkipHalide
            for x = 1:numel(Halide_fields)

                % Current halide calculation type
                Halide_Calculation = regexprep(Halide_fields{x},'^(.+?)_','');

                % Match calculation types
                if ~strcmp(Salt_Calculation,Halide_Calculation)
                    continue
                end

                % Number of different basis sets in Halide data
                N_Hal = size(Halide_Data.(Halide_fields{x}),1);

                for y = 1:N_Hal % Loop through halide Basis sets

                    % Current halide Basis set
                    Halide_Basis_set = Halide_Data.(Halide_fields{x}){y,1};

                    % Check for most accurate basis set for ions (def2-TZVPD)
                    if strcmp(Structure,'Pair')
                        if strcmp(Salt_Basis_set,Halide_Basis_set)
                            matched_Halide = true;
                        elseif strcmp(Salt_H_Basis_set,Halide_Basis_set)
                            matched_Halide = true;
                        else
                            continue
                        end
                    elseif strcmp(Halide_Basis_set,General_Reference_Basis) && ...
                            any(strcmp(Salt_Basis_set,General_Crystal_Basis)) && ~gCP_Corrected
                        matched_Halide = true;
                    elseif strcmp(Halide_Basis_set,gCP_Reference_Basis) && ...
                            any(strcmp(Salt_Basis_set,General_Crystal_Basis)) && gCP_Corrected
                        matched_Halide = true;
                    elseif strcmp(Halide_Basis_set,I_Atom_Aux_Basis) && any(strcmp(Salt_Basis_set,LiI_Crystal_Aux_Basis)) && strcmp(Halide,'I')
                        matched_Halide = true;
                    elseif strcmpi(Halide_Basis_set,Tert_Basis) && any(strcmpi(Salt_Basis_set,Crystal_Tert_Basis))
                        matched_Halide = true;
                    else
                        continue
                    end

                    Data_Types = Salt_Lattice.(Salt_fields{i}){j,2}{1,10};
                    P1_Index = (Data_Types == 4);

                    if contained_in_cell(Structure,{'CsCl' 'Pair'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            QN.*Halide_Data.(Halide_fields{x}){y,2}{1,2});
                    elseif contained_in_cell(Structure,{'Rocksalt' 'Sphalerite'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(~P1_Index) = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(~P1_Index) - ...
                            Halide_Data.(Halide_fields{x}){y,2}{1,2});
                        
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(P1_Index) = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(P1_Index) - ...
                            4.*Halide_Data.(Halide_fields{x}){y,2}{1,2});
                    elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            2*Halide_Data.(Halide_fields{x}){y,2}{1,2};
                    elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            4*Halide_Data.(Halide_fields{x}){y,2}{1,2});
                    else
                        error(['Unknown Crystal Type: ' Structure])
                    end
                    
                    if ~isempty(Salt_Lattice.(Salt_fields{i}){j,5})
                        for z = 1:size(Salt_Lattice.(Salt_fields{i}){j,5},1)
                            Salt_Lattice.(Salt_fields{i}){j,5}{z,16} = ...
                                Salt_Lattice.(Salt_fields{i}){j,5}{z,16} - ...
                                Halide_Data.(Halide_fields{x}){y,2}{1,2};
                        end
                    end
                    
                    break
                end
                if matched_Halide
                    break
                end
            end
        end
        %% Remove data that did not match with halide
        if ~matched_Halide
            Salt_Lattice.(Salt_fields{i})(j,:) = [];
            continue
        end
        
        %% Looping through Metal data
        if ~SkipMetal
            for x = 1:numel(Metal_fields)

                % Current Metal calculation type
                Metal_Calculation = regexprep(Metal_fields{x},'^(.+?)_','');

                % Match calculation types
                if ~strcmp(Salt_Calculation,Metal_Calculation)
                    continue
                end

                % Number of different basis sets in Metal data
                N_Met = size(Metal_Data.(Metal_fields{x}),1);

                for y = 1:N_Met % Loop through metal Basis sets

                    % Current Metal Basis set
                    Metal_Basis_set = Metal_Data.(Metal_fields{x}){y,1};

                    % Check for basis set match (pob-TZVP with def2-TZVPD
                    % considered a match for crystals)
                    if strcmp(Structure,'Pair')
                        if strcmp(Salt_Basis_set,Metal_Basis_set)
                            matched_Metal = true;
                        elseif strcmp(Salt_M_Basis_set,Metal_Basis_set)
                            matched_Metal = true;
                        else
                            continue
                        end
                    elseif strcmp(Metal_Basis_set,General_Reference_Basis) && ...
                            any(strcmp(Salt_Basis_set,General_Crystal_Basis)) && ~gCP_Corrected
                        matched_Metal = true;
                    elseif strcmp(Metal_Basis_set,gCP_Reference_Basis) && ...
                            any(strcmp(Salt_Basis_set,General_Crystal_Basis)) && gCP_Corrected
                        matched_Metal = true;
                    elseif strcmp(Metal_Basis_set,General_Reference_Basis) && any(strcmp(Salt_Basis_set,LiI_Crystal_Aux_Basis)) && strcmp(Halide,'I')
                        matched_Metal = true;
                    elseif strcmpi(Metal_Basis_set,General_Reference_Basis) && any(strcmpi(Salt_Basis_set,Crystal_Tert_Basis))
                        matched_Metal = true;
                    else
                        continue
                    end

                    Data_Types = Salt_Lattice.(Salt_fields{i}){j,2}{1,10};
                    P1_Index = (Data_Types == 4);
                    
                    if contained_in_cell(Structure,{'CsCl' 'Pair'})
                         Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            QN.*Metal_Data.(Metal_fields{x}){y,2}{1,2});
                    elseif contained_in_cell(Structure,{'Rocksalt' 'Sphalerite'})
                         Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(~P1_Index) = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(~P1_Index) - ...
                            Metal_Data.(Metal_fields{x}){y,2}{1,2});
                        
                         Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(P1_Index) = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8}(P1_Index) - ...
                            4.*Metal_Data.(Metal_fields{x}){y,2}{1,2})./4;
                    elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            2*Metal_Data.(Metal_fields{x}){y,2}{1,2})./2;
                    elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
                        Salt_Lattice.(Salt_fields{i}){j,2}{1,8} = ...
                            (Salt_Lattice.(Salt_fields{i}){j,2}{1,8} - ...
                            4*Metal_Data.(Metal_fields{x}){y,2}{1,2})./4;
                    else
                        error(['Unknown Crystal Type: ' Structure])
                    end
                    
                    if ~isempty(Salt_Lattice.(Salt_fields{i}){j,5})
                        for z = 1:size(Salt_Lattice.(Salt_fields{i}){j,5},1)                           
                            Salt_Lattice.(Salt_fields{i}){j,5}{z,16} = ...
                                Salt_Lattice.(Salt_fields{i}){j,5}{z,16} - ...
                                Metal_Data.(Metal_fields{x}){y,2}{1,2};
                        end
                    end
                end
            end
        end
        
        %% Remove data that did not match with metal
        if ~matched_Metal
            Salt_Lattice.(Salt_fields{i})(j,:) = [];
        end
        
    end
    
    %% Prune empty data
    if isempty(Salt_Lattice.(Salt_fields{i}))
    	Salt_Lattice = rmfield(Salt_Lattice,Salt_fields{i});
    end
    
    if (i/K) < 0.9999
        fprintf('%s',[num2str(round(i*100/K)) '%... '])
    else
        fprintf('%s',['100%.' newline])
    end
end

end

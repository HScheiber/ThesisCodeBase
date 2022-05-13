%% Script to check if lowest energy structure is rocksalt
Salts = {'LiI'};
Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' 'BetaBeO' 'FiveFive'};
% = {'HSE06' 'HSEsol' 'M06' 'PBE' 'PW1PW' 'B3LYP'};
Basis_Set = 'pob-TZVP';
%Models = {'HF' 'HSE06' 'HSEsol' 'LSDA' 'M06' 'M062X' 'M06L' 'PBE' 'PW1PW' 'B3LYP'};
Models = {'PBE' 'PW1PW' 'B3LYP'};
D3 = true;
gCP = false;
Data_Type = 2;
Data_Dir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\CRYSTAL';

%% Lengths and stuff
N_Salt = length(Salts);
N_Struc = length(Structures);
N_Model = length(Models);
Ha_kJ = 4.35974e-21; % kiloJoules per Hartree
NA = 6.0221409e23; % Molecules per mole
Ha_kJ_NA = Ha_kJ*NA; % kiloJoule per mole per Hartree
Num_Rocksalt = 0;
Num_Other = 0;
num_skipped = 0;
Skipped_Data = cell(0,1);

% Loop salts
for ii = 1:N_Salt
    Salt = Salts{ii};

    % Loop models
    for jj = 1:N_Model
        Model = Models{jj};

        if D3
            Model = [Model '_D3'];
        end
        if gCP
            Model = [Model '_gCP'];
        end

        % Pre-create array of energies for structures
        Energies = nan(1,N_Struc);
        Temp_Strc_List = cell(1,N_Struc);
        
        % Loop structures
        for kk = 1:N_Struc           
            Structure = Structures{kk};
            Structure_Tag = Structure(1);
            
            Data_Flag = [Salt '_' Structure_Tag '_' Model];

            % Load Data
            DATstrc = load(fullfile(Data_Dir,[Salt '_' Structure '_Lattice_Energies.mat']),'Data');
            Salt_Strc_Data = DATstrc.Data;
            
            if ~isfield(Salt_Strc_Data,Data_Flag)
                num_skipped = num_skipped + 1;
                Skipped_Data{end+1} = [Salt ' ' Model ' ' Structure];
                continue
            end
            Data_Model = Salt_Strc_Data.(Data_Flag);
            
            N_Basis = size(Data_Model,1);
            
            % Find basis set data
            Salt_Model_Struct_Dat = [];
            for ll = 1:N_Basis
                if strcmp(Data_Model{ll,1},Basis_Set)
                    Salt_Model_Struct_Dat = Data_Model{ll,2};
                    break
                end
            end
            if isempty(Salt_Model_Struct_Dat)
                num_skipped = num_skipped + 1;
                Skipped_Data{end+1} = [Salt ' ' Model ' ' Structure];
                continue
            end
            
            % Get Data Type
            DT_Index = Salt_Model_Struct_Dat{10} == Data_Type;
            
            Energy = Salt_Model_Struct_Dat{8}(DT_Index)*Ha_kJ_NA;
            
            if isempty(Energy)
                num_skipped = num_skipped + 1;
                Skipped_Data{end+1} = [Salt ' ' Model ' ' Structure];
                continue
            end
            Energies(kk) = Energy;
            
            % Check if structure converged to different one
            if strcmp(Structure,'Wurtzite')
                FC = Salt_Model_Struct_Dat{7}{DT_Index};
                DZ = abs(mod(FC{1,4},1) - mod(FC{2,4},1));
                
                if abs(DZ - 0.5) < 1e-2
                    Temp_Strc_List{kk} = 'FiveFive';
                else
                    Temp_Strc_List{kk} = Structure;
                end
                
            elseif strcmp(Structure,'BetaBeO')
                a = Salt_Model_Struct_Dat{2}(DT_Index);
                c = Salt_Model_Struct_Dat{4}(DT_Index);
                if abs(a - c) < 1e-2
                    Temp_Strc_List{kk} = 'Rocksalt';
                else
                    Temp_Strc_List{kk} = Structure;
                end
            else
                Temp_Strc_List{kk} = Structure;
            end
        end
        
        % Find lowest energy structure
        [~,Min_Idx] = min(Energies);
        Min_Struc = Temp_Strc_List{Min_Idx};
        
        if strcmp(Min_Struc,'Rocksalt')
            Num_Rocksalt = Num_Rocksalt+1;
        else
            Num_Other = Num_Other+1;
        end
    end
end


disp(['Correct structure: ' num2str(Num_Rocksalt) '/' num2str(Num_Rocksalt+Num_Other) ' ( ' num2str(100*Num_Rocksalt/(Num_Rocksalt+Num_Other),'%2.2f') '% )']);
disp(['Skipped ' num2str(num_skipped) ' Data points:']);
for ii = 1:length(Skipped_Data)
    disp(Skipped_Data{ii});
end
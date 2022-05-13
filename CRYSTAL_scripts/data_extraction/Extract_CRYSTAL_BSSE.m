% SCRIPT FOR EXTRACTING BASSIS SET SUPERPOSITION ERROR ESTIMATES FROM CRYSTAL OUTPUT
% (*.out) FILES.

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'Li', 'F', etc.

function Data = Extract_CRYSTAL_BSSE(Label,Lone_Ion_Data)
%% Level 0 (analysis folder)

% Move to BSSE Estimation folder
BSSE_Folder = find_directories(dir('BSSE*'));
cur_dir = [BSSE_Folder{1} filesep Label];

% Find the folder that matches the label
if ~exist(cur_dir,'dir') 
    Data = [];
    disp(['Warning: No BSSE Data for ' Label '. Excluding from Dataset.'])
    return
end

% Determine if Label is Metal or Halide
if contained_in_cell(Label,{'Li' 'Na' 'K' 'Rb'})
    IonType = 'Metal';
elseif contained_in_cell(Label,{'F' 'Cl' 'Br' 'I'})
    IonType = 'Halide';
else
    error(['Unknown ion type: ' Label])
end

%% Level 1 (outer)

% Find directories in current folder, except . and ..
Folders = find_directories(dir(cur_dir));

% Find unique labels
matches = regexp(Folders,'^.+_.+(?=_)','match');
M = length(matches);
labels = cell(M,1);
for i=1:M
    labels{i} = matches{i}{1};
end
labels = unique(labels);

%% Extended Basis Set System
Data = struct();
R = length(labels);
for i = 1:R
    Structure = regexprep(char(labels(i)),'.+_R','Rocksalt');
    Structure = regexprep(Structure,'.+_W','Wurtzite');
    Structure = regexprep(Structure,'.+_S','Sphalerite');
    Structure = regexprep(Structure,'.+_C','CsCl');
    Structure = regexprep(Structure,'.+_N','NiAs');
    Structure = regexprep(Structure,'.+_B','BetaBeO');
    Structure = regexprep(Structure,'.+_F','FiveFive');
    Structure = regexprep(Structure,'.+_P','Pair');
    labelText = [Structure ' Structure'];
    
    disp(['BSSE correction energy for ' Label ' in ' labelText '...'])
    Data_part = Extract_CRYSTAL_PES(labels{i},cur_dir,Structure,true);
    
    f = fieldnames(Data_part);
    for j = 1:length(f)
    	Data.(f{j}) = Data_part.(f{j});
    end
    disp(['Extracted (' num2str(i) '/' num2str(R) ').'])
end

%% Subtract Original Basis Set Energy
fields = fieldnames(Data);
R = length(fields);
disp(['Subtracting energy from single-ion basis set for ' Label ' BSSE Correction...'])
for i = 1:R % Loop through Salt and calculation types
    
    ct = textscan(fields{i},'%*s %*s %s','Delimiter','_');
    calc_type = ct{1}{1};
    
    % Match with lone ion
    Lone_Ion_Fields = fieldnames(Lone_Ion_Data);
    matched = false;
    for j = 1:length(Lone_Ion_Fields)
        Lone_Ion_ct = textscan(Lone_Ion_Fields{j},'%*s %s','Delimiter','_');
        Lone_Ion_calc_type = Lone_Ion_ct{1}{1};
        
        if strcmp(calc_type,Lone_Ion_calc_type)
            Matched_Lone_Ion_Data = Lone_Ion_Data.(Lone_Ion_Fields{j});
            matched = true;
            break
        end
    end
    if ~matched
        disp(['Warning: No BSSE match for ' calc_type '. Pruning from dataset.'])
        Data = rmfield(Data,fields{i});
        continue
    end
    
    % Loop through basis sets
    N = size(Data.(fields{i}),1);
    for j = N:-1:1
        basis_set = Data.(fields{i}){j,1};
        if strcmp(IonType,'Metal')
            basis_set = regexprep(regexprep(basis_set,'_H.+$',''),'^M','');
        elseif strcmp(IonType,'Halide')
            basis_set = regexprep(basis_set,'^M.+_H','');
        end
        
        % Match with lone ion
        Lone_Ion_Bases = Matched_Lone_Ion_Data(:,1);
        matched = false;
        for k = 1:length(Lone_Ion_Bases)

            if strcmp(basis_set,Lone_Ion_Bases{k})
                Matched_Lone_Ion_Energy = Matched_Lone_Ion_Data{k,2}{2};
                matched = true;
                break
            end
        end
        
        if ~matched
            disp(['Warning: No BSSE match for ' calc_type ' + ' basis_set '. Pruning from dataset.'])
            Data.(fields{i}){j,1} = [];
            if isempty(Data.(fields{i}))
                Data = rmfield(Data,fields{i});
            end
            continue
        end
        
         Data.(fields{i}){j,2}{:,8} = (Data.(fields{i}){j,2}{:,8} - Matched_Lone_Ion_Energy);
    end
end

end
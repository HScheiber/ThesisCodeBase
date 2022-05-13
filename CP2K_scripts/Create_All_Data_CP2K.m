% Load CP2K Dataset
Data_Directory = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data';
CP2K_Data_Obj = load(fullfile(Data_Directory,'CP2K_Data.mat'));
Data = CP2K_Data_Obj.Data;

fileID = fopen('CP2K_Data.csv','wt');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
    'Salt',...
    'Structure',...
    'Theory',...
    'Metal Basis Set',...
    'Halide Basis Set',...
    'Lattice Energy',...
    'a(0 K) [Angstrom]',...
    'b(0 K) [Angstrom]',...
    'c(0 K) [Angstrom]',...
    'a(298 K) [Angstrom]',...
    'b(298 K) [Angstrom]',...
    'c(298 K) [Angstrom]',...
    'Metal x',...
    'Metal y',...
    'Metal z',...
    'Halide x',...
    'Halide y',...
    'Halide z'...
    );

Basis_Sets = {'Sapporo_QZP' 'pob_TZVP'};
Salts_Reg = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Salts_DKH = {'LiBr' 'LiI'};
Structures = {'FiveFive' 'NiAs' 'Rocksalt' 'Sphalerite' 'Wurtzite'};

Frac_Coords.Wurtzite.Metal = [1/3 2/3 0];
Frac_Coords.Wurtzite.Halide = [1/3 2/3 0.375];

Frac_Coords.NiAs.Metal = [0 0 0];
Frac_Coords.NiAs.Halide = [1/3 2/3 1/4];

Frac_Coords.Rocksalt.Metal = [0 0 0];
Frac_Coords.Rocksalt.Halide = [1/2 1/2 1/2];

Frac_Coords.Sphalerite.Metal = [0 0 0];
Frac_Coords.Sphalerite.Halide = [1/4 1/4 1/4];

Frac_Coords.FiveFive.Metal = [1/4 1/6 1/2];
Frac_Coords.FiveFive.Halide = [1/4 1/3 0];

for BS_idx = 1:length(Basis_Sets)
    Basis = Basis_Sets{BS_idx};
    Metal_Basis = 'pob_TZVP';
    
    Theories = fieldnames(Data.(Basis));
    
    for TH_idx = 1:length(Theories)
        Theory = Theories{TH_idx};
        
        is_DKH = ~isempty(regexp(Theory,'DKH','once'));
        
        if is_DKH
            Salts = Salts_DKH;
        else
            Salts = Salts_Reg;
        end
        
        for SL_idx = 1:length(Salts)
            
            Salt = Salts{SL_idx};
            switch Basis
                case 'Sapporo_QZP'
                    if is_DKH
                        Halide_Basis = 'Sapporo-QZP-DKH3';
                    else
                        Halide_Basis = 'Sapporo-QZP';
                    end
                case 'pob_TZVP'
                    if is_DKH
                        Halide_Basis = 'Sapporo-TZP-DKH3';
                    else
                        switch Salt
                            case 'LiI'
                                Halide_Basis = 'Sapporo-TZP';
                            otherwise
                                Halide_Basis = 'pob_TZVP';
                        end
                    end
            end
            
            for SR_idx = 1:length(Structures)
                Structure = Structures{SR_idx};
                
                Lattice_E = Data.(Basis).(Theory).(Salt).(Structure).LE;
                a = Data.(Basis).(Theory).(Salt).(Structure).a;
                b = Data.(Basis).(Theory).(Salt).(Structure).b;
                c = Data.(Basis).(Theory).(Salt).(Structure).c;
                FC_M = Frac_Coords.(Structure).Metal;
                FC_X = Frac_Coords.(Structure).Halide;
                
                fprintf(fileID, '%s,%s,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',...
                    Salt, Structure, strrep(strrep(Theory,'_','-'),'DKH','DKH3'), Metal_Basis, Halide_Basis, Lattice_E,...
                    a,b,c,[],[],[],FC_M(1),FC_M(2),FC_M(3),FC_X(1),FC_X(2),FC_X(3));
            end
        end
    end 
end
fclose(fileID);



%% Build energy table
Structure_List = {'Rocksalt' 'Wurtzite' 'BetaBeO' 'CsCl' 'FiveFive' 'NiAs' 'Sphalerite'};
N_SL = length(Structure_List);
Note = '';
BS = '$^{\ddagger}$';

% Reload the data table
T = readtable('CP2K_Data.csv');

% Index pob-TZVP basis sets
pob_idx = strcmp(T.HalideBasisSet,'pob_TZVP');

% Prune data
T = T(~pob_idx,[1:3 5:6]);

% Sapporo-TZP basis sets
TZP_idx = strcmp(T.HalideBasisSet,'Sapporo-TZP');

% Prune data
T = T(~TZP_idx,:);

% Sapporo-TZP-DKH3 basis sets
TZP_idx = strcmp(T.HalideBasisSet,'Sapporo-TZP-DKH3');

% Prune data
Data_Table = T(~TZP_idx,[1:3 5]);
Data_Table = sortrows(Data_Table,{'Salt' 'Theory' 'Structure'},{'ascend','ascend','ascend'});   

% Find the unique list of species
Q = size(unique(Data_Table(:,1)),1);
TableAdj = ['{l|l' repmat('l', 1,Q) '}'];
    
N = size(Data_Table,1);
    
% Initialize previous salt empty string
Previous_Salt = '';
Energy_Table = {};
    
% Loop through all data entries
for Idx = 1:N
    Current_Salt = Data_Table{Idx,1}{1};
    Current_Structure = Data_Table{Idx,2}{1};
    Current_Theory = Data_Table{Idx,3}{1};
    Current_E = num2str(Data_Table{Idx,4},'%6.2f');
    Current_E_num = Data_Table{Idx,4};

    [~,Str_Index] = ismember(Current_Structure,Structure_List);
    if Str_Index == 0 
        continue
    end

    % Check if salt has changed
    if ~strcmp(Current_Salt,Previous_Salt)

        % Previous Energy Table Complete
        if Idx ~= 1

            % Finish up last row
            [~,MinIdx] = min(LineEnergies);
            LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];

            Favoured = Structure_List{MinIdx};
            Favoured = strrep(Favoured,'Rocksalt','Rock.');
            Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
            Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
            Favoured = strrep(Favoured,'FiveFive','5-5');
            Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
            LinesText{x} = [LinesText{x} '& ' Favoured];

            for i = 1:length(LinesData)
                LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
                LinesText{x} = [LinesText{x} '&' LinesData{i}];
            end

            % Finish title
            for i = 1:N_SL
                LineTitle = [LineTitle pad(['& \textbf{' Shorten_Structure_Name(Structure_List{i}) '} '],23,'right')];
            end

            % Finish up all rows
            LengthT = length(LineTitle);
            LineTitle = [LineTitle '\\ \hline'];
            for i = 1:length(LinesText)
                if i < length(LinesText)
                    LinesText{i} = ['            ' LinesText{i} '\\'];
                else
                    LinesText{i} = ['            ' LinesText{i}];
                end
            end

            % Build table
            Energy_Table{end+1} = [...
                Line1 newline ...
                Line2 newline ...
                Line2b newline ...
                Line3 newline ...
                Line4 newline ...
                LineTitle];

            for i = 1:length(LinesText)
                Energy_Table{end} = [Energy_Table{end} newline ...
                    LinesText{i}];
            end

            Energy_Table{end} = [Energy_Table{end} newline ...
                Line5 newline ...
                Line6 newline ...
                Line7 newline ...
                Line8];

        end

        x = 0; % keep track of theory number

        % Start new table
        Line1 = ['%Table: ' Current_Salt ' Energy Results'];
        Line2 = '\begin{table}[H]';
        Line2b = '    \centering';
        Line3 = '    \begin{adjustwidth}{-0.10in}{-0.10in}% adjust the L and R margins by 0.1 inch';
        Line4 = ['        \begin{tabular}' TableAdj];
        LineTitle = [pad('',46) '& \textbf{Fav.} '];
        Line5 = '        \end{tabular}';
        Line6 = '    \end{adjustwidth}';
        Line7 = ['    \caption{Lattice energies of optimized ' Current_Salt ' crystals with different calculation and crystal structure types. The lowest lattice energy for each calculation is emboldened.}'];
        Line8 = '\end{table}';

        Theory_List = {};
        LinesText = {};
    end

    if ~contained_in_cell(Current_Theory,Theory_List) % New theory

        if strcmp(Current_Salt,Previous_Salt) % Same salt
            [~,MinIdx] = min(LineEnergies);
            LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];

            Favoured = Structure_List{MinIdx};
            Favoured = strrep(Favoured,'Rocksalt','Rock.');
            Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
            Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
            Favoured = strrep(Favoured,'FiveFive','5-5');
            Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
            LinesText{x} = [LinesText{x} '& ' Favoured];

            for i = 1:length(LinesData)
                LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
                LinesText{x} = [LinesText{x} '&' LinesData{i}];
            end
        end

        x = x + 1;
        Theory_List{end+1} = Current_Theory;
        LinesText{x} = pad(['\textbf{' Current_Theory BS '}'],34,'right');
        LinesData = cell(1,N_SL);
        LineEnergies = nan(1,N_SL);
    end

    LinesData{Str_Index} = Current_E;
    LineEnergies(Str_Index) = Current_E_num;
    Previous_Salt = Current_Salt;
end
    
% Build last energy table
% Finish up last row
[~,MinIdx] = min(LineEnergies);
LinesData{MinIdx} = ['\textbf{' LinesData{MinIdx} '}'];

Favoured = Structure_List{MinIdx};
Favoured = strrep(Favoured,'Rocksalt','Rock.');
Favoured = strrep(Favoured,'Wurtzite','Wurtz.');
Favoured = strrep(Favoured,'BetaBeO','$\beta$-BeO');
Favoured = strrep(Favoured,'FiveFive','5-5');
Favoured = pad(strrep(Favoured,'Sphalerite','Sphal.'),14,'right');
LinesText{x} = [LinesText{x} '& ' Favoured];

for i = 1:length(LinesData)
    LinesData{i} = pad([' ' LinesData{i} Note],22,'right');
    LinesText{x} = [LinesText{x} '&' LinesData{i}];
end

% Finish title
for i = 1:N_SL
    LineTitle = [LineTitle pad(['& \textbf{' Shorten_Structure_Name(Structure_List{i}) '} '],23,'right')];
end

% Finish up all rows
LengthT = length(LineTitle);
LineTitle = [LineTitle '\\ \hline'];
for i = 1:length(LinesText)
    if i < length(LinesText)
        LinesText{i} = ['            ' LinesText{i} '\\'];
    else
        LinesText{i} = ['            ' LinesText{i}];
    end
end


% Build table
Energy_Table{end+1} = [...
    Line1 newline ...
    Line2 newline ...
    Line2b newline ...
    Line3 newline ...
    Line4 newline ...
    LineTitle];

for i = 1:length(LinesText)
    Energy_Table{end} = [Energy_Table{end} newline ...
        LinesText{i}];
end

Energy_Table{end} = [Energy_Table{end} newline ...
    Line5 newline ...
    Line6 newline ...
    Line7 newline ...
    Line8];

   
   

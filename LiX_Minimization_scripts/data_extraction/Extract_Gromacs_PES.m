% Script for extracting informtion from gromacs output files.

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'LiF_R', 'NaCl_S', etc.

function Data = Extract_Gromacs_PES(Salt,Structure,Model,Directory)
% Conversion Factor
Ang_per_nm = 10;

%% Level 1: Salts

% Find directories in current folder, except . and ..
Folders = find_directories(dir(Directory));

% Find folders that match the input salt
matches = strcmp(Folders,Salt);

% No matches? return
if sum(matches) == 0
    Data = [];
    return
end

%% Level 2: Structures
% Enter directory i
Current_Dir = fullfile(Directory,Salt);

% Find directories in current folder, except . and ..
Structures = find_directories(dir(Current_Dir));

% Find folders that match the input Structure
matches = strcmp(Structures,Structure);

% No matches? return
if sum(matches) == 0
    Data = [];
    return
end

% Get the label (comes in '_X' form)
Current_Label = Label_replace(Structure);

%% Level 3: Models

% Update to new directory for model
Current_Dir = fullfile(Directory,Salt,Structure,Model);

% Find directories in current folder, except . and ..
LatticeParams = find_directories(dir(Current_Dir));

% No data found? Next model
if isempty(LatticeParams)
    Data = [];
    return
end

O = length(LatticeParams);
rmd = 0;
Data = [];
%% Level 5: Lattice Parameters
for k=1:O

    Current_LatParTag = LatticeParams{k};

    % Update to new directory with current lattice parameter
    Current_Dir = fullfile(Directory,Salt,Structure,...
        Model,Current_LatParTag);

    Current_LatPar = strrep(Current_LatParTag,'LPAR_','');

    if isempty(Current_LatPar)
        rmd = rmd+1;
        continue
    elseif strcmp(Current_LatPar,'CELLOPT') 
        Data_type = 1;
    elseif strcmp(Current_LatPar,'FULLOPT')
        Data_type = 2;
    elseif strcmp(Current_LatPar,'ATOMOPT')
        Data_type = 3;
    elseif strcmp(Current_LatPar,'FULLOPT_SG1')
        Data_type = 4;
    else
        % 0 = normal run type
        Data_type = 0;
    end

    % Find xvg file
    enfile_info = dir([Current_Dir filesep Salt Current_Label '_' Model ...
        '_' Current_LatPar 'energies.xvg']);

    % Find total number of atoms .mat file
    totatoms_info = dir([Current_Dir filesep Salt Current_Label '_' Model ...
        '_' Current_LatPar '.mat']);

    % Find unit cell file
    unitcell_info = dir([Current_Dir filesep Salt Current_Label '_' Model ...
        '_' Current_LatPar '_UnitCell.g96']);

    % If no .xvg, supercell, or unit cell .gro file exists, remove from list and continue
    if isempty(enfile_info) || isempty(unitcell_info) || isempty(totatoms_info)
        rmd = rmd+1;
        continue
    end

    % Get total atoms filename
    Atoms_file = fullfile(Current_Dir,totatoms_info.name);

    % Open the mat file to get total number of atoms/molecules
    try
        N_info = load(Atoms_file,'-mat','N_Cell','N_total');
    catch
        rmd = rmd+1;
        continue
    end
    N_atoms = N_info.N_Cell;
    N_total_mols = N_info.N_total/2;

    % Get unit cell filename
    UnitCell_file = fullfile(Current_Dir,unitcell_info.name);

    % Open the unit cell structure file
    UnitCellText = fileread(UnitCell_file);
    UnitCellText = regexprep(UnitCellText, {'\r', '\n\n+'}, {'', '\n'});
    
    if isempty(UnitCellText)
        rmd = rmd+1;
        continue
    end
    
    %% Get unit cell parameters
    boxtemp = regexp(UnitCellText,'(?<=BOX\n)(.+?)(?=(\n)END)','match');
    boxcoords = textscan(boxtemp{1},'%f %f %f %f %f %f %f %f %f',...
        'Delimiter',' ','MultipleDelimsAsOne',true);

    % Lattice vectors
    a_vec = [boxcoords{1} boxcoords{4} boxcoords{5}];
    b_vec = [boxcoords{6} boxcoords{2} boxcoords{7}];
    c_vec = [boxcoords{8} boxcoords{9} boxcoords{3}];

    % Volume of unit cell in Angstrom^3
    Volume = (dot(cross(a_vec,b_vec),c_vec))*Ang_per_nm^3;

    % Box vector lengths in Angstroms
    a = norm(a_vec)*Ang_per_nm;
    b = norm(b_vec)*Ang_per_nm;
    c = norm(c_vec)*Ang_per_nm;

    % Transformation Matrix
    TM = [a_vec ; b_vec ; c_vec];
    
    %% Get position Data
    Positions = regexp(UnitCellText,'(?<=POSITION\n)(.+?)(?=(\n)END)','match','ONCE');

    % g96 file format
    Coords_Data = textscan(Positions,'%*6c%6c%*6c%*6c%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    
    % Get xyz coordinates of the asymmetric unit
    idx = (N_atoms)/2 + 1;
    Met_XYZ = [Coords_Data{2}(1:idx-1) Coords_Data{3}(1:idx-1) Coords_Data{4}(1:idx-1)];
    Hal_XYZ = [Coords_Data{2}(idx:end) Coords_Data{3}(idx:end) Coords_Data{4}(idx:end)];

    % Convert to Fractional coordinates
    Met_Frac = Met_XYZ/TM;
    Hal_Frac = Hal_XYZ/TM;

    % Save Fractional Coordinates of asymmetric unit in specific format
    if Data_type == 4
        FractionalCoords = cell(N_atoms,4);
        for i=1:(N_atoms)/2
            FractionalCoords{i,1} = strtrim(Coords_Data{1}(i,:));
            for j=2:4
                FractionalCoords{i,j} =  Met_Frac(i,j-1);
            end
        end
        for i=idx:N_atoms
            FractionalCoords{i,1} = strtrim(Coords_Data{1}(i,:));
            for j=2:4
                FractionalCoords{i,j} =  Hal_Frac(i-idx+1,j-1);
            end
        end
        
    else % Only save asymmetric unit
        FractionalCoords = { strtrim(Coords_Data{1}(1,:)) Met_Frac(1,1) Met_Frac(1,2) Met_Frac(1,3) ; ...
            strtrim(Coords_Data{1}(idx,:)) Hal_Frac(1,1) Hal_Frac(1,2) Hal_Frac(1,3) };
    end
    
    % Get cell angles
    alpha = VecAngle(b_vec,c_vec);
    beta = VecAngle(a_vec,c_vec);
    gamma = VecAngle(a_vec,b_vec);

    Cell_angles = [alpha beta gamma];
    
    %% Find minimum Met-Hal bond length
    MH_Bonds = nan(idx-1,idx-1);
    for p = 1:idx-1
        for q = 1:idx-1
            MH_Bonds(p,q) = norm(Met_XYZ(p,:) - Hal_XYZ(q,:));
        end
    end
    % Bond length in angstrom
    MH_BL = min(MH_Bonds(:))*Ang_per_nm;
    
    %% Find mimimum Met-Met bond length
    MM_Bonds = nan(idx-1,idx-1);
    for p = 1:idx-1
        for q = p+1:idx-1
            MM_Bonds(p,q) = norm(Met_XYZ(p,:) - Met_XYZ(q,:));
        end
    end
    % Bond length in angstrom
    if isnan(MM_Bonds)
        % If unit cell only contains 1 metal or halide, then the nearest
        % neighbour is +1 unit cell over;
        MM_BL = min([a b c]);
    else
        MM_BL = min(MM_Bonds(:))*Ang_per_nm;
    end
    
    %% Find minimum Hal-Hal bond length
    HH_Bonds = nan(idx-1,idx-1);
    for p = 1:idx-1
        for q = p+1:idx-1
            HH_Bonds(p,q) = norm(Hal_XYZ(p,:) - Hal_XYZ(q,:));
        end
    end   
    % Bond length in angstrom
    if isnan(HH_Bonds)
        HH_BL = min([a b c]);
    else
        HH_BL = min(HH_Bonds(:))*Ang_per_nm;
    end

    %% Get energy Data
    % filename
    Energy_file = fullfile(Current_Dir,enfile_info.name);

    % Import energy file as text
    Energy_text = fileread(Energy_file);

    % if no text found, end loop
    if isempty(Energy_text)
        rmd = rmd+1;
        continue
    end

    % Find the energy line
    if Data_type == 0
        Energies_Initial = regexp(Energy_text,'    0.000000.+?\n','match');
        Energy_list = textscan(Energies_Initial{1},'%*f %f %f %f %f %*f %*f %f %f %f %f %f %f %*f',...
            'Delimiter',' ','MultipleDelimsAsOne',true);
        Energy_list(cellfun('isempty',Energy_list)) = {0};
    else
        Energies_Final = textscan(Energy_text, '%s', 'delimiter', '\n','MultipleDelimsAsOne',true);

        Energy_list = textscan(Energies_Final{1}{end},'%*f %f %f %f %f %f %f %f %f %f %f',...
            'Delimiter',' ','MultipleDelimsAsOne',true);
        Energy_list(cellfun('isempty',Energy_list)) = {0};
    end
    Energy_array = cell2mat(Energy_list)./N_total_mols;

    LJ_SR = Energy_array(1);
    Coul_SR = Energy_array(2);
    Coul_Recip = Energy_array(3);
    Total_Potential = Energy_array(4);
    Coul_SR_MetMet = Energy_array(5);
    LJ_SR_MetMet = Energy_array(6);
    Coul_SR_Salt = Energy_array(7);
    LJ_SR_Salt = Energy_array(8);
    Coul_SR_HalHal = Energy_array(9);
    LJ_SR_HalHal = Energy_array(10);

    LJ_Recip = 0;

    % Place data in object for output
    if k-rmd == 1
        Data = cell(O,21);
    end
    Data{k-rmd,1} = a; % Angstroms
    Data{k-rmd,2} = b; % Angstroms
    Data{k-rmd,3} = c; % Angstroms
    Data{k-rmd,4} = MH_BL; % (Metal-Halide Bond Length) Angstroms
    Data{k-rmd,5} = Volume; % Angstroms^3
    Data{k-rmd,6} = FractionalCoords; % Frac coords of the asymmetric unit
    Data{k-rmd,7} = Total_Potential; % Potential Energy in kJ/mol of salt pair
    Data{k-rmd,8} = Cell_angles; % alpha beta and gamma
    Data{k-rmd,9} = Data_type; % Data type: 0 is normal run, 1 2 3 are optimization runs
    Data{k-rmd,10} = Coul_SR; % Coulomb energy of direct lattice (short range)
    Data{k-rmd,11} = LJ_SR; % LJ energy of direct lattice (short range)
    Data{k-rmd,12} = Coul_Recip; %  Coulomb energy of reciprocal lattice
    Data{k-rmd,13} = LJ_Recip; % LJ energy of reciprocal lattice
    Data{k-rmd,14} = Coul_SR_Salt; % Coulomb energy of only Metal-Halide interaction (short range)
    Data{k-rmd,15} = LJ_SR_Salt; % LJ energy of only Metal-Halide interaction (short range)
    Data{k-rmd,16} = Coul_SR_MetMet; % Coulomb energy of only Metal-Metal interaction (short range)
    Data{k-rmd,17} = LJ_SR_MetMet; % LJ energy of only Metal-Metal interaction (short range)
    Data{k-rmd,18} = Coul_SR_HalHal; % Coulomb energy of only Halide-Halide interaction (short range)
    Data{k-rmd,19} = LJ_SR_HalHal; % LJ energy of only Halide-Halide interaction (short range)
    Data{k-rmd,20} = MM_BL; % (Metal-Metal Bond Length) Angstroms
    Data{k-rmd,21} = HH_BL; % (Halide-Halide Bond Length) Angstroms
end

if isempty(Data)
    return
end

% Remove any empty rows
while true
    if isempty(Data{end,1})
        Data(end,:) =[];
    else
        break
    end
end

% Re-Sort by lattice parameter a
[~,Idx] = sort([Data{:,1}]);
for i = 1:21
    Data(:,i) = Data(Idx,i);
end

end
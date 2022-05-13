% SCRIPT FOR EXTRACTING ENERGIES FROM CRYSTAL17 OUTPUT (*.out) FILES.

% WARNING: This script does not accept theory types with acronyms longer than
% 10 characters

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'LiF_R', 'NaCl_S', etc.

function Data = Extract_CRYSTAL_PES(Label,directory,Structure,BSSE)
%% DFT method types
DFT_methods = {'SVWN' 'BLYP' 'PBEXC' 'PBESOLXC' 'SOGGAXC' 'B97' 'B3PW' ...
    'B3LYP' 'PBE0' 'PBESOL0' 'WC1LYP' 'B97H' 'PBE0-13' 'PW1PW' 'HSE06' ...
    'HSEsol' 'SC-BLYP' 'HISS' 'RSHXLDA' 'LC-wPBE' 'LC-wPBEsol' 'LC-wBLYP' ...
    'wB97' 'wB97X' 'LC-BLYP' 'CAM-B3LYP' 'M06L' 'M05' 'M052X' 'M06' ...
    'M062X' 'M06HF' 'B2PLYP' 'B2GPPLYP' 'mPW2PLYP'};
Range_cor = {'HSE06' 'HSEsol' 'SC-BLYP' 'HISS' 'RSHXLDA' 'LC-wPBE' ...
    'LC-wPBEsol' 'LC-wBLYP' 'wB97' 'wB97X' 'LC-BLYP' 'CAM-B3LYP'};

%% Level 1 (outer)

% Find directories in current folder, except . and ..
Folders = find_directories(dir(directory));

% Find folders that match the label
matches = regexp(Folders,['^' Label '_']);

% Ensure structure names are good by removing dashes
structnames = regexprep(Folders,'-','_');

Data = struct();

%% Level 2
K = length(Folders);
for i=1:K
    
    % Skip folders that do not match label
    if isempty(matches{i})
        continue
    end
    
    % Enter directory i
    cur_dir = [directory filesep Folders{i}];
    
    % Find directories in current folder, except . and ..
    subFolders = find_directories(dir(cur_dir));
    
    N = length(subFolders);

    % Pre-Create data output for this theory type
    Data.(structnames{i}) = cell(N,6);
    
    % Find type of theory
    Theory_type = regexp(Folders{i},['(' Label '_)' '(.+)'],'tokens');
    Theory_type = Theory_type{1}{2};
    
    %% Level 3
    % Loop through directories
    for j=1:N

        Basis_set = textscan(subFolders{j},'BASS_%s');
        Basis_set = Basis_set{1}{1};
        
        Data.(structnames{i}){j,1} = strrep(Basis_set,'%2A','*');
        
        % Move to subfolder j
        cur_dir2 = [cur_dir filesep subFolders{j}];

        % Find directories in current folder, except . and ..
        subsubFolders = find_directories(dir(cur_dir2));
        
        M = length(subsubFolders);

        %% Level 4
        for k=1:M

            % Move to subsubfolder k
            cur_dir3 = [cur_dir2 filesep subsubFolders{k}];
            
            cell_param = regexp(subsubFolders{k},...
                '(((?!APAR_)[0-9\.]+)|CELLOPT|FULLOPT(.*)|ATOMOPT|VIBANALYSIS(_SS[2-9]){0,1}|QHA(_SS[2-9]){0,1})',...
                'match','ONCE');
            
            if isempty(cell_param)
                continue
            elseif strcmp(cell_param,'CELLOPT') 
                cell_param = 'CellOpt';
                data_type = 1;
            elseif regexp(cell_param,'FULLOPT.*')
                % Pressure or space group P1
                SG1 = regexp(cell_param,'_SG1','match','ONCE');
                Ext_P = regexp(cell_param,'_P([0-9]|\.)+','tokens','ONCE');
                if ~isempty(SG1)
                    cell_param = 'CellOpt';
                    data_type = 4;
                elseif ~isempty(Ext_P)
                    cell_param = 'CellOpt';
                    Ext_P = str2double(Ext_P{1});
                    data_type = 5;
                else
                    cell_param = 'CellOpt';
                    data_type = 2;
                end
            elseif strcmp(cell_param,'ATOMOPT')
                cell_param = 'CellOpt';
                data_type = 3;
            elseif ~isempty(regexp(cell_param,'VIBANALYSIS','ONCE'))
                % How large is supercell
                supercell = regexp(cell_param,'(_SS)([2-9])','tokens');
                cell_param = 'Vib';
                data_type = 2;
            elseif ~isempty(regexp(cell_param,'QHA','ONCE'))
                % How large is supercell
                supercell = regexp(cell_param,'(_SS)([2-9])','tokens');
                cell_param = 'QHA';
                data_type = 2;
            else
                data_type = 0;
            end
            
            % If using D4 or D4TB, collect energy data differently
            if ~isempty(regexp(Theory_type,'_D4','once'))
                D4_Run = true;
                
                % If no OUT.mat file exists, remove from list and continue
                matfile_info = dir([cur_dir3 filesep Label '_' Theory_type ...
                    '_' Basis_set '_' cell_param '_OUT.mat']);
                
                if isempty(matfile_info)
                    continue
                % If multiple log files exist, choose most recent
                elseif size(matfile_info,1) > 1
                    Times = zeros(1,size(matfile_info,1));
                    for y = 1:size(matfile_info,1)
                        Times(y) = datenum(matfile_info(y).date,...
                            'dd-mmm-yyyy HH:MM:SS');
                    end
                    [~,Idx] = max(Times);
                    matfile_info = matfile_info(Idx);
                end
                
                dta = load(fullfile(matfile_info.folder,matfile_info.name));
                if isempty(dta)
                    continue
                end
                
                DOI = dta.Data;
            else
                D4_Run = false;
            end
                
            % If no log file exists, remove from list and continue
            logfile_info = dir([cur_dir3 filesep Label '_' Theory_type ...
                '_' Basis_set '*' cell_param '.out']);

            if isempty(logfile_info)
                continue
            % If multiple log files exist, choose most recent
            elseif size(logfile_info,1) > 1
                Times = zeros(1,size(logfile_info,1));
                for y = 1:size(logfile_info,1)
                    Times(y) = datenum(logfile_info(y).date,...
                        'dd-mmm-yyyy HH:MM:SS');
                end
                [~,Idx] = max(Times);
                logfile_info = logfile_info(Idx);
            end

            % Otherwise, get filename
            logfile = [cur_dir3 filesep logfile_info.name];

            % Open the file
            fid = fopen(logfile,'rt');
            logtext = fread(fid);
            fclose(fid);
            logtext = char(logtext.');

            if isempty(logtext)
                continue
            end
            Toldee = nan;
            
            %% QHA Analysis Runs
            if strcmp(cell_param,'QHA')
                if isempty(supercell)
                    SS_Size = 1;
                else
                    SS_Size = str2double(supercell{1}{2});
                end
                
                % Structural Info
                % Get the lattice parameters
                cell_params = regexp(logtext,...
                    '(PRIMITIVE CELL - CENTRING.+?\n)(.+?\n)(.+?\n)','tokens','ONCE');

                % if the unit cell is already primitive
                if isempty(cell_params)
                	continue
                end

                % output volume and cell parameters
                Vol = regexp(cell_params{1},'(VOLUME= +)(-|\+|\.|[0-9]|E)+','tokens','ONCE');
                if isempty(Vol)
                    continue
                end
                volume = str2double(Vol{2});

                params_num = textscan(cell_params{3},'%f %f %f %f %f %f');
                a = params_num{1};
                b = params_num{2};
                c = params_num{3};
                alpha = params_num{4};
                beta = params_num{5};
                gamma = params_num{6};

                % Get the bond length
                bond_length = regexp(logtext,...
                    '(N = NUMBER OF NEIGHBORS.+?\n)(.+?\n)(.+?\n\n)','tokens','ONCE');
                
                if isempty(bond_length)
                    continue
                end
                BL_Scan = regexprep(bond_length{3},' {10,}.+?\n','');
                
                BL_cells = textscan(BL_Scan,' %*s %*s %*s %6.4f %*[^\n]','MultipleDelimsAsOne',1);
                BL = BL_cells{1}(1);

                % Get the fractional coordinates of atoms in asymmetric unit
                AU_Coords = {Get_Fractional_Coords(logtext)};
                
                % Isolate Energy and QHA Thermo Info from starting text
                endtext = regexp(logtext,'\+\+\+\+\+\+\+ FITTING USING ALL POINTS.+','match','ONCE');
                
                % No QHA analysis found?
                if isempty(endtext) 
                    continue
                end
                
                % Get Equations of state
                MURNAGHAN_1944_Txt = regexp(endtext,'MURNAGHAN 1944 +.+?\n','match','ONCE');
                if isempty(MURNAGHAN_1944_Txt)
                    continue
                end
                MURNAGHAN_1944_Cell = textscan(MURNAGHAN_1944_Txt,'%*s %*f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                MURNAGHAN_1944 = [MURNAGHAN_1944_Cell{:}];
                
                BIRCH_MURNAGHAN_1947_Txt = regexp(endtext,'BIRCH-MURNAGHAN 1947 +.+?\n','match','ONCE');
                if isempty(BIRCH_MURNAGHAN_1947_Txt)
                    continue
                end
                BIRCH_MURNAGHAN_1947_Cell = textscan(BIRCH_MURNAGHAN_1947_Txt,'%*s %*f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                BIRCH_MURNAGHAN_1947 = [BIRCH_MURNAGHAN_1947_Cell{:}];
                
                POIRIER_TARANTOLA_1998_Txt = regexp(endtext,'POIRIER-TARANTOLA 1998 +.+?\n','match','ONCE');
                if isempty(POIRIER_TARANTOLA_1998_Txt)
                    continue
                end
                POIRIER_TARANTOLA_1998_Cell = textscan(POIRIER_TARANTOLA_1998_Txt,'%*s %*f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                POIRIER_TARANTOLA_1998 = [POIRIER_TARANTOLA_1998_Cell{:}];
                
                VINET_1987_Txt = regexp(endtext,'VINET 1987 +.+?\n','match','ONCE');
                if isempty(VINET_1987_Txt)
                    continue
                end
                VINET_1987_Cell = textscan(VINET_1987_Txt,'%*s %*f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                VINET_1987 = [VINET_1987_Cell{:}];
                
                THIRD_ORDER_POLY_Txt = regexp(endtext,'THIRD ORDER POLYNOMIAL +.+?\n','match','ONCE');
                if isempty(THIRD_ORDER_POLY_Txt)
                    continue
                end
                THIRD_ORDER_POLY_Cell = textscan(THIRD_ORDER_POLY_Txt,'%*s %*s %*s %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                THIRD_ORDER_POLY = [THIRD_ORDER_POLY_Cell{:}];
                
                FOURTH_ORDER_POLY_Txt = regexp(endtext,'FOURTH ORDER POLYNOMIAL +.+?\n','match','ONCE');
                if isempty(FOURTH_ORDER_POLY_Txt)
                    continue
                end
                FOURTH_ORDER_POLY_Cell = textscan(FOURTH_ORDER_POLY_Txt,'%*s %*s %*s %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                FOURTH_ORDER_POLY = [FOURTH_ORDER_POLY_Cell{:}];
                
                FIFTH_ORDER_POLY_Txt = regexp(endtext,'FIFTH ORDER POLYNOMIAL +.+?\n','match','ONCE');
                if isempty(FIFTH_ORDER_POLY_Txt)
                    continue
                end
                FIFTH_ORDER_POLY_Cell = textscan(FIFTH_ORDER_POLY_Txt,'%*s %*s %*s %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                FIFTH_ORDER_POLY = [FIFTH_ORDER_POLY_Cell{:}];
                
                % Get Thermodynamic EOS Functions
                MURNAGHAN_1944_Txt = regexp(endtext,'MURNAGHAN 1944 +\n \++.+?\+','match','ONCE');
                if isempty(MURNAGHAN_1944_Txt)
                    continue
                end
                MURNAGHAN_1944_Cell = textscan(MURNAGHAN_1944_Txt,'%f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',7);
                MURNAGHAN_1944_Data = [MURNAGHAN_1944_Cell{:}];
                
                BIRCH_MURNAGHAN_1947_Txt = regexp(endtext,'BIRCH-MURNAGHAN 1947 +\n \++.+?\+','match','ONCE');
                if isempty(BIRCH_MURNAGHAN_1947_Txt)
                    continue
                end
                BIRCH_MURNAGHAN_1947_Cell = textscan(BIRCH_MURNAGHAN_1947_Txt,'%f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',7);
                BIRCH_MURNAGHAN_1947_Data = [BIRCH_MURNAGHAN_1947_Cell{:}];
                
                POIRIER_TARANTOLA_1998_Txt = regexp(endtext,'POIRIER-TARANTOLA 1998 +\n \++.+?\+','match','ONCE');
                if isempty(POIRIER_TARANTOLA_1998_Txt)
                    continue
                end
                POIRIER_TARANTOLA_1998_Cell = textscan(POIRIER_TARANTOLA_1998_Txt,'%f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',7);
                POIRIER_TARANTOLA_1998_Data = [POIRIER_TARANTOLA_1998_Cell{:}];
                
                VINET_1987_Txt = regexp(endtext,'VINET 1987 +\n \++.+?\+','match','ONCE');
                if isempty(VINET_1987_Txt)
                    continue
                end
                VINET_1987_Txt = regexprep(VINET_1987_Txt,' \*{8} ',' NaN ');
                VINET_1987_Cell = textscan(VINET_1987_Txt,'%f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',7);
                VINET_1987_Data = [VINET_1987_Cell{:}];
                
                % THERMAL PROPERTIES FROM GRUNEISEN PARAMETERS
                Thermal_Gruneisen_Txt = regexp(endtext,'GRUNEISEN PARAMETERS\n \*+.+?\*{10}','match','ONCE');
                if isempty(Thermal_Gruneisen_Txt)
                    continue
                end
                Thermal_Gruneisen_Txt = regexp(Thermal_Gruneisen_Txt,'ORDER OF POLYNOMIALS: 2.+?\*{10}','match','ONCE');
                Thermal_Gruneisen_Txt = regexprep(Thermal_Gruneisen_Txt,' \*{6,8} ',' NaN ');
                Thermal_Gruneisen_Cell = textscan(Thermal_Gruneisen_Txt,'%f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',4);
                Thermal_Gruneisen = [Thermal_Gruneisen_Cell{:}];
                
                % THERMAL PROPERTIES FROM HELMHOLTZ FREE ENERGY
                Thermal_Helmholtz_Txt = regexp(endtext,'FROM HELMHOLTZ FREE ENERGY\n \*+.+?\*{10}','match','ONCE');
                if isempty(Thermal_Helmholtz_Txt)
                    continue
                end
                Thermal_Helmholtz_Txt = regexp(Thermal_Helmholtz_Txt,'ORDER OF POLYNOMIALS: 3.+?\*{10}','match','ONCE');
                Thermal_Helmholtz_Txt = regexprep(Thermal_Helmholtz_Txt,' \*{8} ',' NaN ');
                Thermal_Helmholtz_Cell = textscan(Thermal_Helmholtz_Txt,'%f %f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',4);
                Thermal_Helmholtz = [Thermal_Helmholtz_Cell{:}];
                
                % LINEAR THERMAL EXPANSION OF CONVENTIONAL LATTICE PARAMETERS
                Linear_Therm_Exp_Txt = regexp(endtext,'LINEAR THERMAL EXPANSION OF CONVENTIONAL LATTICE PARAMETERS\n \*+.+?\*{10}','match','ONCE');
                if isempty(Linear_Therm_Exp_Txt)
                    continue
                end
                Linear_Therm_Exp_Txt = regexprep(Linear_Therm_Exp_Txt,' \*{8} ',' NaN ');
                Linear_Therm_Exp_Cell = textscan(Linear_Therm_Exp_Txt,'%f %f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',6);
                Linear_Therm_Exp = [Linear_Therm_Exp_Cell{:}];
                
                % PRESSURE-VOLUME-TEMPERATURE RELATION
                PVT_Txt = regexp(endtext,'PRESSURE-VOLUME-TEMPERATURE RELATION\n \*+.+?\*{10}','match','ONCE');
                if isempty(PVT_Txt)
                    continue
                end
                PVT_Txt = regexprep(PVT_Txt,'\n WARNING - PRESSURE OUT OF RANGE - INCREASE VOLUME RANGE\n','');
                PVT_Txt = regexprep(PVT_Txt,' \*{8} ',' NaN ');
                PVT_Cell = textscan(PVT_Txt,'%f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',6);
                PVT = [PVT_Cell{:}];
                
                % THERMAL PROPERTIES FROM HELMHOLTZ FREE ENERGY AT VARIOUS PRESSURES
                Thermal_Helmholtz_P_Txt = regexp(endtext,'THERMAL PROPERTIES FROM HELMHOLTZ FREE ENERGY AT P = +(-|[0-9]|\.)+ GPa\n \*+.+?\*{10}','match');
                if isempty(Thermal_Helmholtz_P_Txt)
                    continue
                end
                Linear_Exp_Therm_Param_Txt = regexp(endtext,'CONVENTIONAL LATTICE PARAMETERS\n \*+.+?(\*{10}|T{10})','match');
                if isempty(Linear_Exp_Therm_Param_Txt)
                    continue
                end
                Linear_Exp_Therm_Param_Txt(1)=[]; % First one already accounted for
                PN = length(Thermal_Helmholtz_P_Txt);
                THP = cell(PN,2);
                for r = 1:PN
                    TH_Pressure = regexp(Thermal_Helmholtz_P_Txt{r},'(P = +)((-|[0-9]|\.)+)','tokens','ONCE');
                    THP{r,1} = str2double(TH_Pressure{2});
                    
                    Thermal_Helmholtz_P_Txt{r} = regexprep(Thermal_Helmholtz_P_Txt{r},' \*{9} ',' NaN ');
                    Thermal_Helmholtz_P_Cell = textscan(Thermal_Helmholtz_P_Txt{r},'%f %f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',6);
                    Thermal_Helmholtz_P = [Thermal_Helmholtz_P_Cell{:}];
                    
                    Linear_Exp_Therm_Param_Txt{r} = regexprep(Linear_Exp_Therm_Param_Txt{r},' \*{8} ',' NaN ');
                    Linear_Exp_Therm_Param_Cell = textscan(Linear_Exp_Therm_Param_Txt{r},'%*f %f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',6);
                    Thermal_Helmholtz_P(:,8:13) = [Linear_Exp_Therm_Param_Cell{:}];
                    THP{r,2} = Thermal_Helmholtz_P;
                end

                % Convert from kJ/mol of unit cells to kJ/mol of ion pairs
                if contained_in_cell(Structure,{'Rocksalt' 'Pair' 'Sphalerite' 'CsCl'})
                    % Do Nothing
                elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
                    % Unit cell contains 2 ion pair
                    volume = volume./2;
                    MURNAGHAN_1944(1:2) = MURNAGHAN_1944(1:2)./2;
                    BIRCH_MURNAGHAN_1947(1:2) = BIRCH_MURNAGHAN_1947(1:2)./2;
                    POIRIER_TARANTOLA_1998(1:2) = POIRIER_TARANTOLA_1998(1:2)./2;
                    VINET_1987(1:2) = VINET_1987(1:2)./2;
                    THIRD_ORDER_POLY(1:2) = THIRD_ORDER_POLY(1:2)./2;
                    FOURTH_ORDER_POLY(1:2) = FOURTH_ORDER_POLY(1:2)./2;
                    FIFTH_ORDER_POLY(1:2) = FIFTH_ORDER_POLY(1:2)./2;
                    MURNAGHAN_1944_Data(:,[1 3]) = MURNAGHAN_1944_Data(:,[1 3])./2;
                    BIRCH_MURNAGHAN_1947_Data(:,[1 3]) = BIRCH_MURNAGHAN_1947_Data(:,[1 3])./2;
                    POIRIER_TARANTOLA_1998_Data(:,[1 3]) = POIRIER_TARANTOLA_1998_Data(:,[1 3])./2;
                    VINET_1987_Data(:,[1 3]) = VINET_1987_Data(:,[1 3])./2;
                    Thermal_Gruneisen(:,4) = Thermal_Gruneisen(:,4)./2;
                    Thermal_Helmholtz(:,[2 5 7]) = Thermal_Helmholtz(:,[2 5 7])./2;
                    PVT(:,[2 4 5]) = PVT(:,[2 4 5])./2;
                    for r = 1:PN
                        THP{r,2}(:,[2 5 6 7]) = THP{r,2}(:,[2 5 6 7])./2;
                    end

                elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
                    % Unit cell contains 4 ion pair
                    volume = volume/4;
                    MURNAGHAN_1944(1:2) = MURNAGHAN_1944(1:2)/4;
                    BIRCH_MURNAGHAN_1947(1:2) = BIRCH_MURNAGHAN_1947(1:2)/4;
                    POIRIER_TARANTOLA_1998(1:2) = POIRIER_TARANTOLA_1998(1:2)/4;
                    VINET_1987(1:2) = VINET_1987(1:2)/4;
                    THIRD_ORDER_POLY(1:2) = THIRD_ORDER_POLY(1:2)/4;
                    FOURTH_ORDER_POLY(1:2) = FOURTH_ORDER_POLY(1:2)/4;
                    FIFTH_ORDER_POLY(1:2) = FIFTH_ORDER_POLY(1:2)/4;
                    MURNAGHAN_1944_Data(:,[1 3]) = MURNAGHAN_1944_Data(:,[1 3])/4;
                    BIRCH_MURNAGHAN_1947_Data(:,[1 3]) = BIRCH_MURNAGHAN_1947_Data(:,[1 3])/4;
                    POIRIER_TARANTOLA_1998_Data(:,[1 3]) = POIRIER_TARANTOLA_1998_Data(:,[1 3])/4;
                    VINET_1987_Data(:,[1 3]) = VINET_1987_Data(:,[1 3])/4;
                    Thermal_Helmholtz(:,[2 5 7]) = Thermal_Helmholtz(:,[2 5 7])/4;
                    Thermal_Gruneisen(:,4) = Thermal_Gruneisen(:,4)/4;
                    PVT(:,[2 4 5]) = PVT(:,[2 4 5])/4;
                    for r = 1:PN
                        THP{r,2}(:,[2 5 6 7]) = THP{r,2}(:,[2 5 6 7])/4;
                    end
                else
                    error(['Unknown Crystal Type: ' Structure])
                end
                
                % Set dispersion key to false, not needed
                dispersion = false;
                D3_gCP = false;
                
            elseif strcmp(cell_param,'Vib')
                if isempty(supercell)
                    SS_Size = 1;
                else
                    SS_Size = str2double(supercell{1}{2});
                end
                
                % Structural Info
                % Get the lattice parameters
                cell_params = regexp(logtext,...
                    '(CRYSTALLOGRAPHIC CELL \(VOLUME=.+?\n)(.+?\n)(.+?\n)','tokens');

                % if the unit cell is already primitive
                if isempty(cell_params)
                    cell_params = regexp(logtext,...
                        '(PRIMITIVE CELL - CENTRING.+?\n)(.+?\n)(.+?\n)','tokens');

                    if isempty(cell_params)
                        continue
                    end
                end

                % output volume and cell parameters
                Vol = regexp(cell_params{1}{1},'(VOLUME= +)(-|\+|\.|[0-9]|E)+','tokens');
                if isempty(Vol)
                    continue
                end
                volume = str2double(Vol{1}{2});

                params_num = textscan(cell_params{1}{3},'%f %f %f %f %f %f');
                a = params_num{1};
                b = params_num{2};
                c = params_num{3};
                alpha = params_num{4};
                beta = params_num{5};
                gamma = params_num{6};

                % Get the bond length
                bond_length = regexp(logtext,...
                    'N = NUMBER OF NEIGHBORS.+?\n.+?\n(.+?\n\n\n)','tokens','ONCE');

                if isempty(bond_length)
                    continue
                end
                BL_Scan = regexprep(bond_length{1},{' {10,}.+?\n' '\n\n'},{'' '\n'});

                BL_cells = textscan(BL_Scan,' %*s %s %f %6.4f %*f %*s %s %*[^\n]','MultipleDelimsAsOne',1);
                BL = BL_cells{3}(1); % nearest neighbour bond length
                
                % Energy and Thermo Info
                endtext = regexp(logtext,'HARMONIC VIBRATIONAL CONTRIBUTIONS.+','match');
                
                % No vib analysis found?
                if isempty(endtext) 
                    continue
                end
                
                endtext = endtext{1};
                
                % Find Total electronic energy in hartrees
                Total_E_txt = regexp(endtext,...
                    '(EL +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                if isempty(Total_E_txt)
                    continue
                end
                Total_E = sscanf(Total_E_txt{1}{2},'%f %*s %*s');
                
                % Find zero-point energy
                Zero_Point_txt = regexp(endtext,...
                    '(E0 +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                if isempty(Zero_Point_txt)
                    continue
                end
                Zero_Point_Energy = sscanf(Zero_Point_txt{1}{2},'%*s %*s %f');
                
                % Find Other thermodynamic properties
                Thermo_txt = regexp(endtext,...
                    'THERMODYNAMIC FUNCTIONS WITH VIBRATIONAL CONTRIBUTIONS.+?\*\*\*\*\*','match');
                if isempty(Thermo_txt)
                    continue
                end
                
                TN = length(Thermo_txt);
                Temp = nan(1,TN);
                Pres = nan(1,TN);
                Vib_E = nan(1,TN);
                PV = nan(1,TN);
                Entropy = nan(1,TN);
                Heat_Capacity = nan(1,TN);
                for idx = 1:TN
                    T_txt = regexp(Thermo_txt{idx},'(\(T = +)(-|\+|\.|[0-9]|E|)+( K)','tokens');
                    Temp(idx) = str2double(T_txt{1}{2});
                    
                    P_txt = regexp(Thermo_txt{idx},'(P = +)(-|\+|\.|[0-9]|E|)+( MPA)','tokens');
                    Pres(idx) = str2double(P_txt{1}{2});
                    
                    Vib_E_txt = regexp(Thermo_txt{idx},'(ET +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                    if ~isempty(Vib_E_txt{1}{2})
                        Vib_E(idx) = sscanf(Vib_E_txt{1}{2},'%*s %*s %f');
                    end
                    
                    PV_txt = regexp(Thermo_txt{idx},'(PV +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                    if ~isempty(PV_txt{1}{2})
                        PV(idx) = sscanf(PV_txt{1}{2},'%*s %*s %f');
                    end
                    
                    Entropy_txt = regexp(Thermo_txt{idx},'(ENTROPY +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                    if ~isempty(Entropy_txt{1}{2})
                        Entropy(idx) = sscanf(Entropy_txt{1}{2},'%*s %*s %f');
                    end
                    
                    Heat_Capacity_txt = regexp(Thermo_txt{idx},'(HEAT CAPACITY +: +)((-|\+|\.|[0-9]|E|)+ *){3}','tokens');
                    if ~isempty(Heat_Capacity_txt{1}{2})
                        Heat_Capacity(idx) = sscanf(Heat_Capacity_txt{1}{2},'%*s %*s %f');
                    end
                end
                
                % Convert from kJ/mol of unit cells to kJ/mol of ion pairs
                if contained_in_cell(Structure,{'Rocksalt' 'Pair' 'Sphalerite' 'CsCl'})
                    % Do Nothing
                elseif contained_in_cell(Structure,{'Wurtzite' 'NiAs'})
                    % Unit cell contains 2 ion pair
                    volume = volume./2;
                    Total_E = Total_E./2;
                    Zero_Point_Energy = Zero_Point_Energy/2;
                    Vib_E = Vib_E./2;
                    PV = PV./2;
                    Entropy = Entropy./2;
                    Heat_Capacity = Heat_Capacity./2;
                elseif contained_in_cell(Structure,{'BetaBeO' 'FiveFive'})
                    % Unit cell contains 4 ion pair
                    volume = volume./4;
                    Total_E = Total_E./4;
                    Zero_Point_Energy = Zero_Point_Energy./4;
                    Vib_E = Vib_E./4;
                    PV = PV./4;
                    Entropy = Entropy./4;
                    Heat_Capacity = Heat_Capacity./4;
                else
                    error(['Unknown Crystal Type: ' Structure])
                end
                
                % Collect phonon density of states data
                Total_PDOS = [];
                Metal_PDOS = [];
                Halide_PDOS = [];
                PDOS_text = regexp(logtext,'TOTAL PHONON DENSITY OF STATES.+','match');
                if ~isempty(PDOS_text)
                    % Total PDOS
                    PDOS_text = PDOS_text{1};
                    Total_PDOS_txt = regexp(PDOS_text,'(FREQUENCY \(CM\*\*-1\)    TOTAL PDOS\n\n)(.+)(\n\n INTEGRAL)','tokens');
                    Total_PDOS_cell = textscan(Total_PDOS_txt{1}{2},'%f %f','Delimiter',' ','MultipleDelimsAsOne',1);
                    Total_PDOS = [Total_PDOS_cell{1} Total_PDOS_cell{2}];
                
                    % PDOS of atomic species
                    PDOS_Atom_text = regexp(PDOS_text,'VALUES PER ATOMIC SPECIES.+','match');
                    PDOS_Atom_text = PDOS_Atom_text{1};
                    PDOS_Atom_text_split = strsplit(PDOS_Atom_text,'ATOMIC NUMBER');
                    
                    for p=2:3
                        current_Atomic_PDOS = PDOS_Atom_text_split{p};
                        atomic_n = regexp(current_Atomic_PDOS,' +([0-9]){1,3}\n','match','ONCE');
                        Z = str2double(strtrim(atomic_n));
                        
                        Atom_PDOS_cell = textscan(current_Atomic_PDOS,'%f %f','Delimiter',' ','HeaderLines',4,'MultipleDelimsAsOne',1);
                        Atom_PDOS = [Atom_PDOS_cell{1} Atom_PDOS_cell{2}];
                        
                        if sum(Z == [3 11 19 37 237]) % Metals
                            Metal_PDOS = Atom_PDOS;
                        elseif sum(Z ==  [9 17 35 235 53 253]) % Halides
                            Halide_PDOS = Atom_PDOS;
                        else
                            error(['Unknown molecular species, atomic number: ' strtrim(atomic_n)])
                        end
                    end
                end
                
                % Set dispersion key to false, not needed
                dispersion = false;
                D3_gCP = false;
                
            %% Optimization Runs
            elseif strcmp(cell_param,'CellOpt')
                
                if D4_Run
                    a = DOI.a;
                    b = DOI.b;
                    c = DOI.c;
                    alpha = DOI.alpha;
                    beta = DOI.beta;
                    gamma = DOI.gamma;
                    disp_E = DOI.E_D4;
                    gCP_E = 0;
                    energy_nocorr =  DOI.E_DFT;
                    energy_tot = DOI.E;
                    dispersion = true;
                    D3_gCP = false;
                    gCP = false;
                    Toldee = DOI.Toldee;
                    
                    % Get the lattice parameters
                    cell_params = regexp(logtext,...
                        '(CRYSTALLOGRAPHIC CELL \(VOLUME=.+?\n)(.+?\n)(.+?\n)','tokens');

                    % if the unit cell is already primitive
                    if isempty(cell_params)
                        cell_params = regexp(logtext,...
                            '(PRIMITIVE CELL - CENTRING.+?\n)(.+?\n)(.+?\n)','tokens');

                        if isempty(cell_params)
                            continue
                        end
                    end

                    % output volume and cell parameters
                    Vol = regexp(cell_params{1}{1},'(VOLUME= +)(-|\+|\.|[0-9]|E)+','tokens');
                    if isempty(Vol)
                        continue
                    end
                    volume = str2double(Vol{1}{2});
                    
                    % Get the bond length
                    bond_length = regexp(logtext,...
                        'N = NUMBER OF NEIGHBORS.+?\n.+?\n(.+?\n\n\n)','tokens','ONCE');

                    if isempty(bond_length)
                        continue
                    end

                    BL_Scan = regexprep(bond_length{1},{' {10,}.+?\n' '\n\n'},{'' '\n'});

                    BL_cells = textscan(BL_Scan,' %*s %s %f %6.4f %*f %*s %s %*[^\n]','MultipleDelimsAsOne',1);
                    BL = BL_cells{3}(1); % nearest neighbour bond length

                    % Separate metal and halide indexes
                    M_Ind = ismember(BL_cells{1},{'LI' 'NA' 'K' 'RB' 'CS'});
                    X_Ind = ~M_Ind;

                    % Pick out MX and MM info
                    M_Neighbours = BL_cells{2}(M_Ind);
                    M_BL = BL_cells{3}(M_Ind);
                    M_Targets = BL_cells{4}(M_Ind);

                    % Pick out XM and XX info
                    X_Neighbours = BL_cells{2}(X_Ind);
                    X_BL = BL_cells{3}(X_Ind);
                    X_Targets = BL_cells{4}(X_Ind);

                    % Get List of MX bond lengths and neighbour numbers
                    MX_Ind = ismember(M_Targets,{'F' 'CL' 'BR' 'I'});
                    MX_BL = {[M_BL(MX_Ind) M_Neighbours(MX_Ind)]};

                    % Get List of MM bond lengths and neighbour numbers
                    MM_Ind = ismember(M_Targets,{'LI' 'NA' 'K' 'RB' 'CS'});
                    MM_BL = {[M_BL(MM_Ind) M_Neighbours(MM_Ind)]};

                    % Get List of XX bond lengths and neighbour numbers
                    XX_Ind = ismember(X_Targets,{'F' 'CL' 'BR' 'I'});
                    XX_BL = {[X_BL(XX_Ind) X_Neighbours(XX_Ind)]};

                    % Get the fractional coordinates of atoms in asymmetric unit
                    AU_Coords = {Get_Fractional_Coords(logtext)};
                    
                else
                    % Grab the convergence criteria
                    Toldee_txt = regexp(logtext,'SCF TOL ON TOTAL ENERGY SET TO +([0-9]+)','tokens','ONCE');
                    if isempty(Toldee_txt)
                        Toldee = nan;
                    else
                        Toldee = str2double(Toldee_txt{1});
                    end

                    % Find text of final optimized parameters
                    endtext = regexp(logtext,'OPT END - CONVERGED.+','match');
                    Eformat = true; % format 1

                    % Check for other format
                    if isempty(endtext) 
                        endtext = regexp(logtext,'CELL/ATOM ITERATION CONVERGED:.+','match');
                        Eformat = false; % format 2

                        % If still no convergence found, do not use this data
                        if isempty(endtext) 
                            continue
                        end

                    end

                    logtext = endtext{1};

                    % Find energy outputs in hartree
                    if Eformat % Format 1
                        energy = regexp(logtext,...
                            '(OPT END - CONVERGED \* E\(AU\): +)(-|\+|\.|[0-9]|E|)+','tokens');
                        energy_num = str2double(energy{1}{2});
                    else % Format 2
                        energytxt = regexp(logtext,'\|E\(CELL\)-E\(ATOM\)\|.+?\n','match');
                        energyvals = textscan(energytxt{1},'%*s %*f %*s %f %*s %f');
                        energy_num = min([energyvals{:}]);
                    end

                    if strcmp(Structure,'Pair') % For pair potentials
                        cell_params = regexp(logtext,...
                            '(N = NUMBER OF NEIGHBORS AT DISTANCE R\n)(.+?\n)(.+?\n)','tokens');

                        if isempty(cell_params)
                            continue
                        end

                        BL_out = textscan(cell_params{1}{3},...
                            '%*s %*s %*s %f %*s %*s %*s %*s %*s %*s');
                        a = nan;
                        b = nan;
                        c = nan;
                        alpha = nan;
                        beta = nan;
                        gamma = nan;
                        volume = 0;
                        BL = BL_out{1};
                        MX_BL = {nan};
                        MM_BL = {nan};
                        XX_BL = {nan};
                        AU_Coords = {Get_Fractional_Coords(logtext)};
                    else % For salts
                        % Get the lattice parameters
                        cell_params = regexp(logtext,...
                            '(CRYSTALLOGRAPHIC CELL \(VOLUME=.+?\n)(.+?\n)(.+?\n)','tokens');

                        % if the unit cell is primitive
                        if isempty(cell_params)
                            cell_params = regexp(logtext,...
                                '(PRIMITIVE CELL - CENTRING.+?\n)(.+?\n)(.+?\n)','tokens');
                        end

                        % output volume and cell parameters
                        Vol = regexp(cell_params{1}{1},'(VOLUME= +)(-|\+|\.|[0-9]|E)+','tokens');
                        volume = str2double(Vol{1}{2});
                        params_num = textscan(cell_params{1}{3},'%f %f %f %f %f %f');
                        a = params_num{1};
                        b = params_num{2};
                        c = params_num{3};
                        alpha = params_num{4};
                        beta = params_num{5};
                        gamma = params_num{6};

                        % Get the bond length
                        bond_length = regexp(logtext,...
                            'N = NUMBER OF NEIGHBORS.+?\n.+?\n(.+?\n\n\n)','tokens','ONCE');

                        if isempty(bond_length)
                            continue
                        end
                        BL_Scan = regexprep(bond_length{1},{' {10,}.+?\n' '\n\n'},{'' '\n'});

                        BL_cells = textscan(BL_Scan,' %*s %s %f %6.4f %*f %*s %s %*[^\n]','MultipleDelimsAsOne',1);
                        BL = BL_cells{3}(1); % nearest neighbour bond length

                        % Separate metal and halide indexes
                        M_Ind = ismember(BL_cells{1},{'LI' 'NA' 'K' 'RB' 'CS'});
                        X_Ind = ~M_Ind;

                        % Pick out MX and MM info
                        M_Neighbours = BL_cells{2}(M_Ind);
                        M_BL = BL_cells{3}(M_Ind);
                        M_Targets = BL_cells{4}(M_Ind);

                        % Pick out XM and XX info
                        X_Neighbours = BL_cells{2}(X_Ind);
                        X_BL = BL_cells{3}(X_Ind);
                        X_Targets = BL_cells{4}(X_Ind);

                        % Get List of MX bond lengths and neighbour numbers
                        MX_Ind = ismember(M_Targets,{'F' 'CL' 'BR' 'I'});
                        MX_BL = {[M_BL(MX_Ind) M_Neighbours(MX_Ind)]};

                        % Get List of MM bond lengths and neighbour numbers
                        MM_Ind = ismember(M_Targets,{'LI' 'NA' 'K' 'RB' 'CS'});
                        MM_BL = {[M_BL(MM_Ind) M_Neighbours(MM_Ind)]};

                        % Get List of XX bond lengths and neighbour numbers
                        XX_Ind = ismember(X_Targets,{'F' 'CL' 'BR' 'I'});
                        XX_BL = {[X_BL(XX_Ind) X_Neighbours(XX_Ind)]};

                        % Get the fractional coordinates of asymmetric unit atoms
                        AU_Coords = {Get_Fractional_Coords(logtext)};
                    end

                    % Get dispersion energy if available
                    dispersiontxt = regexp(logtext,...
                        '(D3 DISPERSION ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+',...
                        'tokens');

                    % Get gCP energy if available
                    gCPtxt = regexp(logtext,...
                        '(GCP ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+',...
                        'tokens');

                    if ~isempty(gCPtxt) && ~isempty(dispersiontxt)
                        disp_E = str2double(dispersiontxt{1}{2});
                        gCP_E = str2double(gCPtxt{1}{2});

                        energy_nocorr = energy_num - disp_E - gCP_E;
                        energy_tot = energy_num;
                        dispersion = false;
                        D3_gCP = true;
                        gCP = false;
                    elseif ~isempty(dispersiontxt)
                        disp_E = str2double(dispersiontxt{1}{2});
                        gCP_E = 0;

                        energy_nocorr = energy_num - disp_E;
                        energy_tot = energy_num;
                        dispersion = true;
                        D3_gCP = false;
                        gCP = false;
                    elseif ~isempty(gCPtxt)
                        disp_E = 0;
                        gCP_E = str2double(gCPtxt{1}{2});

                        energy_nocorr = energy_num - gCP_E;
                        energy_tot = energy_num;
                        dispersion = false;
                        D3_gCP = false;
                        gCP = true;
                    else
                        disp_E = 0;
                        gCP_E = 0;
                        energy_tot = energy_num;
                        dispersion = false;
                        D3_gCP = false;
                        gCP = false;
                    end
                end
                

            %% Non-optimization Runs
            else
                %% Non-optimization Pair runs
                if strcmp(Structure,'Pair') % For pair potentials
                    cell_params = regexp(logtext,...
                        '(N = NUMBER OF NEIGHBORS AT DISTANCE R\n)(.+?\n)(.+?\n)','tokens');
                    
                    if isempty(cell_params)
                        continue
                    end
                    
                    BL_out = textscan(cell_params{1}{3},...
                        '%*s %*s %*s %f %*s %*s %*s %*s %*s %*s');
                    a = nan;
                    b = nan;
                    c = nan;
                    alpha = nan;
                    beta = nan;
                    gamma = nan;
                    volume = 0;
                    BL = BL_out{1};
                    AU_Coords = {Get_Fractional_Coords(logtext)};
                %% Non-optimization crystal runs
                else
                    %% With BSSE
                    if BSSE
                        % Get the lattice parameters
                        cell_params = regexp(logtext,...
                            '(A +B +C +ALPHA +BETA +GAMMA *)(.*?)\n(.+?\n)','tokens');
                        
                        if isempty(cell_params)
                            matfile = [cur_dir3 filesep Label '_' Theory_type ...
                            	'_' Basis_set '-' cell_param '.mat'];
                            if exist(matfile, 'file') == 2
                                ExtData = load(matfile);
                                
                                a = ExtData.a;
                                b = ExtData.b;
                                c = ExtData.c;
                                alpha = 90;
                                beta = 90;
                                gamma = 90;
                                volume = a*b*c;
                                ExtData.FC_Metal{1} = upper(ExtData.FC_Metal{1});
                                ExtData.FC_Halide{1} = upper(ExtData.FC_Halide{1});
                                AU_Coords = {ExtData.FC_Metal{:}; ExtData.FC_Halide{:}};
                            else
                                 % File does not exist.
                                 continue
                            end
                        else
                            
                        for m=1:length(cell_params)
                            if isempty(cell_params{m}{2})
                                params_num = textscan(cell_params{m}{3},'%f %f %f %f %f %f');
                            else
                                % Primitive cell volume
                                vol_prim = textscan(cell_params{m}{3},'%*f %*f %*f %*f %*f %*f %f');
                            end
                        end
                        
                            a = params_num{1};
                            b = params_num{2};
                            c = params_num{3};
                            alpha = params_num{4};
                            beta = params_num{5};
                            gamma = params_num{6};

                            % Unit cell volume
                            if contained_in_cell(Structure,{'CsCl' 'NiAs' 'BetaBeO' 'FiveFive' 'Wurtzite'}) % Primitive unit cells
                                volume = vol_prim{1};
                            elseif contained_in_cell(Structure,{'Rocksalt' 'Sphalerite'}) % Face-centered unit cells
                                volume = vol_prim{1}*4;
                            end

                            % Get the fractional coordinates of atoms in asymmetric unit
                            AU_Coords = Get_Fractional_Coords(logtext);

                            % In units of fractional coordinates
                            for m=1:size(AU_Coords,1)
                                AU_Coords{m,2} = (round(AU_Coords{m,2}*10^5))/((10^5)*a);
                                AU_Coords{m,3} = (round(AU_Coords{m,3}*10^5))/((10^5)*b);
                                AU_Coords{m,4} = (round(AU_Coords{m,4}*10^5))/((10^5)*c);
                            end
                            AU_Coords = {AU_Coords};
                        end
                        
                        % Get the bond length
                        bond_length = regexp(logtext,...
                            'N = NUMBER OF NEIGHBORS.+?\n.+?\n(.+?\n\n\n)','tokens','ONCE');

                        if isempty(bond_length)
                            continue
                        end
                        BL_Scan = regexprep(bond_length{1},{' {10,}.+?\n' '\n\n'},{'' '\n'});

                        BL_cells = textscan(BL_Scan,' %*s %s %f %6.4f %*f %*s %s %*[^\n]','MultipleDelimsAsOne',1);
                        BL = BL_cells{3}(1); % nearest neighbour bond length

                        % Separate metal and halide indexes
                        M_Ind = ismember(BL_cells{1},{'LI' 'NA' 'K' 'RB' 'CS'});
                        X_Ind = ~M_Ind;

                        % Pick out MX and MM info
                        M_Neighbours = BL_cells{2}(M_Ind);
                        M_BL = BL_cells{3}(M_Ind);
                        M_Targets = BL_cells{4}(M_Ind);

                        % Pick out XM and XX info
                        X_Neighbours = BL_cells{2}(X_Ind);
                        X_BL = BL_cells{3}(X_Ind);
                        X_Targets = BL_cells{4}(X_Ind);

                        % Get List of MX bond lengths and neighbour numbers
                        MX_Ind = ismember(M_Targets,{'F' 'CL' 'BR' 'I'});
                        MX_BL = {[M_BL(MX_Ind) M_Neighbours(MX_Ind)]};

                        % Get List of MM bond lengths and neighbour numbers
                        MM_Ind = ismember(M_Targets,{'LI' 'NA' 'K' 'RB' 'CS'});
                        MM_BL = {[M_BL(MM_Ind) M_Neighbours(MM_Ind)]};

                        % Get List of XX bond lengths and neighbour numbers
                        XX_Ind = ismember(X_Targets,{'F' 'CL' 'BR' 'I'});
                        XX_BL = {[X_BL(XX_Ind) X_Neighbours(XX_Ind)]};
                        
                    %% Without BSSE
                    else
                        % Get the lattice parameters
                        cell_params = regexp(logtext,...
                            '(CRYSTALLOGRAPHIC CELL \(VOLUME=.+?\n)(.+?\n)(.+?\n)','tokens');
                        
                        % if the unit cell is already primitive
                        if isempty(cell_params)
                            cell_params = regexp(logtext,...
                                '(PRIMITIVE CELL - CENTRING.+?\n)(.+?\n)(.+?\n)','tokens');
                            
                            if isempty(cell_params)
                                continue
                            end
                        end

                        % output volume and cell parameters
                        Vol = regexp(cell_params{1}{1},'(VOLUME= +)(-|\+|\.|[0-9]|E)+','tokens');
                        if isempty(Vol)
                            continue
                        end
                        volume = str2double(Vol{1}{2});
                        
                        params_num = textscan(cell_params{1}{3},'%f %f %f %f %f %f');
                        a = params_num{1};
                        b = params_num{2};
                        c = params_num{3};
                        alpha = params_num{4};
                        beta = params_num{5};
                        gamma = params_num{6};

                        % Get the bond length
                        bond_length = regexp(logtext,...
                            'N = NUMBER OF NEIGHBORS.+?\n.+?\n(.+?\n\n\n)','tokens','ONCE');

                        if isempty(bond_length)
                            continue
                        end
                        
                        BL_Scan = regexprep(bond_length{1},{' {10,}.+?\n' '\n\n'},{'' '\n'});

                        BL_cells = textscan(BL_Scan,' %*s %s %f %6.4f %*f %*s %s %*[^\n]','MultipleDelimsAsOne',1);
                        BL = BL_cells{3}(1); % nearest neighbour bond length

                        % Separate metal and halide indexes
                        M_Ind = ismember(BL_cells{1},{'LI' 'NA' 'K' 'RB' 'CS'});
                        X_Ind = ~M_Ind;

                        % Pick out MX and MM info
                        M_Neighbours = BL_cells{2}(M_Ind);
                        M_BL = BL_cells{3}(M_Ind);
                        M_Targets = BL_cells{4}(M_Ind);

                        % Pick out XM and XX info
                        X_Neighbours = BL_cells{2}(X_Ind);
                        X_BL = BL_cells{3}(X_Ind);
                        X_Targets = BL_cells{4}(X_Ind);

                        % Get List of MX bond lengths and neighbour numbers
                        MX_Ind = ismember(M_Targets,{'F' 'CL' 'BR' 'I'});
                        MX_BL = {[M_BL(MX_Ind) M_Neighbours(MX_Ind)]};

                        % Get List of MM bond lengths and neighbour numbers
                        MM_Ind = ismember(M_Targets,{'LI' 'NA' 'K' 'RB' 'CS'});
                        MM_BL = {[M_BL(MM_Ind) M_Neighbours(MM_Ind)]};

                        % Get List of XX bond lengths and neighbour numbers
                        XX_Ind = ismember(X_Targets,{'F' 'CL' 'BR' 'I'});
                        XX_BL = {[X_BL(XX_Ind) X_Neighbours(XX_Ind)]};

                        % Get the fractional coordinates of atoms in asymmetric unit
                        AU_Coords = {Get_Fractional_Coords(logtext)};
                    end
                end
                
                % Find energy outputs in hartree, abort if not found        
                energy = regexp(logtext,...
                    '(CONVERGENCE ON (ENERGY|TESTER) +E\(AU\) +)(-|\+|\.|[0-9]|E|)+','tokens');
                
                % Detect energy convergence glitch
                %glitchtext = regexp(logtext,'DETOT  0.00E\+00 tst','ONCE');
                
                if isempty(energy) %|| ~isempty(glitchtext)
                    continue
                end
                
                energy_num = str2double(energy{1}{2});
                
                % Get dispersion energy if available
                dispersiontxt = regexp(logtext,...
                    '(D3 DISPERSION ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+','tokens');
                
                % Get gCP energy if available
                gCPtxt = regexp(logtext,...
                    '(GCP ENERGY \(AU\) +)(-|\+|\.|[0-9]|E|)+',...
                    'tokens');
                
                if ~isempty(gCPtxt) && ~isempty(dispersiontxt)
                    disp_E = str2double(dispersiontxt{1}{2});
                    gCP_E = str2double(gCPtxt{1}{2});
                    
                    energy_nocorr = energy_num;
                    energy_tot = energy_num + disp_E + gCP_E;
                    dispersion = false;
                    gCP = false;
                    D3_gCP = true;
                elseif ~isempty(dispersiontxt)
                    disp_E = str2double(dispersiontxt{1}{2});
                    gCP_E = 0;
                    
                    energy_nocorr = energy_num;
                    energy_tot = energy_num + disp_E;
                    dispersion = true;
                    gCP = false;
                    D3_gCP = false;
                elseif ~isempty(gCPtxt)
                    disp_E = 0;
                    gCP_E = str2double(gCPtxt{1}{2});
                    
                    energy_nocorr = energy_num;
                    energy_tot = energy_num + gCP_E;
                    dispersion = false;
                    gCP = true;
                    D3_gCP = false;
                else
                    disp_E = 0;
                    gCP_E = 0;
                    
                    energy_tot = energy_num;
                    dispersion = false;
                    gCP = false;
                    D3_gCP = false;
                end
            end

            % Extract computation time
            CPU_time = regexp(logtext,...
                '(T+ E[N|R|D]+ +TELAPSE +([0-9]|E|\.)+ +TCPU +)([0-9]|E|\.)+','tokens');
            if isempty(CPU_time)
                Tot_Time = 0;
            else
                Tot_Time = str2double(CPU_time{1}{2});
            end
            
            % For D3_gCP calculations
            if D3_gCP
            
                % Type of theory (no dispersion/gCP)
                Theory_type_no_disp = regexprep(Theory_type,'_D3_gCP','');

                % Input energy (no dispersion/gCP)               
                nocorr_name = [Label '_' Theory_type_no_disp];

                % Add to structure
                if ~isfield(Data,nocorr_name)
                    Data.(nocorr_name) = cell(1,6);
                    Data.(nocorr_name){1,1} = Basis_set;
                    x = 1;
                else
                    % Find the basis set
                    x = NaN;
                    for q = 1:size(Data.(nocorr_name),1)
                        
                        if strcmp(Data.(nocorr_name){q,1},Basis_set)
                            x = q;
                            break
                        end
                    end
                    
                    % If basis set not yet present, add to end
                    if isnan(x)
                        x = q+1;
                        Data.(nocorr_name)(x,:) = cell(1,6);
                        Data.(nocorr_name){x,1} = Basis_set;
                    end
                    
                end                

                if isempty(Data.(nocorr_name){x,2})
                    Data.(nocorr_name){x,3} = Tot_Time;
                    Data.(nocorr_name){x,2} = cell(1,19);
                    Data.(nocorr_name){x,2}{1} = Theory_type_no_disp; % Theory type
                    Data.(nocorr_name){x,2}{2} = a; % A
                    Data.(nocorr_name){x,2}{3} = b;
                    Data.(nocorr_name){x,2}{4} = c;
                    Data.(nocorr_name){x,2}{5} = BL; % Bond length
                    Data.(nocorr_name){x,2}{6} = volume; % Volume
                    Data.(nocorr_name){x,2}{7} = AU_Coords; % Frac coords of asym unit
                    Data.(nocorr_name){x,2}{8} = energy_nocorr; % No-dispersion energy
                    Data.(nocorr_name){x,2}{9} = Tot_Time; % total time of calculation
                    Data.(nocorr_name){x,2}{10} = 0; % Data type
                    Data.(nocorr_name){x,2}{11} = alpha; % alpha angle (degrees)
                    Data.(nocorr_name){x,2}{12} = beta; % beta angle (degrees)
                    Data.(nocorr_name){x,2}{13} = gamma; % gamma angle (degrees)
                    Data.(nocorr_name){x,2}{14} = Toldee; % Convergence criteria
                    Data.(nocorr_name){x,2}{15} = 0; % Dispersion energy
                    Data.(nocorr_name){x,2}{16} = 0; % gCP energy
                    Data.(nocorr_name){x,2}{17} = MX_BL; % MX bond lengths
                    Data.(nocorr_name){x,2}{18} = MM_BL; % MM bond lengths
                    Data.(nocorr_name){x,2}{19} = XX_BL; % XX bond lengths
                else
                    if ~ismemberf(BL,Data.(nocorr_name){x,2}{1,5},'tol',1e-5)
                        Data.(nocorr_name){x,3} = Data.(nocorr_name){x,3} + Tot_Time;
                        Data.(nocorr_name){x,2}{2}(end+1) = a;
                        Data.(nocorr_name){x,2}{3}(end+1) = b;
                        Data.(nocorr_name){x,2}{4}(end+1) = c;
                        Data.(nocorr_name){x,2}{5}(end+1) = BL;
                        Data.(nocorr_name){x,2}{6}(end+1) = volume;
                        Data.(nocorr_name){x,2}{7}(end+1) = AU_Coords;
                        Data.(nocorr_name){x,2}{8}(end+1) = energy_nocorr;
                        Data.(nocorr_name){x,2}{9}(end+1) = Tot_Time;
                        Data.(nocorr_name){x,2}{10}(end+1) = 0;
                        Data.(nocorr_name){x,2}{11}(end+1) = alpha; % alpha angle (degrees)
                        Data.(nocorr_name){x,2}{12}(end+1) = beta; % beta angle (degrees)
                        Data.(nocorr_name){x,2}{13}(end+1) = gamma; % gamma angle (degrees)
                        Data.(nocorr_name){x,2}{14}(end+1) = Toldee; % Convergence criteria
                        Data.(nocorr_name){x,2}{15}(end+1) = 0; % Dispersion energy
                        Data.(nocorr_name){x,2}{16}(end+1) = 0; % gCP energy
                        Data.(nocorr_name){x,2}{17}(end+1) = MX_BL; % MX bond lengths
                        Data.(nocorr_name){x,2}{18}(end+1) = MM_BL; % MM bond lengths
                        Data.(nocorr_name){x,2}{19}(end+1) = XX_BL; % XX bond lengths

                        % Ensure data is sorted correctly by BL
                        [Data.(nocorr_name){x,2}{5},Indx] = sort(Data.(nocorr_name){x,2}{5},'ascend');
                        for y = [2:4 6:length(Data.(nocorr_name){x,2})]
                            Data.(nocorr_name){x,2}{y} = Data.(nocorr_name){x,2}{y}(Indx);
                        end
                    end
                end
            % For only gCP corrected calculations
            elseif gCP
                
                % Type of theory (no dispersion)
                Theory_type_no_gcp = regexprep(Theory_type,'_gCP','');

                % Input energy (no dispersion)               
                nocorr_name = [Label '_' Theory_type_no_gcp];

                % Add to structure
                if ~isfield(Data,nocorr_name)
                    Data.(nocorr_name) = cell(1,6);
                    Data.(nocorr_name){1,1} = Basis_set;
                    x = 1;
                else
                    % Find the basis set
                    x = NaN;
                    for q = 1:size(Data.(nocorr_name),1)
                        
                        if strcmp(Data.(nocorr_name){q,1},Basis_set)
                            x = q;
                            break
                        end
                    end
                    
                    % If basis set not yet present, add to end
                    if isnan(x)
                        x = q+1;
                        Data.(nocorr_name)(x,:) = cell(1,6);
                        Data.(nocorr_name){x,1} = Basis_set;
                    end
                    
                end                

                if isempty(Data.(nocorr_name){x,2})
                    Data.(nocorr_name){x,3} = Tot_Time;
                    Data.(nocorr_name){x,2} = cell(1,19);
                    Data.(nocorr_name){x,2}{1} = Theory_type_no_gcp; % Theory type
                    Data.(nocorr_name){x,2}{2} = a; % A
                    Data.(nocorr_name){x,2}{3} = b;
                    Data.(nocorr_name){x,2}{4} = c;
                    Data.(nocorr_name){x,2}{5} = BL; % Bond length
                    Data.(nocorr_name){x,2}{6} = volume; % Volume
                    Data.(nocorr_name){x,2}{7} = AU_Coords; % Frac coords of asym unit
                    Data.(nocorr_name){x,2}{8} = energy_nocorr; % No-dispersion energy
                    Data.(nocorr_name){x,2}{9} = Tot_Time; % No-dispersion energy
                    Data.(nocorr_name){x,2}{10} = 0; % Data type
                    Data.(nocorr_name){x,2}{11} = alpha; % alpha angle (degrees)
                    Data.(nocorr_name){x,2}{12} = beta; % beta angle (degrees)
                    Data.(nocorr_name){x,2}{13} = gamma; % gamma angle (degrees)
                    Data.(nocorr_name){x,2}{14} = Toldee; % Convergence criteria
                    Data.(nocorr_name){x,2}{15} = 0; % dispersion energy
                    Data.(nocorr_name){x,2}{16} = 0; % gcp energy
                    Data.(nocorr_name){x,2}{17} = MX_BL; % MX bond lengths
                    Data.(nocorr_name){x,2}{18} = MM_BL; % MM bond lengths
                    Data.(nocorr_name){x,2}{19} = XX_BL; % XX bond lengths
                else
                    if ~ismemberf(BL,Data.(nocorr_name){x,2}{1,5},'tol',1e-5)
                        Data.(nocorr_name){x,3} = Data.(nocorr_name){x,3} + Tot_Time;
                        Data.(nocorr_name){x,2}{2}(end+1) = a;
                        Data.(nocorr_name){x,2}{3}(end+1) = b;
                        Data.(nocorr_name){x,2}{4}(end+1) = c;
                        Data.(nocorr_name){x,2}{5}(end+1) = BL;
                        Data.(nocorr_name){x,2}{6}(end+1) = volume;
                        Data.(nocorr_name){x,2}{7}(end+1) = AU_Coords;
                        Data.(nocorr_name){x,2}{8}(end+1) = energy_nocorr;
                        Data.(nocorr_name){x,2}{9}(end+1) = Tot_Time;
                        Data.(nocorr_name){x,2}{10}(end+1) = 0; % Data type
                        Data.(nocorr_name){x,2}{11}(end+1) = alpha; % alpha angle (degrees)
                        Data.(nocorr_name){x,2}{12}(end+1) = beta; % beta angle (degrees)
                        Data.(nocorr_name){x,2}{13}(end+1) = gamma; % gamma angle (degrees)
                        Data.(nocorr_name){x,2}{14}(end+1) = Toldee; % Convergence criteria
                        Data.(nocorr_name){x,2}{15}(end+1) = 0; % dispersion energy
                        Data.(nocorr_name){x,2}{16}(end+1) = 0; % gcp energy
                        Data.(nocorr_name){x,2}{17}(end+1) = MX_BL; % MX bond lengths
                        Data.(nocorr_name){x,2}{18}(end+1) = MM_BL; % MM bond lengths
                        Data.(nocorr_name){x,2}{19}(end+1) = XX_BL; % XX bond lengths

                        % Ensure data is sorted correctly by BL
                        [Data.(nocorr_name){x,2}{5},Indx] = sort(Data.(nocorr_name){x,2}{5},'ascend');
                        for y = [2:4 6:length(Data.(nocorr_name){x,2})]
                            Data.(nocorr_name){x,2}{y} = Data.(nocorr_name){x,2}{y}(Indx);
                        end
                    end
                end
                
                
            % For dispersion calculations
            elseif dispersion

                % Type of theory (no dispersion)
                Theory_type_no_disp = regexprep(Theory_type,'_D[2-4](TB)*','');

                % Input energy (no dispersion)               
                nocorr_name = [Label '_' Theory_type_no_disp];

                % Add to structure
                if ~isfield(Data,nocorr_name)
                    Data.(nocorr_name) = cell(1,6);
                    Data.(nocorr_name){1,1} = Basis_set;
                    x = 1;
                else
                    % Find the basis set
                    x = NaN;
                    for q = 1:size(Data.(nocorr_name),1)
                        
                        if strcmp(Data.(nocorr_name){q,1},Basis_set)
                            x = q;
                            break
                        end
                    end
                    
                    % If basis set not yet present, add to end
                    if isnan(x)
                        x = q+1;
                        Data.(nocorr_name)(x,:) = cell(1,6);
                        Data.(nocorr_name){x,1} = Basis_set;
                    end
                    
                end                

                if isempty(Data.(nocorr_name){x,2})
                    Data.(nocorr_name){x,3} = Tot_Time;
                    Data.(nocorr_name){x,2} = cell(1,19);
                    Data.(nocorr_name){x,2}{1} = Theory_type_no_disp; % Theory type
                    Data.(nocorr_name){x,2}{2} = a; % A
                    Data.(nocorr_name){x,2}{3} = b;
                    Data.(nocorr_name){x,2}{4} = c;
                    Data.(nocorr_name){x,2}{5} = BL; % Bond length
                    Data.(nocorr_name){x,2}{6} = volume; % Volume
                    Data.(nocorr_name){x,2}{7} = AU_Coords; % Frac coords of asym unit
                    Data.(nocorr_name){x,2}{8} = energy_nocorr; % No-dispersion energy
                    Data.(nocorr_name){x,2}{9} = Tot_Time; % No-dispersion energy
                    Data.(nocorr_name){x,2}{10} = 0; % Data type
                    Data.(nocorr_name){x,2}{11} = alpha; % alpha angle (degrees)
                    Data.(nocorr_name){x,2}{12} = beta; % beta angle (degrees)
                    Data.(nocorr_name){x,2}{13} = gamma; % gamma angle (degrees)
                    Data.(nocorr_name){x,2}{14} = Toldee; % Convergence criteria
                    Data.(nocorr_name){x,2}{15} = 0; % dispersion energy
                    Data.(nocorr_name){x,2}{16} = 0; % gcp energy
                    Data.(nocorr_name){x,2}{17} = MX_BL; % MX bond lengths
                    Data.(nocorr_name){x,2}{18} = MM_BL; % MM bond lengths
                    Data.(nocorr_name){x,2}{19} = XX_BL; % XX bond lengths
                else
                    if ~ismemberf(BL,Data.(nocorr_name){x,2}{1,5},'tol',1e-5)
                        Data.(nocorr_name){x,3} = Data.(nocorr_name){x,3} + Tot_Time;
                        Data.(nocorr_name){x,2}{2}(end+1) = a;
                        Data.(nocorr_name){x,2}{3}(end+1) = b;
                        Data.(nocorr_name){x,2}{4}(end+1) = c;
                        Data.(nocorr_name){x,2}{5}(end+1) = BL;
                        Data.(nocorr_name){x,2}{6}(end+1) = volume;
                        Data.(nocorr_name){x,2}{7}(end+1) = AU_Coords;
                        Data.(nocorr_name){x,2}{8}(end+1) = energy_nocorr;
                        Data.(nocorr_name){x,2}{9}(end+1) = Tot_Time;
                        Data.(nocorr_name){x,2}{10}(end+1) = 0;
                        Data.(nocorr_name){x,2}{11}(end+1) = alpha; % alpha angle (degrees)
                        Data.(nocorr_name){x,2}{12}(end+1) = beta; % beta angle (degrees)
                        Data.(nocorr_name){x,2}{13}(end+1) = gamma; % gamma angle (degrees)
                        Data.(nocorr_name){x,2}{14}(end+1) = Toldee; % Convergence criteria
                        Data.(nocorr_name){x,2}{15}(end+1) = 0; % dispersion energy
                        Data.(nocorr_name){x,2}{16}(end+1) = 0; % gcp energy
                        Data.(nocorr_name){x,2}{17}(end+1) = MX_BL; % MX bond lengths
                        Data.(nocorr_name){x,2}{18}(end+1) = MM_BL; % MM bond lengths
                        Data.(nocorr_name){x,2}{19}(end+1) = XX_BL; % XX bond lengths

                        % Ensure data is sorted correctly by BL
                        [Data.(nocorr_name){x,2}{5},Indx] = sort(Data.(nocorr_name){x,2}{5},'ascend');
                        for y = [2:4 6:length(Data.(nocorr_name){x,2})]
                            Data.(nocorr_name){x,2}{y} = Data.(nocorr_name){x,2}{y}(Indx);
                        end
                    end
                end
            end
            
            % For QHA analysis data
            if strcmp(cell_param,'QHA')
               if isempty(Data.(structnames{i}){j,3})
                    Data.(structnames{i}){j,3} = Tot_Time;
                else
                    Data.(structnames{i}){j,3} = Data.(structnames{i}){j,3} + Tot_Time;
                end
                if isempty(Data.(structnames{i}){j,6})
                    Data.(structnames{i}){j,6} = cell(1,23);
                    Data.(structnames{i}){j,6}{1} = SS_Size; % Supercell size
                    Data.(structnames{i}){j,6}{2} = a;
                    Data.(structnames{i}){j,6}{3} = b;
                    Data.(structnames{i}){j,6}{4} = c;
                    Data.(structnames{i}){j,6}{5} = BL; % Bond length (A)
                    Data.(structnames{i}){j,6}{6} = volume; % Volume (A^3) per ion pair
                    Data.(structnames{i}){j,6}{7} = AU_Coords; % Fractional coordinates
                    Data.(structnames{i}){j,6}{8} = MURNAGHAN_1944; % Murnaghan EOS
                    Data.(structnames{i}){j,6}{9} = BIRCH_MURNAGHAN_1947;
                    Data.(structnames{i}){j,6}{10} = POIRIER_TARANTOLA_1998;
                    Data.(structnames{i}){j,6}{11} = VINET_1987;
                    Data.(structnames{i}){j,6}{12} = THIRD_ORDER_POLY;
                    Data.(structnames{i}){j,6}{13} = FOURTH_ORDER_POLY;
                    Data.(structnames{i}){j,6}{14} = FIFTH_ORDER_POLY;
                    Data.(structnames{i}){j,6}{15} = MURNAGHAN_1944_Data;
                    Data.(structnames{i}){j,6}{16} = BIRCH_MURNAGHAN_1947_Data;
                    Data.(structnames{i}){j,6}{17} = POIRIER_TARANTOLA_1998_Data;
                    Data.(structnames{i}){j,6}{18} = VINET_1987_Data;
                    Data.(structnames{i}){j,6}{19} = Thermal_Gruneisen;
                    Data.(structnames{i}){j,6}{20} = Thermal_Helmholtz;
                    Data.(structnames{i}){j,6}{21} = Linear_Therm_Exp;
                    Data.(structnames{i}){j,6}{22} = PVT;
                    Data.(structnames{i}){j,6}{23} = THP; % Thermal Helmholtz data vs pressure
                else
                    % Add new row
                    Data.(structnames{i}){j,6} = [Data.(structnames{i}){j,6}; cell(1,23)];
                    
                    Data.(structnames{i}){j,6}{end,1} = SS_Size; % Supercell size
                    Data.(structnames{i}){j,6}{end,2} = a;
                    Data.(structnames{i}){j,6}{end,3} = b;
                    Data.(structnames{i}){j,6}{end,4} = c;
                    Data.(structnames{i}){j,6}{end,5} = BL; % Bond length (A)
                    Data.(structnames{i}){j,6}{end,6} = volume; % Volume (A^3) per ion pair
                    Data.(structnames{i}){j,6}{end,7} = AU_Coords; % Fractional coordinates
                    Data.(structnames{i}){j,6}{end,8} = MURNAGHAN_1944; % Murnaghan EOS
                    Data.(structnames{i}){j,6}{end,9} = BIRCH_MURNAGHAN_1947;
                    Data.(structnames{i}){j,6}{end,10} = POIRIER_TARANTOLA_1998;
                    Data.(structnames{i}){j,6}{end,11} = VINET_1987;
                    Data.(structnames{i}){j,6}{end,12} = THIRD_ORDER_POLY;
                    Data.(structnames{i}){j,6}{end,13} = FOURTH_ORDER_POLY;
                    Data.(structnames{i}){j,6}{end,14} = FIFTH_ORDER_POLY;
                    Data.(structnames{i}){j,6}{end,15} = MURNAGHAN_1944_Data;
                    Data.(structnames{i}){j,6}{end,16} = BIRCH_MURNAGHAN_1947_Data;
                    Data.(structnames{i}){j,6}{end,17} = POIRIER_TARANTOLA_1998_Data;
                    Data.(structnames{i}){j,6}{end,18} = VINET_1987_Data;
                    Data.(structnames{i}){j,6}{end,19} = Thermal_Gruneisen;
                    Data.(structnames{i}){j,6}{end,20} = Thermal_Helmholtz;
                    Data.(structnames{i}){j,6}{end,21} = Linear_Therm_Exp;
                    Data.(structnames{i}){j,6}{end,22} = PVT;
                    Data.(structnames{i}){j,6}{end,23} = THP; % Thermal Helmholtz data vs pressure
                end
                
            % For vibrational analysis data
            elseif strcmp(cell_param,'Vib')
                if isempty(Data.(structnames{i}){j,3})
                    Data.(structnames{i}){j,3} = Tot_Time;
                else
                    Data.(structnames{i}){j,3} = Data.(structnames{i}){j,3} + Tot_Time;
                end
                if isempty(Data.(structnames{i}){j,5})
                    Data.(structnames{i}){j,5} = cell(1,18);
                    Data.(structnames{i}){j,5}{1} = SS_Size; % Supercell size
                    Data.(structnames{i}){j,5}{2} = Zero_Point_Energy; % Zero point vib energy
                    Data.(structnames{i}){j,5}{3} = Total_PDOS; % Total phonon density of states
                    Data.(structnames{i}){j,5}{4} = Vib_E; % Thermal Vibrational Energy
                    Data.(structnames{i}){j,5}{5} = Pres; % Pressure
                    Data.(structnames{i}){j,5}{6} = Temp; % Temperature
                    Data.(structnames{i}){j,5}{7} = Entropy; % Entropy per mole of ion pair
                    Data.(structnames{i}){j,5}{8} = Heat_Capacity; % Heat capacity per mole of ion pair
                    Data.(structnames{i}){j,5}{9} = PV; % Themodynamic work
                    Data.(structnames{i}){j,5}{10} = a;
                    Data.(structnames{i}){j,5}{11} = b;
                    Data.(structnames{i}){j,5}{12} = c;
                    Data.(structnames{i}){j,5}{13} = BL; % Bond length
                    Data.(structnames{i}){j,5}{14} = volume; % Volume per ion pair
                    Data.(structnames{i}){j,5}{15} = AU_Coords; % Frac coords
                    Data.(structnames{i}){j,5}{16} = Total_E; % Total energy per ion pair
                    Data.(structnames{i}){j,5}{17} = Metal_PDOS; % Metal phonon density of states
                    Data.(structnames{i}){j,5}{18} = Halide_PDOS; % Halide phonon density of states 
                else
                    % Add new row
                    Data.(structnames{i}){j,5} = [Data.(structnames{i}){j,5}; cell(1,18)];
                    Data.(structnames{i}){j,5}{end,1} = SS_Size;
                    Data.(structnames{i}){j,5}{end,2} = Zero_Point_Energy; % Zero point vib energy
                    Data.(structnames{i}){j,5}{end,3} = Total_PDOS; % Total Phonon density of states 
                    Data.(structnames{i}){j,5}{end,4} = Vib_E; % Thermal Vibrational Energy
                    Data.(structnames{i}){j,5}{end,5} = Pres; % Pressure
                    Data.(structnames{i}){j,5}{end,6} = Temp; % Temperature
                    Data.(structnames{i}){j,5}{end,7} = Entropy; % Entropy per mole of ion pair
                    Data.(structnames{i}){j,5}{end,8} = Heat_Capacity; % Heat capacity per mole of ion pair
                    Data.(structnames{i}){j,5}{end,9} = PV; % Themodynamic work
                    Data.(structnames{i}){j,5}{end,10} = a;
                    Data.(structnames{i}){j,5}{end,11} = b;
                    Data.(structnames{i}){j,5}{end,12} = c;
                    Data.(structnames{i}){j,5}{end,13} = BL; % Bond length
                    Data.(structnames{i}){j,5}{end,14} = volume; % Volume per ion pair
                    Data.(structnames{i}){j,5}{end,15} = AU_Coords; % Frac coords
                    Data.(structnames{i}){j,5}{end,16} = Total_E; % Total energy per ion pair
                    Data.(structnames{i}){j,5}{end,17} = Metal_PDOS; % Metal phonon density of states
                    Data.(structnames{i}){j,5}{end,18} = Halide_PDOS; % Halide phonon density of states
                end
            else
                if isempty(Data.(structnames{i}){j,2})
                    Data.(structnames{i}){j,3} = Tot_Time;
                    Data.(structnames{i}){j,2} = cell(1,19);
                    Data.(structnames{i}){j,2}{1} = Theory_type;
                    Data.(structnames{i}){j,2}{2} = a;
                    Data.(structnames{i}){j,2}{3} = b;
                    Data.(structnames{i}){j,2}{4} = c;
                    Data.(structnames{i}){j,2}{5} = BL; % Bond length
                    Data.(structnames{i}){j,2}{6} = volume; % Volume per unit cell
                    Data.(structnames{i}){j,2}{7} = AU_Coords; % Frac coords
                    Data.(structnames{i}){j,2}{8} = energy_tot; % Total energy per unit cell
                    Data.(structnames{i}){j,2}{9} = Tot_Time;
                    Data.(structnames{i}){j,2}{10} = data_type; % Data type
                    Data.(structnames{i}){j,2}{11} = alpha; % alpha angle (degrees)
                    Data.(structnames{i}){j,2}{12} = beta; % beta angle (degrees)
                    Data.(structnames{i}){j,2}{13} = gamma; % gamma angle (degrees)
                    Data.(structnames{i}){j,2}{14} = Toldee; % Convergence criteria
                    Data.(structnames{i}){j,2}{15} = disp_E; % dispersion energy
                    Data.(structnames{i}){j,2}{16} = gCP_E; % gcp energy
                    Data.(structnames{i}){j,2}{17} = MX_BL; % MX bond lengths
                    Data.(structnames{i}){j,2}{18} = MM_BL; % MM bond lengths
                    Data.(structnames{i}){j,2}{19} = XX_BL; % XX bond lengths
                else
                    Data.(structnames{i}){j,3} = Data.(structnames{i}){j,3} + Tot_Time;
                    Data.(structnames{i}){j,2}{2}(end+1) = a;
                    Data.(structnames{i}){j,2}{3}(end+1) = b;
                    Data.(structnames{i}){j,2}{4}(end+1) = c;
                    Data.(structnames{i}){j,2}{5}(end+1) = BL;
                    Data.(structnames{i}){j,2}{6}(end+1) = volume;
                    Data.(structnames{i}){j,2}{7}(end+1) = AU_Coords; % Frac coords
                    Data.(structnames{i}){j,2}{8}(end+1) = energy_tot;
                    Data.(structnames{i}){j,2}{9}(end+1) = Tot_Time;
                    Data.(structnames{i}){j,2}{10}(end+1) = data_type; % Data type
                    Data.(structnames{i}){j,2}{11}(end+1) = alpha; % alpha angle (degrees)
                    Data.(structnames{i}){j,2}{12}(end+1) = beta; % beta angle (degrees)
                    Data.(structnames{i}){j,2}{13}(end+1) = gamma; % gamma angle (degrees)
                    Data.(structnames{i}){j,2}{14}(end+1) = Toldee; % Convergence criteria
                    Data.(structnames{i}){j,2}{15}(end+1) = disp_E; % dispersion energy
                    Data.(structnames{i}){j,2}{16}(end+1) = gCP_E; % gcp energy
                    Data.(structnames{i}){j,2}{17}(end+1) = MX_BL; % MX bond lengths
                    Data.(structnames{i}){j,2}{18}(end+1) = MM_BL; % MM bond lengths
                    Data.(structnames{i}){j,2}{19}(end+1) = XX_BL; % XX bond lengths
                end
            end
        end
    end
    
    % Remove empty data (counting in reverse) and sort
    for j=length(subFolders):-1:1
        if isempty(Data.(structnames{i}){j,2}) && isempty(Data.(structnames{i}){j,5})
            Data.(structnames{i})(j,:) = [];
        elseif isempty(Data.(structnames{i}){j,2})
            % Do nothing
        else
            % Sort the data by increasing BL
            [Data.(structnames{i}){j,2}{5},Index] = sort(Data.(structnames{i}){j,2}{5},'ascend');
            for id = [2:4 6:length(Data.(structnames{i}){j,2})]
                Data.(structnames{i}){j,2}{id} = Data.(structnames{i}){j,2}{id}(Index);
            end
        end
    end
    if isempty(Data.(structnames{i}))
        Data = rmfield(Data,structnames{i});
    end
    
    if (i/K) < 0.9999
        fprintf('%s',[num2str(round(i*100/K)) '%... '])
    end
end
fprintf('%s',['100%.' newline])

%% Add DFT/dispersion markers
fields = fieldnames(Data);
N = length(fields);

for i=1:length(DFT_methods)
    DFT_methods{i} = [Label '_' DFT_methods{i}];
end

for i=1:length(Range_cor)
    Range_cor{i} = [Label '_' Range_cor{i}];
end

for i=1:N
    M = size(Data.(fields{i}),1);
    if ismember(fields{i},Range_cor)
        for j=1:M
            Data.(fields{i}){j,4} = 4; % Range corrected DFT
        end
    elseif ismember(fields{i},DFT_methods) % DFT method
        
        if ~isempty(regexp(fields{i},'_D[2-4].*','ONCE'))
            for j=1:M
                Data.(fields{i}){j,4} = 3; % Dispersion corrected DFT
            end
        else 
            for j=1:M
                Data.(fields{i}){j,4} = 2; % Uncorrected DFT
            end
        end
    else
        for j=1:M
            Data.(fields{i}){j,4} = 1; % HF or post HF
        end
    end
end

end


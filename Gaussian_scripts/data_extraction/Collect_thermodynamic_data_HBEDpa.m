% Collect_thermodynamic_data_HBEDpa
% Define study location and reference parameters
study_loc = '/home/scheiber/project/HBEDpa_study';
ref_water_cluster = 'H2O_Cluster_23';
ref_water_number = 23;
ref_ligand = 'HBEDpa';

% Open a data summary file
cd(study_loc)
fileID = fopen('Data_summary.csv','wt');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','System',...
    'Metal','Spin Mult.','Total G (kJ/mol)','Pot E (kJ/mol)','Thermal E (kJ/mol)',...
    'Thermal G (kJ/mol)', 'ZPVE (kJ/mol)','Vib E (kJ/mol)',...
    'Vib -TS (kJ/mol)','Rot E (kJ/mol)', 'Rot -TS (kJ/mol)',...
    'Trans E (kJ/mol)','Trans -TS (kJ/mol)','Total -TS (kJ/mol)');

%% Gather reference data (ion energies, ref water energy, etc)
% First gather water energies
ref_H2O = [study_loc filesep 'Reference_fragments' filesep ref_water_cluster filesep '*.log*'];
d_H2O = dir(ref_H2O);

H2O_txt = fileread(fullfile(d_H2O.folder,d_H2O.name));
H2O = find_thermal_energies(H2O_txt,ref_water_number);
Multiplicity = 1;

fprintf(fileID, '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', d_H2O.name,'',Multiplicity,H2O.Total_G,...
            H2O.Pot_E,H2O.Thermal_E,H2O.Thermal_G,H2O.ZPVE,H2O.Vib_E,H2O.Vib_TS,...
            H2O.Rot_E,H2O.Rot_TS,H2O.Trans_E,H2O.Trans_TS,H2O.Tot_TS);
Data.H2O = H2O;

% Gather reference ligand energies
ref_ligand = [study_loc filesep 'Reference_fragments' filesep ref_ligand filesep '*.log*'];
d_ligand = dir(ref_ligand);

% Find number of H2O in ref ligand
num_H2O = regexpi(d_ligand.name,'([0-9]+)_H2O','tokens','once');
if isempty(num_H2O)
    num_H2O = 0;
else
    num_H2O = str2double(num_H2O{1});
end

% Read reference ligand log file
ligand_txt = fileread(fullfile(d_ligand.folder,d_ligand.name));

% Find energies and correct for solvent molecules
Free_L = find_thermal_energies(ligand_txt);
Free_L = subtract_solvent(Free_L,H2O,num_H2O);

fprintf(fileID, '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', d_ligand.name,'',Multiplicity,Free_L.Total_G,...
            Free_L.Pot_E,Free_L.Thermal_E,Free_L.Thermal_G,Free_L.ZPVE,Free_L.Vib_E,Free_L.Vib_TS,...
            Free_L.Rot_E,Free_L.Rot_TS,Free_L.Trans_E,Free_L.Trans_TS,Free_L.Tot_TS);
Data.Ligand = Free_L;

% Gather metal ion energies
% List all directories
d = dir([study_loc filesep 'Reference_fragments']);
dfolders = d([d(:).isdir]);
Calc_Directories = dfolders(~ismember({dfolders(:).name},{'.','..'}));

% Find directories labelled [Metal]_Ion
Metal_calcs = regexp({Calc_Directories(:).name},'(.*)_Ion','once','tokens');
Metal_calcs_idx = ~cellfun(@isempty,Metal_calcs);
Ions = Metal_calcs(Metal_calcs_idx);
Ions = [Ions{:}];

% Loop through and gather energies of metal Ion aqua clusters
for idx = 1:length(Ions)
    Ion = Ions{idx};
    
    ref_M = [study_loc filesep 'Reference_fragments' filesep Ion '_Ion'];
    d_M = dir(ref_M);
    d_M = d_M([d_M(:).isdir]);
    d_M = d_M(~ismember({d_M(:).name},{'.','..'}));
    
    % Loop through each subdirectory and find the log file(s)
    M_dat_ref.Total_G = Inf;
    num_H2O_ref = 18;
    for jdx = 1:length(d_M)
        
        % Find number of H2O
        num_H2O = regexpi(d_M(jdx).name,'([0-9]+)_H2O','tokens','once');
        if isempty(num_H2O)
            num_H2O = 0;
        else
            num_H2O = str2double(num_H2O{1});
        end
        
        % Find the spin multiplicity
        spin_m = regexpi(d_M(jdx).name,'HighSpin|LowSpin','match','once');
        if isempty(spin_m)
            Multiplicity = 1;
            Multip_txt = 'Singlet';
        elseif strcmpi('HighSpin',spin_m)
            Multiplicity = 6;
            Multip_txt = 'Sextet';
        else
            Multiplicity = 2;
            Multip_txt = 'Doublet';
        end
        
        % Find the log file(s)
        ref_Mj = fullfile(d_M(jdx).folder,d_M(jdx).name,'*.log*');
        d_Mj = dir(ref_Mj);
        
        % If multiple log files exist, loop through them and find the one
        % with thermal information
        for kdx = 1:length(d_Mj)
            M_txt = fileread( fullfile(d_Mj(kdx).folder,d_Mj(kdx).name) );
            M_dat = find_thermal_energies(M_txt);
            if ~isnan(M_dat.Total_G)
                M_dat = subtract_solvent(M_dat,H2O,num_H2O);
                Mf_Name = d_Mj(kdx).name;
                break
            end
        end
        
        % Print all results to file
        fprintf(fileID, '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', Mf_Name,Ion,Multiplicity,M_dat.Total_G,...
            M_dat.Pot_E,M_dat.Thermal_E,M_dat.Thermal_G,M_dat.ZPVE,M_dat.Vib_E,M_dat.Vib_TS,...
            M_dat.Rot_E,M_dat.Rot_TS,M_dat.Trans_E,M_dat.Trans_TS,M_dat.Tot_TS);
        Data.(Ion).(Multip_txt).Ion = M_dat;

        % Pick out the lowest energy data with the least H2O as the reference
        if(M_dat.Total_G < M_dat_ref.Total_G) && num_H2O == num_H2O_ref %& num_H2O >= num_H2O_ref
            M_dat_ref = M_dat;
            num_H2O_ref = num_H2O;
        end
    end
    
    % Add to ref data structure
    Ref_Data.(Ion) = M_dat_ref;
end

%% Gather complex data
for idx = 1:length(Ions)
    Ion = Ions{idx};
    M_complex_dir = [study_loc filesep Ion '_complexes'];
    
    % Find all calculations within the current metal complex directory
    d_M = dir(M_complex_dir);
    d_M = d_M([d_M(:).isdir]);
    d_M = d_M(~ismember({d_M(:).name},{'.','..'}));
    
    
    % Loop through each subdirectory and find the log file(s)
    for jdx = 1:length(d_M)
        
        % Find number of H2O
        num_H2O = regexpi(d_M(jdx).name,'([0-9]+)_H2O','tokens','once');
        if isempty(num_H2O)
            num_H2O = 0;
        else
            num_H2O = str2double(num_H2O{1});
        end
        
        % Find the spin multiplicity
        spin_m = regexpi(d_M(jdx).name,'HighSpin|LowSpin','match','once');
        if isempty(spin_m)
            Multiplicity = 1;
            Multip_txt = 'Singlet';
        elseif strcmpi('HighSpin',spin_m)
            Multiplicity = 6;
            Multip_txt = 'Sextet';
        else
            Multiplicity = 2;
            Multip_txt = 'Doublet';
        end
        
        % Find the log file(s)
        ref_Mj = fullfile(d_M(jdx).folder,d_M(jdx).name,'*.log*');
        d_Mj = dir(ref_Mj);
        
        % If multiple log files exist, loop through them and find the one
        % with thermal information
        M_dat = find_thermal_energies('');
        Mf_Name = d_M(jdx).name;
        for kdx = 1:length(d_Mj)
            M_txt = fileread( fullfile(d_Mj(kdx).folder,d_Mj(kdx).name) );
            M_dat = find_thermal_energies(M_txt);
            if ~isnan(M_dat.Total_G)
                % Subtract the energy of any solvent molecules, the free
                % ligand, and the free metal
                M_dat = subtract_solvent(M_dat,H2O,num_H2O);
                M_dat = subtract_solvent(M_dat,Free_L,1);
                M_dat = subtract_solvent(M_dat,Ref_Data.(Ion),1);
                Mf_Name = d_Mj(kdx).name;
                break
            end
        end
        
        % Print all results to file
        fprintf(fileID, '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', Mf_Name,Ion,Multiplicity,M_dat.Total_G,...
            M_dat.Pot_E,M_dat.Thermal_E,M_dat.Thermal_G,M_dat.ZPVE,M_dat.Vib_E,M_dat.Vib_TS,...
            M_dat.Rot_E,M_dat.Rot_TS,M_dat.Trans_E,M_dat.Trans_TS,M_dat.Tot_TS);
        
        % Save to data file
        Structure = strrep(d_M(jdx).name,[Ion '_'],'');
        if ~isempty(spin_m)
            Structure = strrep(Structure,['_' spin_m],'');
        end
        Data.(Ion).(Multip_txt).(Structure) = M_dat;
    end
end
fclose(fileID);

% Save data structure as mat file
save(fullfile(study_loc,'Data_Summary.mat'),'Data');
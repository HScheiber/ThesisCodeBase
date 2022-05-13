J_cal = 4.184;
kJ_mol_per_Ha = 2625.4996394799;
T = 298.15;

% List all directories
d = dir;
dfolders = d([d(:).isdir]);
Calc_Categories = dfolders(~ismember({dfolders(:).name},{'.','..'}));

fileID = fopen('Data_summary.csv','wt');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','System','Total G (Ha)','Pot E (Ha)','Thermal E (kJ/mol)',...
    'Thermal G (kJ/mol)', 'ZPE (kJ/mol)','Vib E (kJ/mol)',...
    'Vib -TS (kJ/mol)','Rot E (kJ/mol)', 'Rot -TS (kJ/mol)','Total -TS');

for idx = 1:length(Calc_Categories)
    Calc_Category = [Calc_Categories(idx).folder filesep Calc_Categories(idx).name];
    
    % List all subdirectories
    d = dir(Calc_Category);
    dfolders = d([d(:).isdir]);
    Calc_jobs = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    
    for jdx = 1:length(Calc_jobs)
        Current_job = [Calc_Category filesep Calc_jobs(jdx).name];
        
        % Find all log files in current folder
        d = dir([Current_job filesep '*log*']);
        
        
        % Search each log file for results of interest
        for kdx = 1:length(d)
            jobtxt = fileread(fullfile(Current_job,d(kdx).name));
            
            thermtxt = regexp(jobtxt,'Zero-point correction=(.*?\n){16}','match','once');
            
            if isempty(thermtxt)
                continue
            else
                
                ZPE = regexp(thermtxt,'Zero-point correction= +([0-9]|\.|-)+','tokens','once');
                ZPE = str2double(ZPE{1}); % Ha
                Thermal_E = regexp(thermtxt,'Thermal correction to Energy= +([0-9]|\.|-)+','tokens','once');
                Thermal_E = str2double(Thermal_E{1})*kJ_mol_per_Ha;
                Thermal_G = regexp(thermtxt,'Thermal correction to Gibbs Free Energy= +([0-9]|\.|-)+','tokens','once');
                Thermal_G = str2double(Thermal_G{1})*kJ_mol_per_Ha;
                Total_G = regexp(thermtxt,'Sum of electronic and thermal Free Energies= +([0-9]|\.|-)+','tokens','once');
                Total_G = str2double(Total_G{1});
                Sum_ZPE = regexp(thermtxt,'Sum of electronic and zero-point Energies= +([0-9]|\.|-)+','tokens','once');
                Sum_ZPE = str2double(Sum_ZPE{1});
                Pot_E = Sum_ZPE - ZPE;
                ZPE = ZPE*kJ_mol_per_Ha; % kJ/mol
                Rot_ES = regexp(thermtxt,'Rotational +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
                Rot_E = str2double(Rot_ES{1})*J_cal; % units in kJ/mol
                Rot_TS = -str2double(Rot_ES{3})*J_cal*T/1000; % units in Jk/mol-K
                
                Vib_ES = regexp(thermtxt,'Vibrational +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
                Vib_E = str2double(Vib_ES{1})*J_cal; % units in kJ/mol
                Vib_TS = -str2double(Vib_ES{3})*J_cal*T/1000; % units in kJ/mol-K
                
                Tot_ES = regexp(thermtxt,'Total +([0-9]|\.|-)+ +([0-9]|\.|-)+ +([0-9]|\.|-)+','tokens','once');
                Tot_E = str2double(Tot_ES{1})*J_cal; % units in kJ/mol
                Tot_TS = -str2double(Tot_ES{3})*J_cal*T/1000; % units in kJ/mol
                
                fprintf(fileID, '%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', d(kdx).name,Total_G,...
                    Pot_E,Thermal_E,Thermal_G,ZPE,Vib_E,Vib_TS,Rot_E,Rot_TS,Tot_TS);
                
            end
            
        end
    end
end
fclose(fileID);
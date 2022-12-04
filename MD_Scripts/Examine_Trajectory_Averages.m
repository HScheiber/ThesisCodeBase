DataDir = 'D:\Molten_Salts_MD\NaCl';
ProjectLabel = 'StrucRef';
gmx = 'gmx_d';
wsl = 'wsl source ~/.bashrc; ';
pipe = ' ^| ';

files = dir(fullfile(DataDir,[ProjectLabel '*']));
tempdirloc = tempname;
if ~isfolder(tempdirloc)
    mkdir(tempdirloc)
end

for idx = 1:length(files)
    FileDir = fullfile(DataDir,files(idx).name);
    Ener_file = fullfile(FileDir,[files(idx).name '.edr']);
    Tpr_file = fullfile(FileDir,[files(idx).name '.tpr']);
    Gro_file = fullfile(FileDir,[files(idx).name '.gro']);
    Out_file = fullfile(tempdirloc,[files(idx).name '.xvg']);
    
    gmx_command = [wsl 'echo 0' pipe gmx ' energy -f ' windows2unix(Ener_file)...
        ' -s ' windows2unix(Tpr_file)];
    [err,outpt] = system(gmx_command);
    if err ~= 1
        warndlg('Problem with energy check.')
        return
    end
    en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
    en_opts = en_opts{1};
    En_set = char(regexp(en_opts,'([0-9]{1,2})  Potential','tokens','once'));
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Kinetic-En.','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Total-Energy','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Temperature','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Pressure','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Density','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Enthalpy','tokens','once'))];
    En_set = [En_set ' 0'];
    
    gmx_command = [wsl 'echo ' En_set pipe gmx ' energy -f ' windows2unix(Ener_file)...
        ' -o ' windows2unix(Out_file) ' -s ' windows2unix(Tpr_file) ...
        ' -b 1000 -e 5000'];
    [err,~] = system(gmx_command);
    
    if err ~= 0
        warndlg('Failed to collect data.')
        return
    end

    NSS = load_gro_file(Gro_file);
    NF = NSS.N_atoms/2; % number of ion pairs in system

    Data = import_xvg_as_struc(Out_file); % Gather data
    delete(Out_file) % remove temp output file
    
    PotE =  mean(Data.Potential./NF); % kJ/mol
    KinE =  mean(Data.Kinetic_Energy./NF); % kJ/mol
    TotE =  mean(Data.Total_Energy./NF); % kJ/mol
    Temp =  mean(Data.Temperature); % K
    Press = mean(Data.Pressure); % Bar
    Vol =   mean(Data.Volume.*(1e-7^3).*(6.0221409e+23)./NF); % cm^3/mol
    Dens =  mean(Data.Density); % kg^3/m^3
    Enth =  mean(Data.Enthalpy./NF); % kJ/mol
    
    disp(repmat('*',1,50))
    disp(['System: ' files(idx).name])
    disp(repmat('*',1,50))
    disp(['Average Potential Energy: ' num2str(PotE) ' kJ/mol'])
    disp(['Average Kinetic Energy: ' num2str(KinE) ' kJ/mol'])
    disp(['Average Total Energy: ' num2str(TotE) ' kJ/mol'])
    disp(['Average Temperature: ' num2str(Temp) ' K'])
    disp(['Average Pressure: ' num2str(Press) ' bar'])
    disp(['Average Molar Volume: ' num2str(Vol) ' cm^3/mol'])
    disp(['Average Density: ' num2str(Dens) ' kg^3/m^3'])
    disp(['Average Enthalpy: ' num2str(Enth) ' kJ/mol'])
end

rmdir(tempdirloc)
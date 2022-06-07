Settings = Initialize_MD_Settings;
Settings.Project_Directory_Name = 'Melting_Point_Studies';
DataSetName = 'Melting_Point_Data.mat';
DataKeyword = 'MP_Prod2';
ProjectDir = fullfile(Settings.project,Settings.Project_Directory_Name);
SaveDataDir = fullfile(Settings.home,'data',DataSetName);
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' ...
         'NaF' 'NaCl' 'NaBr' 'NaI' ...
         'KF' 'KCl' 'KBr' 'KI' ...
         'RbF' 'RbCl' 'RbBr' 'RbI' ...
         'CsF' 'CsCl' 'CsBr' 'CsI'};
Theories = {'JC3P' 'JCSD' 'JC4P' 'JC' 'TF'};
Structures = {'Rocksalt' 'Wurtzite' 'CsCl' 'FiveFive'};
Data = load(SaveDataDir,'Data').Data;


for idx = 1:length(Salts)
    Salt = Salts{idx};
    for jdx = 1:length(Theories)
        Theory = Theories{jdx};
        for kdx = 1:length(Structures)
            Structure = Structures{kdx};
            try
                T_dat = Data.(Salt).(Theory).(Structure).(DataKeyword);
                if any(T_dat.Freeze_Trace & T_dat.Melt_Trace)
                    disp([Salt ' ' Theory ' ' Structure ': Tm = ' num2str(mean(T_dat.dT),'%.1f') ])
                end
            catch
            end
        end
    end
end
% Run_Project_Structure_Classifier
Keyword = 'ExpRef';
Project_dir = '/home/haydensc/project/Molten_Salts_MD';
Scratch_dir = '/home/haydensc/scratch/Molten_Salts_MD';
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl' 'CsCl'};

for idx = 1:length(Salts)
    Salt = Salts{idx};
    
    calcs = dir(fullfile(Project_dir,Salt,[Keyword '*']));
    
    for jdx = 1:length(calcs)
        JobName = calcs(jdx).name;
        WorkDir = fullfile(Project_dir,Salt,JobName);
        SaveDir = fullfile(Scratch_dir,Salt,JobName);
        
        if ~isfolder(SaveDir)
            mkdir(SaveDir)
        end
        
        % Check if job already completed
        if ~isfile(fullfile(SaveDir,['V2NN-5-1_' Salt '_' JobName '.png']))
            disp(['Starting Job: ' Salt ' ' JobName '.'])
            py.LiXStructureDetector.Check_Structures(WorkDir, Salt,...
                        pyargs('SystemName',JobName,...
                        'SaveTrajectory',true,...
                        'SavePredictionsImage',true,...
                        'ML_TimeLength',5,...
                        'ML_TimeStep',1,...
                        'TimePerFrame',20,...
                        'FileType','gro',...
                        'Version',2,...
                        'SaveDir',SaveDir,...
                        'Verbose',false));
            disp(['Job Complete: ' Salt ' ' JobName '.'])
        else
            disp(['Job ' Salt ' ' JobName ' already complete, skipping.'])
        end
    end
end




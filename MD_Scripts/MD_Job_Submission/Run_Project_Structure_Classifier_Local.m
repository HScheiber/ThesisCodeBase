% Run_Project_Structure_Classifier_Local
Project_dir = '/mnt/e/Molten_Salts_MD';
ML_Models = {[0 0] [0 1] [3 1] [5 1] [11 1] [11 5] [21 5] [51 5] [51 10]};
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl' 'CsCl'};

for ML_idx = 1:length(ML_Models)
    
    ML_TimeLength = ML_Models{ML_idx}(1);
    ML_TimeStep = ML_Models{ML_idx}(2);
    ML_Model = ['V2NN-' num2str(ML_TimeLength) '-' num2str(ML_TimeStep)];
    
    for idx = 1:length(Salts)
        Salt = Salts{idx};
        
        % Find calculation folders in project and discard hidden folders
        calcs = dir(fullfile(Project_dir,Salt));
        for i = length(calcs):-1:1
        	dirname = calcs(i).name;
        	if strncmp( dirname,'.',1) || ~calcs(i).isdir
            	calcs(i) = [];
        	end
        end

        for jdx = 1:length(calcs)
            JobName = calcs(jdx).name;
            WorkDir = fullfile(Project_dir,Salt,JobName);

            % Check if job already completed
            if ~isfile(fullfile(WorkDir,[ML_Model '_' Salt '_' JobName '.png']))
                disp(['Starting Job: ' Salt ' ' JobName '.'])
                py.LiXStructureDetector.Check_Structures(WorkDir, Salt,...
                            pyargs('SystemName',JobName,...
                            'SaveTrajectory',true,...
                            'SavePredictionsImage',true,...
                            'ML_TimeLength',ML_TimeLength,...
                            'ML_TimeStep',ML_TimeStep,...
                            'TimePerFrame',20,...
                            'FileType','gro',...
                            'Version',2,...
                            'SaveDir',WorkDir,...
                            'Verbose',true,...
                            'InMemory',true));
                disp(['Job Complete: ' Salt ' ' JobName '.'])
            else
                disp(['Job ' Salt ' ' JobName ' already complete, skipping.'])
            end
        end
    end
end
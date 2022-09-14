

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; %  'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Theory = 'BH';
ModelID = 'MC';
BestOnly = false;
Reps = [1:5];

ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';

% Find reps of models
N_Salts = numel(Salts);
Models = cell(1,length(Salts));
for idx = 1:N_Salts
    files = {dir(fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' ModelID '*_data.mat'])).name};
    Models{idx} = unique(regexp(files,[ModelID '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models{idx});
    Models{idx} = Models{idx}(nonempty_idx);
end
Models = unique(horzcat(Models{:}));
N_Reps = numel(Reps);
Models_of_Interest = cell(1,N_Reps);
for idx = 1:N_Reps
    Models_of_Interest{idx} = [ModelID num2str(Reps(idx))];
end
Models = intersect(Models_of_Interest,Models);
N_Models = numel(Models);


Settings = Initialize_MD_Settings;
for idx = 1:N_Salts
    Salt = Salts{idx};
    Settings.Salt = Salt;
    for iidx = 1:N_Models
        Model = Models{iidx};

        % Find the fully optimized model
        dat_file = fullfile(ML_results_dir,Salt,[Salt '_' Theory '_Model_' Model '_data.mat']);

        if ~isfile(dat_file)
            disp(['Could not load results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end

        % Load data
        try
            data = load(dat_file).full_data;
            params = data.full_opt_point;
            Settings = data.Settings;
            Ptypes = bayesopt_params(Settings);
            Param_names = {Ptypes.Name};
            
            Param = table();
            for jdx = 1:numel(Param_names)
                Param.(Param_names{jdx}) = params(jdx);
            end
            Settings = Get_Scaling_Params(Settings,Param);
            
            disp(repmat('*',1,50))
            disp([Salt ', ' Theory ', Model ' Model '.'])
            disp(repmat('*',1,50))
            for jdx = 1:numel(Param_names)
                disp([Param_names{jdx} ' : ' num2str(Param.(Param_names{jdx}))]);
            end
            disp(repmat('*',1,50))
            disp(['C6_MM :' num2str(Settings.S.D.MM)])
            disp(['C6_XX :' num2str(Settings.S.D.XX)])
            disp(['C6_MX :' num2str(Settings.S.D.MX)])
            disp(['B_MM :' num2str(Settings.S.R.MM)])
            disp(['B_XX :' num2str(Settings.S.R.XX)])
            disp(['B_MX :' num2str(Settings.S.R.MX)])
            disp(['Alpha_MM :' num2str(Settings.S.A.MM)])
            disp(['Alpha_XX :' num2str(Settings.S.A.XX)])
            disp(['Alpha_MX :' num2str(Settings.S.A.MX)])
            
        catch
            disp(['Skipping: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end
    end
end


% Pars order: 






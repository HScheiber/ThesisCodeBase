% Fix_RR_Problem
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'}; %  'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Theory = 'BH';
ModelID = 'LA';
%ML_results_dir = 'C:\Users\Hayden\Documents\Patey_Lab\BO_Models';
ML_results_dir = '/home/user/project/BO_Models';


N_Salts = numel(Salts);
% Find reps of models
Models = cell(1,length(Salts));
for idx = 1:N_Salts
    files = {dir(fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' ModelID '*_data.mat'])).name};
    Models{idx} = unique(regexp(files,[ModelID '([0-9])+'],'match','once'));
    nonempty_idx = ~cellfun(@isempty,Models{idx});
    Models{idx} = Models{idx}(nonempty_idx);
end
Models = unique(horzcat(Models{:}));
N_Models = numel(Models);

% Make sure the models are sorted
if ~isempty(Models)
    [~,idx] = sort(cellfun(@str2double,regexp(Models,'[0-9]+','match','once')));
    Models = Models(idx);
end

% Loop through all models
for idx = 1:N_Salts
    for jdx = 1:N_Models
        dat_file = fullfile(ML_results_dir,Salts{idx},[Salts{idx} '_' Theory '_Model_' Models{jdx} '_data.mat']);
        if ~isfile(dat_file)
            disp(['Could not load results found for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end
        
        % Load data
        try
            full_data = load(dat_file).full_data;
            ptable = table;
            params = bayesopt_params(full_data.Settings);
            pnames = {params.Name};
            full_opt_point = full_data.full_opt_point;
            Minimization_Data = full_data.Minimization_Data;
            for kdx = 1:numel(pnames)
                ptable.(pnames{kdx}) = full_opt_point(kdx);
            end
            Structures = cell(1,numel(Minimization_Data));
            for kdx = 1:numel(Minimization_Data)
                Structures{kdx} = Minimization_Data{kdx}.Structure;
            end
        catch
            disp(['Could not obtain crystal minimization data for: ' Salt ', ' Theory ', Model ' Model '.']);
            continue
        end
            
        switch Theory
        case 'JC'
            if ptable.sigma_MM > ptable.sigma_XX
                disp([Salts{idx} ' ' Theory ' ' Models{jdx} ' has reversed radius ratios! Switching...']);
                sigmas_idx = [find(strcmp(pnames,'sigma_MM')) find(strcmp(pnames,'sigma_XX'))];
                epsilons_idx = [find(strcmp(pnames,'epsilon_MM')) find(strcmp(pnames,'epsilon_XX'))];
                NiAs_idx = find(strcmp(Structures,'NiAs'));
                AntiNiAs_idx = find(strcmp(Structures,'AntiNiAs'));
                bBeO_idx = find(strcmp(Structures,'BetaBeO'));
                
                full_opt_point(sigmas_idx) = full_opt_point(fliplr(sigmas_idx));
                full_opt_point(epsilons_idx) = full_opt_point(fliplr(epsilons_idx));
                for kdx = 1:numel(pnames)
                    ptable.(pnames{kdx}) = full_opt_point(kdx);
                end
                
                % Switch NiAs <-> AntiNiAs
                CoulSR_MM = Minimization_Data{NiAs_idx}.CoulSR_MM;
                LJSR_MM   = Minimization_Data{NiAs_idx}.LJSR_MM;
                CoulSR_XX = Minimization_Data{NiAs_idx}.CoulSR_XX;
                LJSR_XX   = Minimization_Data{NiAs_idx}.LJSR_XX;
                FC_Metal  = Minimization_Data{NiAs_idx}.FC_Metal;
                FC_Halide = Minimization_Data{NiAs_idx}.FC_Halide;
                Minimization_Data{NiAs_idx}.CoulSR_MM = CoulSR_XX;
                Minimization_Data{NiAs_idx}.LJSR_MM   = LJSR_XX;
                Minimization_Data{NiAs_idx}.CoulSR_XX = CoulSR_MM;
                Minimization_Data{NiAs_idx}.LJSR_XX   = LJSR_MM;
                Minimization_Data{NiAs_idx}.FC_Metal  = FC_Halide;
                Minimization_Data{NiAs_idx}.FC_Halide = FC_Metal;
                Minimization_Data{NiAs_idx}.Structure = 'AntiNiAs';
                
                CoulSR_MM = Minimization_Data{AntiNiAs_idx}.CoulSR_MM;
                LJSR_MM   = Minimization_Data{AntiNiAs_idx}.LJSR_MM;
                CoulSR_XX = Minimization_Data{AntiNiAs_idx}.CoulSR_XX;
                LJSR_XX   = Minimization_Data{AntiNiAs_idx}.LJSR_XX;
                FC_Metal  = Minimization_Data{AntiNiAs_idx}.FC_Metal;
                FC_Halide = Minimization_Data{AntiNiAs_idx}.FC_Halide;
                Minimization_Data{AntiNiAs_idx}.CoulSR_MM = CoulSR_XX;
                Minimization_Data{AntiNiAs_idx}.LJSR_MM   = LJSR_XX;
                Minimization_Data{AntiNiAs_idx}.CoulSR_XX = CoulSR_MM;
                Minimization_Data{AntiNiAs_idx}.LJSR_XX   = LJSR_MM;
                Minimization_Data{AntiNiAs_idx}.FC_Metal  = FC_Halide;
                Minimization_Data{AntiNiAs_idx}.FC_Halide = FC_Metal;
                Minimization_Data{AntiNiAs_idx}.Structure = 'NiAs';
                
                Minimization_Data([NiAs_idx AntiNiAs_idx]) = Minimization_Data([AntiNiAs_idx NiAs_idx]);
                
                % Recalculate BetaBeO energy
                Settings = Potential_Scaling(full_data.Settings,ptable);
                Settings.Structure = 'BetaBeO';
                Settings.Diary_Loc = '';
                Settings.Parallel_Bayesopt = false;
                Settings.Parallel_LiX_Minimizer = false;
                Settings.MinMDP.Parallel_Min = false;
                [Settings.home,Settings.project,Settings.computer,Settings.slurm,Settings.BO_Models,...
                    Settings.qsub,Settings.passlog,Settings.pipe,Settings.wsl,~] = find_home;
                Settings.scratch_dir = pwd;
                [~,Settings.gmx,Settings.gmx_loc,Settings.mdrun_opts,Settings.MLModelDir] = MD_Batch_Template(Settings.JobSettings);
                bBeO_Energies = Structure_Minimization(Settings,'Extra_Properties',true);
                Minimization_Data{bBeO_idx} = bBeO_Energies;
                
                % Update data and re-save
                full_data.full_opt_point = full_opt_point;
                full_data.Minimization_Data = Minimization_Data;
                full_data.Pars = ptable;
                save(dat_file,'full_data');
            end
        end
    end
end
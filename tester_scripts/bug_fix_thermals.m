Proj_dir = '/home/scheiber/project/Model_Building'; %/LiI/JC_Model_EX5
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};

for idx = 1:length(Salts)
    Salt = Salts{idx};
    
    all_files = dir(fullfile(Proj_dir,Salt));
    dirFlags = [all_files.isdir];
    subFolders = all_files(dirFlags);
    model_dirs = {subFolders(3:end).name};
    
    for jdx = 1:length(model_dirs)
        Model_Name = model_dirs{jdx};
        
        Complete_calc_dir  = fullfile(Proj_dir,Salt,Model_Name);
        Complete_calc_file = fullfile(Complete_calc_dir,[Salt '_' Model_Name '_fullopt.mat']);
        if isfile(Complete_calc_file)
            
            % Find the final MP folder, if it exists
            files = dir(Complete_calc_dir);
            dirFlags = [files.isdir];
            subFolders = files(dirFlags);
            subFolders = subFolders(3:end);
            if isempty(subFolders)
                continue
            end
            [~,kdx] = sort([subFolders.datenum],'descend');
            subFolders = subFolders(kdx);
            
            % Check for BestPoint_Thermal
            BPT_idx = cellfun(@(x) strcmp(x,'BestPoint_Thermal'),{subFolders.name});
            if any(BPT_idx)
                subFolders = subFolders(~BPT_idx);
                if isempty(subFolders)
                    continue
                else
                    % Delete all other folders
                    disp([num2str(idx) ' ' num2str(jdx) ' ' Salt ' ' Model_Name ' BestPoint_Thermal found'])
%                     for mdx = 1:length(subFolders)
%                         rmdir(fullfile(subFolders(mdx).folder,subFolders(mdx).name),'s')
%                     end
                end
            else
                % If no BestPoint_Thermal folder exists, make one and add
                % last 3 thermal calculations to it, delete the rest
                disp([num2str(idx) ' ' num2str(jdx) ' ' Salt ' ' Model_Name ' no BestPoint_Thermal found'])
                try
                    Bestpoint_thermal_folders = subFolders(1:3);
                    for mdx = 1:length(Bestpoint_thermal_folders)
                        if ~isempty(regexp(Bestpoint_thermal_folders(mdx).name,'_LP','once'))
                            newdir = 'Liq_Properties_at_MP';
                        elseif ~isempty(regexp(Bestpoint_thermal_folders(mdx).name,'_SP','once'))
                            newdir = 'Sol_Properties_at_MP';
                        else
                            newdir = 'Melting_Point';
                        end
                        src = fullfile(Complete_calc_dir,Bestpoint_thermal_folders(mdx).name);
                        dst = fullfile(Complete_calc_dir,'BestPoint_Thermal',newdir);
                        if ~isfolder(fullfile(Complete_calc_dir,'BestPoint_Thermal'))
                            mkdir(fullfile(Complete_calc_dir,'BestPoint_Thermal'))
                        end
                        movefile(src,dst)
                    end

                    Other_thermal_folders = subFolders(4:end);
                    for mdx = 1:length(Other_thermal_folders)
                        rmdir(fullfile(Other_thermal_folders(mdx).folder,Other_thermal_folders(mdx).name),'s')
                    end
                catch
                    for mdx = 1:length(subFolders)
                        rmdir(fullfile(subFolders(mdx).folder,subFolders(mdx).name),'s')
                    end
                    delete(Complete_calc_file)
                end
            end
        end
    end
end
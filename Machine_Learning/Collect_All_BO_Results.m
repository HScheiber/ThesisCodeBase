% Collect_All_BO_Results
proj_dir = '/home/scheiber/project/Model_Building';
destination_folder = '/home/scheiber/project/Model_Building/Completed';

if ~isfolder(destination_folder)
    mkdir(destination_folder);
end

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI'};
for idx = 1:length(Salts)
    Salt = Salts{idx};
    
    % Find model directories
    files = dir(fullfile(proj_dir,Salt));
    diridx = [files.isdir];
    files = files(diridx);
    Models = files(~ismember({files.name},{'.','..'}));
    
    for jdx = 1:length(Models)
        Model_Name = Models(jdx).name;
        Current_Model_dir = fullfile(Models(jdx).folder,Model_Name);

        fullopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_fullopt.mat']);
        fullopt_history = dir([Current_Model_dir filesep 'intermediate_secondary_opt.mat']);
        bayesopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_bayesopt.mat']);
        if ~isempty(fullopt)
            fullopt_history_name = strrep(fullopt.name,'fullopt','fullopt_history');
            
            destfile = fullfile(destination_folder,fullopt.name);
            destdir = dir(destfile);
            sourcefile = fullfile(fullopt.folder,fullopt.name);
            if ~isfile(destfile) || ( datenum(fullopt.date) > datenum(destdir.date) )
                copyfile(sourcefile,destfile);
                disp(['Copied: ' fullopt.name])
            end
            
            destfile = fullfile(destination_folder,fullopt_history_name);
            destdir = dir(destfile);
            if isempty(fullopt_history)
                disp([Salt '_' Model_Name ': No fullopt history found!'])
            elseif ~isfile(destfile)
                sourcefile = fullfile(fullopt_history.folder,fullopt_history.name);
                copyfile(sourcefile,destfile);
                disp(['Copied: ' fullopt_history_name])
            elseif datenum(fullopt_history.date) > datenum(destdir.date)
                sourcefile = fullfile(fullopt_history.folder,fullopt_history.name);
                copyfile(sourcefile,destfile);
                disp(['Copied: ' fullopt_history_name])
            end
            
            destfile = fullfile(destination_folder,bayesopt.name);
            destdir = dir(destfile);
            sourcefile = fullfile(bayesopt.folder,bayesopt.name);
            if ~isfile(destfile) || ( datenum(bayesopt.date) > datenum(destdir.date) )
                copyfile(sourcefile,destfile);
                disp(['Copied: ' bayesopt.name])
            end
            
            % Clean up folder
            del_cmd = ['find "' fullopt.folder  '" -type f \( -name *.stde -o -name *.stdo -o -name *.subm \) -delete'];
            system(del_cmd);
        else
            disp([Salt '_' Model_Name ': Job Not Complete!'])
        end
    end
end
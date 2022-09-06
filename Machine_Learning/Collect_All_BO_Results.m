% Collect_All_BO_Results
function Collect_All_BO_Results(varargin)

if nargin > 0
    parallel = varargin{1};
else
    parallel = false;
end

proj_dir = '/home/scheiber/project/Model_Building';
destination_folder = '/home/scheiber/project/Model_Building/Completed';

if ~isfolder(destination_folder)
    mkdir(destination_folder);
end

Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
for idx = 1:length(Salts)
    Salt = Salts{idx};
    
    % Find model directories
    files = dir(fullfile(proj_dir,Salt));
    diridx = [files.isdir];
    files = files(diridx);
    Models = files(~ismember({files.name},{'.','..'}));
    
    if parallel
        parfor jdx = 1:length(Models)
            warning('off','MATLAB:class:InvalidSuperClass')
            warning('off','MATLAB:load:classNotFound')
            warning('off','MATLAB:load:classError')
            Model_Name = Models(jdx).name;
            Current_Model_dir = fullfile(Models(jdx).folder,Model_Name);
            
            inpfile = fullfile(Current_Model_dir,[Salt '_' Model_Name '.inp']);
            fullopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_fullopt.mat']);
            fullopt_history = dir([Current_Model_dir filesep 'intermediate_secondary_opt.mat']);
            bayesopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_bayesopt.mat']);
            destfile = fullfile(destination_folder,Salt,[Salt '_' Model_Name '_data.mat']);
            
            if ~isempty(fullopt) && ~isfile(destfile)
                try
                    fullopt_dat = load(fullfile(Current_Model_dir,fullopt.name));
                    fullopt_hist_dat = load(fullfile(Current_Model_dir,fullopt_history.name));
                    bayesopt_dat = load(fullfile(Current_Model_dir,bayesopt.name));
                    
                    full_data = fullopt_dat;
                    if isfile(inpfile)
                        input_model = load(inpfile,'-mat');
                        full_data.Settings = input_model.Model;
                    end
                    full_data.secondary_result = fullopt_hist_dat.intermediate_data;
                    full_data.bayesopt_results = bayesopt_dat.results;
                    if ~isfield(full_data,'Finite_T_Data')
                        Settings = struct();
                        Settings.Salt = Salt;
                        full_data.Finite_T_Data = Initialize_Finite_T_Data(Settings);
                    end
                    
                    if ~isfolder(fullfile(destination_folder,Salt))
                        mkdir(fullfile(destination_folder,Salt))
                    end
                    mySave(destfile,full_data);
                    disp(['Copied: ' Salt ' ' Model_Name])
                    
                    % Clean up folder
                    del_cmd = ['find "' fullopt.folder  '" -type f \( -name *.stde -o -name *.stdo -o -name *.subm \) -delete'];
                    system(del_cmd);
                catch me
                    disp(['Error extracting data for model: ' Salt ' ' Model_Name])
                    disp(me.message)
                    disp('Skipping model')
                end
            elseif isfile(destfile)
                %disp([Salt ' ' Model_Name ': Job Already Copied.'])
            else
                disp([Salt ' ' Model_Name ': Job Not Complete!'])
            end
        end
    else
        for jdx = 1:length(Models)
            warning('off','MATLAB:class:InvalidSuperClass')
            warning('off','MATLAB:load:classNotFound')
            warning('off','MATLAB:load:classError')
            Model_Name = Models(jdx).name;
            Current_Model_dir = fullfile(Models(jdx).folder,Model_Name);
            
            inpfile = fullfile(Current_Model_dir,[Salt '_' Model_Name '.inp']);
            fullopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_fullopt.mat']);
            fullopt_history = dir([Current_Model_dir filesep 'intermediate_secondary_opt.mat']);
            bayesopt = dir([Current_Model_dir filesep Salt '_' Model_Name '_bayesopt.mat']);
            destfile = fullfile(destination_folder,Salt,[Salt '_' Model_Name '_data.mat']);
            
            if ~isempty(fullopt) && ~isfile(destfile)
                try
                    fullopt_dat = load(fullfile(Current_Model_dir,fullopt.name));
                    fullopt_hist_dat = load(fullfile(Current_Model_dir,fullopt_history.name));
                    bayesopt_dat = load(fullfile(Current_Model_dir,bayesopt.name));
                    
                    full_data = fullopt_dat;
                    if isfile(inpfile)
                        input_model = load(inpfile,'-mat');
                        full_data.Settings = input_model.Model;
                    end
                    full_data.secondary_result = fullopt_hist_dat.intermediate_data;
                    full_data.bayesopt_results = bayesopt_dat.results;
                    if ~isfield(full_data,'Finite_T_Data')
                        Settings = struct();
                        Settings.Salt = Salt;
                        full_data.Finite_T_Data = Initialize_Finite_T_Data(Settings);
                    end
                    
                    if ~isfolder(fullfile(destination_folder,Salt))
                        mkdir(fullfile(destination_folder,Salt))
                    end
                    mySave(destfile,full_data);
                    disp(['Copied: ' Salt ' ' Model_Name])
                    
                    % Clean up folder
                    del_cmd = ['find "' fullopt.folder  '" -type f \( -name *.stde -o -name *.stdo -o -name *.subm \) -delete'];
                    system(del_cmd);
                catch me
                    disp(['Error extracting data for model: ' Salt ' ' Model_Name])
                    disp(me.message)
                    disp('Skipping model')
                end
            elseif isfile(destfile)
                %disp([Salt ' ' Model_Name ': Job Already Copied.'])
            else
                disp([Salt ' ' Model_Name ': Job Not Complete!'])
            end
        end
    end
end
end
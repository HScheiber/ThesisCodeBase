% Outer script to generate Gromacs data
% Generate_Trajectory_Data
%% Input parameters
Salts = {'LiF' 'NaCl'};
Structures = {'Rocksalt' 'Wurtzite' 'FiveFive'};
Parallel = false; % Run in parallel or serial

analysis_folder_PC = 'C:\Users\Hayden\Documents\Patey_Lab\GTesting';
analysis_folder_unix = '/home/scheiber/project/PURE_SALT_MD';

datadir_PC = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\PURE_SALT_MD';
datadir_unix = '/home/scheiber/ThesisCodeBase/data/PURE_SALT_MD';

if ispc
    Directory = analysis_folder_PC;
    datadir = datadir_PC;
    addpath(Directory);
elseif isunix
    Directory = analysis_folder_unix;
    datadir = datadir_unix;
end

% Crystal types
if strcmp('All',Structures)
    Structures = {'Rocksalt' 'Wurtzite' 'Sphalerite' 'CsCl' 'NiAs' ...
        'BetaBeO' 'FiveFive'};
else
    if ~iscell(Structures)
        Structures = {Structures};
    end
end

if strcmp('All',Salts)
    Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
else
    if ~iscell(Salts)
        Salts = {Salts};
    end
end

N = length(Salts);
M = length(Structures);

if Parallel
    if ispc 
        Parcores = 12;
        if isempty(gcp('nocreate'))
            parpool(Parcores);
        end
    elseif isunix
        if isempty(gcp('nocreate'))
            [~,coretxt]= system('echo $SLURM_CPUS_ON_NODE');
            Parcores = str2double(coretxt);
            if isnan(Parcores)
                Parcores = 6;
            end    
            MyProfile = parcluster('local');
            MyProfile.NumWorkers = Parcores;
            MyProfile.NumThreads = Parcores;
            MyProfile.saveProfile
            while true
                try
                    poolobj = parpool(MyProfile,Parcores);
                    break
                catch
                    Parcores = Parcores-1;
                end
            end
        end
    end
end

for i=1:N
    Salt = Salts{i};
    for j=1:M
        Structure = Structures{j};
        Data = struct;

        % Find models
        Current_Dir = fullfile(Directory,Salt,Structure);
        Models = find_directories(dir(Current_Dir));
        O = length(Models);
        for k =1:O
            Model = Models{k};
            Model_Tag = strrep(strrep(Model,'-','N'),'.','P');
            disp(['Extracting ' Salt ' ' Structure ' ' Model ' Data...']);
            tic;

            % Find Thermostats
            Current_Dir = fullfile(Directory,Salt,Structure,Model);
            Thermostats = find_directories(dir(Current_Dir));
            P = length(Thermostats);
            for l = 1:P
                Thermostat = Thermostats{l};
                Thermostat_Tag = strrep(strrep(Thermostat,'-','N'),'.','P');
                disp(['Current Thermostat: ' Thermostat '...']);
                
                % Find Barostats
                Current_Dir = fullfile(Directory,Salt,Structure,Model,Thermostat);
                Barostats = find_directories(dir(Current_Dir));
                Q = length(Barostats);
                for m = 1:Q
                    Barostat = Barostats{m};
                    Barostat_Tag = strrep(strrep(Barostat,'-','N'),'.','P');
                    disp(['Current Barostat: ' Barostat '...']);

                    % Find MD / simulated annealing
                    Current_Dir = fullfile(Directory,Salt,Structure,Model,Thermostat,Barostat);
                    Annealing = find_directories(dir(Current_Dir));
                    R = length(Annealing);
                    for n = 1:R
                        Anneal = Annealing{n};
                        Anneal_Tag = strrep(strrep(Anneal,'-','N'),'.','P');
                        disp(['Sim. Annealing? ' Anneal '...']);

                        Current_Dir = fullfile(Directory,Salt,Structure,...
                            Model,Thermostat,Barostat,Anneal);
                        
                        Data.(Salt).(Structure).(Model_Tag).(Thermostat_Tag).(Barostat_Tag).(Anneal_Tag) = ...
                            Trajectory_and_RDF(Current_Dir,Parallel);
                    end
                end
            end
            timer = toc;
            disp(['Finished ' Salt ' ' Structure ' ' Model ' Data Extraction. Time Elapsed: ' num2str(timer) ' s.']);
        end
        filename = [Salt '_' Structure '_MD_Data.mat'];
        parallelsave(fullfile(datadir,filename),{'Data'},Data)
    end
end

%#ok<*UNRCH>
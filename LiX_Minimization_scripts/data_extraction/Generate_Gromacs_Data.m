% Outer script to generate Gromacs data
% Generate_Gromacs_Data
%% Input parameters
Salts = {'LiCl'}; % 'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'
Structures = {'Rocksalt' 'Wurtzite' 'FiveFive'}; %Rocksalt, Wurtzite, Sphalerite, CsCl, NiAs, BetaBeO, FiveFive
Watermodel = 'SPC/E'; % For pair potential generator
Parallel = true; % Run in parallel or serial

if ispc
    Directory = 'C:\Users\Hayden\Documents\Patey_Lab\GTesting';
    datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\GROMACS';
elseif isunix
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
    
    if strcmp(Server,'ced') || strcmp(Server,'cdr') % Cedar
        Directory = '/home/scheiber/project/GROMACS_LiHalides';
        datadir = '/home/scheiber/ThesisCodeBase/data/GROMACS';
        
    elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log') % Sockeye
        Directory = '/home/haydensc/scratch/GROMACS_LiHalides';
        datadir = '/home/haydensc/ThesisCodeBase/data/GROMACS';
    elseif strcmp(Server,'pat')
        Directory = '/home/user/project/GROMACS_LiHalides';
        datadir = '/home/user/ThesisCodeBase/data/GROMACS';
    end
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
    
    % Set up matlab parallel features
    PrefCores = feature('numcores');    
    if ~isempty(gcp('nocreate'))
        Cur_Pool = gcp;
        Cur_Workers = Cur_Pool.NumWorkers;
        if Cur_Workers ~= PrefCores
            delete(Cur_Pool);
            ppool = parpool(PrefCores);
        else
            ppool = Cur_Pool;
        end
    else
        ppool = parpool(PrefCores); 
    end
    
    for i=1:N
        Salt = Salts{i};
        for j=1:M
        	Structure = Structures{j};
            if strcmp(Structure,'Pair')
                tic;
                disp(['Extracting ' Salt ' ' Structure ' Data...']);
                FileName = [Salt '_' Structure '_Lattice_Energies.mat'];
                Data = Pair_Potential_Saver(0,1,0.001,Salt,Watermodel);
                parallelsave(fullfile(datadir,FileName),{'Data'},Data)
                Timer = toc;
                disp(['Finished ' Salt ' ' Structure ' Data Extraction. Time Elapsed: ' num2str(Timer) ' s.']);
                continue
            end
            
            % Find models
            Current_Dir = fullfile(Directory,Salt,Structure);
            Models = find_directories(dir(Current_Dir));
            O = length(Models);
            Data_parts = cell(1,O);
            Model_Tag = cell(1,O);
            parfor k =1:O
                Model = Models{k};
                disp(['Extracting ' Salt ' ' Structure ' ' Model ' Data...']);
                tic;
                Model_Tag{k} = strrep(strrep(Model,'-','N'),'.','P');
                Data_parts{k} = Extract_Gromacs_PES(Salt,Structure,Model,Directory);
                timer = toc;
                disp(['Finished ' Salt ' ' Structure ' ' Model ' Data Extraction. Time Elapsed: ' num2str(timer) ' s.']);
            end
            
            % Merge data parts
            Data = [];
            for k =1:O
                if ~isempty(Data_parts{k})
                    Data.(Salt).(Structure).(Model_Tag{k}) = Data_parts{k};
                end
            end
            filename = [Salt '_' Structure '_Lattice_Energies.mat'];
            parallelsave(fullfile(datadir,filename),{'Data'},Data)
        end
    end
else
    for i=1:N
        Salt = Salts{i};
        for j=1:M
        	Structure = Structures{j};
            if strcmp(Structure,'Pair')
                tic;
                disp(['Extracting ' Salt ' ' Structure ' Data...']);
                FileName = [Salt '_' Structure '_Lattice_Energies.mat'];
                Data = Pair_Potential_Saver(0,1,0.001,Salt,Watermodel);
                parallelsave(fullfile(datadir,FileName),{'Data'},Data)
                Timer = toc;
                disp(['Finished ' Salt ' ' Structure ' Data Extraction. Time Elapsed: ' num2str(Timer) ' s.']);
                continue
            else
                Data = struct;
            end
            
            % Find models
            Current_Dir = fullfile(Directory,Salt,Structure);
            Models = find_directories(dir(Current_Dir));
            O = length(Models);
            Data_parts = cell(1,O);
            Model_Tag = cell(1,O);
            for k = 1:O
                Model = Models{k};
                disp(['Extracting ' Salt ' ' Structure ' ' Model ' Data...']);
                tic;
                Model_Tag{k} = strrep(strrep(Model,'-','N'),'.','P');
                Data_parts{k} = Extract_Gromacs_PES(Salt,Structure,Model,Directory);
                timer = toc;
                disp(['Finished ' Salt ' ' Structure ' ' Model ' Data Extraction. Time Elapsed: ' num2str(timer) ' s.']);
            end
            % Merge data parts
            Data = [];
            for k =1:O
                if ~isempty(Data_parts{k})
                    Data.(Salt).(Structure).(Model_Tag{k}) = Data_parts{k};
                end
            end
            filename = [Salt '_' Structure '_Lattice_Energies.mat'];
            parallelsave(fullfile(datadir,filename),{'Data'},Data)
        end
    end
end

%#ok<*UNRCH>
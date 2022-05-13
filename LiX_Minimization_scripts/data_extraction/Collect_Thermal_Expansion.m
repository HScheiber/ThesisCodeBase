% Collect_Thermal_Expansion
Salts = {'LiF' 'LiCl' 'LiBr' 'LiI' 'NaCl'};
Structures = {'Rocksalt'};
Begin_frame = 5000; % In ps
nm_Ang = 10;

if ispc % for testing
    gmx = 'wsl source ~/.bashrc; gmx_d';
    projdir = 'C:\Users\Hayden\Documents\Patey_Lab\GTesting';
    homedir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase'; % PC
    datadir = 'C:\Users\Hayden\Documents\Patey_Lab\ThesisCodeBase\data\PURE_SALT_MD';
    sys = @(inp) system(inp); 
elseif isunix
    [~,Servertxt] = system('hostname -s | cut -c 1-3');
    Server = strtrim(Servertxt);
    if strcmpi(Server,'ced') || strcmpi(Server,'cdr') || strcmpi(Server,'sea')
        gmx = 'gmx_d';
        projdir = '/home/scheiber/project/PURE_SALT_MD';
        homedir = '/home/scheiber/ThesisCodeBase'; % Cedar/Graham/orcinus
        datadir = '/home/scheiber/ThesisCodeBase/data/PURE_SALT_MD';
        sys = @(inp) system(inp);
    elseif ~isempty(regexp(Server,'se[0-9]','ONCE')) || strcmpi(Server,'log')
        gmx = 'gmx_mpi_d';
        projdir = '/home/haydensc/scratch/PURE_SALT_MD';
        homedir = '/home/haydensc/ThesisCodeBase'; % Sockeye
        datadir = '/home/haydensc/ThesisCodeBase/data/PURE_SALT_MD';
        sys = @(inp) system(inp);
    elseif strcmpi(Server,'bel')
        gmx = 'gmx_d';
        projdir = '/home/scheiber/project/PURE_SALT_MD';
        homedir = '/home/scheiber/ThesisCodeBase'; % Beluga
        datadir = '/home/scheiber/ThesisCodeBase/data/PURE_SALT_MD';
        sys = @(inp) system_def(inp); % Needed to circumvent error
    elseif strcmpi(Server,'pat')
        gmx = 'gmx_d';
        projdir = '/media/user/LaCie_2TB/project/PURE_SALT_MD';
        homedir = '/home/user/Documents/Thesis_work/ThesisCodeBase'; % PC
        datadir = '/home/user/ThesisCodeBase/data/PURE_SALT_MD';
        sys = @(inp) system(inp); 
    end
else
    error('Unknown machine type.')
end

N_Salt = length(Salts);
N_Struc = length(Structures);
Data = struct;
for idx1=1:N_Salt % Loop through Salts
    Salt = Salts{idx1};
    [Metal,Halide] = Separate_Metal_Halide(Salt);
    for idx2=1:N_Struc % Loop through Structures
        Structure = Structures{idx2};
        

        % Find models
        Current_Dir = fullfile(projdir,Salt,Structure);
        Models = find_directories(dir(Current_Dir));
        N_Model = length(Models);
        for idx3 =1:N_Model % Loop through models
            Model = Models{idx3};
            Model_Tag = strrep(strrep(Model,'-','N'),'.','P');
            disp(['Extracting ' Salt ' ' Structure ' ' Model ' Data...']);
            tic;

            % Find Thermostats
            Current_Dir = fullfile(projdir,Salt,Structure,Model);
            Thermostats = find_directories(dir(Current_Dir));
            N_Therm = length(Thermostats);
            for idx4 = 1:N_Therm
                Thermostat = Thermostats{idx4};
                Thermostat_Tag = strrep(strrep(Thermostat,'-','N'),'.','P');
                disp(['Current Thermostat: ' Thermostat '...']);
                
                % Find Barostats
                Current_Dir = fullfile(projdir,Salt,Structure,Model,Thermostat);
                Barostats = find_directories(dir(Current_Dir));
                N_Baro = length(Barostats);
                for idx5 = 1:N_Baro
                    Barostat = Barostats{idx5};
                    Barostat_Tag = strrep(strrep(Barostat,'-','N'),'.','P');
                    disp(['Current Barostat: ' Barostat '...']);

                    % Find MD / simulated annealing
                    Current_Dir = fullfile(projdir,Salt,Structure,Model,Thermostat,Barostat);
                    Annealing = find_directories(dir(Current_Dir));
                    N_MD = length(Annealing);
                    for idx6 = 1:N_MD
                        Anneal = Annealing{idx6};
                        Anneal_Tag = strrep(strrep(Anneal,'-','N'),'.','P');
                        disp(['Sim. Annealing? ' Anneal '...']);

                        Current_Dir = fullfile(projdir,Salt,Structure,...
                            Model,Thermostat,Barostat,Anneal);
                        
                        % Find the edr (energy) file
                        Files = dir([Current_Dir filesep '*.edr']);
                        
                        % Get filebase from trr name
                        [~,Filebase,~] = fileparts(Files.name);
                        
                        % Binary energy file
                        EDR_File = [Current_Dir filesep Filebase '.edr'];
                        
                        % XVG energy file for data analysis
                        XVG_Energy_File = [Current_Dir filesep Filebase '_ThermExp.xvg']; 
                        
                        % Contains info about number of particles
                        MAT_File = [Current_Dir filesep Filebase '.mat']; 

                        % Generate total average box size
                        Energy_Command = ['echo 9 10 11 0 | ' gmx ' energy -f ' EDR_File ' -o ' XVG_Energy_File ' -b ' num2str(Begin_frame)];
                        [~,Eng_Output] = sys(Energy_Command);
                        
                        Logend = regexp(Eng_Output,'Energy +Average +Err.+?\n\n','match','once');
                        BoxDims_cell = textscan(Logend,'%*s %f %*f %*f %*f %*s','HeaderLines',2,'Delimiter',' ',...
                            'MultipleDelimsAsOne',true);
                        BoxDims = BoxDims_cell{1};
                        
                        MatData = load(MAT_File);
                        SuperCellSize = (MatData.N_total/MatData.N_Cell)^(1/3);
                        
                        a = (BoxDims(1)/SuperCellSize)*nm_Ang;
                        b = (BoxDims(2)/SuperCellSize)*nm_Ang;
                        c = (BoxDims(3)/SuperCellSize)*nm_Ang;
                        
                        
                        Data.(Salt).(Structure).(Model_Tag) = [a b c];
                    end
                end
            end
            timer = toc;
            disp(['Finished ' Salt ' ' Structure ' ' Model ' Data Extraction. Time Elapsed: ' num2str(timer) ' s.']);
        end
    end
end
filename = 'T298K_Thermal_Expansion_Data.mat';
parallelsave(fullfile(datadir,filename),{'Data'},Data);


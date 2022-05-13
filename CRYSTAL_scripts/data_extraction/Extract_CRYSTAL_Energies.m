% SCRIPT FOR EXTRACTING SINGLE POINT ENERGIES FROM CRYSTAL17 OUTPUT
% (*.log) FILES.

% WARNING: This script does not accept job types with acronyms longer than
% 10 characters

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'Li', 'F', etc.

function Data = Extract_CRYSTAL_Energies(Label,Labend)
%% Level 0 (analysis folder)

% Find directories in current folder, except . and ..
Directories = find_directories(dir);

% Find folders that match the label
matches0 = regexp(Directories,['^' Label '_' Labend]);

Data = struct();
for x=1:length(Directories)
    %% Level 1 (first level such as "LiF")

    % Skip folders that do not match label
    if isempty(matches0{x})
        continue
    end
    
    cur_dir = Directories{x};
    
    % Find directories in current folder, except . and ..
    Folders = find_directories(dir(cur_dir));

    % Find folders that match the label
    matches = regexp(Folders,['^' Label '_']);

    % Ensure structure names are good by removing dashes
    structnames = regexprep(Folders,'-','_');

    %% Level 2 (Theory Type folder)
    R = length(Folders);
    for i=1:R

        % Skip folders that do not match label
        if isempty(matches{i})
            continue
        end

        % Enter directory i
        cur_dir2 = [cur_dir filesep Folders{i}];

        % Find directories in current folder, except . and ..
        subFolders = find_directories(dir(cur_dir2));

        N = length(subFolders);

        % Pre-Create data output for this theory type
        Data.(structnames{i}) = cell(N,3);
        
        % Find type of theory
        Theory_type = regexp(Folders{i},['(?<=' Label '_)' '.+'],'match');
        Theory_type = Theory_type{1};

        %% Level 3 (Basis set folder)
        % Loop through directories
        for j=1:N

            Basis_set = textscan(subFolders{j},'BASS_%s');
            Data.(structnames{i}){j,1} = regexprep(Basis_set{1}{1},'%2A','*');

            % Move to subfolder j
            cur_dir3 = [cur_dir2 filesep subFolders{j}];

            logfile = dir([cur_dir3 filesep '*' Basis_set{1}{1} '.out']);

            % If no log file exists, remove from list and continue
            if isempty(logfile)
                continue
            end

            % Open the file
            fid = fopen([cur_dir3 filesep logfile.name],'rt');
            logtext = fread(fid);
            fclose(fid);
            logtext = char(logtext.');

            % Find energy outputs
            energy = regexp(logtext,...
                '(?<=SCF ENDED - CONVERGENCE ON (ENERGY|TESTER) +E\(AU\) )(-|\+|\.|[0-9]|E|)+','match');

            % Find CPU time ouputs in seconds
            if ~isempty(energy)
                energy = str2double(energy{1});
                CPU_time = regexp(logtext,...
                    '(?<=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT END         TELAPSE +([0-9]|E|\.)+ +TCPU +)([0-9]|E|\.)+','match');
                Tot_Time = str2double(CPU_time{1});
            else
                continue
            end

            Data_Array = cell(1,2);

            % Input type of theory
            Data_Array{1,1} = Theory_type;

            % Input energy
            Data_Array{1,2} = energy;

            % Add data to output data structure
            Data.(structnames{i}){j,2} = Data_Array;
            Data.(structnames{i}){j,3} = Tot_Time;
        end

        % Remove empty data (counting in reverse)
        for j=length(subFolders):-1:1
            if isempty(Data.(structnames{i}){j,2})
                Data.(structnames{i})(j,:) = [];
            end
        end

        if (i/R) < 0.9999
            fprintf('%s',[num2str(round(i*100/R)) '%... '])
        end
    end
    fprintf('%s',['100%.' newline])
end
end


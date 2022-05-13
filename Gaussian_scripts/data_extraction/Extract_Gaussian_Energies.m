% SCRIPT FOR EXTRACTING SINGLE POINT ENERGIES FROM GAUSSIAN OUTPUT
% (*.log) FILES.

% WARNING: This script does not accept job types with acronyms longer than
% 10 characters

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'Li', 'F', etc.

function Data = Extract_Gaussian_Energies(Label)
    %% Level 0 (analysis folder)

    % Find directories in current folder, except . and ..
    Directories = find_directories(dir);

    % Find folders that match the label
    matches0 = regexp(Directories,['^' Label '_']);

    Data = struct();

for x=1:length(Directories)
    %% Level 1 (first level such as "LiF")

    % Skip folders that do not match label
    if isempty(matches0{x})
        continue
    else
        cd(Directories{x})
    end
    
    % Find directories in current folder, except . and ..
    Folders = find_directories(dir);

    % Find folders that match the label
    matches = regexp(Folders,['^' Label '_']);

    % Ensure structure names are good by removing dashes
    structnames = regexprep(Folders,'-','_');

    %% Level 2 (Theory Type folder)
    for i=1:length(Folders)

        % Skip folders that do not match label
        if isempty(matches{i})
            continue
        end

        % Enter directory i
        cd(Folders{i})   

        % Find directories in current folder, except . and ..
        subFolders = find_directories(dir);

        N = length(subFolders);

        % Pre-Create data output for this theory type
        Data.(structnames{i}) = cell(N,3);

        %% Level 3 (Basis set folder)
        % Loop through directories
        for j=1:N

            Basis_set = textscan(subFolders{j},'BASS_%s');
            Data.(structnames{i}){j,1} = regexprep(Basis_set{1}{1},'%2A','*');

            % Move to subfolder j
            cd(subFolders{j})

            logfile = dir('*.log');

            % If no log file exists, remove from list and continue
            if isempty(logfile)
                cd('..')
                continue
            end

            % Open the file
            fid = fopen(logfile.name,'rt');
            logtext = fread(fid);
            fclose(fid);
            logtext = char(logtext.');

            % Find energy outputs
            energies = regexp(logtext,...
                '\\(.){1,10}?([0-9]|-|,|\.| |\n|e){9,}?(?=\\)','match');

            % Find CPU time ouputs
            if ~isempty(energies)
                CPU_time = regexp(logtext,...
                    '(?<=Elapsed time:)( )*[0-9]+(\.)*[0-9]*( )*days( )*[0-9]+(\.)*[0-9]*( )*hours( )*[0-9]+(\.)*[0-9]*( )*minutes( )*[0-9]+(\.)*[0-9]*( )*seconds.','match');
                NFtime = textscan(CPU_time{1},'%f %*s','delimiter',' ','MultipleDelimsAsOne',1);

                % Convert total time to seconds
                Tot_Time = NFtime{1}(1)*86400 + NFtime{1}(2)*3600 + ...
                    NFtime{1}(3)*60 + NFtime{1}(4);
            else
                Tot_Time = [];
            end

            Data_Array = cell(length(energies),2);
            for k=1:length(energies)

                % Remove newlines and spaces
                energies{k} = regexprep(energies{k},'(\n| )','');

                % Find type of theory
                Theory_type = regexp(energies{k},'\\(.+)=','match');
                Data_Array{k,1} = regexprep(Theory_type{1},'(\\|=)','');

                % Remove the header
                energies{k} = regexprep(energies{k},'\\(.+)=','');

                % Split string at commas and convert to double
                Data_Array{k,2} = str2double(energies{k});
            end

            % Add data to output data structure
            Data.(structnames{i}){j,2} = Data_Array;
            Data.(structnames{i}){j,3} = Tot_Time;
            clearvars Data_Array Tot_Time

            % Return to level 2
            cd('..')
        end

        % Remove empty data (counting in reverse)
        for j=length(subFolders):-1:1
            if isempty(Data.(structnames{i}){j,2})
                Data.(structnames{i})(j,:) = [];
            end
        end

        % Return to level 1
        cd('..')
    end
    % Return to level 0
    cd('..')
end
end


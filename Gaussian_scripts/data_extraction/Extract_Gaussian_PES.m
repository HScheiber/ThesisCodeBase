% SCRIPT FOR EXTRACTING POTENTIAL ENERGY SURFACES FROM GAUSSIAN OUTPUT
% (*.log) FILES.
% WARNING: This script does not accept job types with acronyms longer than
% 10 characters

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'LiF', 'FF', etc.

function Data = Extract_Gaussian_PES(Label)
% Independent variables and parameters
x1 = 1:0.1:10;
x2 = 1:0.1:4;

DFT_methods = {'B3LYP','BPBE','B2PLYP','B2PLYPD','B2PLYPD3','B97D','B97D3','LSDA','PBE','wB97XD','BPBE_D3BJ','B3LYP_D3BJ','DSDPBEP86'};
Dispersion = {'B2PLYPD','B2PLYPD3','B97D','B97D3','wB97XD','BPBE_D3BJ','B3LYP_D3BJ','DSDPBEP86_D3BJ'};

%% Level 1 (outer)

% Find directories in current folder, except . and ..
Folders = find_directories(dir);

% Find folders that match the label
matches = regexp(Folders,['^' Label '_']);

% Ensure structure names are good by removing dashes
structnames = regexprep(Folders,'-','_');

Data = struct();

%% Level 2
for i=1:length(Folders)
    
    % Skip folders that do not match label
    if isempty(matches{i})
        continue
    end
    
    % Enter directory i
    cd(Folders{i})   
    
    % Find directories in current folder, except . and ..
    subFolders = find_directories(dir);

    % Pre-Create data output for this theory type
    Data.(structnames{i}) = cell(length(subFolders),3);
    
    %% Level 3
    % Loop through directories
    for j=1:length(subFolders)

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
            '\\(.){1,10}?([0-9]|-|,|\.| |\n|e){300,}?(?=\\)','match');
        
        % Find CPU time ouputs
        if ~isempty(energies)
            CPU_time = regexp(logtext,...
                '(?<=cpu time:)( )*[0-9]+(\.)*[0-9]*( )*days( )*[0-9]+(\.)*[0-9]*( )*hours( )*[0-9]+(\.)*[0-9]*( )*minutes( )*[0-9]+(\.)*[0-9]*( )*seconds.','match');
            NFtime = textscan(CPU_time{1},'%f %*s','delimiter',' ','MultipleDelimsAsOne',1);
            
            % Convert total time to seconds
            Tot_Time = NFtime{1}(1)*86400 + NFtime{1}(2)*3600 + ...
                NFtime{1}(3)*60 + NFtime{1}(4);
        else
            Tot_Time = [];
        end

        Data_Array = cell(length(energies),3);
        for k=1:length(energies)
            
            % Dependent variable
            if strcmp(Folders{i},'LiF_CCSD_T')
                Data_Array{k,2} = x2;
            else
                Data_Array{k,2} = x1;
            end
            
            % Remove newlines and spaces
            energies{k} = regexprep(energies{k},'(\n| )','');
            
            % Find type of theory
            Theory_type = regexp(energies{k},'\\(.+)=','match');
            Data_Array{k,1} = regexprep(Theory_type{1},'(\\|=)','');
            
            % Remove the header
            energies{k} = regexprep(energies{k},'\\(.+)=','');
            
            % Split string at commas and convert to double
            Data_Array{k,3} = str2double( strsplit(energies{k},',') );
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


%% Add DFT/dispersion markers
fields = fieldnames(Data);
N = length(fields);

for i=1:length(DFT_methods)
    DFT_methods{i} = [Label '_' DFT_methods{i}];
end

for i=1:length(Dispersion)
    Dispersion{i} = [Label '_' Dispersion{i}];
end

for i=1:N
    M = size(Data.(fields{i}),1);
    
    if ismember(fields{i},DFT_methods) && ismember(fields{i},Dispersion)
        for j=1:M
            Data.(fields{i}){j,4} = 3;
        end
    elseif ismember(fields{i},DFT_methods)
        for j=1:M
            Data.(fields{i}){j,4} = 2;
        end
    else
        for j=1:M
            Data.(fields{i}){j,4} = 1;
        end
    end
end


end


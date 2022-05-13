% SCRIPT FOR EXTRACTING SINGLE POINT ENERGIES FROM GAUSSIAN OUTPUT
% (*.log) FILES.

% WARNING: This script does not accept job types with acronyms longer than
% 10 characters

% the input 'label' should be the label at the beginning of each directory
% of interest such as 'Li', 'F', etc.

function Data = Extract_Gaussian_DFT_PES(Label)
%% DFT method types
DFT_methods = {'LSDA','LC_LSDA','LC-PBEPBE','PBEPBE','DSDPBEP86','B2GPPLYP','HSEH1PBE','wB97X','wB97','B3LYP',...
    'CAM_B3LYP','LC_wHPBE','B2PLYP','PW91PW91','mPW1PW91','PW1PW'};
Range_cor = {'LC_LSDA','LC_PBEPBE','wB97X','wB97','CAM_B3LYP','LC_wHPBE'};
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
    
    N = length(subFolders);

    % Pre-Create data output for this theory type
    Data.(structnames{i}) = cell(N,4);
    
    %% Level 3
    % Loop through directories
    for j=1:N

        Basis_set = textscan(subFolders{j},'BASS_%s');
        Data.(structnames{i}){j,1} = regexprep(Basis_set{1}{1},'%2A','*');
        
        % Move to subfolder j
        cd(subFolders{j})

        % Find directories in current folder, except . and ..
        subsubFolders = find_directories(dir);
        
        M = length(subsubFolders);

        %% Level 4
        % Loop through directories
        a = zeros(1,M);
        for k=1:M
            cell_param = textscan(subsubFolders{k},'APAR_%s');
            a(k) = str2double(cell_param{1}{1});
            
            % Move to subsubfolder k
            cd(subsubFolders{k})
            
            logfile = dir('*.out');
            
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
                '\\((?!D\s*i\s*p\s*o\s*l\s*e).){1,10}?=([0-9]|-|,|\.| |\n|e){9,}?(?=\\)','match');
            
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
            
            for l=1:length(energies)

                % Remove newlines and spaces
                energies{l} = regexprep(energies{l},'(\n| )','');

                % Find type of theory
                Theory_type = regexp(energies{l},'\\(.+)=','match');
                Data_Array{l,1} = regexprep(Theory_type{1},'(\\|=)','');

                % Remove the header
                energies{l} = regexprep(energies{l},'\\(.+)=','');

                % Convert to double
                Data_Array{l,2} = str2double(energies{l});
            end
            
            if isempty(Data.(structnames{i}){j,2})
                if ~isempty(Data_Array)
                    Data.(structnames{i}){j,2} = cell(size(Data_Array,1),4);
                    Data.(structnames{i}){j,2}(:,1) = Data_Array(:,1);
                    Data.(structnames{i}){j,2}(:,2) = {a(k)};
                    Data.(structnames{i}){j,2}(:,3) = Data_Array(:,2);
                    Data.(structnames{i}){j,2}(:,4) = {Tot_Time};
                    Data.(structnames{i}){j,3} = Tot_Time;
                else
                    Data.(structnames{i}){j,3} = 0;
                end
            else
                if ~isempty(Data_Array)
                    for x=1:size(Data_Array,1)
                        Data.(structnames{i}){j,2}{x,2}(end+1) = a(k);
                        Data.(structnames{i}){j,2}{x,3}(end+1) = Data_Array{x,2};
                        Data.(structnames{i}){j,2}{x,4}(end+1) = Tot_Time;
                    end
                    Data.(structnames{i}){j,3} = Data.(structnames{i}){j,3} + Tot_Time;
                end
            end
            clearvars Data_Array Tot_Time
            cd('..')
        end
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

for i=1:length(Range_cor)
    Range_cor{i} = [Label '_' Range_cor{i}];
end


for i=1:N
    M = size(Data.(fields{i}),1);
    if ismember(fields{i},Range_cor)
        for j=1:M
            Data.(fields{i}){j,4} = 4; % Range corrected DFT
        end
    elseif ismember(fields{i},Dispersion)
        for j=1:M
            Data.(fields{i}){j,4} = 3; % Dispersion corrected DFT
        end
    elseif ismember(fields{i},DFT_methods)
        for j=1:M
            Data.(fields{i}){j,4} = 2; % Uncorrected DFT
        end
    else
        for j=1:M
            Data.(fields{i}){j,4} = 1; % HF or post HF
        end
    end
end

end


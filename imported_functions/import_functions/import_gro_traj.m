%% import_gro_traj.m
% * This function imports a .gro trajectory
function [XYZ_labels,Traj,Box_dim] = import_gro_traj(filename)

% Get the number of atoms and frames
fileID = fopen(filename,'r');
N_Atoms = cell2mat(textscan(fileID,'%f',1,'delimiter','\n','HeaderLines',1));
frewind(fileID) % Reset to start of file
Lines = textscan(fileID,'%s',inf,'delimiter','\n'); % Number of lines (does not count empty lines)
N_Frames = length(Lines{1})/(N_Atoms+3); % number of time points

% Initialize arrays
Box_dim = zeros(N_Frames,9);
XYZ_labels = cell(N_Atoms,1);
Traj = zeros(N_Atoms,3,N_Frames);
%Velo = zeros(N_Atoms,3,N_Frames);
for t = 1:N_Frames
    
    frewind(fileID)
    StartRow = 2*t + (N_Atoms+1)*(t-1);
    
    % Scan the block
    formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%*[^\n\r]';
    dataArray = textscan(fileID, formatSpec, N_Atoms, 'Delimiter','',...
        'WhiteSpace','','EmptyValue',NaN,'HeaderLines',StartRow,'ReturnOnError',false);
    
    % Assign position coordinates
    Traj(:,1,t) = dataArray{5}; % x
    Traj(:,2,t) = dataArray{6}; % y
    Traj(:,3,t) = dataArray{7}; % z

    % Assign velocity coordinates
%     Velo(:,1,t) = dataArray{8}; % x
%     Velo(:,2,t) = dataArray{9}; % y
%     Velo(:,3,t) = dataArray{10}; % z
    
    % Assign Box dimensions
    Box_cell = textscan(fileID, '%f %f %f %f %f %f %f %f %f',1,'delimiter', '\n','HeaderLines', 1);
    Box_cell(cellfun(@isempty,Box_cell)) = {0};
    Box_dim(t,:) = cell2mat(Box_cell);
    
    % Assign Box labels
    if t == 1
        XYZ_labels = strtrim(dataArray{3});
    end
end
fclose(fileID);
end


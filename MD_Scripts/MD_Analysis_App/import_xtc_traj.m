%% import_gro_traj.m
% * This function imports an xtc trajectory and initial gro structure and
% ouputs the coordinates, labels, and box vectors (units of length = nm)
function [XYZ_labels,Traj,Box_a,Box_b,Box_c] = import_xtc_traj(trajfile,grofile)

% Import trajectory with compiled c++ gromacsMex 
[Traj,~,~,~,Box] = gromacsMex(1,0,trajfile);

% Import atom info from grofile
N_atoms = size(Traj,1);
fileID = fopen(grofile,'r');
formatSpec = '%5s%5s%5s%5.0f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f%*[^\n\r]';
dataArray = textscan(fileID, formatSpec,N_atoms,'Delimiter','',...
    'WhiteSpace','','EmptyValue',NaN,'HeaderLines',2,'ReturnOnError',false);
XYZ_labels = strtrim(dataArray{3});

% Assign Box dimensions
Box_cell = textscan(fileID, '%f %f %f %f %f %f %f %f %f',1,'delimiter', '\n','HeaderLines', 1);
Box_cell(cellfun(@isempty,Box_cell)) = {0};
Box_dim = cell2mat(Box_cell);
a_hat = Box_dim([1 4 5])./norm(Box_dim([1 4 5]));
b_hat = Box_dim([6 2 7])./norm(Box_dim([6 2 7]));
c_hat = Box_dim([8 9 3])./norm(Box_dim([8 9 3]));

Box_a = Box(:,1).*a_hat;
Box_b = Box(:,2).*b_hat;
Box_c = Box(:,3).*c_hat;

fclose(fileID);
end


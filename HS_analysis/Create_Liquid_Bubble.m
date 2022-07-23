function Geometry = Create_Liquid_Bubble(Settings)

% IP = load(fullfile(WorkDir,'TempJobInfo.mat'));
% Settings = IP.Settings;

warning('off','MATLAB:UndefinedFunction')
Ref_Density = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
warning('on','MATLAB:UndefinedFunction')

pos_filename = 'positions.dat';
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
Ref_M = fullfile(Setings.home,'templates','GRO_Templates',[Metal '_Box.gro']);
Ref_X = fullfile(Setings.home,'templates','GRO_Templates',[Halide '_Box.gro']);
nmol_cluster = Settings.Cluster_N/2;
nmol = Settings.N_atoms/2;

% Calculate the radius needed for the given density
Vol_cluster = nmol_cluster/Ref_Density; % Volume in cubic nm
R_cluster = ( (3/(4*pi))*Vol_cluster )^(1/3); % Radius in nm
Inner_Box_Length = 2*(R_cluster);

Vol_total = nmol/Ref_Density; % Volume in cubic nm
Full_Box_length = Vol_total^(1/3); % Side length of box

Box_Center = Full_Box_length.*[1/2 1/2 1/2]; % Box center in nm

%% How many atoms should fit into a box with the same side length as the sphere?
Box_Vol = Inner_Box_Length^3; % nm^3 volume
nmol_box = ceil(Ref_Density*Box_Vol);

% Create a positions.dat file
fid = fopen(pos_filename,'wt');
fprintf(fid, '# Box\n');
for i = 1:nmol_box
    fprintf(fid, num2str(Box_Center));
    fprintf(fid,'\n');
end
fclose(fid);


%% Add atoms to the file
% Add metal
tmp_file = 'tmp.gro';
L = num2str(Full_Box_length);
cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_M) ' -ip ' windows2unix(pos_filename)...
    ' -o ' windows2unix(tmp_file) ' -nmol ' num2str(nmol_box) ' -box ' L ' ' L ' ' L ...
    ' -dr ' num2str(R_cluster)];
system(cmd,'-echo');

% Add Halide
MX_out_File = 'full_box.gro';
cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(tmp_file) ...
    ' -o ' windows2unix(MX_out_File) ' -nmol ' num2str(nmol_box) ' -ip ' windows2unix(pos_filename) ...
    ' -dr ' num2str(R_cluster) ' -try 200'];
system(cmd,'-echo');


%% Carve box into a sphere
[atom,data] = import_atom_gro(MX_out_File);
delete(pos_filename)
delete(tmp_file)
delete(MX_out_File)

% split atom into metal and halide data
m_idx = 1;
x_idx = 1;
for idx = 1:data.nAtoms
    if strcmp(atom(idx).resname{1},Metal)
        atom_M(m_idx) = atom(idx);
        m_idx = m_idx+1;
    else
        atom_X(x_idx) = atom(idx);
        x_idx = x_idx+1;
    end
end


xyz_coords_M = [atom_M.x; atom_M.y; atom_M.z]';
xyz_coords_X = [atom_X.x; atom_X.y; atom_X.z]';
R0_M = vecnorm(xyz_coords_M-10.*Box_Center,2,2); % Distances in Angstrom from box center
R0_X = vecnorm(xyz_coords_X-10.*Box_Center,2,2); % Distances in Angstrom from box center

[~,idb_M] = sort(R0_M);
[~,idb_X] = sort(R0_X);
idx_Sorted_M = idb_M(1:nmol_cluster);
idx_Sorted_X = idb_X(1:nmol_cluster);
    
i = 1;
for idx = 1:length(idx_Sorted_M)
    idq = idx_Sorted_M(idx);
    atom_fin(i) = atom_M(idq);
    atom_fin(i).index = i;
    atom_fin(i).molid = 0;
    i = i+1;
end
for idx = 1:length(idx_Sorted_X)
    idq = idx_Sorted_X(idx);
    atom_fin(i) = atom_X(idq);
    atom_fin(i).index = i;
    atom_fin(i).molid = 0;
    i = i+1;
end

% Save
write_atom_gro(atom_fin,data.Box_dim,Settings.Coordinate_File)

Geometry.a = Full_Box_length*10;
Geometry.b = Full_Box_length*10;
Geometry.c = Full_Box_length*10;
Geometry.alpha = 90;
Geometry.beta  = 90;
Geometry.gamma = 90;
Geometry.Transform = eye(3);
Geometry.FC_Metal = [];
Geometry.FC_Halide = [];
Geometry.N = Settings.N_Atoms;
Geometry.Label = 'L';

end
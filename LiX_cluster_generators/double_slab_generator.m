Salt = 'FeO';
Structure = 'Rocksalt';
Theory = 'TMTPSS-rVV10L';
Conventional = true;
DFT = Load_Best_DFT_Data('Theory',Theory,'Conventional',Conventional);
N_Cluster = 1026; % Number of atoms in each slab. must be an even number to ensure charge neutrality
Major_axis = [1 1 1]; % The major axis of the cylinder in terms of [a b c]
Minor_axis = [1 0 -1]*eul2rotm([deg2rad(0) 0 0]); % The minor axis defines the direction of the first point (A) of the Prism. Must be orthogonal to Major axis
N_points = 4; % number of sides on the prism (6 = hexagon)
Cluster_Center = [1/8 1/8 1/8]; % Location of cluster center in units of a, b, and c vectors
Prism_height = 1; % Maximum cylinder height in units of the chosen axis
Cutoff = 50; % Angstroms, maximum search distance
slab_distance = 45; % Slab separation in Angstrom

% Angles to rotate the final coordinate system by
% Leave as Rot = {} for no rotation
% Rotations are performed left to right.
% Give axis of rotation followed by the angle in radians
Rot = {[0 1 0 deg2rad(45)] [1 0 0 deg2rad(-35.264389682754661)] [1 0 0 deg2rad(0)]}; 



Settings = Initialize_MD_Settings;
Settings.Salt = Salt;
Settings.Structure = Structure;
Settings.Use_Conv_cell = Conventional;
% Crystal = Default_Crystal(Settings);
% Crystal.a = DFT.(Salt).(Structure).a;
% Crystal.b = DFT.(Salt).(Structure).b;
% Crystal.c = DFT.(Salt).(Structure).c;
% Ratios 1:3:1:8 (Al:Si:K:O)
Crystal = Default_Crystal(Settings);
Crystal.a = 6.35800;
Crystal.b = 6.35800;
Crystal.c = 6.35800;
[Metal,Halide] = Separate_Metal_Halide(Salt);
N = Crystal.N;

Frac_Coordinates = [Crystal.FC_Metal; Crystal.FC_Halide];

Identities_UC = [repmat(string(Metal),N/2,1);
    repmat(string(Halide),N/2,1)];

%% Generate transformation matrix from crystal basis to cartesian basis
abc = [Crystal.a 0 0; 0 Crystal.b 0; 0 0 Crystal.c];
Transform_Matrix = abc*Crystal.Transform;
a_vec = Transform_Matrix(1,:);
b_vec = Transform_Matrix(2,:);
c_vec = Transform_Matrix(3,:);

Prism_center = Cluster_Center*Transform_Matrix;
Ref_axis_major = Major_axis*Transform_Matrix./norm(Major_axis*Transform_Matrix);
Ref_axis_minor = (Minor_axis*Transform_Matrix./norm(Minor_axis*Transform_Matrix))*axang2rotm([Ref_axis_major 2*pi/8]);
Prism_height = norm(Major_axis*Transform_Matrix)*Prism_height;

%% Check orthogonality of minor and major axes
if abs(dot(Ref_axis_major,Ref_axis_minor)) > eps
    error('Major and minor axes must be orthogonal')
end

%% Calculate size of search space in each direction based on given cutoff and shortest distance between planes
Ma = ceil(max(Cutoff/(Crystal.a*cosd(Crystal.gamma-90)),...
    Cutoff/(Crystal.a*cosd(Crystal.beta-90))));
Mb = ceil(max(Cutoff/(Crystal.b*cosd(Crystal.gamma-90)),...
    Cutoff/(Crystal.b*cosd(Crystal.alpha-90))));
Mc = ceil(max(Cutoff/(Crystal.c*cosd(Crystal.alpha-90)),...
    Cutoff/(Crystal.c*cosd(Crystal.beta-90))));

%% Expand unit cell into 2Ma+1 x 2Mb+1 x 2Mc+1 Supercell of cartesian coordinates       
Search_N = Frac_Coordinates*Transform_Matrix;
Search_a = [0 -Ma:-1 1:Ma]' * a_vec;
Search_b = [0 -Mb:-1 1:Mb]' * b_vec;
Search_c = [0 -Mc:-1 1:Mc]' * c_vec;

Search_Combs = combvec(1:N,1:(2*Ma+1),1:(2*Mb+1),1:(2*Mc+1));
Identities = Search_Combs(1,:);
SCa = Search_Combs(2,:);
SCb = Search_Combs(3,:);
SCc = Search_Combs(4,:);
Cartesian_Coordinates = Search_N(Identities,:) + Search_a(SCa,:) + Search_b(SCb,:) + Search_c(SCc,:) - Prism_center;

%% Calculate projection of each atom along the chosen major axis
N_atoms = size(Cartesian_Coordinates,1);
d_Major = dot(Cartesian_Coordinates,repmat(Ref_axis_major,N_atoms,1),2);

%% Grab all atoms within the chosen height range
Cylinder_min = -Prism_height/2;
Cylinder_max = Prism_height/2;

d_Major_sel_idx = (d_Major >= Cylinder_min) & (d_Major <= Cylinder_max);

Identities_d_Major_sel = Identities(d_Major_sel_idx);
Cartesian_Coordinates_d_Major_sel = Cartesian_Coordinates(d_Major_sel_idx,:);

%% Calculate the unit vectors of the Prism
Prism_point_unit_vectors = zeros(N_points,3);

for i = 0:(N_points-1)
    Angle = 2*pi*i/N_points; % Radians
    Prism_point_unit_vectors(i+1,:) =  Ref_axis_minor*axang2rotm([Ref_axis_major Angle]); % unit vectors
end

%% Project onto the plane orthogonal to the major axis
Cartesian_Coords_projected = (Cartesian_Coordinates_d_Major_sel) - d_Major(d_Major_sel_idx).*Ref_axis_major;
N_atoms_proj = size(Cartesian_Coords_projected,1);

%% Calculate angle between atom vector and prism vectors, get the prism distance
Angle_to_unit_vec = zeros(N_atoms_proj,N_points);
for i = 1:N_points
     Angle_to_unit_vec(:,i) = VecAngle(Cartesian_Coords_projected(:,:),repmat(Prism_point_unit_vectors(i,:),N_atoms_proj,1),2);
end

% Sort by Angle
[Angle_to_point_unit_vectors_sorted,idx] = sort(Angle_to_unit_vec,2);

% Get the corresponding unit vectors for the Prism in each case
fvec_a = Prism_point_unit_vectors(idx(:,1),:);
fvec_b = Prism_point_unit_vectors(idx(:,2),:);

% Find the unit vector that splits the two Prism vectors
x_hat = (fvec_a + fvec_b) ./ vecnorm(fvec_a + fvec_b,2,2);

% project onto this unit vector to get the Prismal distance
x_distance = dot(Cartesian_Coords_projected,x_hat,2);

% Sort by distance
[~,Sort_idx] = sort(x_distance,1);
Sorted_Identities = Identities_UC(Identities_d_Major_sel(Sort_idx));
Sorted_Coordinates = Cartesian_Coordinates_d_Major_sel(Sort_idx,:);

%% Rotate system
for idx = 1:length(Rot)
    Sorted_Coordinates = Sorted_Coordinates*axang2rotm(Rot{idx});
end


%% Zero out the z-axis
min_z = min(Sorted_Coordinates(:,3));
Sorted_Coordinates(:,3) = Sorted_Coordinates(:,3) - min_z;
max_z = max(Sorted_Coordinates(:,3));

%% Separate metal and halide indexes
M_Ind = ismember(Sorted_Identities,["Li" "Na" "K" "Rb" "Cs" "Fe"]);
X_Ind = ~M_Ind;

% Pick out M vs X info
M_Sorted_Coordinates = Sorted_Coordinates(M_Ind,:);
X_Sorted_Coordinates = Sorted_Coordinates(X_Ind,:);

% Pick out closest N_Cluster/2 of both and combine
M_Sel_Coordinates = M_Sorted_Coordinates(1:N_Cluster/2,:);
X_Sel_Coordinates = X_Sorted_Coordinates(1:N_Cluster/2,:);

%% Add in second slab
Slab_2_M = M_Sel_Coordinates;
Slab_2_X = X_Sel_Coordinates;

% Flip z-coordinates
Slab_2_M(:,3) = -Slab_2_M(:,3);
Slab_2_X(:,3) = -Slab_2_X(:,3);

% Place second slab higher by given distance along z coordinate
Slab_2_M(:,3) = Slab_2_M(:,3) + 2*max_z + slab_distance;
Slab_2_X(:,3) = Slab_2_X(:,3) + 2*max_z + slab_distance;

% Recombine
M_Sel_Coordinates = [M_Sel_Coordinates; Slab_2_M];
X_Sel_Coordinates = [X_Sel_Coordinates; Slab_2_X];
N_Cluster = N_Cluster*2;

%% Print into xyz format
filename = fullfile(pwd,[Salt '_' Structure '_N' num2str(N_Cluster) '_Cylinder.xyz']);
fid = fopen(filename, 'a');
% First two lines
fprintf(fid,[num2str(N_Cluster) newline ...
    Salt ' ' Structure ' Cluster of size ' num2str(N_Cluster) newline]);

for jdx = 1:N_Cluster/2
    line = [Metal char(9) num2str(M_Sel_Coordinates(jdx,:),'\t%.10f') newline];
    fprintf(fid,line);
end
for jdx = 1:N_Cluster/2
    line = [Halide char(9) num2str(X_Sel_Coordinates(jdx,:),'\t%.10f') newline];
    fprintf(fid,line);
end
fclose(fid);
system(['"C:\Program Files\VESTA-win64\VESTA.exe" "' filename '"']);
delete(filename)
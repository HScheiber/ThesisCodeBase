
Settings.Salt = 'LiI';
Settings.Structure = 'Rocksalt';
Settings.Use_Conv_cell = true;
Settings.Target_T = 1000;
Settings.Theory = 'JC';
Settings.Model = '';

Settings.Expand_a = 1;
Settings.Expand_b = 1;
Settings.Expand_c = 1;

lattice_plane = [0 0 1]; % The lattice plane to cut along
Minor_axis = [1 0 0]*eul2rotm([deg2rad(0) deg2rad(0) deg2rad(0)]); % The minor axis defines the direction of the first point (A) of the Prism. Must be orthogonal to Major axis
nmol = 512;
Prism_height = 2; % Maximum cylinder height in units of the chosen axis

% This is from the input
Geometry = Default_Crystal(Settings,'Center_Coordinates',true);

% Expand the volume for the given temperature
Exp_Coeff = LiX_Thermal_Expansion(Settings);

Geometry.a = Geometry.a*Settings.Expand_a*Exp_Coeff;
Geometry.b = Geometry.b*Settings.Expand_b*Exp_Coeff;
Geometry.c = Geometry.c*Settings.Expand_c*Exp_Coeff;


N_Cluster = nmol/2; % Number of atoms in the cluster. 
N_points = 4; % number of sides on the prism (6 = hexagon)
Cluster_Center = [1/8 1/8 1/8]; % Location of cluster center in units of a, b, and c vectors
Cutoff = 50; % Angstroms, maximum search distance

% Angles to rotate the final coordinate system by
% Leave as Rot = {} for no rotation
% Rotations are performed left to right.
% Give axis of rotation followed by the angle in radians

Rot = {[0 0 1 deg2rad(45)]};%{[0 1 0 deg2rad(0)] [1 0 0 deg2rad(0)] [1 0 0 deg2rad(0)]}; 

% Crystal = Default_Crystal(Salt,Structure,Conventional);
% Crystal.a = DFT.(Salt).(Structure).a;
% Crystal.b = DFT.(Salt).(Structure).b;
% Crystal.c = DFT.(Salt).(Structure).c;
% Ratios 1:3:1:8 (Al:Si:K:O)


% Get crystal properties
Frac_Coordinates = Geometry.FC;
Identities_UC = Geometry.AtomNames;
N = Geometry.N;

%% Generate transformation matrix from crystal basis to cartesian basis
abc = [Geometry.a 0 0; 0 Geometry.b 0; 0 0 Geometry.c];
Transform_Matrix = abc*Geometry.Transform;
a1_vec = Transform_Matrix(1,:);
a2_vec = Transform_Matrix(2,:);
a3_vec = Transform_Matrix(3,:);


%% Generate reciprocal lattice vectors
V = dot(a1_vec,cross(a2_vec,a3_vec));
b1_vec = (1/V)*cross(a2_vec,a3_vec);
b2_vec = (1/V)*cross(a3_vec,a1_vec);
b3_vec = (1/V)*cross(a1_vec,a2_vec);
Recip_Transform_Matrix = [b1_vec; b2_vec; b3_vec];


Prism_center = Cluster_Center*Transform_Matrix;
Ref_axis_major = lattice_plane*Recip_Transform_Matrix./norm(lattice_plane*Recip_Transform_Matrix);
Ref_axis_minor = (Minor_axis*Transform_Matrix./norm(Minor_axis*Transform_Matrix));%*axang2rotm([Ref_axis_major 2*pi/8]);
Prism_height = norm(( (Ref_axis_major/Transform_Matrix)./norm(Ref_axis_major/Transform_Matrix) )*Transform_Matrix)*Prism_height;

%% Check orthogonality of minor and major axes
if abs(dot(Ref_axis_major,Ref_axis_minor)) > eps
    warning('Major and minor axes are not orthogonal, projecting minor axis onto the plane orthogonal to the major axis.')
    
    proj_nb = dot(Ref_axis_minor,Ref_axis_major).*Ref_axis_major; %
    Ref_axis_minor = (Ref_axis_minor - proj_nb)./norm(Ref_axis_minor - proj_nb);
end

%% Calculate size of search space in each direction based on given cutoff and shortest distance between planes
Ma = ceil(max(Cutoff/(Geometry.a*cosd(Geometry.gamma-90)),...
    Cutoff/(Geometry.a*cosd(Geometry.beta-90))));
Mb = ceil(max(Cutoff/(Geometry.b*cosd(Geometry.gamma-90)),...
    Cutoff/(Geometry.b*cosd(Geometry.alpha-90))));
Mc = ceil(max(Cutoff/(Geometry.c*cosd(Geometry.alpha-90)),...
    Cutoff/(Geometry.c*cosd(Geometry.beta-90))));

%% Expand unit cell into 2Ma+1 x 2Mb+1 x 2Mc+1 Supercell of cartesian coordinates       
Search_N = Frac_Coordinates*Transform_Matrix;
Search_a = [0 -Ma:-1 1:Ma]' * a1_vec;
Search_b = [0 -Mb:-1 1:Mb]' * a2_vec;
Search_c = [0 -Mc:-1 1:Mc]' * a3_vec;

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


%% Rotate system
Sorted_Coordinates = Cartesian_Coordinates_d_Major_sel(Sort_idx,:);
% Define a signed rotation function
vec_angle_signed = @(Va,Vb,Vn) rad2deg(atan2( dot( cross(Vb,Va), Vn), dot(Va,Vb) ));

a1_vec = Geometry.Transform(1,:);
a2_vec = Geometry.Transform(2,:);
a3_vec = Geometry.Transform(3,:);

% vector normal to the plane defined by the z_old and z_new axes
V_normal = cross(a3_vec,Ref_axis_major)./norm(cross(a3_vec,Ref_axis_major));
if ~isnan(V_normal)
    % Calculate the signed angle between the z_old and z_new axes
    z_angle = vec_angle_signed(Ref_axis_major,a3_vec,V_normal);

    % Rotate everything along the normal axis
    Sorted_Coordinates = Sorted_Coordinates*axang2rotm([V_normal deg2rad(z_angle)]);
end
%         ax = axes;
%         hold(ax,'on');
%         vectarrow([0 0 0],[1 0 0])
%         hold on
%         vectarrow([0 0 0],[0 1 0])
%         hold on
%         vectarrow([0 0 0],[0 0 1])
%         daspect(gca,[1 1 1])
%         hold on
%         vectarrow([0 0 0],Ref_axis_major)
%         hold on
%         vectarrow([0 0 0],a3_vec)
%         hold on
%         vectarrow([0 0 0],V_normal)
%         daspect(gca,[1 1 1])

%% Align the new x axis with the minor axis
z_angle = vec_angle_signed(a1_vec,Ref_axis_minor,a3_vec);
Sorted_Coordinates = Sorted_Coordinates*axang2rotm([a3_vec deg2rad(z_angle)]);

for idx = 1:length(Rot)
    Sorted_Coordinates = Sorted_Coordinates*axang2rotm(Rot{idx});
end

%% Zero out the z-axis
min_z = min(Sorted_Coordinates(:,3));
Sorted_Coordinates(:,3) = Sorted_Coordinates(:,3) - min_z;

%% Separate metal and halide indexes
M_Ind = ismember(Sorted_Identities,["Li" "Na" "K" "Rb" "Cs" "Fe"]);
X_Ind = ~M_Ind;

% Pick out M vs X info
M_Sorted_Coordinates = Sorted_Coordinates(M_Ind,:);
X_Sorted_Coordinates = Sorted_Coordinates(X_Ind,:);

% Pick out closest N_Cluster/2 of both and combine
M_Sel_Coordinates = M_Sorted_Coordinates(1:N_Cluster/2,:);
X_Sel_Coordinates = X_Sorted_Coordinates(1:N_Cluster/2,:);

% Recombine back into single cluster
Final_Coordinates = [M_Sel_Coordinates; X_Sel_Coordinates];
Final_Identities = [repmat({Geometry.Metal},N_Cluster/2,1); repmat({Geometry.Halide},N_Cluster/2,1)];


%% Print into xyz format
filename = fullfile(pwd,[Settings.Salt '_' Settings.Structure '_N' num2str(N_Cluster) '_Cylinder.xyz']);
fid = fopen(filename, 'a');
% First two lines
fprintf(fid,[num2str(N_Cluster) newline ...
    Settings.Salt ' ' Settings.Structure ' Cluster of size ' num2str(N_Cluster) newline]);

for jdx = 1:N_Cluster
    line = [char(9) Final_Identities{jdx} char(9) num2str(Final_Coordinates(jdx,:),'\t%.10f') newline];
    fprintf(fid,line);
end
fclose(fid);
system(['"C:\Program Files\VESTA-win64\VESTA.exe" "' filename '"']);
delete(filename)


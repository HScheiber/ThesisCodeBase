Settings.Salt = 'LiI';
Settings.Structure = 'Rocksalt';
Settings.Use_Conv_cell = true;
Settings.Target_T = 1000;
Settings.Theory = 'JC';
Settings.Model = '';

Settings.Expand_a = 1;
Settings.Expand_b = 1;
Settings.Expand_c = 1;

Settings.lattice_plane = [0 0 1]; % The lattice plane to cut along in the z-dimension
Settings.Minor_axis = [1 0 0]*eul2rotm([deg2rad(0) deg2rad(0) deg2rad(0)]); % The minor axis defines the direction of the first point (A) of the Prism. Must be orthogonal to Major axis
Settings.N_Cluster = 512; % number of atoms in the cluster
Settings.Prism_height = 2; % Maximum cylinder height in units of the chosen axis
Settings.Cluster_Center = [1/8 1/8 1/8]; % Location of cluster center in units of a, b, and c vectors
Settings.ShowCluster = true; % Opens the cluster for visualization when true

% This is from the input
Geometry = Default_Crystal(Settings,'Center_Coordinates',true);

% Expand the volume for the given temperature
Exp_Coeff = LiX_Thermal_Expansion(Settings);

Geometry.a = Geometry.a*Settings.Expand_a*Exp_Coeff;
Geometry.b = Geometry.b*Settings.Expand_b*Exp_Coeff;
Geometry.c = Geometry.c*Settings.Expand_c*Exp_Coeff;

[Final_Coordinates,Final_Identities] = GenPrismCluster(Settings,Geometry);


% Calculate solid volume per atom
abc = [Geometry.a 0 0; 0 Geometry.b 0; 0 0 Geometry.c];
Transform_Matrix = abc*Geometry.Transform;
a1_vec = Transform_Matrix(1,:);
a2_vec = Transform_Matrix(2,:);
a3_vec = Transform_Matrix(3,:);
Vol = abs(dot(cross(a1_vec,a2_vec),a3_vec))/Geometry.N; % [A^3 / atom]
Occupied_Vol = Settings.N_Cluster*Vol; % [A^3]




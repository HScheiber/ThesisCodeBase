Salt = 'LiI';
Structure = 'NiAs';
Theory = 'TMTPSS-rVV10L';
Conventional = true;
DFT = Load_Best_DFT_Data('Theory',Theory,'Conventional',Conventional);
N_Cluster = 200; % Must be an even number to ensure charge neutrality
Major_axis = [0 0 1]; % The major axis of the cylinder in terms of [a b c]
Cluster_Center = [1/2 1/2 1/2]; % Location of cluster center in units of a, b, and c vectors
Cylinder_height = 0.99; % Maximum cylinder height in units of the chosen axis
Cutoff = 50; % Angstroms, maximum search distance

% Angles to rotate the final coordinate system by
% Leave as Rot = {} for no rotation
% Rotations are performed left to right.
% Give axis of rotation followed by the angle in radians
Rot = {[0 1 0 deg2rad(45)] [1 0 0 deg2rad(-35.264389682754661)]}; 


Settings = Initialize_MD_Settings;
Settings.Salt = Salt;
Settings.Structure = Structure;
Settings.Use_Conv_cell = Conventional;
Crystal = Default_Crystal(Settings);
Crystal.a = DFT.(Salt).(Structure).a;
Crystal.b = DFT.(Salt).(Structure).b;
Crystal.c = DFT.(Salt).(Structure).c;
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

Cylinder_center = Cluster_Center*Transform_Matrix;
Ref_axis = Major_axis*Transform_Matrix;
Cylinder_height = norm(Ref_axis)*Cylinder_height;

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
Cartesian_Coordinates = Search_N(Identities,:) + Search_a(SCa,:) + Search_b(SCb,:) + Search_c(SCc,:);

%% Calculate projection of each atom along the chosen reference axis
N_atoms = size(Cartesian_Coordinates,1);
d_Major = dot(Cartesian_Coordinates-Cylinder_center,repmat(Ref_axis,N_atoms,1),2)./norm(Ref_axis);

%% Grab all atoms within the chosen height range
Cylinder_min = -Cylinder_height/2;
Cylinder_max = Cylinder_height/2;

d_Major_sel_idx = (d_Major >= Cylinder_min) & (d_Major <= Cylinder_max);

Identities_d_Major_sel = Identities(d_Major_sel_idx);
Cartesian_Coordinates_d_Major_sel = Cartesian_Coordinates(d_Major_sel_idx,:);

%% Calculate all atoms distance radially from the chosen reference axis
major_axis_projections = d_Major(d_Major_sel_idx).*Ref_axis./norm(Ref_axis);
distance_from_major_axis = vecnorm(Cartesian_Coordinates_d_Major_sel - major_axis_projections,2,2);

%% Sort atoms by distance from the chosen reference axis
[Sorted_rAB,Sort_idx] = sort(distance_from_major_axis);
Sorted_Identities = Identities_UC(Identities_d_Major_sel(Sort_idx));
Sorted_Coordinates = Cartesian_Coordinates_d_Major_sel(Sort_idx,:);

%% Rotate system
for idx = 1:length(Rot)
    Sorted_Coordinates = Sorted_Coordinates*axang2rotm(Rot{idx});
end

%% Separate metal and halide indexes
M_Ind = ismember(Sorted_Identities,["Li" "Na" "K" "Rb" "Cs"]);
X_Ind = ~M_Ind;

% Pick out M vs X info
M_Sorted_rAB = Sorted_rAB(M_Ind);
M_Sorted_Coordinates = Sorted_Coordinates(M_Ind,:);
X_Sored_rAB = Sorted_rAB(X_Ind);
X_Sorted_Coordinates = Sorted_Coordinates(X_Ind,:);

% Pick out closest N_Cluster/2 of both and combine
M_Sel_Coordinates = M_Sorted_Coordinates(1:N_Cluster/2,:);
X_Sel_Coordinates = X_Sorted_Coordinates(1:N_Cluster/2,:);

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
system(['"C:\Program Files\VESTA-win64\VESTA.exe" "' filename '"'])
delete(filename)
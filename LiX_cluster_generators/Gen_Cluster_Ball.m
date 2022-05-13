Salt = 'LiI';
Structure = 'NiAs';
Theory = 'TMTPSS-rVV10L';
Conventional = true;
DFT = Load_Best_DFT_Data('Theory',Theory,'Conventional',Conventional);
N_Cluster = 200; % Must be an even number to ensure charge neutrality
Cutoff = 50; % Angstroms
Cluster_Center = [1/2 1/2 1/2]; % Location of cluster center in units of a, b, and c vectors

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

%% Initialize reference
% Coordinates of reference atom
Reference_loc = Cluster_Center*Transform_Matrix;

%% Calculate distances from reference atom
rAB = vecnorm(Reference_loc - Cartesian_Coordinates,2,2);


%% Sort atoms by distance from reference atom
[Sorted_rAB,Sort_idx] = sort(rAB);
Sorted_Identities = Identities_UC(Identities(Sort_idx));
Sorted_Coordinates = Cartesian_Coordinates(Sort_idx,:);

% Separate metal and halide indexes
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
filename = fullfile(pwd,[Salt '_' Structure '_N' num2str(N_Cluster) '_Ball.xyz']);
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
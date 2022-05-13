% A function to generate stars of coordination shells given input structure
function [MA_Clus,MX_Clus,MM_Clus,XA_Clus,XM_Clus,XX_Clus] = Gen_Stars_Cartesian(CryStruc,Structure,Metal,Halide,StartPoint,Cutoff,MaxNum)
Met_Ind = strcmpi(CryStruc.FC(:,1),Metal);
Hal_Ind = strcmpi(CryStruc.FC(:,1),Halide);
FC_Metal_AU = mod([CryStruc.FC{Met_Ind,2:4}],1);
FC_Halide_AU = mod([CryStruc.FC{Hal_Ind,2:4}],1);

% number of atoms in unit cell
switch lower(Structure)
    case 'rocksalt'
        N = 8; 
    case 'wurtzite'
        N = 4;
    case 'betabeo'
        N = 8;
    case 'cscl'
        N = 2;
    case 'fivefive'
        N = 8;
    case 'nias'
        N = 4;
    case 'sphalerite'
        N = 8;
    otherwise
        error(['Unknown structure: ' Structure])
end

[CryStruc.FC_Metal,CryStruc.FC_Halide] = UnitCell_FractionalCoords(FC_Metal_AU,FC_Halide_AU,Structure);

Frac_Coordinates = [CryStruc.FC_Metal; CryStruc.FC_Halide];

Identities_UC = [repmat(string(Metal),N/2,1);
    repmat(string(Halide),N/2,1)];

Idx_M = 1;
Idx_X = size(CryStruc.FC_Metal,1)+1;

%% Generate transformation matrix from crystal basis to cartesian basis
abc = [CryStruc.a 0 0; 0 CryStruc.b 0; 0 0 CryStruc.c];
Transform_Matrix = abc*GenTransformMatrix(CryStruc);
a_vec = Transform_Matrix(1,:);
b_vec = Transform_Matrix(2,:);
c_vec = Transform_Matrix(3,:);

%% Calculate size of search space in each direction based on given cutoff and shortest distance between planes
Ma = ceil(max(Cutoff/(CryStruc.a*cosd(CryStruc.gamma-90)),...
    Cutoff/(CryStruc.a*cosd(CryStruc.beta-90))));
Mb = ceil(max(Cutoff/(CryStruc.b*cosd(CryStruc.gamma-90)),...
    Cutoff/(CryStruc.b*cosd(CryStruc.alpha-90))));
Mc = ceil(max(Cutoff/(CryStruc.c*cosd(CryStruc.alpha-90)),...
    Cutoff/(CryStruc.c*cosd(CryStruc.beta-90))));

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

%% Generate List of stars
for idx = [Idx_M Idx_X] % For each atom in asymmetric unit
    
    %% Initialize reference
    % Coordinates of reference atom
    Reference_Atom = Cartesian_Coordinates(idx,:);

    %% Pick out atoms closer than Cutoff from reference atom
    rAB = vecnorm(Reference_Atom - Cartesian_Coordinates,2,2);

    % Which atoms are closer than cutoff and further than first search point?
    Cutoff_idx = (rAB <= Cutoff & rAB >= StartPoint);
    % Skip reference atom
    Cutoff_idx(idx) = false;

    % Take items within cutoff
    Cutoff_Identities = Identities(Cutoff_idx);
    Cutoff_rAB = rAB(Cutoff_idx);
    Cutoff_Cartesian_Coordinates = Cartesian_Coordinates(Cutoff_idx,:);
    
    % Sort by distance
    [~,Sort_idx] = sort(Cutoff_rAB);
    Sorted_Identities = Identities_UC(Cutoff_Identities(Sort_idx));
    Sorted_Cartesian_Coordinates = Cutoff_Cartesian_Coordinates(Sort_idx,:);
    
    % Separate metal and halide indexes
    M_Ind = ismember(Sorted_Identities,["Li" "Na" "K" "Rb" "Cs"]);
    X_Ind = ~M_Ind;

    % Pick out M vs X info
    M_xyz = Sorted_Cartesian_Coordinates(M_Ind,:);
    X_xyz = Sorted_Cartesian_Coordinates(X_Ind,:);
    
    % Limit cluster sizes to max number
    if size(Sorted_Cartesian_Coordinates,1) > MaxNum
        Sorted_Cartesian_Coordinates = Sorted_Cartesian_Coordinates(1:MaxNum,:);
        M_xyz = M_xyz(1:MaxNum,:);
        X_xyz = X_xyz(1:MaxNum,:);
    end
    if size(M_xyz,1) > MaxNum
        M_xyz = M_xyz(1:MaxNum,:);
    end
    if size(X_xyz,1) > MaxNum
        X_xyz = X_xyz(1:MaxNum,:);
    end
    
    % Transform to bond vectors
    if idx == Idx_M
        MA_Clus = Sorted_Cartesian_Coordinates - Reference_Atom;
        MX_Clus = X_xyz - Reference_Atom;
        MM_Clus = M_xyz - Reference_Atom;
    elseif idx == Idx_X
        XA_Clus = Sorted_Cartesian_Coordinates - Reference_Atom;
        XM_Clus = M_xyz - Reference_Atom;
        XX_Clus = X_xyz - Reference_Atom;
    end
end

end
% A function to generate stars of coordination shells given input structure
function [MX_BL,MM_BL,XX_BL,XM_BL] = Gen_First_Stars(CryStruc,Structure,Metal,Halide,Cutoff)
tol = 1e-8; % Two radii are equivalent when there difference is smaller than this (Angstroms)
Met_Ind = strcmpi(CryStruc.FC(:,1),Metal);
Hal_Ind = strcmpi(CryStruc.FC(:,1),Halide);
FC_Metal_AU = mod([CryStruc.FC{Met_Ind,2:4}],1);
FC_Halide_AU = mod([CryStruc.FC{Hal_Ind,2:4}],1);

% number of atoms in unit cell
switch lower(Structure)
    case 'rocksalt'
        N = 8;
        MX_CN = 6;
        MM_CN = 12;
        XX_CN = 12;
    case 'wurtzite'
        N = 4;
        MX_CN = 4;
        MM_CN = 12;
        XX_CN = 12;
    case 'betabeo'
        N = 8;
        MX_CN = 4;
        MM_CN = 11;
        XX_CN = 11;
    case 'cscl'
        N = 2;
        MX_CN = 8;
        MM_CN = 6;
        XX_CN = 6;
    case 'fivefive'
        N = 8;
        MX_CN = 5;
        MM_CN = 12;
        XX_CN = 12;
    case 'nias'
        N = 4;
        MX_CN = 6;
        MM_CN = 8;
        XX_CN = 12;
    case 'sphalerite'
        N = 8;
        MX_CN = 4;
        MM_CN = 12;
        XX_CN = 12;
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

    % Which atoms are closer than cutoff?
    Cutoff_idx = (rAB <= Cutoff);
    % Skip reference atom
    Cutoff_idx(idx) = false;

    % Take items within cutoff
    Cutoff_Identities = Identities(Cutoff_idx);
    Cutoff_rAB = rAB(Cutoff_idx);
    
    % Sort by distance
    [Sorted_rAB,Sort_idx] = sort(Cutoff_rAB);
    Sorted_Identities = Identities_UC(Cutoff_Identities(Sort_idx));

    % Separate metal and halide indexes
    M_Ind = ismember(Sorted_Identities,["Li" "Na" "K" "Rb" "Cs"]);
    X_Ind = ~M_Ind;

    % Pick out M vs X info
    M_BL = Sorted_rAB(M_Ind);
    X_BL = Sorted_rAB(X_Ind);

    % Collapse into stars
    [M_BL_U,IA,~] = uniquetol(M_BL,tol);
    M_StarSize = diff([IA; length(M_BL)+1]);
    [X_BL_U,IA,~] = uniquetol(X_BL,tol);
    X_StarSize = diff([IA; length(X_BL)+1]);
    
    if idx == Idx_M
        MX_BL = [X_BL_U X_StarSize];
        
        MM_BL = [M_BL_U M_StarSize];
    elseif idx == Idx_X
        XM_BL = [M_BL_U M_StarSize];
        
        XX_BL = [X_BL_U X_StarSize];
    end
end

% Collapse coordination shells to specified size
while MX_BL(1,2) < MX_CN
    MX_BL(1,1) = (MX_BL(1,1)*MX_BL(1,2) + MX_BL(2,1)*MX_BL(2,2))/...
        (MX_BL(1,2) + MX_BL(2,2)); % new distance = weighted mean
    MX_BL(1,2) = MX_BL(1,2) + MX_BL(2,2); % New CN is sum of old
    MX_BL(2,:) = []; % Erase second shell
end
while MM_BL(1,2) < MM_CN
    MM_BL(1,1) = (MM_BL(1,1)*MM_BL(1,2) + MM_BL(2,1)*MM_BL(2,2))/...
        (MM_BL(1,2) + MM_BL(2,2)); % new distance = weighted mean
    MM_BL(1,2) = MM_BL(1,2) + MM_BL(2,2); % New CN is sum of old
    MM_BL(2,:) = []; % Erase second shell
end
while XX_BL(1,2) < XX_CN
    XX_BL(1,1) = (XX_BL(1,1)*XX_BL(1,2) + XX_BL(2,1)*XX_BL(2,2))/...
        (XX_BL(1,2) + XX_BL(2,2)); % new distance = weighted mean
    XX_BL(1,2) = XX_BL(1,2) + XX_BL(2,2); % New CN is sum of old
    XX_BL(2,:) = []; % Erase second shell
end

end
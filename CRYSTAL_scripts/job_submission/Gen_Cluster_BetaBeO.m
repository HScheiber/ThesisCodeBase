function Coordinates_Out = Gen_Cluster_BetaBeO(Salt,Coordinates_In,BSSE_id,...
	NStars,Rmax,a,b,c)

[Metal,Halide] = Separate_Metal_Halide(Salt);
Metal_Z = elements('Symbol',Metal,'atomic_number');
Halide_Z = elements('Symbol',Halide,'atomic_number');

N = 3; % Search space (i.e. N x N x N supercell centered around a central cell)
M = 8; % Atoms per unit cell

% Grab coordinates from input as string
Frac_Coordinates = regexp(Coordinates_In,['[0-9]+\n(' num2str(Metal_Z) ...
    ') (([0-9]|-|\+|\.|E)+) (([0-9]|-|\+|\.|E)+) (([0-9]|-|\+|\.|E)+)\n'...
    '([0-9]{0,2}' num2str(Halide_Z) ...
    ') (([0-9]|-|\+|\.|E)+) (([0-9]|-|\+|\.|E)+) (([0-9]|-|\+|\.|E)+)'],'tokens');

% Convert to floats
Frac_Coordinates = cellfun(@str2num,Frac_Coordinates{1});

FC_Metal = Frac_Coordinates(2:4);
FC_Halide = Frac_Coordinates(6:8);

Metal_ZCard = num2str(Frac_Coordinates(1));
Halide_ZCard = num2str(Frac_Coordinates(5));

% Ensure coordinates are in reference unit cell
FC_Metal  = arrayfun(@(a) mod(a,1),FC_Metal);
FC_Halide = arrayfun(@(a) mod(a,1),FC_Halide);

FC_Metal  = mod([FC_Metal; ...
    -FC_Metal(1) -FC_Metal(2) FC_Metal(3); ...
    (1/2)-FC_Metal(2) (1/2)+FC_Metal(1) (1/2)+FC_Metal(3); ...
    (1/2)+FC_Metal(2) (1/2)-FC_Metal(1) (1/2)-FC_Metal(3)],1);
FC_Halide  = mod([FC_Halide; ...
    -FC_Halide(1) -FC_Halide(2) FC_Halide(3); ...
    (1/2)-FC_Halide(2) (1/2)+FC_Halide(1) (1/2)+FC_Halide(3); ...
    (1/2)+FC_Halide(2) (1/2)-FC_Halide(1) (1/2)-FC_Halide(3)],1);

a_vec = [a 0 0];
b_vec = [0 b 0];
c_vec = [0 0 c];

Transform_Mat = [a_vec; b_vec; c_vec];

Metal_Cart_Coordinates = nan(4,3);
Halide_Cart_Coordinates = nan(4,3);
for i = 1:4
    Metal_Cart_Coordinates(i,:) = FC_Metal(i,:)*Transform_Mat;
    Halide_Cart_Coordinates(i,:) = FC_Halide(i,:)*Transform_Mat;
end

Metal_Cart_Coordinates = nan((M/2)*(2*N+1)^3,3);
Halide_Cart_Coordinates = nan((M/2)*(2*N+1)^3,3);
x = 0;
for i = [0 -N:-1 1:N]
    for j = [0 -N:-1 1:N]
        for k = [0 -N:-1 1:N]
            for idx = 1:4
                x = x+1;
                Metal_Cart_Coordinates(x,:) = FC_Metal(idx,:)*Transform_Mat + (i.*a_vec + j.*b_vec + k.*c_vec);
                Halide_Cart_Coordinates(x,:) = FC_Halide(idx,:)*Transform_Mat + (i.*a_vec + j.*b_vec + k.*c_vec);
            end
        end
    end
end

% Remove reference atom from coordinates array
if strcmpi(BSSE_id,'Metal')
    Reference_Atom = Metal_Cart_Coordinates(1,:);
    Metal_Cart_Coordinates = Metal_Cart_Coordinates(2:end,:);
elseif strcmpi(BSSE_id,'Halide')
    Reference_Atom = Halide_Cart_Coordinates(1,:);
    Halide_Cart_Coordinates = Halide_Cart_Coordinates(2:end,:);
else
    error(['Unknown argument for BSSE_id: ' BSSE_id])
end

% Pick out atoms closer than Search_Distance from reference atom
Metal_Coordinates_Of_Interest = nan((M/2)*(2*N+1)^3,3);
Halide_Coordinates_Of_Interest = nan((M/2)*(2*N+1)^3,3);
x = 0;
for i = 1:size(Metal_Cart_Coordinates,1)
    if norm(Reference_Atom - Metal_Cart_Coordinates(i,:)) <= Rmax
        x = x+1;
        Metal_Coordinates_Of_Interest(x,:) = Metal_Cart_Coordinates(i,:);
    end
end
Metal_Coordinates_Of_Interest = Metal_Coordinates_Of_Interest(all(~isnan(Metal_Coordinates_Of_Interest),2),:);
x = 0;
for i = 1:size(Halide_Cart_Coordinates,1)
    if norm(Reference_Atom - Halide_Cart_Coordinates(i,:)) <= Rmax
        x = x + 1;
        Halide_Coordinates_Of_Interest(x,:) = Halide_Cart_Coordinates(i,:);
    end
end
Halide_Coordinates_Of_Interest = Halide_Coordinates_Of_Interest(all(~isnan(Halide_Coordinates_Of_Interest),2),:);

Coordinates_Of_interest = [Metal_Coordinates_Of_Interest; Halide_Coordinates_Of_Interest];

Vectors_From_Reference = Coordinates_Of_interest - repmat(Reference_Atom,size(Coordinates_Of_interest,1),1);


Distance_From_Reference = nan(size(Vectors_From_Reference,1),1);
for i = 1:size(Vectors_From_Reference,1)
    Distance_From_Reference(i) = norm(Vectors_From_Reference(i,:));
end
% Array of distances of shells from reference atom, ordered
Star_Distances = unique(Distance_From_Reference);

% Pick out distance given my input Stars
Max_Dist = Star_Distances(NStars);

NM = size(Metal_Coordinates_Of_Interest,1);
NH = size(Halide_Coordinates_Of_Interest,1);

Metal_Coordinates_Final = nan(NM,3);
Halide_Coordinates_Final = nan(NH,3);
x = 0;
for i = 1:size(Metal_Coordinates_Of_Interest,1)
    if norm(Reference_Atom - Metal_Coordinates_Of_Interest(i,:)) <= Max_Dist
        x = x+1;
        Metal_Coordinates_Final(x,:) = Metal_Coordinates_Of_Interest(i,:);
    end
end
Metal_Coordinates_Final = Metal_Coordinates_Final(all(~isnan(Metal_Coordinates_Final),2),:);
x = 0;
for i = 1:size(Halide_Coordinates_Of_Interest,1)
    if norm(Reference_Atom - Halide_Coordinates_Of_Interest(i,:)) <= Max_Dist
        x = x + 1;
        Halide_Coordinates_Final(x,:) = Halide_Coordinates_Of_Interest(i,:);
    end
end
Halide_Coordinates_Final = Halide_Coordinates_Final(all(~isnan(Halide_Coordinates_Final),2),:);

if strcmp(BSSE_id,'Metal')
    Ref_Cart = sprintf([Metal_ZCard ' %18.15f %18.15f %18.15f\n'],Reference_Atom);
elseif strcmp(BSSE_id,'Halide')
    Ref_Cart = sprintf([Halide_ZCard ' %18.15f %18.15f %18.15f\n'],Reference_Atom);
end

Metal_Cart = sprintf('100 %18.15f %18.15f %18.15f\n',Metal_Coordinates_Final');
Halide_Cart = sprintf('200 %18.15f %18.15f %18.15f\n',Halide_Coordinates_Final');

Tot_Num = 1 + size(Metal_Coordinates_Final,1) + size(Halide_Coordinates_Final,1);

Coordinates_Out = [num2str(Tot_Num) newline Ref_Cart Metal_Cart regexprep(Halide_Cart,'\n$','')];
end

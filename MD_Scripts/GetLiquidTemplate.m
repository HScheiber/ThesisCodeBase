function [Filename,Geometry,CoordType] = GetLiquidTemplate(home,Salt,Target_T,Target_P,N_atoms,CoordType,Geometry)
Ang_per_nm = 10;
[Metal,Halide] = Separate_Metal_Halide(Salt);

% Grab liquid templates
Templatedir = [home filesep 'templates' filesep upper(CoordType) '_Templates'];
Templates = dir([Templatedir filesep 'Liquid*.' CoordType]);
Template_Names = {Templates.name};

if isempty(Template_Names) && strcmpi(CoordType,'g96')
    CoordType = 'gro';
    disp('Unable to find g96 template, using gro template')
    Templatedir = [home filesep 'templates' filesep upper(CoordType) '_Templates'];
    Templates = dir([Templatedir filesep 'Liquid*.' CoordType]);
    Template_Names = {Templates.name};
end

% Grab templates of the correct salt
SaltIDX = cellfun(@(x) contains(x,Salt),Template_Names);
Template_Salts = Template_Names(SaltIDX);

% Get ref T
RefTemps_c = cellfun(@(x) regexp(x,'_T([0-9]+)_','tokens','once'),Template_Salts);
RefTemps = cellfun(@str2double,RefTemps_c);

% Get ref N
RefNum_c = cellfun(@(x) regexp(x,'_N([0-9]+)_','tokens','once'),Template_Salts);
RefNum = cellfun(@str2double,RefNum_c);

% Get ref P
RefPres_c = cellfun(@(x) regexp(x,['_P([0-9]+)\.' CoordType],'tokens','once'),Template_Salts);
RefPres = cellfun(@str2double,RefPres_c);

% Grab only those with N >= N_atoms
N_Idx = RefNum >= N_atoms;
RefTemps = RefTemps(N_Idx);
RefNum = RefNum(N_Idx);
RefPres = RefPres(N_Idx);
Template_Salts = Template_Salts(N_Idx);

% Grab those remaining closest to T of interest
if length(Template_Salts) > 1
    Min_dT = min(abs(RefTemps-Target_T));
    T_Idx = (Min_dT == abs(RefTemps-Target_T));
    RefNum = RefNum(T_Idx);
    RefPres = RefPres(T_Idx);
    Template_Salts = Template_Salts(T_Idx);
elseif isempty(Template_Salts)
    error('No template available for requested initial conditions');
end

% Grab those remaining closest to P of interest
if length(Template_Salts) > 1
    Min_dP = min(abs(RefPres-Target_P(1)));
    P_Idx = (Min_dP == abs(RefPres-Target_P(1)));
    RefNum = RefNum(P_Idx);
    Template_Salts = Template_Salts(P_Idx);
end

% Grab those remaining closest to N of interest
if length(Template_Salts) > 1
    Min_dN = min(abs(RefNum-N_atoms));
    N_Idx = (Min_dN == abs(RefNum-N_atoms));
    RefNum = RefNum(N_Idx);
    Template_Salts = Template_Salts(N_Idx);
end

% If still multiple outputs, choose first one
Filename = [Templatedir filesep Template_Salts{1}];
Geometry.N = RefNum(1);
Geometry.NF = Geometry.N/2;
Geometry.Salt = Salt;

% Get structure info
GeomText = fileread(Filename);
GeomText = regexprep(GeomText,'(\r|\n)+','\n');
if strcmpi(CoordType,'gro')
    
    grodat = load_gro_file(Filename);
    Num_Mols = textscan(GeomText,'%f',1,'Delimiter','','whitespace','','headerlines',1);
    N = Num_Mols{1};
    
    Coord_Data = textscan(GeomText,'%5f%-5s%5s%5f%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',N,...
        'Delimiter','','whitespace','','headerlines',2);
    
    
    M_Ind = cellfun(@(x) strcmpi(x,pad(Metal,5,'right')),Coord_Data{2});
    X_Ind = cellfun(@(x) strcmpi(x,pad(Halide,5,'right')),Coord_Data{2});
    
    Coords = [Coord_Data{5} Coord_Data{6} Coord_Data{7}];
    Met_XYZ = Coords(M_Ind,:);
    Hal_XYZ = Coords(X_Ind,:);
    
    boxcoords = textscan(GeomText,'%f %f %f %f %f %f %f %f %f',...
        'Delimiter',' ','MultipleDelimsAsOne',true,'headerlines',N+2);

    % Check for empty coords
    indx = cellfun(@isnan,boxcoords); % true for empty cells
    boxcoords(indx) = {0}; % replace by a cell with a zero
    
    Geometry.xyz = grodat.xyz;
    Geometry.boxcoords = boxcoords;
    Geometry.res_number = grodat.res_number;
    Geometry.res_name = strtrim(grodat.res_name);
    Geometry.atom_name = strtrim(grodat.atom_name);
    Geometry.atom_number = grodat.atom_number;
    
elseif strcmpi(CoordType,'g96')
    Positions = regexp(GeomText,'(?<=POSITION\n)(.+?)(?=(\n)END)','match','ONCE');
    Box = regexp(GeomText,'(?<=BOX\n)(.+?)(?=(\n)END)','match','ONCE');
    
    if isempty(Positions)
        Positions = regexp(GeomText,'(?<=POSITIONRED\n)(.+?)(?=(\n)END)','match','ONCE');
        coord_format = '%15.9f%15.9f%15.9f\n';
        Coord_Data = textscan(Positions,coord_format,...
            'Delimiter','','whitespace','');
        Coords = [Coord_Data{1} Coord_Data{2} Coord_Data{3}];
    else
        coord_format = '%6f%6s%6s%6f%15.9f%15.9f%15.9f\n';
        % g96 file format
        Coord_Data = textscan(Positions,coord_format,...
            'Delimiter','','whitespace','');
        Coords = [Coord_Data{5} Coord_Data{6} Coord_Data{7}];
    end
    
    M_Ind = cellfun(@(x) strcmpi(x,pad(Metal,6,'right')),Coord_Data{1});
    X_Ind = cellfun(@(x) strcmpi(x,pad(Halide,6,'right')),Coord_Data{1});
    
    Met_XYZ = Coords(M_Ind,:);
    Hal_XYZ = Coords(X_Ind,:);
    
    boxcoords = textscan(Box,'%f %f %f %f %f %f %f %f %f',...
        'Delimiter',' ','MultipleDelimsAsOne',true);

    % Check for empty coords
    indx = cellfun(@isempty,boxcoords); % true for empty cells
    boxcoords(indx) = {0}; % replace by a cell with a zero
    
    Geometry.xyz = Coords;
    Geometry.boxcoords = boxcoords;
    Geometry.res_number = Coord_Data{1};
    Geometry.res_name = strtrim(Coord_Data{2});
    Geometry.atom_name = Coord_Data{3};
    Geometry.atom_number = Coord_Data{4};
end

% Lattice vectors
a_vec = [boxcoords{1} boxcoords{4} boxcoords{5}]; % units nm
b_vec = [boxcoords{6} boxcoords{2} boxcoords{7}]; % units nm
c_vec = [boxcoords{8} boxcoords{9} boxcoords{3}]; % units nm

% Transformation matrix (unit size)
Geometry.Transform = [a_vec./norm(a_vec) ; b_vec./norm(b_vec) ; c_vec./norm(c_vec)];
TM = [a_vec ; b_vec ; c_vec];

% Box vector lengths in Angstroms
Geometry.a = norm(a_vec)*Ang_per_nm;
Geometry.b = norm(b_vec)*Ang_per_nm;
Geometry.c = norm(c_vec)*Ang_per_nm;

Geometry.alpha = VecAngle(b_vec,c_vec);
Geometry.beta = VecAngle(a_vec,c_vec);
Geometry.gamma = VecAngle(a_vec,b_vec);

% Convert to Fractional coordinates
Geometry.FC_Metal = Met_XYZ/TM;
Geometry.FC_Halide = Hal_XYZ/TM;
Geometry.FC = Coords/TM;

end
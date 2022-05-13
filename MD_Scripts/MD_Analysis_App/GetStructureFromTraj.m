function [Filename,Geometry] = GetStructureFromTraj(Prev_geom_loc,Salt,WorkDir,...
    Prev_geom_time,Geometry,JobName,gmx,CoordType)

Ang_per_nm = 10;
[Metal,Halide] = Separate_Metal_Halide(Salt);
t0 = num2str(Prev_geom_time*1000); % ps

% Grab trajectory and topology files
traj_file_dir = dir([Prev_geom_loc filesep '*.trr']);
if ~isfield(traj_file_dir,'folder')
    traj_file_dir.folder = Prev_geom_loc;
end

traj_file = windows2unix(fullfile(traj_file_dir.folder,traj_file_dir.name));
tpr_file_dir  = dir([Prev_geom_loc filesep strrep(traj_file_dir.name,'.trr','.gro')]);
if ~isfield(tpr_file_dir,'folder')
    tpr_file_dir.folder = Prev_geom_loc;
end
tpr_file = windows2unix(fullfile(tpr_file_dir.folder,tpr_file_dir.name));

Filename = fullfile(WorkDir,[JobName '.' CoordType]);

if ispc
    gmx = strrep(gmx,'gmx','echo 0 ^| gmx');
else
    gmx = strrep(gmx,'gmx','echo 0 | gmx');
end


% Grab the structure in question using gromacs trjconv
gmx_input = [gmx ' trjconv -b ' t0 ' -e ' t0 ...
    ' -f ' traj_file ' -s ' tpr_file ' -o ' windows2unix(Filename)];

disp('Extracting initial geometry from previous trajectory with trjconv...')
[err,outp] = system(gmx_input);
if err ~= 0 
    error(['Problem with gmx trjconv command:' newline outp]);
end
disp('Trajectory configuration extracted!')


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
    % Find the gro file for positions reference
    grodat = load_gro_file(strrep(fullfile(tpr_file_dir.folder,tpr_file_dir.name),'.tpr','.gro'));
    
    Postxt = regexp(GeomText,'POSITIONR*E*D*\n(.+?)\nEND','tokens','once');
    Position_data = textscan(Postxt{1},'%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    Coords = [Position_data{1} Position_data{2} Position_data{3}];
    
    N = grodat.N_atoms;
    
    Boxtxt = regexp(GeomText,'BOX\n(.+?)\nEND','tokens','once');
    M_Ind = cellfun(@(x) strcmpi(x,Metal),strtrim(grodat.atom_name));
    X_Ind = cellfun(@(x) strcmpi(x,Halide),strtrim(grodat.atom_name));
    
    Met_XYZ = Coords(M_Ind,:);
    Hal_XYZ = Coords(X_Ind,:);
    
    boxcoords = textscan(Boxtxt{1},'%f %f %f %f %f %f %f %f %f',...
        'Delimiter',' ','MultipleDelimsAsOne',true);

    % Check for empty coords
    indx = cellfun(@isempty,boxcoords); % true for empty cells
    boxcoords(indx) = {0}; % replace by a cell with a zero
    
    Geometry.xyz = Coords;
    Geometry.boxcoords = boxcoords;
    Geometry.res_number = grodat.res_number;
    Geometry.res_name = strtrim(grodat.res_name);
    Geometry.atom_name = strtrim(grodat.atom_name);
    Geometry.atom_number = grodat.atom_number;
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

Geometry.N = N;
Geometry.NF = Geometry.N/2;
Geometry.FC = Coords/TM;
Geometry.Salt = Salt;

% Re-save the g96 file to get rid of the POSITIONRED format
if strcmpi(CoordType,'g96')
    SaveG96File(Filename,Geometry)
end

end
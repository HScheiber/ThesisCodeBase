function ndx_add = add_polarization_shells(Settings,filename,varargin)

p = inputParser;
p.FunctionName = 'add_polarization_shells';
addOptional(p,'ndx_filename','')
addOptional(p,'make_index',true)
addOptional(p,'add_shells',true)
parse(p,varargin{:});
ndx_filename = p.Results.ndx_filename; % Gathers components of energy at the end in addition to the total energy.
make_index = p.Results.make_index;
add_shells = p.Results.add_shells;

if ~Settings.Polarization
    ndx_add = '';
    return
end

% Load gro file
Geometry = load_gro_file(filename);

if add_shells
    % Add in shells
    Geometry.res_number = repelem(Geometry.res_number,2);
    Geometry.res_name = repelem(Geometry.res_name,2);
    Geometry.atom_name = repelem(Geometry.atom_name,2);
    Geometry.atom_name(2:2:end) = cellfun(@(x) strtrim([x '_s']),Geometry.atom_name(2:2:end),'UniformOutput',false);
    Geometry.N_atoms = Geometry.N_atoms*2;
    Geometry.atom_number = int32(transpose(1:Geometry.N_atoms));
    Geometry.xyz = repelem(Geometry.xyz,2,1);
    Geometry.vel = repelem(Geometry.vel,2,1);
    Geometry.Salt = Settings.Salt;
    Geometry.N = Geometry.N_atoms;

    % Save file
    SaveGroFile(filename,Geometry,true);
end

% Generate index file
% Create an index file to keep track of frozen atom numbers
system_atnums = double(Geometry.atom_number);
[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
energy_groups = {Metal [Metal '_s'] Halide [Halide '_s']};

% Create index file
addnan = 15 - mod(length(system_atnums),15);
system_atnums(end+1:end+addnan) = nan;
system_atnums = reshape(system_atnums,15,[])';
system_atnum_txt = char(regexprep(strjoin(string(num2str(system_atnums)),newline),{' *NaN' '\n\n'},{'' '\n'}));
ndx_text = ['[ System ]' newline system_atnum_txt newline];
ndx_text = [ndx_text '[ Other ]' newline system_atnum_txt newline];

% Get indeces for atom types
for idx = 1:numel(energy_groups)
    Grp_idx = find(strcmp(strtrim(Geometry.atom_name),energy_groups{idx}));
    Grp_atnums = double(Geometry.atom_number(Grp_idx));
    
    addnan = 15 - mod(length(Grp_atnums),15);
    Grp_atnums(end+1:end+addnan) = nan;
    Grp_atnums = reshape(Grp_atnums,15,[])';
    Grp_atnums_txt = char(regexprep(strjoin(string(num2str(Grp_atnums)),newline),{' *NaN' '\n\n'},{'' '\n'}));
    ndx_text = [ndx_text '[ ' energy_groups{idx} ' ]' newline Grp_atnums_txt newline]; %#ok<AGROW>
end

% Save index file
if make_index
    if isempty(ndx_filename)
        ndx_filename = fullfile(Settings.WorkDir,[Settings.FileBase '.ndx']);
    end
    fidNDX = fopen(ndx_filename,'wt');
    fwrite(fidNDX,regexprep(ndx_text,'\r',''));
    fclose(fidNDX);

    ndx_add = [' -n ' windows2unix(ndx_filename)];
else
    ndx_add = '';
end

end
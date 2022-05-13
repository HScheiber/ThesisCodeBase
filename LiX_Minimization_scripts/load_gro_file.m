function output = load_gro_file(filename)
[~,~,ext] = fileparts(filename);

if strcmp(ext,'.gro')
    gro_text = fileread(filename);
    gro_text = regexprep(gro_text,'(\r|\n)+','\n');
    comment_line = regexp(gro_text,'(.+?)\n','once','tokens');
    
    % Get number of molecules line
    Num_Mols = textscan(gro_text,'%f',1,'Delimiter','','whitespace','','headerlines',1);
    N = Num_Mols{1};
    
    % Grab position data
    Position_data = textscan(gro_text,'%5f%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',N,...
        'Delimiter','','whitespace','','headerlines',2);

    % Get box data
    latpar = textscan(gro_text, '%f %f %f %f %f %f %f %f %f',1, 'delimiter', ' ',...
        'whitespace','','MultipledelimsAsOne',true,'headerlines',2+N);
    
    output.res_number = Position_data{1};
    output.res_name   = Position_data{2};
    output.atom_name  = Position_data{3};
    output.atom_number = Position_data{4};
    output.xyz = [Position_data{5} Position_data{6} Position_data{7}];
    output.vel = [Position_data{8} Position_data{9} Position_data{10}];
    output.N_atoms = N;
    output.comment = comment_line{1};
    
    % v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    latpar(cellfun(@isempty,latpar)) = {0};
    latpar(cellfun(@isnan,latpar)) = {0};
    output.a_vec = [latpar{1} latpar{4} latpar{5}];
    output.b_vec = [latpar{6} latpar{2} latpar{7}];
    output.c_vec = [latpar{8} latpar{9} latpar{3}];
    output.Volume = Cell_Volume(output.a_vec,output.b_vec,output.c_vec);
    output.boxcoords = latpar;
    
    % Fix atom numbers for atoms over 99,999
    N_add = int32(floor((1:length(output.atom_number))./1e5)');
    output.atom_number = output.atom_number + N_add.*1e5;
    
elseif strcmp(ext,'.g96')
    g96_text = fileread(filename);
    g96_text = regexprep(g96_text,'(\r|\n)+','\n');
    
    Postxt = regexp(g96_text,'POSITIONR*E*D*\n(.+?)\nEND','tokens','once');
    Veltxt = regexp(g96_text,'VELOCITYR*E*D*\n(.+?)\nEND','tokens','once');
    Boxtxt = regexp(g96_text,'BOX\r*\n(.+?)\nEND','tokens','once');
    
    Position_data = textscan(Postxt{1},'%6f%6s%6s%6f%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    N = size(Position_data{1},1);
    
    if ~isempty(Veltxt)
        Velocity_data = textscan(Veltxt{1},'%*6f%*6s%*6s%*6f%15.9f%15.9f%15.9f\n',...
            'Delimiter','','whitespace','');
    else
        Velocity_data{1} = zeros(N,1);
        Velocity_data{2} = zeros(N,1);
        Velocity_data{3} = zeros(N,1);
    end
    
    output.res_number = Position_data{1};
    output.res_name   = strtrim(Position_data{2});
    output.atom_name  = strtrim(Position_data{3});
    output.atom_number = Position_data{4};
    output.xyz = [Position_data{5} Position_data{6} Position_data{7}];
    output.vel = [Velocity_data{1} Velocity_data{2} Velocity_data{3}];
    output.N_atoms = N;
    
    
    latpar = textscan(strtrim(Boxtxt{1}), '%f %f %f %f %f %f %f %f %f', 'delimiter', ' ',...
        'whitespace','','MultipledelimsAsOne',true);
    latpar(cellfun(@isempty,latpar)) = {0};
    latpar(cellfun(@isnan,latpar)) = {0};
    output.a_vec = [latpar{1}(end) latpar{4}(end) latpar{5}(end)];
    output.b_vec = [latpar{6}(end) latpar{2}(end) latpar{7}(end)];
    output.c_vec = [latpar{8}(end) latpar{9}(end) latpar{3}(end)];
    output.Volume = Cell_Volume(output.a_vec,output.b_vec,output.c_vec);
    output.boxcoords = latpar;
    output.comment = '';
end
end
% Input a CoordType = 'gro' or 'g96' structure file (Settings.SuperCellFile) and build a cluster of Cluster_N atoms
function Build_Cluster(Settings)

% Structural info
[~,~,CoordType] = fileparts(Settings.SuperCellFile);
GeomText = fileread(Settings.SuperCellFile);
GeomText = regexprep(GeomText,'(\r|\n)+','\n');

if strcmpi(CoordType,'.gro')
    % Get number of molecules line
    Num_Atoms = textscan(GeomText,'%f',1,'Delimiter','','whitespace','','headerlines',1);
    N = Num_Atoms{1};
    
    Coord_Data = textscan(GeomText,'%5f%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n',N,...
        'Delimiter','','whitespace','','headerlines',2);
    
    M_Ind = cellfun(@(x) strcmpi(x,pad(Settings.Metal,5,'right')),Coord_Data{2});
    X_Ind = cellfun(@(x) strcmpi(x,pad(Settings.Halide,5,'right')),Coord_Data{2});
    
    Coords = [Coord_Data{5}(1:N) Coord_Data{6}(1:N) Coord_Data{7}(1:N)];
    Met_XYZ = Coords(M_Ind,:);
    Hal_XYZ = Coords(X_Ind,:);
    
    boxcoords = textscan(GeomText, '%f %f %f %f %f %f %f %f %f',1, 'delimiter', ' ',...
        'whitespace','','MultipledelimsAsOne',true,'headerlines',2+N);

    % Check for empty coords
    boxcoords(cellfun(@isempty,boxcoords)) = {0};
    boxcoords(cellfun(@isnan,boxcoords)) = {0};
    
    % Lattice vectors
    a_vec = [boxcoords{1} boxcoords{4} boxcoords{5}]; % units nm
    b_vec = [boxcoords{6} boxcoords{2} boxcoords{7}]; % units nm
    c_vec = [boxcoords{8} boxcoords{9} boxcoords{3}]; % units nm
    
elseif strcmpi(CoordType,'.g96')
    
    Postxt = regexp(GeomText,'POSITIONR*E*D*\n(.+?)\nEND','tokens','once');
    Boxtxt = regexp(GeomText,'BOX\r*\n(.+?)\nEND','tokens','once');
    
    Position_data = textscan(Postxt{1},'%6f%6s%6s%6f%15.9f%15.9f%15.9f\n',...
        'Delimiter','','whitespace','');
    
    Coords = [Position_data{5} Position_data{6} Position_data{7}];
    
    M_Ind = cellfun(@(x) strcmpi(x,pad(Settings.Metal,6,'right')),Position_data{3});
    X_Ind = cellfun(@(x) strcmpi(x,pad(Settings.Halide,6,'right')),Position_data{3});
    Met_XYZ = Coords(M_Ind,:);
    Hal_XYZ = Coords(X_Ind,:);
    
    boxcoords = textscan(strtrim(Boxtxt{1}), '%f %f %f %f %f %f %f %f %f', 'delimiter', ' ',...
        'whitespace','','MultipledelimsAsOne',true);
    boxcoords(cellfun(@isempty,boxcoords)) = {0};
    boxcoords(cellfun(@isnan,boxcoords)) = {0};
    a_vec = [boxcoords{1}(end) boxcoords{4}(end) boxcoords{5}(end)]; % units nm
    b_vec = [boxcoords{6}(end) boxcoords{2}(end) boxcoords{7}(end)]; % units nm
    c_vec = [boxcoords{8}(end) boxcoords{9}(end) boxcoords{3}(end)]; % units nm
    
end

Box_Center = (1/2).*(a_vec + b_vec + c_vec);

% Go through all coordinates and calculate distances from box center
Met_r = vecnorm(Box_Center - Met_XYZ,2,2);
Hal_r = vecnorm(Box_Center - Hal_XYZ,2,2);

% Sort coordinates by distance
[~,M_idx] = sort(Met_r);
Met_XYZ_Sorted = Met_XYZ(M_idx,:);
[~,X_idx] = sort(Hal_r);
Hal_XYZ_Sorted = Hal_XYZ(X_idx,:);

% Take first N/2 atoms from each and recenter
Met_XYZ_Clus = Met_XYZ_Sorted(1:Settings.Cluster_N/2,:) - Box_Center;
Hal_XYZ_Clus = Hal_XYZ_Sorted(1:Settings.Cluster_N/2,:) - Box_Center;

% Create new square box with the same volume and recenter the coordinates
if strcmpi(Settings.pbc,'on')
    if Settings.Manual_Box
        Volume = Settings.Vol_Total; % Update the total volume
    else
        Volume = Cell_Volume(a_vec,b_vec,c_vec);
    end
    Newbox_L = Volume^(1/3); % nm
    Newbox_Center = [Newbox_L Newbox_L Newbox_L]./2;
    Met_XYZ_Clus = Met_XYZ_Clus + Newbox_Center;
    Hal_XYZ_Clus = Hal_XYZ_Clus + Newbox_Center;
    
    boxcoords = [Newbox_L Newbox_L Newbox_L 0 0 0 0 0 0];  % units nm
else
    boxcoords = [0 0 0 0 0 0 0 0 0];
end

if strcmpi(CoordType,'.gro')

    Out_text = [Settings.Metal Settings.Halide ' ' num2str(Settings.Cluster_N) ' atom salt cluster' newline ...
                num2str(Settings.Cluster_N)];
    idx = 1;
    for i = 1:Settings.Cluster_N/2
        Out_text = [Out_text newline sprintf('%5d%-5s%5s%5d%8.3f%8.3f%8.3f',mod(idx,1e5),Settings.Metal,...
            Settings.Metal,mod(idx,1e5),Met_XYZ_Clus(i,:)) newline ...
            sprintf('%5d%-5s%5s%5d%8.3f%8.3f%8.3f',mod(idx+1,1e5),Settings.Halide,...
            Settings.Halide,mod(idx+1,1e5),Hal_XYZ_Clus(i,:))];
        idx = idx+2;
    end
    Out_text = [Out_text newline sprintf('%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f',boxcoords) newline];
elseif strcmpi(CoordType,'.g96')
    Out_text = ['TITLE' newline Settings.Metal Settings.Halide ' ' num2str(Settings.Cluster_N) ...
        ' atom salt cluster' newline 'END' newline 'POSITION'];
    
    idx = 1;
    for i = 1:Settings.Cluster_N/2
        Out_text = [Out_text newline sprintf('%6d%6s%6s%6d%15.9f%15.9f%15.9f',idx,Settings.Metal,...
            Settings.Metal,idx,Met_XYZ_Clus(i,:)) newline ...
            sprintf('%6d%6s%6s%6d%15.9f%15.9f%15.9f',idx+1,Settings.Halide,...
            Settings.Halide,idx+1,Hal_XYZ_Clus(i,:))];
        idx = idx+2;
    end
    if strcmp(Settings.pbc,'on')
        Out_text = [Out_text newline 'END' newline 'BOX' ...
            sprintf('%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f',boxcoords) ...
            newline 'END'];
    else
        Out_text = [Out_text newline 'END'];
    end
    
end

% Remove the old supercell file and replace with a new one
delete(Settings.SuperCellFile)
fid = fopen(Settings.SuperCellFile,'wt');
fprintf(fid,Out_text);
fclose(fid);

end
function Create_Liquid_Box(WorkDir)

Dat = load(fullfile(WorkDir,'TempJobInfo.mat'));
Settings = Dat.Settings;

% Sanity checks
if Settings.GenCluster && abs(Settings.c_over_a - 1) > sqrt(eps)
    Settings.c_over_a = 1;
    disp('C/A ratio reset to 1 for liquid cluster generation')
end

[Metal,Halide] = Separate_Metal_Halide(Settings.Salt);
Ref_M = fullfile(Settings.home,'templates','GRO_Templates',[Metal '_Box.gro']);
Ref_X = fullfile(Settings.home,'templates','GRO_Templates',[Halide '_Box.gro']);
nmol = Settings.N_atoms/2; % number of atoms needed

% Find the target density based on the temperature and pressure
warning('off','MATLAB:UndefinedFunction')
Ref_Density = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
warning('on','MATLAB:UndefinedFunction')

% Calculate the dimensions for the box at the given density
Vol = nmol/Ref_Density; % Volume in cubic nm

Settings.Geometry.a = ((Vol/Settings.c_over_a)^(1/3)); % nm
Settings.Geometry.b = Settings.Geometry.a; % nm
Settings.Geometry.c = Settings.c_over_a*Settings.Geometry.a; % nm

a_vec = Settings.Geometry.Transform(1,:).*Settings.Geometry.a.*Settings.Expand_a_SC; % nm
b_vec = Settings.Geometry.Transform(2,:).*Settings.Geometry.b.*Settings.Expand_b_SC; % nm
c_vec = Settings.Geometry.Transform(3,:).*Settings.Geometry.c.*Settings.Expand_c_SC; % nm

% Check the box vectors to make sure the cutoff isn't too long
LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                    Settings.Geometry.Skew_b*norm(b_vec) ...
                    Settings.Geometry.Skew_c*norm(c_vec)]);

% Calculate the largest cutoff distance
Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]); % nm

if LatticeLength/2 <= Longest_Cutoff*Settings.Cutoff_Buffer
    old_atnum = Settings.N_atoms;
    
    New_a = 2*Longest_Cutoff*Settings.Cutoff_Buffer;
    New_c_tot = 2*Longest_Cutoff*Settings.Cutoff_Buffer;
    if Settings.c_over_a > New_c_tot/New_a
        % Expand c to maintain the user-selected c/a ratio
        New_c_tot = New_a*Settings.c_over_a;
    elseif Settings.c_over_a < New_c_tot/New_a
        % Expand a to maintain the user-selected c/a ratio
        New_a = New_c_tot/Settings.c_over_a;
    end
    New_b = New_a;
    
    if Settings.Liquid_Interface
        New_c = New_c_tot*(1-Settings.Liquid_Fraction);
    else
        New_c = New_c_tot;
    end
    
    % Update number of atoms to match the new box size
    Vol = New_a*New_b*New_c_tot;
    nmol = ceil(Ref_Density*Vol);
    Settings.N_atoms = nmol*2;
    
    % Update lattice vectors
    a_vec = Settings.Geometry.Transform(1,:).*New_a.*Settings.Expand_a_SC; % nm
    b_vec = Settings.Geometry.Transform(2,:).*New_b.*Settings.Expand_b_SC; % nm
    c_vec = Settings.Geometry.Transform(3,:).*New_c.*Settings.Expand_c_SC; % nm
    
    disp(['Warning: With ' num2str(old_atnum) ...
        ' atoms, the cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element.'])
    disp(['Expanding the box to ' num2str(Settings.N_atoms) ' atoms.'])
end

% Create an empty box
tmp_empty_box = fullfile(WorkDir,['tmp_empty_box.' Settings.CoordType]);
Settings.Geometry.N = 0;
Settings.Geometry.boxcoords = {a_vec(1) b_vec(2) c_vec(3) a_vec(2) a_vec(3) b_vec(1) b_vec(3) c_vec(1) c_vec(2)};
SaveGroFile(tmp_empty_box,Settings.Geometry,true);

%% Add atoms to the file
switch Settings.Salt
    case 'LiF'
        r0 = 0.1;
    case 'LiCl'
        r0 = 0.11;
    case 'LiBr'
        r0 = 0.12;
    case 'LiI'
        r0 = 0.13;
    otherwise
        r0 = 0.12;
end

% If making a liquid cluster (bubble)
if Settings.GenCluster
    
    Settings.Cluster_N = ceil(Settings.N_atoms*Settings.Cluster_Frac);
    if mod(Settings.Cluster_N,2) > 0
        Settings.Cluster_N = Settings.Cluster_N + 1; % Ensure it's even
    end
    
    nmol_cluster = Settings.Cluster_N/2;
    Box_Center = norm(a_vec).*[1/2 1/2 1/2];  
    
    % Calculate the radius needed for the given density
    Vol_cluster = nmol_cluster/Ref_Density; % Volume in cubic nm
    R_cluster = ( (3/(4*pi))*Vol_cluster )^(1/3); % Radius in nm
    Inner_Box_Length = 2*(R_cluster);

    %% How many atoms should fit into a box with the same side length as the sphere?
    Box_Vol = Inner_Box_Length^3; % nm^3 volume
    nmol_box = ceil(Ref_Density*Box_Vol);

    % Create a positions.dat file
    pos_filename = fullfile(WorkDir,'positions.dat');
    fid = fopen(pos_filename,'wt');
    fprintf(fid, '# Box\n');
    for i = 1:nmol_box
        fprintf(fid, num2str(Box_Center));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %% Add atoms to the file
    % Add metal
    disp(['Randomly adding ' num2str(nmol_box) ' ' Metal ' ions to liquid box...'])
    mtimer = tic;
    tmp_metal_only_file = fullfile(WorkDir,['tmp_metal_only.' Settings.CoordType]);
    cmd = [Settings.gmx_loc ' insert-molecules -f ' windows2unix(tmp_empty_box) ' -ci ' windows2unix(Ref_M) ...
        ' -o ' windows2unix(tmp_metal_only_file) ' -nmol ' num2str(nmol_box) ' -radius ' num2str(r0) ' -try 200' ...
        ' -dr ' num2str(R_cluster) ' -ip ' windows2unix(pos_filename)];
    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Metal ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Metal ' atoms added. Epalsed Time: ' datestr(seconds(toc(mtimer)),'HH:MM:SS')])
    
    % Add Halide
    tmp_liquid_file = fullfile(WorkDir,['tmp_liquid.' Settings.CoordType]);
    disp(['Randomly adding ' num2str(nmol_box) ' ' Halide ' ions to liquid box...'])
    htimer = tic;
    cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(tmp_metal_only_file) ...
        ' -o ' windows2unix(tmp_liquid_file) ' -nmol ' num2str(nmol_box) ' -radius ' num2str(r0) ' -dr ' num2str(R_cluster) ...
        ' -try 400  -ip ' windows2unix(pos_filename)];

    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Halide ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Halide ' atoms added. Epalsed Time: ' datestr(seconds(toc(htimer)),'HH:MM:SS')])


    %% Carve box into a sphere
    [atom,data] = import_atom_gro(tmp_liquid_file);

    % split atom into metal and halide data
    m_idx = 1;
    x_idx = 1;
    for idx = 1:data.nAtoms
        if strcmp(atom(idx).resname{1},Metal)
            atom_M(m_idx) = atom(idx);
            m_idx = m_idx+1;
        else
            atom_X(x_idx) = atom(idx);
            x_idx = x_idx+1;
        end
    end
    
    xyz_coords_M = [atom_M.x; atom_M.y; atom_M.z]';
    xyz_coords_X = [atom_X.x; atom_X.y; atom_X.z]';
    R0_M = vecnorm(xyz_coords_M-10.*Box_Center,2,2); % Distances in Angstrom from box center
    R0_X = vecnorm(xyz_coords_X-10.*Box_Center,2,2); % Distances in Angstrom from box center

    [~,idb_M] = sort(R0_M);
    [~,idb_X] = sort(R0_X);
    idx_Sorted_M = idb_M(1:nmol_cluster);
    idx_Sorted_X = idb_X(1:nmol_cluster);

    i = 1;
    for idx = 1:length(idx_Sorted_M)
        idq = idx_Sorted_M(idx);
        atom_fin(i) = atom_M(idq);
        atom_fin(i).index = i;
        atom_fin(i).molid = 0;
        i = i+1;
    end
    for idx = 1:length(idx_Sorted_X)
        idq = idx_Sorted_X(idx);
        atom_fin(i) = atom_X(idq);
        atom_fin(i).index = i;
        atom_fin(i).molid = 0;
        i = i+1;
    end

    % Save
    write_atom_gro(atom_fin,data.Box_dim,tmp_liquid_file)
    
else

    % Add metal
    disp(['Randomly adding ' num2str(nmol) ' ' Metal ' ions to liquid box...'])
    mtimer = tic;
    tmp_metal_only_file = fullfile(WorkDir,['tmp_metal_only.' Settings.CoordType]);
    cmd = [Settings.gmx_loc ' insert-molecules -f ' windows2unix(tmp_empty_box) ' -ci ' windows2unix(Ref_M) ...
        ' -o ' windows2unix(tmp_metal_only_file) ' -nmol ' num2str(nmol) ' -radius ' num2str(r0) ' -try 200'];
    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Metal ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Metal ' atoms added. Epalsed Time: ' datestr(seconds(toc(mtimer)),'HH:MM:SS')])

    % Add Halide
    tmp_liquid_file = fullfile(WorkDir,['tmp_liquid.' Settings.CoordType]);
    disp(['Randomly adding ' num2str(nmol) ' ' Halide ' ions to liquid box...'])
    htimer = tic;
    cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(tmp_metal_only_file) ...
        ' -o ' windows2unix(tmp_liquid_file) ' -nmol ' num2str(nmol) ' -radius ' num2str(r0) ' -try 400'];

    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Halide ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Halide ' atoms added. Epalsed Time: ' datestr(seconds(toc(htimer)),'HH:MM:SS')])

    % Load current randomly-generated liquid file data
    Liquid_file_data = load_gro_file(tmp_liquid_file);

    % Check to make sure all atoms were added
    if ~Liquid_file_data.N_atoms == nmol*2
        error('Not all requested liquid atoms were added!')
    end

end

%% Minimize the randomly-generated liquid: build MDP file
MDP_Minimization_txt = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
'MDP.template'));

MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##NSTEPS##',pad(num2str(Settings.MinMDP.nsteps_min),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##INTEGR##',pad(Settings.MinMDP.min_integrator,18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'dt                       = ##TIMEST##; Time step (ps)',...
    ['emtol                    = ' num2str(Settings.MinMDP.emtol)]);
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##MET##',Metal);
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##HAL##',Halide);
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##COULOMB##',pad(num2str(Settings.MDP.CoulombType),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##FOURIER##',pad(num2str(Settings.MDP.Fourier_Spacing),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##PMEORDER##',pad(num2str(Settings.MDP.PME_Order),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##EWALDTOL##',pad(num2str(Settings.MDP.Ewald_rtol),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));

if Settings.Table_Req || strcmp(Settings.Theory,'BH')
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##VDWTYPE##',pad('user',18));
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##CUTOFF##',pad('group',18));
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'ewald-rtol-lj.+?\n','');
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'lj-pme-comb-rule.+?\n','');
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'verlet-buffer-tolerance.+?\n','');
    
    % For minimization, add in a close-range repulsive wall to the
    % potential with the following function
    Settings.Geometry = Settings.Geometry;
    Settings.WorkDir = WorkDir;
    Settings.TableFile_MX = MakeTables(Settings);
else
    % Modify the MDP file
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
    MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##CUTOFF##',pad(Settings.MDP.CutOffScheme,18));
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'energygrp-table.+?\n','');
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'ewald-rtol-lj.+?\n','');
    MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'lj-pme-comb-rule.+?\n','');

    % Add in Verlet Settings
    if strcmp(Settings.MDP.CutOffScheme,'Verlet')
        MDP_Minimization_txt = strrep(MDP_Minimization_txt,'##VerletBT##',pad(num2str(Settings.MDP.VerletBT),18));
    else
        MDP_Minimization_txt = regexprep(MDP_Minimization_txt,'verlet-buffer-tolerance.+?\n','');
    end
end

% Add in dispersion corrections
if Settings.MDP.Disp_Correction && ~(Settings.Table_Req || strcmp(Settings.Theory,'BH'))
    MDP_Minimization_txt = [MDP_Minimization_txt newline newline...
        '; Long-range dispersion correction' newline ...
        'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
elseif Settings.MDP.Disp_Correction && Settings.MDP.Disp_Correction_Tables
    MDP_Minimization_txt = [MDP_Minimization_txt newline newline ...
        '; Long-range dispersion correction' newline ...
        'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
end

% Save MDP file
MDP_Filename = fullfile(WorkDir,'tmp_liquid.mdp');
fidMDP = fopen(MDP_Filename,'wt');
fwrite(fidMDP,regexprep(MDP_Minimization_txt,'\r',''));
fclose(fidMDP);

% Complete a topology file for the liquid box to be minimized
Atomlist = copy_atom_order(tmp_liquid_file);
Top_Filename = fullfile(WorkDir,'tmp_liquid.top');

% Buckingham potentials require recreating the topology text
if strcmp(Settings.Theory,'BH')

    % Load Topology template
    Topology_Text = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
    'Topology.template'));

    % Add in global parameters to Topology template
    Topology_Text = strrep(Topology_Text,'##GENPAIRS##',Settings.Top_gen_pairs);
    Topology_Text = strrep(Topology_Text,'##FUDGELJ##',num2str(Settings.Top_fudgeLJ));
    Topology_Text = strrep(Topology_Text,'##FUDGEQQ##',num2str(Settings.Top_fudgeQQ));

    Metal_Info = elements('Sym',Settings.Metal);
    Halide_Info = elements('Sym',Settings.Halide);

    % Insert element info into topology template
    Topology_Text = strrep(Topology_Text,'##MET##',pad(Settings.Metal,2));
    Topology_Text = strrep(Topology_Text,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
    Topology_Text = strrep(Topology_Text,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));
    Topology_Text = strrep(Topology_Text,'##MCHRG##',pad(num2str(Settings.S.Q),2));

    Topology_Text = strrep(Topology_Text,'##HAL##',pad(Settings.Halide,2));
    Topology_Text = strrep(Topology_Text,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
    Topology_Text = strrep(Topology_Text,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));
    Topology_Text = strrep(Topology_Text,'##XCHRG##',pad(num2str(-Settings.S.Q),2));

    % Add number of unit cells to topology file
    Topology_Text = strrep(Topology_Text,'##N##x##N##x##N##',num2str(Settings.nmol_liquid));
    Topology_Text = strrep(Topology_Text,'##GEOM##','molecule liquid');

    % Define the function type as 1 (needed for tabulated functions)
    Topology_Text = strrep(Topology_Text,'##NBFUNC##','1');

    % Define the combination rules (Lorenz-berthelot)
    Topology_Text = strrep(Topology_Text,'##COMBR##','1');

    % Define all the parameters as 1.0 (already included in potentials)
    Topology_Text = strrep(Topology_Text,'##METMETC##',pad('1.0',10));
    Topology_Text = strrep(Topology_Text,'##HALHALC##',pad('1.0',10));
    Topology_Text = strrep(Topology_Text,'##METHALC##',pad('1.0',10));
    Topology_Text = strrep(Topology_Text,'##METMETA##','1.0');
    Topology_Text = strrep(Topology_Text,'##HALHALA##','1.0');
    Topology_Text = strrep(Topology_Text,'##METHALA##','1.0');
end
Topology_Text = strrep(Topology_Text,'##LATOMS##',Atomlist);
Topology_Text = strrep(Topology_Text,'##N##x##N##x##N## ##GEOM##',...
    ['Liquid with ' num2str(Settings.N_atoms) ' atoms.']);

% Save topology file
fidTOP = fopen(Top_Filename,'wt');
fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
fclose(fidTOP);

TPR_File = fullfile(WorkDir,'tmp_liquid.tpr');
MDPout_File = fullfile(WorkDir,'tmp_liquid_out.mdp');
GrompLog_File = fullfile(WorkDir,'tmp_liquid_Grompplog.log');

FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(tmp_liquid_file) ...
    ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Top_Filename) ...
    ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
    ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GrompLog_File)];
[state,~] = system(FMin_Grompp);
% Catch error in grompp
if state ~= 0
    error(['Error running GROMPP. Problem command: ' newline FMin_Grompp]);
else
    delete(GrompLog_File)
end

% Prepare minimization mdrun command
Log_File = fullfile(WorkDir,'tmp_liquid.log');
Energy_file = fullfile(WorkDir,'tmp_liquid.edr');
TRR_File = fullfile(WorkDir,'tmp_liquid.trr');

mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
    ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
    ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Settings.SuperCellFile) ...
    Settings.mdrun_opts];

if Settings.Table_Req || strcmp(Settings.Theory,'BH')
    mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
end

% Final minimization
disp('Begining liquid minimization...')
mintimer = tic;
[state,mdrun_output] = system(mdrun_command);
if state == 0
    disp(['System Successfully Minimized! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
else
    disp('Minimization failed. Stopping.')
    disp(mdrun_output);
    error(['Error running mdrun for system minimization. Problem command: ' newline mdrun_command]);
end

% Generate final topology file for molecular dynamics
Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
Settings.Topology_Text = strrep(Settings.Topology_Text,'##N##x##N##x##N## ##GEOM##',...
    ['Liquid with ' num2str(Settings.N_atoms) ' atoms.']);
fidTOP = fopen(Settings.Topology_File,'wt');
fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
fclose(fidTOP);

GROMPP_command = [Settings.gmx_loc ' grompp -c ' windows2unix(Settings.SuperCellFile) ...
    ' -f ' windows2unix(Settings.MDP_in_File) ' -p ' windows2unix(Settings.Topology_File) ...
    ' -o ' windows2unix(Settings.Traj_Conf_File) ' -po ' windows2unix(Settings.MDP_out_File) ...
    ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(Settings.GrompLog_File)];
[errcode,~] = system(GROMPP_command);

% Catch error in grompp
if errcode ~= 0
    error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
end

% Remove minimization temporary folder
if Settings.Delete_Backups
    system([Settings.wsl 'find ' windows2unix(WorkDir) ...
        ' -iname "#*#" -delete']);
end
if Settings.Delete_Minimize
    system([Settings.wsl 'rm -r ' windows2unix(WorkDir) '/*']);
    rmdir(WorkDir);
end

end
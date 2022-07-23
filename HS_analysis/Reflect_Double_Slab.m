filename_in = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\slab_doubling\Sapphire_slab_001_trimmed_tight_boxed.gro';
filename_out = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\slab_doubling\Sapphire_slab_001_doubled_boxed_processed.gro';
gmx = 'wsl source ~/.bashrc; gmx_d';
nmol = 1250; % number of LiI molecules to insert
tmp_file = 'tmp.gro';

%% Import data
[atom,data] = import_atom_gro(filename_in);

% Copy data
atom_copy = atom;

%% Reflect structure in the z dimension
N_A = 6.0221409e23; % molecules/(mol of ion pair)
cm3_per_Ang3 = (1e-8)^3; % cm^3/A^3
Ref_Density_SI = 1/58.3; % mol/cm^3. Reference comes from equilibrated LiI system with modelA, at 2000 K
Ref_Density = Ref_Density_SI*N_A*cm3_per_Ang3; % molecules/A^3
box_area = data.Box_dim(1)*data.Box_dim(2); % A^2

vaccuum_height = nmol/(box_area*Ref_Density);

for idx = 1:data.nAtoms
    atom_copy(idx).z = -atom(idx).z + data.Box_dim(3)*2 + vaccuum_height; % Angstroms
end

%% Combine
atom = [atom atom_copy];

%% Modify box size
data.Box_dim(3) = data.Box_dim(3)*2 + vaccuum_height;

%% Save
write_atom_gro(atom,data.Box_dim,tmp_file)

%% Process with gromacs and remove the temp file
cmd = [gmx ' editconf -f ' windows2unix(tmp_file) ' -o ' windows2unix(filename_out)];
system(cmd);
delete(tmp_file) 

%% Add molecules with gromacs
M_Source_File = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\slab_doubling\Li_Box.gro';
X_Source_File = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\slab_doubling\I_Box.gro';
MX_out_File = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\slab_doubling\LiI_L_Sapphire_slab_001_No_Vaccuum.gro';

% Add metal
cmd = [gmx ' insert-molecules -ci ' windows2unix(M_Source_File) ' -f ' windows2unix(filename_out) ...
    ' -o ' windows2unix(tmp_file) ' -nmol ' num2str(nmol) ' -scale 1'];
system(cmd);

% Add Halide
cmd = [gmx ' insert-molecules -ci ' windows2unix(X_Source_File) ' -f ' windows2unix(tmp_file) ...
    ' -o ' windows2unix(MX_out_File) ' -nmol ' num2str(nmol) ' -scale 1'];
system(cmd,'-echo');
delete(tmp_file)


%% Add vaccuum spaces
[atom,data] = import_atom_gro(MX_out_File);
atom_copy = atom;
vacuum_space = 15; % Angstrom

for idx = 1:data.nAtoms
    element = regexp(atom(idx).resname{1},'(Al|O|Li|I)[0-9]*','once','tokens');
    if ~isempty( regexp(atom(idx).resname{1},'(Al|O)[0-9]+','once') ) % Al2O3 atoms
        if atom_copy(idx).z >= data.Box_dim(3)/2 % Al2O3 atoms in upper half of box
            
            atom_copy(idx).z = atom(idx).z + 2.5*vacuum_space; % Angstroms
        else % Al2O3 atoms in lower half of box
            
            atom_copy(idx).z = atom(idx).z + 0.5*vacuum_space; % Angstroms
        end
        atom_copy(idx).resname = {'AlO'};
        atom_copy(idx).molid = 0;
    else % LiI atoms
        atom_copy(idx).resname = {'LiI'};
        atom_copy(idx).molid = 1;
        atom_copy(idx).z = atom(idx).z + 1.5*vacuum_space; % Angstroms
    end
    atom_copy(idx).fftype = element;
    atom_copy(idx).type = element;
end

% Modify box size
data.Box_dim(3) = data.Box_dim(3) + 3*vacuum_space;

% Save
write_atom_gro(atom_copy,data.Box_dim,tmp_file)

% Process with gromacs and remove the temp file
Completed_Filename = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\LiI_L_Sapphire_slab_001.gro';
cmd = [gmx ' editconf -f ' windows2unix(tmp_file) ' -o ' windows2unix(Completed_Filename)];
system(cmd);
delete(tmp_file)

%% Make index file
index_file = 'C:\Users\Hayden\Documents\Patey_Lab\LiI_L_Sapphire_2_Slabs_001_JCmodelA\LiI_L_Sapphire_slab_001.ndx';
cmd = ['wsl source ~/.bashrc; echo -ne "case\na Al\na O\na Li\na I\nq\n" ^| gmx_d make_ndx -f ' ...
    windows2unix(Completed_Filename) ' -o ' windows2unix(index_file)];
system(cmd);

output = copy_atom_order(Completed_Filename);
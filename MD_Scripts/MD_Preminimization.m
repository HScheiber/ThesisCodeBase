function MD_Preminimization(Directory)


%% Move to input directory and load input variables
Settings = load(fullfile(Directory,'TempJobInfo.mat'));

% Place "Settings" sub-structure into the main input structure
if isfield(Settings,'Settings')
    f = fieldnames(Settings.Settings);
    for i = 1:length(f)
        Settings.(f{i}) = Settings.Settings.(f{i});
    end
    Settings = rmfield(Settings,'Settings');
end
Settings.WorkDir = Directory;

% Select symmetry settings
if strcmpi(Settings.Structure,'Liquid') || Settings.Found_DataMatch
    Run_Min = false;
else
    Run_Min = true;
end

if Run_Min
    OutputFile = fullfile(Settings.WorkDir,[Settings.JobName '_MinOut.mat']);
    
    if isfile(OutputFile)
        dat = load(OutputFile,'MinOut');
        MinOut = dat.MinOut;
    else
        % Save Minimization Diary
        DiaryFile = fullfile(Settings.WorkDir,[Settings.JobName '.minlog']);
        diary(DiaryFile);

        % Initial Structure Minimization
        disp(['*************** ' Settings.Salt ' ' Settings.Structure ' ' Settings.Full_Model_Name ' Optimization ****************'])
        disp('********************** Convergence Requirements ***********************')
        disp(['Max ' char(916) 'E between cycles: ' num2str(Settings.MinMDP.Energy_Tol,'%4.4E') ' a.u.']);
        disp(['Max RMS ' char(8711) 'E: ' num2str(Settings.MinMDP.Gradient_Tol_RMS,'%4.4E') ' a.u. / ' char(0197)]);
        disp(['Max component ' char(8711) 'E: ' num2str(Settings.MinMDP.Gradient_Tol_Max,'%4.4E') ' a.u. / ' char(0197)]);
        disp('***********************************************************************')
        MinOut = Structure_Minimization(Settings);
        diary off
        save(OutputFile,'MinOut')
    end
    
    % Gather output
    Settings.Geometry.a = MinOut.a;
    Settings.Geometry.b = MinOut.b;
    Settings.Geometry.c = MinOut.c;
    Settings.Geometry.FC_Metal = MinOut.FC_Metal;
    Settings.Geometry.FC_Halide = MinOut.FC_Halide;
    
    % Expand the volume
    if Settings.Expand_LP && Settings.Thermal_Solid
        Exp_Coeff = LiX_Thermal_Expansion(Settings);
    else
        Exp_Coeff = 1;
    end

    Settings.Geometry.a = Settings.Geometry.a*Settings.Expand_a*Exp_Coeff;
    Settings.Geometry.b = Settings.Geometry.b*Settings.Expand_b*Exp_Coeff;
    Settings.Geometry.c = Settings.Geometry.c*Settings.Expand_c*Exp_Coeff;
end

Longest_Cutoff = max([Settings.MDP.RList_Cutoff Settings.MDP.RCoulomb_Cutoff Settings.MDP.RVDW_Cutoff]);

% Load mdp template location
Settings.MDP.Temp_Model = fullfile(Settings.home,'templates','Gromacs_Templates',...
'MDP.template');

% Load MDP template
Settings.MDP.Temp_Model = fileread(Settings.MDP.Temp_Model);

% Add in global parameters to MDP template
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##NSTEPS##',pad(num2str(Settings.MinMDP.nsteps_point),18));
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##INTEGR##',pad(Settings.MinMDP.point_integrator,18));
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##TIMEST##',pad(num2str(Settings.MinMDP.dt),18));

% Determine type of optimization
if Settings.MinMDP.OptPos
    Settings.OptTxt = 'FULLOPT';
else
    Settings.OptTxt = 'CELLOPT';
end

% Insert salt components into MDP template
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##MET##',Settings.Metal);
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##HAL##',Settings.Halide);

% Update File Base Name
FileBase = [Settings.Salt '_' Settings.Geometry.Label '_' Settings.Full_Model_Name '_' Settings.OptTxt];

% Update Topology and MDP files
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##COULOMB##',pad(num2str(Settings.MDP.CoulombType),18));
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##FOURIER##',pad(num2str(Settings.MDP.Fourier_Spacing),18));
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##PMEORDER##',pad(num2str(Settings.MDP.PME_Order),18));
Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##EWALDTOL##',pad(num2str(Settings.MDP.Ewald_rtol),18));

if Settings.Table_Req || strncmp(Settings.Theory,'BH',2)
    
    % Modify the MDP file
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##VDWTYPE##',pad('user',18));
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##CUTOFF##',pad('group',18));
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'ewald-rtol-lj.+?\n','');
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'lj-pme-comb-rule.+?\n','');
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'verlet-buffer-tolerance.+?\n','');
    
else

    % Modify the MDP file
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##CUTOFF##',pad(Settings.MDP.CutOffScheme,18));
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'energygrp-table.+?\n','');
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'ewald-rtol-lj.+?\n','');
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'lj-pme-comb-rule.+?\n','');

    % Add in Verlet Settings
    if strcmp(Settings.MDP.CutOffScheme,'Verlet')
        Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##VerletBT##',pad(num2str(Settings.MDP.VerletBT),18));
    else
        Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'verlet-buffer-tolerance.+?\n','');
    end
end

if strcmp(Settings.pbc,'off')
    comm_txt = 'comm_mode                = Angular           ; Remove center of mass translational and rotational velocity around the center of mass';
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'(pbc += )xyz(.+?)\n',['$1no $2' newline comm_txt]);
end

ReRun_OptPos = false;
if Run_Min
    
    % Check to make sure supercell is large enough for the cutoffs by finding the shortest lattice parameter
    a_vec = Settings.Geometry.Transform(1,:).*Settings.Expand_a_SC*Settings.Geometry.a*Settings.N_Supercell_a/10; % nm
    b_vec = Settings.Geometry.Transform(2,:).*Settings.Expand_b_SC*Settings.Geometry.b*Settings.N_Supercell_b/10; % nm
    c_vec = Settings.Geometry.Transform(3,:).*Settings.Expand_c_SC*Settings.Geometry.c*Settings.N_Supercell_c_tot/10; % nm

    if Settings.GenCluster % When generating clusters, reset the box vectors to cubic
        L = Cell_Volume(a_vec,b_vec,c_vec)^(1/3); % nm^3
        a_vec = L.*[1 0 0];
        b_vec = L.*[0 1 0];
        c_vec = L.*[0 0 1];
    end
    LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                        Settings.Geometry.Skew_b*norm(b_vec) ...
                        Settings.Geometry.Skew_c*norm(c_vec)]);
    
    if LatticeLength/2 <= Longest_Cutoff*Settings.Cutoff_Buffer
        
        old_atnum = (Settings.N_Supercell_a*Settings.N_Supercell_b*Settings.N_Supercell_c_tot)*(Settings.Geometry.N);
        
        N_Supercell_a = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_a*Settings.Geometry.a/10));
        N_Supercell_b = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_b*Settings.Geometry.b/10));
        N_Supercell_c_tot = ceil(2*(Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_c*Settings.Geometry.c/10));
        
        a_current = N_Supercell_a*Settings.Geometry.a/10;
        c_current = N_Supercell_c_tot*Settings.Geometry.c/10;
        
        
        if Settings.c_over_a > c_current/a_current
            % Expand c to maintain the user-selected c/a ratio
            N_Supercell_c_tot = ceil((Settings.Geometry.a/Settings.Geometry.c)*N_Supercell_a*Settings.c_over_a);
        elseif Settings.c_over_a < c_current/a_current
            % Expand a to maintain the user-selected c/a ratio
            N_Supercell_a = ceil((Settings.Geometry.c/Settings.Geometry.a)*N_Supercell_c_tot/Settings.c_over_a);
            N_Supercell_b = ceil((Settings.Geometry.c/Settings.Geometry.b)*N_Supercell_c_tot/Settings.c_over_a);
        end
        
        if Settings.Liquid_Interface && Settings.GenCluster
            N_Supercell_c = round(N_Supercell_c_tot*Settings.Sol_fraction);
        else
            N_Supercell_c = N_Supercell_c_tot;
        end
        disp(['Warning: With ' num2str(old_atnum) ...
            ' atoms, the cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element.'])
        disp(['Expanding the box to ' num2str((N_Supercell_a*N_Supercell_b*N_Supercell_c_tot)*(Settings.Geometry.N)) ' atoms.'])
    else
        N_Supercell_a = Settings.N_Supercell_a;
        N_Supercell_b = Settings.N_Supercell_b;
        N_Supercell_c = Settings.N_Supercell_c;
        N_Supercell_c_tot = Settings.N_Supercell_c_tot;
    end
    
    % Add number of unit cells to topology file (this will not change)
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##N##x##N##x##N##',...
        [num2str(N_Supercell_a) 'x' num2str(N_Supercell_b) 'x' num2str(N_Supercell_c)]);

    % Save number of atoms into .mat file (this wont change)
    NumberFile = fullfile(Settings.WorkDir,[Settings.JobName '.mat']);
    N_Cell = Settings.Geometry.N;
    N_total = (N_Supercell_a*N_Supercell_b*N_Supercell_c_tot)*N_Cell;
    save(NumberFile,'N_total','N_Cell');
    Na = num2str(N_Supercell_a);
    Nb = num2str(N_Supercell_b);
    Nc = num2str(N_Supercell_c);

    %% Geometry Editing of minimized cell

    % Add unit cell coordinates
    Settings.Coordinate_Text = AddCartesianCoord(Settings.Coordinate_Text,Settings.Geometry,1,false,Settings.CoordType);

    % Save unit cell .gro file into main directory
    fid = fopen(Settings.UnitCellFile,'wt');
    fwrite(fid,regexprep(Settings.Coordinate_Text,'\r',''));
    fclose(fid);

    % Create supercell
    temp_SuperCellFile = fullfile(Settings.WorkDir,[Settings.JobName '.' Settings.CoordType]);
    Supercell_command = [Settings.gmx_loc ' genconf -f ' windows2unix(Settings.UnitCellFile) ...
         ' -o ' windows2unix(temp_SuperCellFile) ' -nbox ' Na ' ' Nb ' ' Nc];
    [errcode,output] = system(Supercell_command);

    if errcode ~= 0
        disp(output);
        error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
    end
    
    % build a cluster if requested
    if Settings.GenCluster
        Build_Cluster(Settings)
    end
    
    if Settings.Expand_a_SC ~= 1 || Settings.Expand_b_SC ~= 1 || Settings.Expand_c_SC ~= 1
        % Geometry Editing: expand supercell by requested amount
        a_sc = num2str(Settings.Expand_a_SC*Settings.Geometry.a*N_Supercell_a/10,'%10.8e'); % supercell a length in nm
        b_sc = num2str(Settings.Expand_b_SC*Settings.Geometry.b*N_Supercell_b/10,'%10.8e'); % supercell b length in nm
        c_sc = num2str(Settings.Expand_c_SC*Settings.Geometry.c*N_Supercell_c/10,'%10.8e'); % supercell c length in nm

        % Cell angles
        bc = num2str(Settings.Geometry.alpha,'%10.4e');
        ac = num2str(Settings.Geometry.beta,'%10.4e');
        ab = num2str(Settings.Geometry.gamma,'%10.4e');

        Expand_command = [Settings.gmx_loc ' editconf -f ' windows2unix(temp_SuperCellFile) ...
             ' -o ' windows2unix(temp_SuperCellFile) ' -box ' a_sc ' ' b_sc ' ' c_sc ' ' ...
             '-noc -angles ' bc ' ' ac ' ' ab];
        [errcode,output] = system(Expand_command);

        if errcode ~= 0
            disp(output);
            error(['Error expanding supercell with editconf. Problem command: ' newline Expand_command]);
        end
        if Settings.MinMDP.OptPos
            ReRun_OptPos = true;
        end
    end
    copyfile(temp_SuperCellFile,Settings.SuperCellFile);
else
    disp('Minimization of lattice parameters skipped.')
end

%% Gromacs final minimize after geometry editing
if Run_Min && ReRun_OptPos
    
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'= md                ; What type of calculation is run',...
        ['= ' Settings.MinMDP.min_integrator newline 'emtol                    = ' num2str(Settings.MinMDP.emtol)]);
    Settings.MDP.Temp_Model = regexprep(Settings.MDP.Temp_Model,'nsteps                   = [0-9|\.|\-]+',...
        ['nsteps                   = ' num2str(Settings.MinMDP.nsteps_min)]);

    % Determine cutoff length
    R_List_Cutoff = pad(num2str(Settings.MDP.RList_Cutoff),18);
    R_Coulomb_Cutoff = pad(num2str(Settings.MDP.RCoulomb_Cutoff),18);
    R_VDW_Cutoff = pad(num2str(Settings.MDP.RVDW_Cutoff),18);
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##RLIST##',R_List_Cutoff);
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##RCOULOMB##',R_Coulomb_Cutoff);
    Settings.MDP.Temp_Model = strrep(Settings.MDP.Temp_Model,'##RVDW##',R_VDW_Cutoff);

    % Save MDP file in current directory
    MDP_in_file = fullfile(Settings.WorkDir,[FileBase '.mdp']);
    fidMDP = fopen(MDP_in_file,'wt');
    fwrite(fidMDP,regexprep(Settings.MDP.Temp_Model,'\r',''));
    fclose(fidMDP);

    % Generate topology file
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##',Settings.Structure);
    Topology_File = fullfile(Settings.WorkDir,[FileBase '.top']);
    Atomlist = copy_atom_order(Settings.SuperCellFile);
    Topology_text_new = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
    fidTOP = fopen(Topology_File,'wt');
    fwrite(fidTOP,regexprep(Topology_text_new,'\r',''));
    fclose(fidTOP);

    Trajectory_File = fullfile(Settings.WorkDir,[FileBase '.tpr']);
    MDP_out_File = fullfile(Settings.WorkDir,[FileBase '_out.mdp']);
    GrompLog_File = fullfile(Settings.WorkDir,[FileBase '_Grompplog.log']);

    FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Settings.SuperCellFile) ...
        ' -f ' windows2unix(MDP_in_file) ' -p ' windows2unix(Topology_File) ...
        ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDP_out_File) ...
        ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GrompLog_File)];
    [state,~] = system(FMin_Grompp);
    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline FMin_Grompp]);
    else
        delete(GrompLog_File)
    end

    % Prepare minimization mdrun command
    Log_File = fullfile(Settings.WorkDir,[FileBase '.log']);

    Energy_file = fullfile(Settings.WorkDir,[FileBase '.edr']);

    TRR_File = fullfile(Settings.WorkDir,[FileBase '.trr']);
    
    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(Trajectory_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Settings.SuperCellFile) ...
        Settings.mdrun_opts];
    
    if Settings.Table_Req
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end
    mdrun_command = strrep(Settings.mdrun,'##MDRUN##',mdrun_command);
    
    % Final minimization
    [state,mdrun_output] = system(mdrun_command);
    if state ~= 0
        disp(mdrun_output);
        error(['Error running mdrun for final minimization. Problem command: ' newline mdrun_command]);
    end
end

% Add a liquid interface if requested        
if Settings.Liquid_Interface
    Create_Liquid_Solid_Interface(Settings);
end

% Generate final topology file for molecular dynamics
Atomlist = copy_atom_order(Settings.SuperCellFile);
Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
fidTOP = fopen(Settings.Topology_File,'wt');
fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
fclose(fidTOP);

if isfile(Settings.Traj_Conf_File)
    delete(Settings.Traj_Conf_File)
end

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
if Settings.Delete_Minimize
    system([Settings.wsl 'rm -r ' windows2unix(Settings.WorkDir) '/*']);
    rmdir(Settings.WorkDir);
end

end

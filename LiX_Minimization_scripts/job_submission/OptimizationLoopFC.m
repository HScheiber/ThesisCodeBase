function [N_Supercell_a,N_Supercell_b,N_Supercell_c] = OptimizationLoopFC(Settings)
    
    % Assign lattice parameter (a)
    % Find shortest lattice parameter
    a_vec = Settings.Geometry.Transform(1,:).*Settings.Geometry.a;
    b_vec = Settings.Geometry.Transform(2,:).*Settings.Geometry.b;
    c_vec = Settings.Geometry.Transform(3,:).*Settings.Geometry.c;
    LatticeLength = min([Settings.Geometry.Skew_a*norm(a_vec) ...
                        Settings.Geometry.Skew_b*norm(b_vec) ...
                        Settings.Geometry.Skew_c*norm(c_vec)]);
    
    % Create directory if it does not exist
    if ~exist(Settings.WorkDir,'dir')
        mkdir(Settings.WorkDir)
    end
    
    % Save topology file into directory
    Topology_File = fullfile(Settings.WorkDir,...
        [Settings.FileBase '.top']);
    
    % Determine cutoff length
    if Settings.MinMDP.Auto_Cutoff
        R_List_Cutoff = pad(num2str(LatticeLength/2-0.01),18);
        R_Coulomb_Cutoff = R_List_Cutoff;
        R_VDW_Cutoff = R_List_Cutoff;
    else
        R_List_Cutoff = pad(num2str(Settings.MinMDP.RList_Cutoff),18);
        R_Coulomb_Cutoff = pad(num2str(Settings.MinMDP.RCoulomb_Cutoff),18);
        R_VDW_Cutoff = pad(num2str(Settings.MinMDP.RVDW_Cutoff),18);

        % If chosen cutoffs are too large, increase supercell size to at least
        % large enough that the shortest supercell dimension is greater
        % than double the longest cutoff
        if LatticeLength/2 <= Settings.Longest_Cutoff
            N_Supercell_a = ceil(2*(Settings.Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_a*Settings.Geometry.a/10));
            N_Supercell_b = ceil(2*(Settings.Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_b*Settings.Geometry.b/10));
            N_Supercell_c = ceil(2*(Settings.Longest_Cutoff*Settings.Cutoff_Buffer)/(Settings.Geometry.Skew_c*Settings.Geometry.c/10));
        else
            N_Supercell_a = Settings.N_Supercell;
            N_Supercell_b = Settings.N_Supercell;
            N_Supercell_c = Settings.N_Supercell;
        end
    end
    
    % Update MDP file with cutoffs
    if strcmpi(Settings.MinMDP.VDWType,'Switch') || strcmpi(Settings.MinMDP.VDWType,'Shift')
        rvdw_line = 'rvdw                     = ##RVDW##; Van der Waals cutoff radius';
        vdw_switch_line = ['rvdw-switch              = ' ...
            pad(num2str(Settings.MinMDP.rvdw_switch),18) ...
            '; Where to start switching the LJ potential'];
        Settings.MDP_Template = regexprep(Settings.MDP_Template,rvdw_line,[rvdw_line newline vdw_switch_line]);
    end
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##RLIST##',R_List_Cutoff);
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##RCOULOMB##',R_Coulomb_Cutoff);
    Settings.MDP_Template = strrep(Settings.MDP_Template,'##RVDW##',R_VDW_Cutoff);
    
    % Save MDP file in current directory
    MDP_File = fullfile(Settings.WorkDir,[Settings.FileBase '.mdp']);
    fidMDP = fopen(MDP_File,'wt');
    fwrite(fidMDP,regexprep(Settings.MDP_Template,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fidMDP);
    
    % Add coordinates in xyz space (lattice parameter-dependent)
    Template_text_Lat = AddCartesianCoord(Settings.Coordinate_Template,...
        Settings.Geometry,1,false,Settings.CoordType);
    
    % Unit Cell Filename
    UnitCellFile = fullfile(Settings.WorkDir,[Settings.FileBase '_UnitCell.' Settings.CoordType]);
    
    % Supercell Filename
    SuperCellFile = fullfile(Settings.WorkDir,[Settings.FileBase '.' Settings.CoordType]);
    
    % Save unit cell .gro file into current directory
    fid = fopen(UnitCellFile,'wt');
    fwrite(fid,regexprep(Template_text_Lat,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fid);
    
    % Convert to Na x Nb x Nc supercell
    N_Cell = Settings.Geometry.N;
    N_a = num2str(N_Supercell_a);
    N_b = num2str(N_Supercell_b);
    N_c = num2str(N_Supercell_c);
    N_total = (N_Supercell_a*N_Supercell_b*N_Supercell_c)*N_Cell;
    
    Supercell_command = [Settings.gmx_loc Settings.genconf ' -f ' windows2unix(UnitCellFile) ...
         ' -o ' windows2unix(SuperCellFile) ' -nbox ' N_a ' ' N_b ' ' N_c];
    [err,output] = system(Supercell_command);
    
    if err ~= 0
        disp(output);
        error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
    end
    
    % Save number of atoms into .mat file
    NumberFile = fullfile(Settings.WorkDir,[Settings.FileBase '.mat']);
    save(NumberFile,'N_total','N_Cell')
    
    % Generate topology file
    Atomlist = copy_atom_order(SuperCellFile);
    Topology_text_new = strrep(Settings.Topology_Template,'##LATOMS##',Atomlist);
    fidTOP = fopen(Topology_File,'wt');
    fwrite(fidTOP,regexprep(Topology_text_new,'\r',''));
    fclose(fidTOP);
    
    % Create name for mdpout file
    MDPout_File = fullfile(Settings.WorkDir,...
        [Settings.FileBase '_out.mdp']);
    
    % Grompp log file
    GrompLog_File = fullfile(Settings.WorkDir,...
        [Settings.FileBase '_Grompplog.log']);
    
    % Prepare trajectory file
    Trajectory_File = fullfile(Settings.WorkDir,...
        [Settings.FileBase '.tpr']);
    
    if ispc
        passlog = ' ^&^> ';
    else
        passlog = ' &> ';
    end
    
    GROMPP_command = [Settings.gmx_loc Settings.grompp ' -c ' windows2unix(SuperCellFile) ...
        ' -f ' windows2unix(MDP_File) ' -p ' windows2unix(Topology_File) ...
        ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDPout_File) ...
        passlog windows2unix(GrompLog_File)];
    [state,~] = system(GROMPP_command);
    
    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
    end
    
    % Clean up MDP out file if set
    if Settings.Delete_MDPout
        delete(MDPout_File)
    end
    
    % Prepare mdrun command
    Log_File = fullfile(Settings.WorkDir,[Settings.FileBase '.log']);
    Energy_file = fullfile(Settings.WorkDir,[Settings.FileBase '.edr']);
    TRR_File = fullfile(Settings.WorkDir,[Settings.FileBase '.trr']);
    ConfOut_File = fullfile(Settings.WorkDir,[Settings.FileBase 'OutConf.' Settings.CoordType]);
    
    mdrun_command = [Settings.gmx Settings.mdrun ' -s ' windows2unix(Trajectory_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(ConfOut_File) ...
        ' -rerun ' windows2unix(SuperCellFile) Settings.mdrun_opts];
    
    if ~isempty(Settings.TableFile_MX) % If potential requires table
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end
    
    % Run it
    [state,~] = system(mdrun_command);
    
    if state ~= 0
        error(['Error running mdrun. Problem command: ' newline mdrun_command]);
    end
    
    % Remove files if requested
    if Settings.Delete_Supercell
        delete(SuperCellFile)
    end
    
    if Settings.Delete_MDlog
        delete(Log_File) %#ok<*UNRCH>
    end
    if Settings.Delete_ConfOut
        delete(ConfOut_File)
    end
    if Settings.Delete_TRR
        delete(TRR_File)
    end
    if Settings.Delete_TPR
        delete(Trajectory_File)
    end
    
    % Check energy options
    gmx_command = [Settings.wsl 'echo "0" ' Settings.pipe ...
        ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
        ' -f ' windows2unix(Energy_file)];
    
    [~,outpt] = system(gmx_command);
    
    en_opts = regexp(outpt,'End your selection with an empty line or a zero.\n-+(.+?)\n\n','tokens','once');
    en_opts = en_opts{1};
    En_set = char(regexp(en_opts,'([0-9]{1,2})  Potential','tokens','once'));
    En_set = regexprep([En_set ' 0'],' +',' ');
    
    % Convert the energy file into readable form
    Energy_output = fullfile(Settings.WorkDir,[Settings.FileBase 'energies.xvg']);
    eneconv_cmd = [Settings.wsl 'echo ' En_set ' ' Settings.pipe ...
        ' ' strrep(Settings.gmx_loc,Settings.wsl,'') Settings.g_energy ...
        ' -f ' windows2unix(Energy_file) ' -o ' windows2unix(Energy_output)];
    
    % Run it
    [state,eneconv_output] = system(eneconv_cmd);
    
    if state ~= 0
        disp(eneconv_output);
        error(['Error running eneconv. Problem command: ' newline eneconv_cmd]);
    end
    
    % Remove the energy input file if selected
    if Settings.Delete_EnergyIn
        delete(Energy_file)
    end
end
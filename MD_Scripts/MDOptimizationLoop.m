function N_Supercell_Out = MDOptimizationLoop(IP,Geometry,Directory,MDP_Template,Longest_Cutoff,Coordinate_text,EnergySetting)

    FileBase = [IP.Salt '_' Geometry.Label '_' IP.Model_Scaled '_' IP.OptTxt];
    
    % Assign lattice parameter (a)
    % Find shortest lattice parameter
    aLatpar = Geometry.a;
    bLatpar = Geometry.b;
    cLatpar = Geometry.c;
    LatticeLength = min([aLatpar bLatpar cLatpar])*IP.N_Supercell/10; % in nm

    % Save topology file into directory
    Topology_File = fullfile(Directory,[FileBase '.top']);

    % Determine cutoff length
    R_List_Cutoff = pad(num2str(IP.MDP_RList_Cutoff),18);
    R_Coulomb_Cutoff = pad(num2str(IP.MDP_RCoulomb_Cutoff),18);
    R_VDW_Cutoff = pad(num2str(IP.MDP_RVDW_Cutoff),18);

    % If chosen cutoffs are too large, increase supercell size to at least
    % large enough that the shortest supercell dimension is greater
    % than double the longest cutoff
    if LatticeLength/2 <= Longest_Cutoff
        IP.N_Supercell = ceil(2*(Longest_Cutoff + 0.01)/(min([Geometry.a ...
        Geometry.b Geometry.c]./10)));
    end

    % Update MDP file with cutoffs
    if strcmpi(IP.MDP_VDWType,'Switch') || strcmpi(IP.MDP_VDWType,'Shift')
        rvdw_line = 'rvdw                     = ##RVDW##; Van der Waals cutoff radius';
        vdw_switch_line = ['rvdw-switch              = ' ...
            pad(num2str(IP.MDP_rvdw_switch),18) ...
            '; Where to start switching the LJ potential'];
        MDP_Template = regexprep(MDP_Template,rvdw_line,[rvdw_line newline vdw_switch_line]);
    end
    MDP_Template = strrep(MDP_Template,'##RLIST##',R_List_Cutoff);
    MDP_Template = strrep(MDP_Template,'##RCOULOMB##',R_Coulomb_Cutoff);
    MDP_Template = strrep(MDP_Template,'##RVDW##',R_VDW_Cutoff);

    % Save MDP file in current directory
    MDP_File = fullfile(Directory,[FileBase '.mdp']);
    fidMDP = fopen(MDP_File,'wt');
    fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
    fclose(fidMDP);

    % Add coordinates in xyz space (lattice parameter-dependent)
    Template_text_Lat = AddCartesianCoord(Coordinate_text,Geometry,1,...
        false,IP.CoordType);

    % Unit Cell Filename
    UnitCellFile = fullfile(Directory,[FileBase '_UnitCell.' IP.CoordType]);

    % Supercell Filename
    SuperCellFile = fullfile(Directory,[FileBase '.' IP.CoordType]);

    % Save unit cell file into current directory
    fid = fopen(UnitCellFile,'wt');
    fwrite(fid,regexprep(Template_text_Lat,'\r',''));
    fclose(fid);

    % Convert to N x N x N supercell
    N_Cell = Geometry.N;
    N = num2str(IP.N_Supercell);
    N_total = (IP.N_Supercell^3)*N_Cell;

    Supercell_command = [IP.gmx ' genconf -f ' windows2unix(UnitCellFile) ...
         ' -o ' windows2unix(SuperCellFile) ' -nbox ' N ' ' N ' ' N];
    [state,output] = system(Supercell_command);

    if state ~= 0
        disp(output);
        error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
    end

    % Save number of atoms into .mat file
    NumberFile = fullfile(Directory,[FileBase '.mat']);
    save(NumberFile,'N_total','N_Cell')

    % Generate topology file
    Atomlist = copy_atom_order(SuperCellFile);
    Topology_text_new = strrep(IP.Topology_Text,'##LATOMS##',Atomlist);
    
    
    fidTOP = fopen(Topology_File,'wt');
    fwrite(fidTOP,regexprep(Topology_text_new,'\r',''));
    fclose(fidTOP);

    % Create name for mdpout file
    MDPout_File = fullfile(Directory,[FileBase '_out.mdp']);
    
    % Grompp log file
    GrompLog_File = fullfile(Directory,[FileBase '_Grompplog.log']);

    % Prepare trajectory file
    Trajectory_File = fullfile(Directory,[FileBase '.tpr']);
    
    GROMPP_command = [IP.gmx ' grompp -c ' windows2unix(SuperCellFile) ...
        ' -f ' windows2unix(MDP_File) ' -p ' windows2unix(Topology_File) ...
        ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDPout_File) ...
        IP.passlog windows2unix(GrompLog_File)];
    [state,~] = system(GROMPP_command);

    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
    else
        delete(GrompLog_File)
    end

    % Prepare mdrun command
    Log_File = fullfile(Directory,[FileBase '.log']);

    Energy_file = fullfile(Directory,[FileBase '.edr']);

    TRR_File = fullfile(Directory,[FileBase '.trr']);

    ConfOut_File = fullfile(Directory,[FileBase 'OutConf.' IP.CoordType]);
    
    mdrun_command = [IP.gmx ' mdrun -s ' windows2unix(Trajectory_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(ConfOut_File) ...
        ' -rerun ' windows2unix(SuperCellFile)];

    if ~isempty(IP.TableFile_MX) % TF potential requires table
        mdrun_command = [mdrun_command ' -table ' windows2unix(IP.TableFile_MX)];
    end
    
    % Run it
    [state,mdrun_output] = system(mdrun_command);

    if state ~= 0
        q=1;
        while state ~= 0 || q < 5
            [state,mdrun_output] = system(mdrun_command);
            q = q+1;
        end
        if q >= 5
            disp(mdrun_output);
            v = ['Failed to run mdrun at ' num2str(aLatpar,'%5.3f')];
            fid = fopen(fullfile(Directory,[FileBase '_mdrun_FailedRun.log']), 'at');
            fprintf(fid, '%s\n', v);
            fclose(fid);
        end
    end

    % Convert the energy file into readable form
    Energy_output = fullfile(Directory,[FileBase 'energies.xvg']);
    if ispc
        pipe = '^|';
    else
        pipe = '|';
    end

    eneconv_cmd = regexprep(IP.gmx,'(gmx[^\s]*)',['echo ' EnergySetting ' ' pipe ' $1 energy -f ' ...
        windows2unix(Energy_file) ' -o ' windows2unix(Energy_output)]);

    % Run it
    [state,eneconv_output] = system(eneconv_cmd);

    if state ~= 0
        disp(eneconv_output);
        error(['Error running eneconv. Problem command: ' newline eneconv_cmd]);
    end
    N_Supercell_Out = IP.N_Supercell;
end
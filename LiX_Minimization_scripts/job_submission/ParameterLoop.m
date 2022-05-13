function ParameterLoop(gmx,Cry,Salt,Structure,Model_Scaled,Label,N_Supercell,Maindir,...
    MDP,Longest_Cutoff,Coordinate_text,Settings,Topology_text,TableFile,...
    EnergySetting,pin,m)

    N_Expanded = [];
    eneconv_cmd = '';
    LatParTimer = tic;
    % Assign lattice parameter (a)
    % Find shortest lattice parameter
    aLatpars = Crystal_Info(Cry,Structure,'a'); %#ok<*PFBNS>
    bLatpars = Crystal_Info(Cry,Structure,'b');
    cLatpars = Crystal_Info(Cry,Structure,'c');
    aLatpar = aLatpars(m);
    bLatpar = bLatpars(m);
    cLatpar = cLatpars(m);
    LatticeLength = min([aLatpar bLatpar cLatpar])*N_Supercell/10;

    % Directory name
    LatparDir = pad(num2str(aLatpar,'%5.3f'),6,'left','0');

    % Update Directory
    Current_Directory = fullfile(Maindir,Salt,...
        Structure,Model_Scaled,['LPAR_' LatparDir]);

    % Create directory if it does not exist
    if ~exist(Current_Directory,'dir')
        mkdir(Current_Directory)
    end

    % Save topology file into directory
    Topology_File = fullfile(Current_Directory,...
        [Salt '_' Label '_' Model_Scaled '_' LatparDir '.top']);

    % Copy over MDP file
    MDP.text = MDP.Temp_Model;

    % Determine cutoff length
    if MDP.Auto_Cutoff
        R_List_Cutoff = pad(num2str(LatticeLength/2-0.01),18);
        R_Coulomb_Cutoff = R_List_Cutoff;
        R_VDW_Cutoff = R_List_Cutoff;

        Use_Expanded_Cell = false;
    else
        R_List_Cutoff = pad(num2str(MDP.RList_Cutoff),18);
        R_Coulomb_Cutoff = pad(num2str(MDP.RCoulomb_Cutoff),18);
        R_VDW_Cutoff = pad(num2str(MDP.RVDW_Cutoff),18);

        % If chosen cutoffs are too large, increase box size
        if LatticeLength/2 <= Longest_Cutoff
            N_Expanded = ceil(20*(Longest_Cutoff + 0.01)/(min([Cry.(Structure).a(m) ...
            Cry.(Structure).b(m) Cry.(Structure).c(m)])));
            Use_Expanded_Cell = true;
        else
            Use_Expanded_Cell = false;
        end
    end

    % Update MDP file with cutoffs
    MDP.text = strrep(MDP.text,'##RLIST##',R_List_Cutoff);
    MDP.text = strrep(MDP.text,'##RCOULOMB##',R_Coulomb_Cutoff);
    MDP.text = strrep(MDP.text,'##RVDW##',R_VDW_Cutoff);

    % Save MDP file in current directory
    MDP.File = fullfile(Current_Directory,...
        [Salt '_' Label '_' Model_Scaled '_' LatparDir '.mdp']);
    fidMDP = fopen(MDP.File,'wt');
    fwrite(fidMDP,MDP.text);
    fclose(fidMDP);

    % Add coordinates in xyz space (lattice parameter-dependent)
    Template_text_Lat = AddCartesianCoord(Coordinate_text,...
        Crystal_Info(Cry,Structure,''),m,Settings.FindParams,Settings.CoordType);

    % Unit Cell Filename
    UnitCellFile = fullfile(Current_Directory,[Salt '_' Label '_' ...
        Model_Scaled '_' LatparDir '_UnitCell.' Settings.CoordType]);

    % Supercell Filename
    SuperCellFile = fullfile(Current_Directory,[Salt '_' Label '_' ...
        Model_Scaled '_' LatparDir '.' Settings.CoordType]);

    % Save unit cell .gro file into current directory
    fid = fopen(UnitCellFile,'wt');
    fwrite(fid,regexprep(Template_text_Lat,{'\r', '\n\n+'}',{'', '\n'}));
    fclose(fid);

    % Convert to N x N x N supercell
    N_Cell = Crystal_Info(Cry,Structure,'N');
    if Use_Expanded_Cell
        N = num2str(N_Expanded);
        N_total = (N_Expanded^3)*N_Cell;
    else
        N = num2str(N_Supercell);
        N_total = (N_Supercell^3)*N_Cell;
    end
    Supercell_command = [gmx ' genconf -f ' windows2unix(UnitCellFile) ...
         ' -o ' windows2unix(SuperCellFile) ' -nbox ' N ' ' N ' ' N];
    [state,output] = system(Supercell_command);

    if state ~= 0
        disp(output);
        error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
    end

    % Save number of atoms into .mat file
    NumberFile = fullfile(Current_Directory,[Salt '_' Label '_' ...
        Model_Scaled '_' LatparDir '.mat']);
    parallelsave(NumberFile,{'N_total' 'N_Cell'},N_total,N_Cell)

    % Generate topology file
    Atomlist = copy_atom_order(SuperCellFile);
    Topology_text_new = strrep(Topology_text,'##LATOMS##',Atomlist);
    fidTOP = fopen(Topology_File,'wt');
    fwrite(fidTOP,Topology_text_new);
    fclose(fidTOP);

    % Create name for mdpout file
    MDPout_File = fullfile(Current_Directory,...
        [Salt '_' Label '_' Model_Scaled '_' LatparDir '_out.mdp']);
    
    % Grompp log file
    GrompLog_File = fullfile(Current_Directory,...
        [Salt '_' Label '_' Model_Scaled '_' LatparDir '_Grompplog.log']);

    % Prepare trajectory file
    Trajectory_File = fullfile(Current_Directory,...
        [Salt '_' Label '_' Model_Scaled '_' LatparDir '.tpr']);
    
    if ispc
        passlog = ' ^&^> ';
    else
        passlog = ' &> ';
    end
    
    GROMPP_command = [gmx ' grompp -c ' windows2unix(SuperCellFile) ...
        ' -f ' windows2unix(MDP.File) ' -p ' windows2unix(Topology_File) ...
        ' -o ' windows2unix(Trajectory_File) ' -po ' windows2unix(MDPout_File) ...
        passlog windows2unix(GrompLog_File)];
    [state,~] = system(GROMPP_command);

    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline GROMPP_command]);
    else
        delete(GrompLog_File)
    end

    % Clean up MDP out file if set
    if Settings.Delete_MDPout
        delete(MDPout_File)
        MDPout_File = '';
    end

    % Prepare mdrun command
    FileBase = [Salt '_' Label '_' Model_Scaled '_' LatparDir];

    Log_File = fullfile(Current_Directory,[FileBase '.log']);

    Energy_file = fullfile(Current_Directory,[FileBase '.edr']);

    TRR_File = fullfile(Current_Directory,[FileBase '.trr']);

    ConfOut_File = fullfile(Current_Directory,[FileBase 'OutConf.' Settings.CoordType]);
    
    mdrun_command = [gmx ' mdrun -s ' windows2unix(Trajectory_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(ConfOut_File) ...
        ' -rerun ' windows2unix(SuperCellFile) pin];

    if contains(Model_Scaled,'TF') % TF potential requires table
        mdrun_command = [mdrun_command ' -table ' windows2unix(TableFile)];
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
            fid = fopen(fullfile(Current_Directory,[FileBase '_mdrun_FailedRun.log']), 'at');
            fprintf(fid, '%s\n', v);
            fclose(fid);
        end
    end
    
    % Remove files if requested
    if Settings.Delete_Supercell
        delete(SuperCellFile)
        SuperCellFile = '';
    end

    if Settings.Delete_MDlog
        delete(Log_File) %#ok<*UNRCH>
        Log_File = '';
    end
    if Settings.Delete_ConfOut
        delete(ConfOut_File)
        ConfOut_File = '';
    end
    if Settings.Delete_TRR
        delete(TRR_File)
        TRR_File = '';
    end
    if Settings.Delete_TPR
        delete(Trajectory_File)
        Trajectory_File = '';
    end

    % Convert the energy file into readable form
    Energy_output = fullfile(Current_Directory,[FileBase 'energies.xvg']);
    if ispc
        eneconv_cmd = ['wsl source ~/.bashrc; echo ' EnergySetting ' ^| gmx_d energy -f ' ...
            windows2unix(Energy_file) ' -o ' windows2unix(Energy_output)];
    elseif isunix
        eneconv_cmd = ['echo ' EnergySetting ' | ' gmx ' energy -f ' ...
            Energy_file ' -o ' Energy_output];
    end

    % Run it
    [state,eneconv_output] = system(eneconv_cmd);

    if state ~= 0
        disp(eneconv_output);
        error(['Error running eneconv. Problem command: ' newline eneconv_cmd]);
    end

    % Remove the energy input file if selected
    if Settings.Delete_EnergyIn
        delete(Energy_file)
        Energy_file = '';
    end
    LatParElapsed = num2str(toc(LatParTimer));
    disp(['Lattice Parameter ' num2str(aLatpar) ' Finished. (' LatParElapsed ' s.)']);
end
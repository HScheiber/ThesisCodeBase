function Minimize_Liquid_Cluster(Settings)
    
    if Settings.Verbose
        disp('*** Combined Minimization of Liquid Selected ***')
    end
    Supercell_file_data = load_gro_file(Settings.SuperCellFile);
    
    % Calculate density of cluster
    a_vec = Settings.Geometry.a/10.*Settings.Geometry.Transform(1,:); % nm
    b_vec = Settings.Geometry.b/10.*Settings.Geometry.Transform(2,:); % nm
    c_vec = Settings.Geometry.c/10.*Settings.Geometry.Transform(3,:); % nm
    
    % Volume of unit cell in nm^3 and corresponding density of solid
    Volume_UC = Cell_Volume(a_vec,b_vec,c_vec); % [nm^3] Volume of unit cell 
    solid_density = Settings.Geometry.NF/Volume_UC; % molecules/nm^3
    Cluster_Vol = (Settings.Cluster_N/2)/solid_density; % nm^3
    
    % Volume required for liquid
    Liquid_Vol =  Settings.nmol_liquid/Settings.Ref_Density;
    
    Total_Vol = Cluster_Vol + Liquid_Vol;
    L = Total_Vol^(1/3); % box length in nm
    Box_center_old = (norm(Supercell_file_data.a_vec)/2).*[1 1 1];
    Box_center_new = (L/2).*[1 1 1];
    
    Supercell_file_data.a_vec  = [L 0 0];
    Supercell_file_data.b_vec  = [0 L 0];
    Supercell_file_data.c_vec  = [0 0 L];
    Supercell_file_data.boxcoords = {L L L 0 0 0 0 0 0};
    Supercell_file_data.Volume = Total_Vol;
    Supercell_file_data.Salt = Settings.Salt;
    Supercell_file_data.N = Settings.Cluster_N;
    
    % Update box geometry, recenter
    Supercell_file_data.xyz = Supercell_file_data.xyz + (Box_center_new - Box_center_old);
    
    Solid_Geometry.Liquid_Vol = Liquid_Vol; % nm^3
    Solid_Geometry.Solid_Vol = Cluster_Vol; % nm^3
    Solid_Geometry.Solid_R = ( 3*Cluster_Vol/(4*pi) )^(1/3); % cluster radius in nm
    
    % Save solid-only supercell file info
    SS_file_info = dir(Settings.SuperCellFile);
    Solid_Geometry.a_vec = Supercell_file_data.a_vec;
    Solid_Geometry.b_vec = Supercell_file_data.b_vec;
    Solid_Geometry.c_vec = Supercell_file_data.c_vec;
    Solid_Geometry.c_vec_tot = c_vec;
    Solid_Geometry.L = L; % Length of liquid box
    Solid_Geometry.N_atoms = Supercell_file_data.N_atoms;
    Solid_Geometry.Volume = Supercell_file_data.Volume;
    Solid_Geometry_file = fullfile(SS_file_info.folder,[Settings.JobName '_SolInfo.mat']);
    save(Solid_Geometry_file,'Solid_Geometry');

    % Rename solids resname to SOL and save
    Supercell_file_data.res_name(:) = {'SOL  '};
    Comb_Equil_Sol_File = fullfile(Settings.WorkDir,['Comb_Equil_Sol.' Settings.CoordType]);
    SaveGroFile(Comb_Equil_Sol_File,Supercell_file_data,true);
    
    % Generate a vdwradii.dat file such that no liquid particles end up inside the solid
    r0 = 0.8*min(pdist(Supercell_file_data.xyz(1:min(Supercell_file_data.N,5e4),:)))/2;
    Halide_IDX = cellfun(@(x) strcmp(x,pad(Settings.Halide,5,'left')),Supercell_file_data.atom_name);
    r0_sol = max(min(pdist(Supercell_file_data.xyz(~Halide_IDX,:))),...
        min(pdist(Supercell_file_data.xyz(Halide_IDX,:))))/2; % minimum atomic radius in nm

    vdwradii_txt = ['SOL  '             pad(Settings.Metal,2)   '    ' num2str(r0_sol) newline ...
                    'SOL  '             pad(Settings.Halide,2)  '    ' num2str(r0_sol) newline ...
                    pad(Settings.Metal,2)  '   ' pad(Settings.Metal,2)   '    ' num2str(r0)     newline ...
                    pad(Settings.Halide,2) '   ' pad(Settings.Halide,2)  '    ' num2str(r0)];
    vdwDatFile = fullfile(Settings.WorkDir,'vdwradii.dat');
    fidMDP = fopen(vdwDatFile,'wt');
    fwrite(fidMDP,regexprep(vdwradii_txt,'\r',''));
    fclose(fidMDP);

    %% Add liquid atoms to the file
    OldDir = pwd;
    cd(Settings.WorkDir);
    
    Ref_M = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Metal '_Box.gro']);
    Ref_X = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Halide '_Box.gro']);
    
    % Add metal
    if Settings.Verbose
        disp(['Randomly adding ' num2str(Settings.nmol_liquid) ' ' Settings.Metal ' ions to liquid box...'])
    end
    mtimer = tic;
    Comb_Equil_Sol_Metal_File = fullfile(Settings.WorkDir,['Comb_Equil_Sol_Metal.' Settings.CoordType]);
    cmd = [Settings.gmx_loc ' insert-molecules -f ' windows2unix(Comb_Equil_Sol_File) ' -ci ' windows2unix(Ref_M) ...
        ' -o ' windows2unix(Comb_Equil_Sol_Metal_File) ' -nmol ' num2str(Settings.nmol_liquid) ' -radius ' num2str(r0) ...
        ' -scale 1 -try 200'];
    [errcode,output] = system(cmd);

    if errcode ~= 0
        if Settings.Verbose
            disp(output);
        end
        error(['Error adding ' Settings.Metal ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    if Settings.Verbose
        disp([Settings.Metal ' atoms added. Epalsed Time: ' datestr(seconds(toc(mtimer)),'HH:MM:SS')])
    end

    % Add Halide
    Comb_Equil_File = fullfile(Settings.WorkDir,['Comb_Equil.' Settings.CoordType]);
    if Settings.Verbose
        disp(['Randomly adding ' num2str(Settings.nmol_liquid) ' ' Settings.Halide ' ions to liquid box...'])
    end
    htimer = tic;
    cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(Comb_Equil_Sol_Metal_File) ...
        ' -o ' windows2unix(Comb_Equil_File) ' -nmol ' num2str(Settings.nmol_liquid) ...
        ' -scale 1 -radius ' num2str(r0) ' -try 400'];
    [errcode,output] = system(cmd);

    if errcode ~= 0
        if Settings.Verbose
            disp(output);
        end
        error(['Error adding ' Settings.Halide ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    if Settings.Verbose
        disp([Settings.Halide ' atoms added. Epalsed Time: ' datestr(seconds(toc(htimer)),'HH:MM:SS')])
    end

    % Load current randomly-generated liquid file data
    Comb_Equil_data = load_gro_file(Comb_Equil_File);

    % Check to make sure all atoms were added
    if ~(Comb_Equil_data.N_atoms == Settings.nmol_liquid*2 + Supercell_file_data.N)
        error('Not all requested liquid atoms were added!')
    end
    cd(OldDir); % Return to previous directory

    %% Minimize the randomly-generated liquid: build MDP file
    MDP = Settings.MDP;
    MDP.Minimization_txt = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
    'MDP.template'));

    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##NSTEPS##',pad(num2str(Settings.MinMDP.nsteps_min),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##INTEGR##',pad(Settings.MinMDP.min_integrator,18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'dt                       = ##TIMEST##; Time step (ps)',...
        ['emtol                    = ' num2str(Settings.MinMDP.emtol)]);
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##MET##',Settings.Metal);
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##HAL##',Settings.Halide);
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##COULOMB##',pad(num2str(MDP.CoulombType),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##FOURIER##',pad(num2str(MDP.Fourier_Spacing),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##PMEORDER##',pad(num2str(MDP.PME_Order),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##EWALDTOL##',pad(num2str(MDP.Ewald_rtol),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##RLIST##',pad(num2str(MDP.RList_Cutoff),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##RCOULOMB##',pad(num2str(MDP.RCoulomb_Cutoff),18));
    MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##RVDW##',pad(num2str(MDP.RVDW_Cutoff),18)); 

    if Settings.Table_Req || strcmp(Settings.Theory,'BH')
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##VDWTYPE##',pad('user',18));
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##VDWMOD##',pad(MDP.vdw_modifier,18));
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##CUTOFF##',pad('group',18));
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'ewald-rtol-lj.+?\n','');
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'lj-pme-comb-rule.+?\n','');
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'verlet-buffer-tolerance.+?\n','');

        % For minimization, add in a close-range repulsive wall to the
        % potential with the following function
        Settings.TableFile_MX = MakeTables(Settings);
    else
        % Modify the MDP file
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##VDWTYPE##',pad(MDP.VDWType,18));
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##VDWMOD##',pad(MDP.vdw_modifier,18));
        MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##CUTOFF##',pad(MDP.CutOffScheme,18));
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'energygrp-table.+?\n','');
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'ewald-rtol-lj.+?\n','');
        MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'lj-pme-comb-rule.+?\n','');

        % Add in Verlet Settings
        if strcmp(MDP.CutOffScheme,'Verlet')
            MDP.Minimization_txt = strrep(MDP.Minimization_txt,'##VerletBT##',pad(num2str(MDP.VerletBT),18));
        else
            MDP.Minimization_txt = regexprep(MDP.Minimization_txt,'verlet-buffer-tolerance.+?\n','');
        end
    end

    % Create an index file to keep track of solid and liquid atom numbers
    sol_atnums = double(Supercell_file_data.atom_number);
    liq_atnums = double(Comb_Equil_data.atom_number(Supercell_file_data.N_atoms+1:end));
    system_atnums = [sol_atnums; liq_atnums];

    % Get indeces for the metal and halide
    M_idx = find(contains(Comb_Equil_data.atom_name,Settings.Metal));
    M_atnums = double(Comb_Equil_data.atom_number(M_idx));
    X_idx = find(contains(Comb_Equil_data.atom_name,Settings.Halide));
    X_atnums = double(Comb_Equil_data.atom_number(X_idx));

    % Create index file
    addnan = 15 - mod(length(system_atnums),15);
    system_atnums(end+1:end+addnan) = nan;
    system_atnums = reshape(system_atnums,15,[])';
    system_atnum_txt = strtrim(char(regexprep(strjoin(string(num2str(system_atnums)),newline),' +NaN','')));

    addnan = 15 - mod(length(sol_atnums),15);
    sol_atnums(end+1:end+addnan) = nan;
    sol_atnums = reshape(sol_atnums,15,[])';
    sol_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(sol_atnums)),newline),' +NaN','')));

    addnan = 15 - mod(length(M_atnums),15);
    M_atnums(end+1:end+addnan) = nan;
    M_atnums = reshape(M_atnums,15,[])';
    M_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(M_atnums)),newline),' +NaN','')));

    addnan = 15 - mod(length(X_atnums),15);
    X_atnums(end+1:end+addnan) = nan;
    X_atnums = reshape(X_atnums,15,[])';
    X_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(X_atnums)),newline),' +NaN','')));

    ndx_text = ['[ System ]'     newline system_atnum_txt newline ...
                '[ Solid ]'      newline sol_atnums_txt newline ...
                '[ ' Settings.Metal ' ]'  newline M_atnums_txt newline ...
                '[ ' Settings.Halide ' ]' newline X_atnums_txt newline];

    % Save index file
    NDX_Filename = fullfile(Settings.WorkDir,'Comb_Equil.ndx');
    fidNDX = fopen(NDX_Filename,'wt');
    fwrite(fidNDX,regexprep(ndx_text,'\r',''));
    fclose(fidNDX);

    % Freeze the solid by adding freeze groups to the mdp file
    Freezegrps_text = [newline newline '; Freezing Groups' newline ...
    'freezegrps               = Solid' newline ...
    'freezedim                = Y Y Y'];
    MDP.Minimization_txt = [MDP.Minimization_txt Freezegrps_text];

    % Add in dispersion corrections
    if Settings.MDP.Disp_Correction && ~(Settings.Table_Req || strcmp(Settings.Theory,'BH'))
        MDP.Minimization_txt = [MDP.Minimization_txt newline newline...
            '; Long-range dispersion correction' newline ...
            'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
    elseif Settings.MDP.Disp_Correction && Settings.MDP.Disp_Correction_Tables
        if Settings.Verbose
            disp('Warning: enabling long-range dispersion correction for tabulated potential!')
        end
        MDP.Minimization_txt = [MDP.Minimization_txt newline newline ...
            '; Long-range dispersion correction' newline ...
            'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
    end

    % Save MDP file
    MDP_Filename = fullfile(Settings.WorkDir,'Comb_Equil.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP.Minimization_txt,'\r',''));
    fclose(fidMDP);

    % Complete a topology file for the liquid box to be minimized
    Atomlist = copy_atom_order(Comb_Equil_File);
    Top_Filename = fullfile(Settings.WorkDir,'Comb_Equil.top');
    
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
        Topology_Text = strrep(Topology_Text,'##LATOMS##',Atomlist);
    else
        Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
    end
    
    % Save topology file
    fidTOP = fopen(Top_Filename,'wt');
    fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
    fclose(fidTOP);

    TPR_File = fullfile(Settings.WorkDir,'Comb_Equil.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Comb_Equil_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Comb_Equil_Grompplog.log');

    FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Comb_Equil_File) ...
        ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Top_Filename) ...
        ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
        ' -n ' windows2unix(NDX_Filename) ' -maxwarn ' num2str(Settings.MaxWarn) ...
         Settings.passlog windows2unix(GrompLog_File)];
    [state,~] = system(FMin_Grompp);

    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline FMin_Grompp]);
    else
        delete(GrompLog_File)
    end

    % Prepare minimization mdrun command
    Log_File = fullfile(Settings.WorkDir,'Comb_Equil.log');
    Energy_file = fullfile(Settings.WorkDir,'Comb_Equil.edr');
    TRR_File = fullfile(Settings.WorkDir,'Comb_Equil.trr');

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Settings.SuperCellFile) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Comb_Equil')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req || strcmp(Settings.Theory,'BH')
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Remove previous supercell
    delete(Settings.SuperCellFile)

    % Final minimization, output geometry is the new supercell file
    if Settings.Verbose
        disp('Begining Combined System minimization...')
    end
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        if Settings.Verbose
            disp(['System Successfully Minimized! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
        end
    else
        if Settings.Verbose
            disp(mdrun_output);
        end
        error(['Error running mdrun for system minimization. Problem command: ' newline mdrun_command]);
    end
    
    if Settings.Verbose
        disp('*** Minimization of Liquid Complete ***')
    end
    
    % Remove minimization temporary folder
    if Settings.Delete_Backups
        system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
            ' -iname "#*#" -delete']);
    end
    
end
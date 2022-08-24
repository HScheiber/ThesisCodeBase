function Output = Minimize_Equilibrate_Liquid_Interface(Settings)

    % What this function does:
    % Initialize density based on experiment, and the required number of atoms in the liquid phase
    % Create a box that matches the x-y dimensions of the solid box
    % Set the z-length of the box so that it will fit the required number of atoms in the liquid phase at the experimental density
    % Randomly fill box with the required number of atoms for the full simulation
    % Use Make_Tables function if necessary, run a [Settings.MinMDP.nsteps_min] step minimization
    % Set up for fast NPT equilibration: Use semiisotropic barostat. Set box compressibility in X-Y to zero and Z to appropriate value. Use Berendsen Barostat and Berendsen Thermostat, both with time constants equal to the time step. 
    % Equilibrate for [Settings.Equilibrate_Liquid] amount of time with fast equilibration settings (Berendsen baro and Velocity-Rescale thermo)
    % Attach equilibrated liquid box to solid such that X-Y dimensions align, Displace all liquid atoms by Z length of solid box
    % Fix all atoms except atoms within -+0.5 Angstroms of the solid-liquid interface interface. Re-minimize interfacial atoms using Make-Tables    
    Output.StructureChange = false;
    Output.SolidMelted = false;
    Output.LiquidFroze = false;
    Output.LiquidAmorphous = false;
    Output.Aborted = false;
    
    if Settings.Verbose
        disp('*** Separate Equilibration of Liquid Selected ***')
    end
    Inp_Settings = Settings;
    MinDir = Settings.WorkDir;
    Settings.WorkDir = fullfile(Settings.WorkDir,'Equil_Liq');
    if ~isfolder(Settings.WorkDir)
        mkdir(Settings.WorkDir)
    end
    
    % Load solid box
    Supercell_file_data = load_gro_file(Settings.SuperCellFile);
    
    % Calculate the Z box dimension needed for the given density
    XY_area = norm(cross(Supercell_file_data.a_vec,Supercell_file_data.b_vec)); % nm^2
    Liq_Vol = Settings.nmol_liquid/(Settings.Ref_Density*Settings.ScaleInitialLiqDensity); % Volume of Liq box in nm^3
    Liq_Z = Liq_Vol/XY_area; % Liquid Z box length in nm
    R0 = num2str(min(0.5*((3/(4*pi))*(Liq_Vol/Settings.nmol_liquid))^(1/3),0.57),'%0.3f');
    
    if Liq_Z < norm(Supercell_file_data.c_vec)
        Liq_Z = norm(Supercell_file_data.c_vec);
    end
    
    % Box vectors of the Liquid box
    a_vec = Supercell_file_data.a_vec;
    b_vec = Supercell_file_data.b_vec;
    c_vec = (Supercell_file_data.c_vec/norm(Supercell_file_data.c_vec)).*Liq_Z;
    
    % Initialize empty box with correct dimensions
    Liq_Box.a_vec  = a_vec;
    Liq_Box.b_vec  = b_vec;
    Liq_Box.c_vec  = c_vec;
    Liq_Box.boxcoords = {a_vec(1) b_vec(2) c_vec(3) a_vec(2) a_vec(3) b_vec(1) b_vec(3) c_vec(1) c_vec(2)};
    Liq_Box.Volume = abs(dot(cross(a_vec,b_vec),c_vec));
    Liq_Box.Salt = Settings.Salt;
    Liq_Box.N = 0;
    Liq_Box.N_atoms = 0;
    Liq_Box.comment = [Settings.Salt ' liquid cell'];
    Liq_Box.res_number = double.empty(0,1);
    Liq_Box.res_name = cell.empty(0,1);
    Liq_Box.atom_name = cell.empty(0,1);
    Liq_Box.atom_number = int32.empty(0,1);
    Liq_Box.xyz = double.empty(0,3);
    Liq_Box.vel = double.empty(0,3);
    
    % Save empty box
    Liq_Box_File = fullfile(Settings.WorkDir,['Prep_Liq_Empty.' Settings.CoordType]);
    SaveGroFile(Liq_Box_File,Liq_Box,true);

    Ref_M = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Metal '_Box.gro']);
    Ref_X = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Halide '_Box.gro']);

    % Add metal
    if Settings.Verbose
        disp(['Randomly adding ' num2str(Settings.nmol_liquid) ' ' Settings.Metal ' ions to liquid box...'])
    end
    mtimer = tic;
    Prep_Liq_Metal_Only = fullfile(Settings.WorkDir,['Prep_Liq_Metal_Only.' Settings.CoordType]);
    cmd = [Settings.gmx_loc ' insert-molecules -f ' windows2unix(Liq_Box_File) ' -ci ' windows2unix(Ref_M) ...
        ' -o ' windows2unix(Prep_Liq_Metal_Only) ' -nmol ' num2str(Settings.nmol_liquid) ' -try 200 -scale ' R0 ' -radius ' R0];
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
    Prep_Liq_Random_Liq = fullfile(Settings.WorkDir,['Prep_Liq_Random_Liq.' Settings.CoordType]);
    if Settings.Verbose
        disp(['Randomly adding ' num2str(Settings.nmol_liquid) ' ' Settings.Halide ' ions to liquid box...'])
    end
    htimer = tic;
    cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(Prep_Liq_Metal_Only) ...
        ' -o ' windows2unix(Prep_Liq_Random_Liq) ' -nmol ' num2str(Settings.nmol_liquid) ' -try 400 -scale ' R0 ' -radius ' R0];

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
    Random_Liquid_data = load_gro_file(Prep_Liq_Random_Liq);

    % Check to make sure all atoms were added
    if Random_Liquid_data.N_atoms ~= Settings.nmol_liquid*2
        error('Not all requested liquid atoms were added!')
    end
    
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
        if isfield(Settings,'TableFile_MX')
            TableFile_MX_old = Settings.TableFile_MX;
        else
            TableFile_MX_old = '';
        end
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
    MDP_Filename = fullfile(Settings.WorkDir,'Prep_Liq.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP.Minimization_txt,'\r',''));
    fclose(fidMDP);

    % Complete a topology file for the liquid box to be minimized
    Atomlist = copy_atom_order(Prep_Liq_Random_Liq);
    Top_Filename = fullfile(Settings.WorkDir,'Prep_Liq.top');
    
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
        Topology_Text_Min = Topology_Text;
        Topology_Text = strrep(Topology_Text,'##LATOMS##',Atomlist);
    else
        Topology_Text_Min = Settings.Topology_Text;
        Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
    end
    
    % Save topology file
    fidTOP = fopen(Top_Filename,'wt');
    fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
    fclose(fidTOP);
    
    TPR_File = fullfile(Settings.WorkDir,'Prep_Liq.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Prep_Liq_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Prep_Liq_Grompplog.log');

    FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Prep_Liq_Random_Liq) ...
        ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Top_Filename) ...
        ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
        ' -maxwarn ' num2str(Settings.MaxWarn) ...
         Settings.passlog windows2unix(GrompLog_File)];
    [state,~] = system(FMin_Grompp);
    % Catch error in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline FMin_Grompp]);
    else
        delete(GrompLog_File)
    end

    % Prepare minimization mdrun command
    Log_File = fullfile(Settings.WorkDir,'Prep_Liq.log');
    Energy_file = fullfile(Settings.WorkDir,'Prep_Liq.edr');
    TRR_File = fullfile(Settings.WorkDir,'Prep_Liq.trr');
    Minimized_Geom_File = fullfile(Settings.WorkDir,['Prep_Liq_Min.' Settings.CoordType]);

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Minimized_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Prep_Liq')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req || strcmp(Settings.Theory,'BH')
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Liquid Minimization
    if Settings.Verbose
        disp('Beginning Liquid Minimization...')
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
        error(['Error running mdrun for liquid system minimization. Problem command: ' newline mdrun_command]);
    end
    
    %% System is now minimized, run a fast equilibration, allowing box to equilibrate in Z dimension only
    if Settings.Table_Req
    	Settings.TableFile_MX = TableFile_MX_old;
    end
    % Buckingham potentials require recreating the topology text
    if strcmp(Settings.Theory,'BH')
        % Save topology file
        Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
        fidTOP = fopen(Top_Filename,'wt');
        fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
        fclose(fidTOP);
    end
    
    % Set the number of steps
    timesteps = Settings.Equilibrate_Liquid/Settings.MDP.dt;
    %Compressibility = Get_Alkali_Halide_Compressibility(Settings.Salt,'Isotropy','isotropic','Molten',true);
    Compressibility = ['0 ' num2str(Settings.QECompressibility)]; % bar^(-1)
    tau_p = Settings.MDP.dt; % ps
    tau_t = Settings.MDP.dt; % ps
    
    nstpcouple = Get_nstcouple(tau_p,Settings.MDP.dt);
    nsttcouple = Get_nstcouple(tau_t,Settings.MDP.dt);
    
    Target_P = [num2str(Settings.Target_P(1)) ' ' num2str(Settings.Target_P(1))];
    
    % Ensure fast equilibration with Berendsen barostat + small time constant
    xyz_out = num2str(0.1 / Settings.MDP.dt); % Output coords every 0.1 ps
    MDP_Template = regexprep(Settings.MDP_Template,'(nsteps += *)(.+?)( *);',['$1' num2str(timesteps) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = regexprep(MDP_Template,'(nstcalcenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = regexprep(MDP_Template,'(nstxout += *)(.+?)( *);',['$1' xyz_out '$3;']);
    MDP_Template = regexprep(MDP_Template,'(pcoupl += *)(.+?)( *);','$1Berendsen$3;');
    MDP_Template = regexprep(MDP_Template,'(pcoupltype += *)(.+?)( *);','$1semiisotropic$3;');
    MDP_Template = regexprep(MDP_Template,'(tau-p += *)(.+?)( *);',['$1 ' num2str(tau_p) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstpcouple += *)(.+?)( *);',['$1 ' num2str(nstpcouple) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(compressibility += *)(.+?)( *);',['$1 ' Compressibility '$3;']);
    MDP_Template = regexprep(MDP_Template,'(ref-p += *)(.+?)( *);',['$1 ' Target_P '$3;']);
    MDP_Template = regexprep(MDP_Template,'(dt += *)(.+?)( *);',['$1' num2str(Settings.MDP.dt) '$3;']);
    
    % Pair it with velocity rescale thermostat + small time constant
    MDP_Template = regexprep(MDP_Template,'(tcoupl += *)(.+?)( +);','$1v-rescale$3;');
    MDP_Template = regexprep(MDP_Template,'(tau-t += *)(.+?)( +);',['$1 ' num2str(tau_t) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nsttcouple += *)(.+?)( +);',['$1 ' num2str(nsttcouple) '$3;']);
    
    % Save MDP file, topology file can be reused
    MDP_Filename = fullfile(Settings.WorkDir,'Equil_Liq.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
    fclose(fidMDP);
    
    TPR_File = fullfile(Settings.WorkDir,'Equil_Liq.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Equil_Liq_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Equil_Liq_Grompplog.log');
    
    FEquil_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Minimized_Geom_File) ...
        ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Top_Filename) ...
        ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
        ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GrompLog_File)];
    [state,~] = system(FEquil_Grompp);
    
    % Catch errors in grompp
    if state ~= 0
        error(['Error running GROMPP. Problem command: ' newline FEquil_Grompp]);
    else
        delete(GrompLog_File)
    end
    
    % Prepare Equilibration mdrun command
    Log_File = fullfile(Settings.WorkDir,'Equil_Liq.log');
    Energy_file = fullfile(Settings.WorkDir,'Equil_Liq.edr');
    TRR_File = fullfile(Settings.WorkDir,'Equil_Liq.trr');
    Equilibrated_Geom_File = fullfile(MinDir,['Equil_Liq.' Settings.CoordType]);
    
    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...   
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Equilibrated_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Equil_Liq')) ...
        Settings.mdrun_opts];
    
    if Settings.Table_Req
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end
    
    % Run Liquid Equilibration
    if Settings.Verbose
        disp(['Beginning Liquid Equilibration for ' num2str(Settings.Equilibrate_Liquid) ' ps...'] )
    end
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        if Settings.Verbose
            disp(['Liquid Successfully Equilibrated! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
        end
    else
        try % Clean up
            [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
        catch me
            disp(me.message)
        end
        
        Settings = Inp_Settings;
%         if ~isfield(Settings,'QECompressibility_init')
%             Settings.QECompressibility_init = Settings.QECompressibility;
%         end
%         if Settings.QECompressibility > 1e-8 % Retry until compressibility is very tight
%             if Settings.Verbose
%                 disp('Equilibration failed. Retrying with stiffer compressibility.')
%             end
%             Settings.QECompressibility = Settings.QECompressibility/2;
%             Output = Minimize_Equilibrate_Liquid_Interface(Settings);
%             return
        if Settings.MDP.dt > 1e-4
            if Settings.Verbose
                disp('Equilibration failed. Reducing time step.')
            end
            %Settings.QECompressibility = Settings.QECompressibility_init;
            Settings.MDP.dt = Settings.MDP.dt/2;
            Settings.Output_Coords = Settings.Output_Coords*2;
            Output = Minimize_Equilibrate_Liquid_Interface(Settings);
            return
        else
            if Settings.Verbose
                disp('Equilibration failed. Reducing time step did not resolve.')
                disp(mdrun_output);
                disp(['Error running mdrun for liquid equilibration. Problem command: ' newline mdrun_command]);
            end
            Output.Aborted = true;
            TDir = fullfile(strrep(MinDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
            [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
            copyfile(Settings.WorkDir,TDir)
            try
                if Settings.Delete_Equil
                    rmdir(Settings.WorkDir,'s')
                end
            catch
                disp(['Unable to remove directory: ' Settings.WorkDir])
            end
            return
        end
    end
    
%     system(['wsl source ~/.bashrc; echo "5 8 15 0" ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o ' windows2unix(strrep(Energy_file,'.edr','.xvg'))])
%     En_xvg_file = fullfile(Settings.WorkDir,'Equil_Liq.xvg');
%     Data = import_xvg(En_xvg_file);
%     plot(Data(:,1),Data(:,2)./Settings.nmol_liquid) % Potential
%     plot(Data(:,1),Data(:,3)./Settings.nmol_liquid) % Conversved Energy
%     plot(Data(:,1),(10^3).*Data(:,4)./Settings.nmol_liquid) % Volume
    
    % Check to ensure system remained liquid
    copyfile(Minimized_Geom_File,fullfile(Settings.WorkDir,['Equil_Liq.' Settings.CoordType]));
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(Settings.WorkDir, Settings.Salt, ...
        pyargs('SystemName','Equil_Liq',...
        'RefStructure',Settings.Structure,...
        'CheckFullTrajectory',true,...
        'FileType',Settings.CoordType,...
        'ML_TimeLength',0,...
        'ML_TimeStep',0,...
        'SaveTrajectory',true,...
        'SavePredictionsImage',true));
    Sol_Fraction = PyOut{4};
    Liq_Fraction = PyOut{5};
    
    if Liq_Fraction < (1 - Settings.MeltFreezeThreshold)
        if Settings.Verbose
            disp('Detected Liquid Phase change.')
        end
        if (1-Liq_Fraction-Sol_Fraction) >= Settings.MeltFreezeThreshold
            Output.StructureChange = true;
        else
            Output.LiquidFroze = true;
        end
        Output.Aborted = true;
        TDir = fullfile(strrep(MinDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
        [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
        copyfile(Settings.WorkDir,TDir)
        return
    end
    
    %% Optional: Check if liquid is properly mobile (i.e. not amorphous solid)
    if Settings.CheckAmorphousLiquid
        
        % Set the number of steps
        MD_nsteps = Settings.Liquid_Test_Time/Settings.MDP.dt;

        % Ensure fast equilibration with Berendsen barostat + small time constant
        xyz_out = num2str(0.1 / Settings.MDP.dt); % Output coords every 0.1 ps
        MDP_Template = regexprep(Settings.MDP_Template,'(nsteps += *)(.+?)( *);',['$1' num2str(MD_nsteps) '$3;']);
        MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$11$3;');
        MDP_Template = regexprep(MDP_Template,'(nstcalcenergy += *)(.+?)( *);','$11$3;');
        MDP_Template = regexprep(MDP_Template,'(nstxout += *)(.+?)( *);',['$1' xyz_out '$3;']);
        MDP_Template = regexprep(MDP_Template,'(dt += *)(.+?)( *);',['$1' num2str(Settings.MDP.dt) '$3;']);
        MDP_Template = regexprep(MDP_Template,'(gen-vel += *)(.+?)( *);','$1no$3;');
        MDP_Template = regexprep(MDP_Template,'gen-temp.+?\n','');
        
        % Save MDP file
        MDP_Filename = fullfile(Settings.WorkDir,'MSD_Liq.mdp');
        fidMDP = fopen(MDP_Filename,'wt');
        fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
        fclose(fidMDP);
        
        % Generate TPR with grompp
        TPR_File = fullfile(Settings.WorkDir,'MSD_Liq.tpr');
        MDPout_File = fullfile(Settings.WorkDir,'MSD_Liq_out.mdp');
        GrompLog_File = fullfile(Settings.WorkDir,'MSD_Liq_Grompplog.log');

        FEquil_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Equilibrated_Geom_File) ...
            ' -f ' windows2unix(MDP_Filename) ' -p ' windows2unix(Top_Filename) ...
            ' -o ' windows2unix(TPR_File) ' -po ' windows2unix(MDPout_File) ...
            ' -maxwarn ' num2str(Settings.MaxWarn) Settings.passlog windows2unix(GrompLog_File)];
        [state,~] = system(FEquil_Grompp);

        % Catch errors in grompp
        if state ~= 0
            error(['Error running GROMPP. Problem command: ' newline FEquil_Grompp]);
        else
            delete(GrompLog_File)
        end
        
        % Prepare Equilibration mdrun command
        Log_File = fullfile(Settings.WorkDir,'MSD_Liq.log');
        Energy_file = fullfile(Settings.WorkDir,'MSD_Liq.edr');
        TRR_File = fullfile(Settings.WorkDir,'MSD_Liq.trr');
        Dynamics_Geom_File = fullfile(Settings.WorkDir,['MSD_Liq_out.' Settings.CoordType]);

        mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
            ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
            ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Dynamics_Geom_File) ...
            ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'MSD_Liq')) ...
            Settings.mdrun_opts];

        if Settings.Table_Req
            mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
        end

        % Run Liquid Equilibration
        if Settings.Verbose
            disp(['(2/2) Running liquid with realistic dynamics for ' num2str(Settings.Liquid_Test_Time) ' ps...'] )
        end
        mintimer = tic;
        [state,~] = system(mdrun_command);
        if state == 0
            if Settings.Verbose
                disp(['(2/2) Liquid dynamics simulation complete. Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
            end
            
            % Check to ensure system remained liquid
            copyfile(Equilibrated_Geom_File,fullfile(Settings.WorkDir,['MSD_Liq.' Settings.CoordType]));
            PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(Settings.WorkDir, Settings.Salt, ...
                pyargs('SystemName','MSD_Liq',...
                'RefStructure',Settings.Structure,...
                'CheckFullTrajectory',true,...
                'FileType',Settings.CoordType,...
                'ML_TimeLength',0,...
                'ML_TimeStep',0,...
                'SaveTrajectory',true,...
                'SavePredictionsImage',true));
            Sol_Fraction = PyOut{4};
            Liq_Fraction = PyOut{5};
            
            if Liq_Fraction < (1 - Settings.MeltFreezeThreshold)
                if Settings.Verbose
                    disp('Detected Liquid Phase change.')
                end
                if (1-Liq_Fraction-Sol_Fraction) >= Settings.MeltFreezeThreshold
                    Output.StructureChange = true;
                else
                    Output.LiquidFroze = true;
                end
                Output.Aborted = true;
                TDir = fullfile(strrep(MinDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
                [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
                [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.OuterDir) ' -iname "*core*" ' Settings.pipe ' xargs rm']);
                copyfile(Settings.WorkDir,TDir)
                return
            end

            MSD_File = fullfile(Settings.WorkDir,'Equil_Liq_MSD.xvg');
            MSD_Log_File = fullfile(Settings.WorkDir,'Equil_Liq_MSD.log');
            msd_command = [Settings.wsl 'echo ' Settings.Metal ' ' Settings.pipe ' '  strrep(Settings.gmx_loc,Settings.wsl,'') ' msd -f ' windows2unix(TRR_File) ...   
                ' -s ' windows2unix(TPR_File) ' -o ' windows2unix(MSD_File) ' -b 0 -e ' num2str(Settings.Liquid_Test_Time) ...
                ' -trestart 0.1 -beginfit 1 -endfit ' num2str(0.75*Settings.Liquid_Test_Time) Settings.passlog windows2unix(MSD_Log_File)];
            [~,~] = system(msd_command);
            outp = fileread(MSD_Log_File);
            Diff_txt = regexp(outp,['D\[ *' Settings.Metal '] *([0-9]|\.|e|-)+ *(\(.+?\)) *([0-9]|\.|e|-)+'],'tokens','once');
            Diff_const = str2double(Diff_txt{1})*str2double(Diff_txt{3}); % cm^2 / s

            Exp = Load_Experimental_Data;
            if Diff_const <= Exp.(Settings.Salt).Liquid.DM_mp/100
                if Settings.Verbose
                    disp('Detected liquid has hardened to amorphous solid.')
                end
                Output.LiquidAmorphous = true;
                Output.Aborted = true;
                TDir = fullfile(strrep(MinDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
                [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
                [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.OuterDir) ' -iname "*core*" ' Settings.pipe ' xargs rm']);
                copyfile(Settings.WorkDir,TDir);
                try
                    if Settings.Delete_Equil
                        rmdir(Settings.WorkDir,'s')
                    end
                catch
                    if Settings.Verbose
                        disp(['Unable to remove directory: ' Settings.WorkDir])
                    end
                end
                return
            end
        else
            if Settings.Verbose
                disp(['(2/2) Liquid Dynamics Failed! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
            end
            try
                % Attempt to check to ensure system remained liquid
                copyfile(Equilibrated_Geom_File,fullfile(Settings.WorkDir,['MSD_Liq.' Settings.CoordType]));
                PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(Settings.WorkDir, Settings.Salt, ...
                    pyargs('SystemName','MSD_Liq',...
                    'RefStructure',Settings.Structure,...
                    'CheckFullTrajectory',true,...
                    'FileType',Settings.CoordType,...
                    'ML_TimeLength',0,...
                    'ML_TimeStep',0,...
                    'SaveTrajectory',true,...
                    'SavePredictionsImage',true));
                Sol_Fraction = PyOut{4};
                Liq_Fraction = PyOut{5};

                if Liq_Fraction < (1 - Settings.MeltFreezeThreshold)
                    if Settings.Verbose
                        disp('Detected Liquid Phase change.')
                    end
                    if (1-Liq_Fraction-Sol_Fraction) >= Settings.MeltFreezeThreshold
                        Output.StructureChange = true;
                    else
                        Output.LiquidFroze = true;
                    end
                    Output.Aborted = true;
                    TDir = fullfile(strrep(MinDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
                    [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
                    [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.OuterDir) ' -iname "*core*" ' Settings.pipe ' xargs rm']);
                    copyfile(Settings.WorkDir,TDir)
                    return
                end
            catch
                if Settings.Verbose
                    disp('Unable to check liquid dynamics trajectory for phase change!');
                end
            end
        end
    end
    
    try
        if Settings.Delete_Equil
            rmdir(Settings.WorkDir,'s')
        end
    catch
        if Settings.Verbose
            disp(['Unable to remove directory: ' Settings.WorkDir])
        end
    end
    
    if Settings.Verbose
        disp('*** Separate Equilibration of Liquid Complete ***')
    end
    
    %% Liquid is now equilibrated, attach solid + liquid and minimize interface
    
    if Settings.Verbose
        disp('*** Attaching Equilibrated Liquid to Solid and Minimizing Interface ***')
    end
    Settings.WorkDir = MinDir;
    
    Solid_file_data = load_gro_file(Settings.SuperCellFile);
    Liquid_file_data = load_gro_file(Equilibrated_Geom_File);
    
    MinInterfaceWidth = Settings.Geometry.c/10; % Convert to nm. Minimize the liquid one unit cell above/below the solid
    
    Liq_Z = norm(Liquid_file_data.c_vec);
    Liq_atoms_at_edge = Liquid_file_data.xyz(:,3) < MinInterfaceWidth | Liquid_file_data.xyz(:,3) > (Liq_Z - MinInterfaceWidth); % index of atoms within +-1 unit cell of of liquid-solid interfaces
%     Sol_Z = norm(Solid_file_data.c_vec);
%     Sol_atoms_at_edge = Solid_file_data.xyz(:,3) < MinInterfaceWidth | Solid_file_data.xyz(:,3) > (Sol_Z - MinInterfaceWidth); % index of atoms within +-1 angstrom of liquid-solid interfaces
    
    buffer_vec = 0.1.*(Liquid_file_data.c_vec)/norm(Liquid_file_data.c_vec); % 1 Angstrom buffer
    
    a_vec = Solid_file_data.a_vec;
    b_vec = Solid_file_data.b_vec;
    c_vec = Solid_file_data.c_vec + Liquid_file_data.c_vec + buffer_vec;
    
    Combined_file_data = Solid_file_data;
    Combined_file_data.N_atoms = Solid_file_data.N_atoms + Liquid_file_data.N_atoms;
    Combined_file_data.c_vec = c_vec;
    Combined_file_data.xyz = [Solid_file_data.xyz; (Liquid_file_data.xyz + Solid_file_data.c_vec + 0.5.*buffer_vec)];
    Combined_file_data.vel = [Liquid_file_data.vel; Liquid_file_data.vel];
    Combined_file_data.res_number = [Solid_file_data.res_number; (Liquid_file_data.res_number + Solid_file_data.res_number(end))];
    Combined_file_data.res_name = [Solid_file_data.res_name; Liquid_file_data.res_name];
    Combined_file_data.atom_name = [Solid_file_data.atom_name; Liquid_file_data.atom_name];
    Combined_file_data.atom_number = [Solid_file_data.atom_number; (Liquid_file_data.atom_number + Solid_file_data.atom_number(end))];
    Combined_file_data.boxcoords = {a_vec(1) b_vec(2) c_vec(3) a_vec(2) a_vec(3) b_vec(1) b_vec(3) c_vec(1) c_vec(2)};
    Combined_file_data.Salt = 'NaCl';
    Combined_file_data.N = Combined_file_data.N_atoms;
    
    Comb_Equil_Geom_File = fullfile(Settings.WorkDir,['Comb_Equil.' Settings.CoordType]);
    SaveGroFile(Comb_Equil_Geom_File,Combined_file_data,true);
    
    Solid_Geometry.Liquid_Vol = Cell_Volume(Supercell_file_data.a_vec,...
        Supercell_file_data.b_vec,Liq_Z.*(c_vec./norm(c_vec))); % nm^3
    Solid_Geometry.Solid_Vol = Cell_Volume(Supercell_file_data.a_vec,...
        Supercell_file_data.b_vec,Supercell_file_data.c_vec); % nm^3
    Solid_Geometry.Solid_R = ( 3*Solid_Geometry.Solid_Vol/(4*pi) )^(1/3); % cluster radius in nm
    
    % Save solid-only supercell file info
    SS_file_info = dir(Settings.SuperCellFile);
    Solid_Geometry.a_vec = Supercell_file_data.a_vec;
    Solid_Geometry.b_vec = Supercell_file_data.b_vec;
    Solid_Geometry.c_vec = Supercell_file_data.c_vec;
    Solid_Geometry.c_vec_tot = c_vec;
    Solid_Geometry.L = Liq_Z; % Length of liquid box
    Solid_Geometry.N_atoms = Supercell_file_data.N_atoms;
    Solid_Geometry.Volume = Supercell_file_data.Volume;
    Solid_Geometry_file = fullfile(SS_file_info.folder,[Settings.JobName '_SolInfo.mat']);
    save(Solid_Geometry_file,'Solid_Geometry');
    
    % Freeze all atoms except those near the edge
    Atoms_at_Edge = Liq_atoms_at_edge; %[Sol_atoms_at_edge; Liq_atoms_at_edge];
    Freeze_atnums = double(Combined_file_data.atom_number(~Atoms_at_Edge));
    
    % Create an index file to keep track of frozen atom numbers
    system_atnums = double(Combined_file_data.atom_number);
    
    % Get indeces for the metal and halide
    M_idx = find(contains(Combined_file_data.atom_name,Settings.Metal));
    M_atnums = double(Combined_file_data.atom_number(M_idx));
    X_idx = find(contains(Combined_file_data.atom_name,Settings.Halide));
    X_atnums = double(Combined_file_data.atom_number(X_idx));
    
    % Create index file
    addnan = 15 - mod(length(system_atnums),15);
    system_atnums(end+1:end+addnan) = nan;
    system_atnums = reshape(system_atnums,15,[])';
    system_atnum_txt = strtrim(char(regexprep(strjoin(string(num2str(system_atnums)),newline),' +NaN','')));
    
    addnan = 15 - mod(length(Freeze_atnums),15);
    Freeze_atnums(end+1:end+addnan) = nan;
    Freeze_atnums = reshape(Freeze_atnums,15,[])';
    sol_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(Freeze_atnums)),newline),' +NaN','')));
    
    addnan = 15 - mod(length(M_atnums),15);
    M_atnums(end+1:end+addnan) = nan;
    M_atnums = reshape(M_atnums,15,[])';
    M_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(M_atnums)),newline),' +NaN','')));

    addnan = 15 - mod(length(X_atnums),15);
    X_atnums(end+1:end+addnan) = nan;
    X_atnums = reshape(X_atnums,15,[])';
    X_atnums_txt = strtrim(char(regexprep(strjoin(string(num2str(X_atnums)),newline),' +NaN','')));

    ndx_text = ['[ System ]'     newline system_atnum_txt newline ...
    '[ Freeze ]'      newline sol_atnums_txt newline ...
    '[ ' Settings.Metal ' ]'  newline M_atnums_txt newline ...
    '[ ' Settings.Halide ' ]' newline X_atnums_txt newline];

    % Save index file
    NDX_Filename = fullfile(Settings.WorkDir,'Comb_Equil.ndx');
    fidNDX = fopen(NDX_Filename,'wt');
    fwrite(fidNDX,regexprep(ndx_text,'\r',''));
    fclose(fidNDX);
    
    % Reuse the previous minimization file, freeze the solid and liquid by adding freeze groups to the mdp file
    Freezegrps_text = [newline newline '; Freezing Groups' newline ...
    'freezegrps               = Freeze' newline ...
    'freezedim                = Y Y Y'];
    MDP.Minimization_txt = [MDP.Minimization_txt Freezegrps_text];
    
    % Save MDP file
    MDP_Filename = fullfile(Settings.WorkDir,'Comb_Equil.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP.Minimization_txt,'\r',''));
    fclose(fidMDP);
    
    % Complete a topology file for the liquid box to be minimized
    Atomlist = copy_atom_order(Comb_Equil_Geom_File);
    Topology_Text = strrep(Topology_Text_Min,'##LATOMS##',Atomlist);
    Top_Filename = fullfile(Settings.WorkDir,'Comb_Equil.top');
    
    % Save topology file
    fidTOP = fopen(Top_Filename,'wt');
    fwrite(fidTOP,regexprep(Topology_Text,'\r',''));
    fclose(fidTOP);
    
    % Minimization preprocessing
    TPR_File = fullfile(Settings.WorkDir,'Comb_Equil.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Comb_Equil_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Comb_Equil_Grompplog.log');
    
    FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Comb_Equil_Geom_File) ...
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
        Settings.TableFile_MX = MakeTables(Settings);
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Remove previous supercell
    delete(Settings.SuperCellFile)
    
    % Final minimization, saving output to supercell file in outer directory
    if Settings.Verbose
        disp('Beginning Combined System minimization...')
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

%     system(['wsl source ~/.bashrc; echo "4 0" ^| gmx_d energy -f ' windows2unix(Energy_file) ' -o ' windows2unix(strrep(Energy_file,'.edr','.xvg'))])
%     system(['echo "4 0" | gmx_d energy -f ' Energy_file ' -o ' strrep(Energy_file,'.edr','.xvg')])
%     En_xvg_file = fullfile(Settings.WorkDir,'Comb_Equil.xvg');
%     Data = import_xvg(En_xvg_file);
%     plot(Data(:,1),Data(:,2)./(2*Settings.nmol_liquid)) % Potential
%     ylim([-1000 1000])
    if Settings.Verbose
        disp('*** Minimization of Solid-Liquid Interface Complete ***')
    end
    
    % Remove backups
    if Settings.Delete_Backups
        system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
            ' -iname "#*#" -delete']);
    end
    
end
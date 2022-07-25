function Equil_density = Equilibrate_Liquid(Settings)

    %% What this function does:
    % Initialize density based on experiment
    % Create a cubic box with the smallest possible system size for a given cutoff
    % Increase that box size by Settings.Cutoff_Buffer to account for possibility of shrinkage
    % Randomly fill box with appropriate number of atoms to reach experimental density
    % Use Make_Tables function if necessary, run a brief ~1000 step minimization
    % Run NPT simulation for [Settings.Equilibrate_Liquid] amount of time using fast equilibration settings (Berendsen baro and Velocity-Rescale thermo)
    % Calculate average density of equilibrated box based on last 25% of simulation
    % Give new density as output

    %% Possibility 1 Implemented
    disp('*** Separate Equilibration of Liquid Selected ***')
    
    Settings.WorkDir = fullfile(Settings.WorkDir,'Equil_Liq');
    if ~isfolder(Settings.WorkDir)
        mkdir(Settings.WorkDir)
    end

    % Grab reference density, cutoff, and corresponding box size
    Ref_Density = Get_LiX_Liquid_Density(Settings); % molecules/nm^3
    L = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer; % nm, the box dimension
    Volume = L^3;
    nmol_liquid = round(Volume*Ref_Density);

    % Initialize empty box
    Box.a_vec  = [L 0 0];
    Box.b_vec  = [0 L 0];
    Box.c_vec  = [0 0 L];
    Box.boxcoords = {L L L 0 0 0 0 0 0};
    Box.Volume = Volume;
    Box.Salt = Settings.Salt;
    Box.N = 0;
    Box.N_atoms = 0;
    Box.comment = [Settings.Salt ' liquid cell'];
    Box.res_number = double.empty(0,1);
    Box.res_name = cell.empty(0,1);
    Box.atom_name = cell.empty(0,1);
    Box.atom_number = int32.empty(0,1);
    Box.xyz = double.empty(0,3);
    Box.vel = double.empty(0,3);

    % Save empty box
    Liq_Box_File = fullfile(Settings.WorkDir,['Prep_Liq.' Settings.CoordType]);
    SaveGroFile(Liq_Box_File,Box,true);

    Ref_M = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Metal '_Box.gro']);
    Ref_X = fullfile(Settings.home,'templates','GRO_Templates',[Settings.Halide '_Box.gro']);

    % Add metal
    disp(['Randomly adding ' num2str(nmol_liquid) ' ' Settings.Metal ' ions to liquid box...'])
    mtimer = tic;
    Prep_Liq_Metal_Only = fullfile(Settings.WorkDir,['Prep_Liq_Metal_Only.' Settings.CoordType]);
    cmd = [Settings.gmx_loc ' insert-molecules -f ' windows2unix(Liq_Box_File) ' -ci ' windows2unix(Ref_M) ...
        ' -o ' windows2unix(Prep_Liq_Metal_Only) ' -nmol ' num2str(nmol_liquid) ' -try 200'];
    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Settings.Metal ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Settings.Metal ' atoms added. Epalsed Time: ' datestr(seconds(toc(mtimer)),'HH:MM:SS')])

    % Add Halide
    Prep_Liq_Random_Liq = fullfile(Settings.WorkDir,['Prep_Liq_Random_Liq.' Settings.CoordType]);
    disp(['Randomly adding ' num2str(nmol_liquid) ' ' Settings.Halide ' ions to liquid box...'])
    htimer = tic;
    cmd = [Settings.gmx_loc ' insert-molecules -ci ' windows2unix(Ref_X) ' -f ' windows2unix(Prep_Liq_Metal_Only) ...
        ' -o ' windows2unix(Prep_Liq_Random_Liq) ' -nmol ' num2str(nmol_liquid) ' -try 400'];

    [errcode,output] = system(cmd);

    if errcode ~= 0
        disp(output);
        error(['Error adding ' Settings.Halide ' atoms with insert-molecules. Problem command: ' newline cmd]);
    end
    disp([Settings.Halide ' atoms added. Epalsed Time: ' datestr(seconds(toc(htimer)),'HH:MM:SS')])


    % Load current randomly-generated liquid file data
    Random_Liquid_data = load_gro_file(Prep_Liq_Random_Liq);

    % Check to make sure all atoms were added
    if Random_Liquid_data.N_atoms ~= nmol_liquid*2
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
        disp('Warning: enabling long-range dispersion correction for tabulated potential!')
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
        Topology_Text = strrep(Topology_Text,'##LATOMS##',Atomlist);
    else
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
    Minimized_Geom_File = fullfile(Settings.WorkDir,['Minimized_Liq.' Settings.CoordType]);

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Minimized_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Prep_Liq')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req || strcmp(Settings.Theory,'BH')
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end
    
    % Liquid Minimization
    disp('Begining Liquid Minimization...')
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        disp(['System Successfully Minimized! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
    else
        disp(mdrun_output);
        error(['Error running mdrun for liquid system minimization. Problem command: ' newline mdrun_command]);
    end
    
    %% System is now minimized, run a fast equilibration
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
    Compressibility = Settings.QECompressibility; % bar^(-1)
    tau_p = Settings.MDP.dt; % ps
    tau_t = Settings.MDP.dt; % ps

    nstpcouple = max(round(tau_p/(20*Settings.MDP.dt)),1);
    nsttcouple = max(round(tau_t/(20*Settings.MDP.dt)),1);

    % Ensure fast equilibration with Berendsen barostat + small time constant
    MDP_Template = regexprep(Settings.MDP_Template,'(nsteps += *)(.+?)( *);',['$1' num2str(timesteps) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = regexprep(MDP_Template,'(pcoupl += *)(.+?)( *);','$1Berendsen$3;');
    MDP_Template = regexprep(MDP_Template,'(pcoupltype += *)(.+?)( *);','$1isotropic$3;');
    MDP_Template = regexprep(MDP_Template,'(tau-p += *)(.+?)( *);',['$1 ' num2str(tau_p) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstpcouple += *)(.+?)( *);',['$1 ' num2str(nstpcouple) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(compressibility += *)(.+?)( *);',['$1 ' num2str(Compressibility) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(ref-p += *)(.+?)( *);',['$1 ' num2str(Settings.Target_P(1)) '$3;']);

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
    Equilibrated_Geom_File = fullfile(Settings.WorkDir,['Equil_Liq_out.' Settings.CoordType]);

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Equilibrated_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Equil_Liq')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Run Liquid Equilibration
    disp(['Begining Liquid Equilibration for ' num2str(Settings.Equilibrate_Liquid) ' ps...'] )
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        disp(['Liquid Successfully Equilibrated! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
    else
        disp('Equilibration failed. Stopping.')
        disp(mdrun_output);
        error(['Error running mdrun for solid equilibration. Problem command: ' newline mdrun_command]);
    end

    En_xvg_file = fullfile(Settings.WorkDir,'Equil_Liq_Energy.xvg');

    % Check energy options
    gmx_command = [strrep(Settings.gmx_loc,'gmx',['echo 0 ' Settings.pipe ' gmx']) ...
        ' energy -f ' windows2unix(Energy_file) ...
        ' -s ' windows2unix(TPR_File)];
    [~,outpt] = system(gmx_command);

    en_opts = regexp(outpt,'-+\n.+?-+\n','match','once');
    En_set = '';
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
    %En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Potential','tokens','once'))];
    En_set = [En_set ' 0'];
    En_set = regexprep(En_set,' +',' ');

    % Grab data from results
    startpoint = Settings.Equilibrate_Liquid*0.5; % ps
    gmx_command = [strrep(Settings.gmx_loc,'gmx',['echo' En_set ' ' Settings.pipe ' gmx']) ...
    ' energy -f ' windows2unix(Energy_file)...
    ' -o ' windows2unix(En_xvg_file) ' -s ' windows2unix(TPR_File) ...
    ' -b ' num2str(startpoint) ' -e ' num2str(Settings.Equilibrate_Liquid)];

    [err,~] = system(gmx_command);
    if err ~= 0
        warndlg('Failed to collect data.')
        return
    end

    Data = import_xvg(En_xvg_file); % Gather X,Y,Z lengths
    delete(En_xvg_file) % remove temp output file

    Equil_Volume = mean(Data(:,2)); % nm^3

    % plot(Data(:,1),(10^3).*Data(:,2)./nmol_liquid)
    % V = mean((10^3).*Data(timesteps/2:end,2)./nmol_liquid) % A^3/molecule
    % stdevV = std((10^3).*Data(timesteps/2:end,2)./nmol_liquid) % A^3/molecule

    Equil_density = nmol_liquid/Equil_Volume;

    if Settings.Delete_Equil
        rmdir(Settings.WorkDir,'s')
    end
    
    disp('*** Separate Equilibration of Liquid Complete ***')

end
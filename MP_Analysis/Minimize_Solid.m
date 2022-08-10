function Minimize_Solid(Settings)
    
    if Settings.Verbose
        disp('*** Full Minimization of Solid ***')
    end
    Supercell_file_data = load_gro_file(Settings.SuperCellFile);
    Settings.JobName = 'Min_Sol';
    
    %% Minimize the solid: build MDP file
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
    
    Settings.Table_Req = IsGmxTableRequired(Settings);
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
    MDP_Filename = fullfile(Settings.WorkDir,'Min_Sol.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP.Minimization_txt,'\r',''));
    fclose(fidMDP);
    
    % Complete a topology file for the solid box to be minimized
    Atomlist = copy_atom_order(Settings.SuperCellFile);
    Top_Filename = fullfile(Settings.WorkDir,'Min_Sol.top');
    
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
        Topology_Text = strrep(Topology_Text,'##N##x##N##x##N##',num2str(Supercell_file_data.N_atoms));
        Topology_Text = strrep(Topology_Text,'##GEOM##','atom solid');
        
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

    TPR_File = fullfile(Settings.WorkDir,'Min_Sol.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Min_Sol_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Min_Sol_Grompplog.log');
    
    FMin_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Settings.SuperCellFile) ...
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
    Log_File = fullfile(Settings.WorkDir,'Min_Sol.log');
    Energy_file = fullfile(Settings.WorkDir,'Min_Sol.edr');
    TRR_File = fullfile(Settings.WorkDir,'Min_Sol.trr');

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Settings.SuperCellFile) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Min_Sol')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req || strcmp(Settings.Theory,'BH')
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Remove previous supercell
    delete(Settings.SuperCellFile)

    % Final minimization, output geometry is the new supercell file
    if Settings.Verbose
        disp('Begining Solid Minimization...')
    end
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        if Settings.Verbose
            disp(['System Successfully Minimized! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
        end
    else
        disp(mdrun_output);
        error(['Error running mdrun for system minimization. Problem command: ' newline mdrun_command]);
    end
    
    if Settings.Verbose
        disp('*** Minimization of Solid Complete ***')
    end
    
    % Remove minimization temporary folder
    if Settings.Delete_Backups
        system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ...
            ' -iname "#*#" -delete']);
    end
    
end
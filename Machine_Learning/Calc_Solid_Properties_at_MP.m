function Output = Calc_Solid_Properties_at_MP(Settings,varargin)
    p = inputParser;
    p.FunctionName = 'Calc_Solid_Properties_at_MP';
    addOptional(p,'Verbose',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
    addOptional(p,'Skip_Cell_Construction',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
    
    parse(p,varargin{:});
    Verbose = p.Results.Verbose;
    Settings.Skip_Cell_Construction =  p.Results.Skip_Cell_Construction;
	Inp_Settings = Settings;
    
    if ~isfield(Settings,'WorkDir')
        Settings.WorkDir = GetMDWorkdir(Settings);
        Settings.WorkDir = [Settings.WorkDir '_SP'];
    end
    if ~isfolder(Settings.WorkDir)
        mkdir(Settings.WorkDir)
    end
    save(fullfile(Settings.WorkDir,'Calc_Settings.mat'),'Settings')
    
    if Verbose
        disp('*** Separate Equilibration of Solid Selected ***')
    end
    
    % Generate unit cell: Unit Cell Filename
    Settings.JobName = 'Equil_Sol';
    
    if ~Settings.Skip_Cell_Construction
        Settings.UnitCellFile = fullfile(Settings.WorkDir,[Settings.JobName '_UnitCell.' Settings.CoordType]);
        Settings.Geometry = Convert_to_Primitive(Settings.Geometry);
        if contained_in_cell(Settings.Structure,{'Rocksalt' 'Sphalerite' 'FiveFive'})
            Settings.Coordinate_File = fullfile(Settings.home,'templates',[upper(Settings.CoordType) '_Templates'],...
                [Settings.Structure '_Primitive.' Settings.CoordType]);
        else
            Settings.Coordinate_File = fullfile(Settings.home,'templates',[upper(Settings.CoordType) '_Templates'],...
                [Settings.Structure '.' Settings.CoordType]);
        end

        % Generate unit cell coord file text
        Coordinate_Text = fileread(Settings.Coordinate_File);

        % Add Metal and Halide symbols
        if strcmp(Settings.CoordType,'gro')
            Met = pad(Settings.Metal,2,'left');
            Hal = pad(Settings.Halide,2,'left');
        elseif strcmp(Settings.CoordType,'g96')
            Met = pad(Settings.Metal,2,'right');
            Hal = pad(Settings.Halide,2,'right');
        end
        Coordinate_Text = strrep(Coordinate_Text,'##MET##',Met);
        Coordinate_Text = strrep(Coordinate_Text,'##HAL##',Hal);
        Coordinate_Text = AddCartesianCoord(Coordinate_Text,Settings.Geometry,1,false,Settings.CoordType); %input coordinates

        % Save unit cell .gro file
        fid = fopen(Settings.UnitCellFile,'wt');
        fwrite(fid,regexprep(Coordinate_Text,'\r',''));
        fclose(fid);

        % Generate a box with the structure of interest containing the smallest
        % possible number of atoms for the given cutoff, plus a buffer in case of contraction
        La = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_a; % nm, the minimum box dimension
        Lb = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_b; % nm, the minimum box dimension
        Lc = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_c; % nm, the minimum box dimension

        Na = ceil(La/(Settings.Geometry.a/10));
        Nb = ceil(Lb/(Settings.Geometry.b/10));
        Nc = ceil(Lc/(Settings.Geometry.c/10));

        % Calculate number of formula units
        nmol_solid = Na*Nb*Nc*Settings.Geometry.NF;
        
        while nmol_solid < 500 % Enforce a minimum of 1000 atoms
            La = La*1.1;
            Lb = Lb*1.1;
            Lc = Lc*1.1;
            Na = ceil(La/(Settings.Geometry.a/10));
            Nb = ceil(Lb/(Settings.Geometry.b/10));
            Nc = ceil(Lc/(Settings.Geometry.c/10));
            nmol_solid = Na*Nb*Nc*Settings.Geometry.NF;
        end

        Settings.SuperCellFile = fullfile(Settings.WorkDir,['Equil_Sol.' Settings.CoordType]);
        Supercell_command = [Settings.gmx_loc ' genconf -f ' windows2unix(Settings.UnitCellFile) ...
             ' -o ' windows2unix(Settings.SuperCellFile) ' -nbox ' num2str(Na) ' ' num2str(Nb) ' ' num2str(Nc)];
        [errcode,output] = system(Supercell_command);

        if errcode ~= 0
            disp(output);
            error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
        end
    else
        Supercell_file_data = load_gro_file(Settings.SuperCellFile);
        nmol_solid = Supercell_file_data.N_atoms/2;
    end
    
    % Set the number of steps
    MD_nsteps = Settings.Solid_Test_Time/Settings.MDP.dt;
    %Compressibility = Get_Alkali_Halide_Compressibility(Settings.Salt);
    Compressibility = Settings.QECompressibility;
    tau_p = Settings.MDP.dt; % ps
    tau_t = Settings.MDP.dt; % ps

    nstpcouple = max(round(tau_p/(20*Settings.MDP.dt)),1);
    nsttcouple = max(round(tau_t/(20*Settings.MDP.dt)),1);
    
    if abs(Settings.Geometry.a - Settings.Geometry.b) < sqrt(eps) && abs(Settings.Geometry.a - Settings.Geometry.c) < sqrt(eps) % a = b = c
        isotropy = 'isotropic';
        ref_p = num2str(Settings.Target_P(1));
        Compresstxt = num2str(Compressibility); % bar^(-1)
    elseif abs(Settings.Geometry.a - Settings.Geometry.b) < sqrt(eps) % a = b != c
        isotropy = 'semiisotropic';
        ref_p = regexprep(num2str(ones(1,2).*Settings.Target_P(1)),' +',' ');
        Compresstxt = regexprep(num2str(ones(1,2).*Compressibility),' +',' '); % bar^(-1)
    else % a != b != c
        isotropy = 'anisotropic';
        ref_p = regexprep(num2str(ones(1,6).*Settings.Target_P(1)),' +',' ');
        Compresstxt = regexprep(num2str([ones(1,3) zeros(1,3)].*Compressibility),' +',' '); % bar^(-1)
    end
    
    Table_Req = IsGmxTableRequired(Settings);
    Metal_Info = elements('Sym',Settings.Metal);
    Halide_Info = elements('Sym',Settings.Halide);
    
    % Load new topology and MPD templates
    Topology_Template_file = fullfile(Settings.home,'templates','Gromacs_Templates',...
    'Topology.template');
    MDP_Template = fileread(fullfile(Settings.home,'templates','Gromacs_Templates',...
    'MDP_MD.template'));
    
    % Topology filename and directory for output
    Settings.Topology_File = fullfile(Settings.WorkDir,[Settings.JobName '.top']);
    Settings.Topology_Text = fileread(Topology_Template_file);
    
    % Add in global parameters to Topology template
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##GENPAIRS##',Settings.Top_gen_pairs);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##FUDGELJ##',num2str(Settings.Top_fudgeLJ));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##FUDGEQQ##',num2str(Settings.Top_fudgeQQ));
    
    % Insert element info into topology template
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##MET##',pad(Settings.Metal,2));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METZ##',pad(num2str(Metal_Info.atomic_number),3));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMASS##',pad(num2str(Metal_Info.atomic_mass),7));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##MCHRG##',pad(num2str(Settings.S.Q),2));
    
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HAL##',pad(Settings.Halide,2));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALZ##',pad(num2str(Halide_Info.atomic_number),3));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALMASS##',pad(num2str(Halide_Info.atomic_mass),7));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##XCHRG##',pad(num2str(-Settings.S.Q),2));
    
    if Table_Req
        % Define the function type as 1 (required for custom functions)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');
        
        % Define the combination rules (Lorenz-berthelot)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','1');
        
        % Define all the parameters as 1.0 (already included in potentials)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad('1.0',10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad('1.0',10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad('1.0',10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##','1.0');
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##','1.0');
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##','1.0');
        
        % Generate tables of the potential
        if strcmp(Settings.Theory,'TF')
            [U_MX, U_MM, U_XX] = TF_Potential_Generator(Settings);
        elseif strcmp(Settings.Theory,'BH')
            [U_MX, U_MM, U_XX] = BH_Potential_Generator(Settings);
        elseif contains(Settings.Theory,'JC')
            switch Settings.Theory
                case 'JC'
                    Settings.WaterModel = 'SPC/E';
                case 'JC3P'
                    Settings.WaterModel = 'TIP3P';
                case 'JC4P'
                    Settings.WaterModel = 'TIP4PEW';
                case 'JCSD'
                    Settings.WaterModel = 'SD';
            end
            [U_MX, U_MM, U_XX] = JC_Potential_Generator(Settings);
        else
            error(['Warning: Unknown theory type: "' Settings.Theory '".'])
        end
        
        TableName = [Settings.JobName '_Table'];
        Settings.TableFile_MX = fullfile(Settings.WorkDir,[TableName '.xvg']);
        Settings.TableFile_MM = fullfile(Settings.WorkDir,[TableName '_' Settings.Metal '_' Settings.Metal '.xvg']);
        Settings.TableFile_XX = fullfile(Settings.WorkDir,[TableName '_' Settings.Halide '_' Settings.Halide '.xvg']);
        
        % Save tables into current directory
        fidMX = fopen(Settings.TableFile_MX,'wt');
        fwrite(fidMX,regexprep(U_MX,'\r',''));
        fclose(fidMX);
        
        fidMM = fopen(Settings.TableFile_MM,'wt');
        fwrite(fidMM,regexprep(U_MM,'\r',''));
        fclose(fidMM);
        
        fidXX = fopen(Settings.TableFile_XX,'wt');
        fwrite(fidXX,regexprep(U_XX,'\r',''));
        fclose(fidXX);
        
        % Modify the MDP file
        MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad('user',18));
        MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
        MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
        MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
        MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
        MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad('none',18));
        
    elseif contains(Settings.Theory,'JC')
        switch Settings.Theory
            case 'JC'
                Settings.WaterModel = 'SPC/E';
            case 'JC3P'
                Settings.WaterModel = 'TIP3P';
            case 'JC4P'
                Settings.WaterModel = 'TIP4PEW';
            case 'JCSD'
                Settings.WaterModel = 'SD';
        end
        
        % Definte the function type as 1 (LJ)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','1');
        
        % Define the combination rules (Lorenz-berthelot in sigma-epsilon form)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','2');
        
        % Get JC parameters
        [MX_JC_Param,MM_JC_Param,XX_JC_Param] = JC_Potential_Parameters(Settings);
        
        % Add parameters to topology text
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',pad(num2str(MM_JC_Param.sigma,'%10.8e'),10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',pad(num2str(XX_JC_Param.sigma,'%10.8e'),10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',pad(num2str(MX_JC_Param.sigma,'%10.8e'),10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',num2str(MM_JC_Param.epsilon,'%10.8e'));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',num2str(XX_JC_Param.epsilon,'%10.8e'));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',num2str(MX_JC_Param.epsilon,'%10.8e'));
        
        % Modify the MDP file
        MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
        MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad(Settings.MDP.CutOffScheme,18));
        MDP_Template = regexprep(MDP_Template,'energygrp-table.+?\n','');
        MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
        MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
        MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));

        % Add in Verlet Settings
        if strcmp(Settings.MDP.CutOffScheme,'Verlet')
            MDP_Template = strrep(MDP_Template,'##VerletBT##',pad(num2str(Settings.MDP.VerletBT),18));
        else
            MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
        end
    elseif contains(Settings.Theory,'BH')
        
        Settings.TableFile_MX = '';
        
        % Definte the function type as 2 (Buckingham)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##NBFUNC##','2');
        
        % Define the combination rule as 1 (Buckingham only has 1 comb rule)
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##COMBR##','1');
        
        % Get BH parameters (cross terms are pre-computed using my combining rules)
        [U_MX,U_MM,U_XX] = BH_Potential_Parameters(Settings);
        
        % Add parameters to topology text
        % For BH potentials, parameter are B*exp(-alpha*r) + C/r^6
        % Parameter order is B alpha C
        Settings.Topology_Text = strrep(Settings.Topology_Text,'ptype  C          A','ptype   a              b           c6');
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETC##',[num2str(U_MM.B,'%10.8e') ' ' num2str(U_MM.alpha,'%10.8e')]);
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METMETA##',pad(num2str(U_MM.C,'%10.8e'),10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALC##',[num2str(U_XX.B,'%10.8e') ' ' num2str(U_XX.alpha,'%10.8e')]);
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##HALHALA##',pad(num2str(U_XX.C,'%10.8e'),10));
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALC##',[num2str(U_MX.B,'%10.8e') ' ' num2str(U_MX.alpha,'%10.8e')]);
        Settings.Topology_Text = strrep(Settings.Topology_Text,'##METHALA##',pad(num2str(U_MX.C,'%10.8e'),10));
        
        % Modify the MDP file
        MDP_Template = strrep(MDP_Template,'##VDWTYPE##',pad(Settings.MDP.VDWType,18));
        MDP_Template = strrep(MDP_Template,'##CUTOFF##',pad('group',18));
        MDP_Template = regexprep(MDP_Template,'energygrp-table.+?\n','');
        MDP_Template = regexprep(MDP_Template,'ewald-rtol-lj.+?\n','');
        MDP_Template = regexprep(MDP_Template,'lj-pme-comb-rule.+?\n','');
        MDP_Template = strrep(MDP_Template,'##RLIST##',pad(num2str(Settings.MDP.RList_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RCOULOMB##',pad(num2str(Settings.MDP.RCoulomb_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##RVDW##',pad(num2str(Settings.MDP.RVDW_Cutoff),18));
        MDP_Template = strrep(MDP_Template,'##VDWMOD##',pad(Settings.MDP.vdw_modifier,18));
        MDP_Template = regexprep(MDP_Template,'verlet-buffer-tolerance.+?\n','');
    else
        error(['Warning: Unknown theory type: "' Settings.Theory '".'])
    end
    
    if Settings.MDP.Disp_Correction && ~Table_Req
        MDP_Template = [MDP_Template newline newline ...
            '; Long-range dispersion correction' newline ...
            'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
    elseif Settings.MDP.Disp_Correction && Settings.MDP.Disp_Correction_Tables
        disp('Warning: enabling long-range dispersion correction for tabulated potential!')
        MDP_Template = [MDP_Template newline newline ...
            '; Long-range dispersion correction' newline ...
            'DispCorr                 = EnerPres          ; apply long range dispersion corrections for Energy and pressure'];
    end
    
    % Add in parameters to MDP template
    MDP_Template = strrep(MDP_Template,'##NSTEPS##',pad(num2str(MD_nsteps),18));
    MDP_Template = strrep(MDP_Template,'##INTEGR##',pad(Settings.MDP.integrator,18));
    MDP_Template = strrep(MDP_Template,'##TIMEST##',pad(num2str(Settings.MDP.dt),18));
    MDP_Template = strrep(MDP_Template,'##CONTINUE##',pad(Settings.MDP.continuation,18));
    MDP_Template = strrep(MDP_Template,'##REFTINIT##',pad(num2str(Settings.T0),18));
    MDP_Template = strrep(MDP_Template,'##COULOMB##',pad(Settings.MDP.CoulombType,18));
    MDP_Template = strrep(MDP_Template,'##FOURIER##',pad(num2str(Settings.MDP.Fourier_Spacing),18));
    MDP_Template = strrep(MDP_Template,'##PMEORDER##',pad(num2str(Settings.MDP.PME_Order),18));
    MDP_Template = strrep(MDP_Template,'##EWALDTOL##',pad(num2str(Settings.MDP.Ewald_rtol),18));
    MDP_Template = strrep(MDP_Template,'##LISTUPDATE##',pad(num2str(Settings.Update_NeighbourList),18));
    MDP_Template = strrep(MDP_Template,'##POSOUT##',pad(num2str(Settings.Output_Coords),18));
    MDP_Template = strrep(MDP_Template,'##POSOUTCOMP##',pad(num2str(Settings.Output_Coords_Compressed),18));
    MDP_Template = strrep(MDP_Template,'##VELOUT##',pad(num2str(Settings.Output_Velocity),18));
    MDP_Template = strrep(MDP_Template,'##FORCEOUT##',pad(num2str(Settings.Output_Forces),18));
    MDP_Template = strrep(MDP_Template,'##ENOUT##',pad(num2str(Settings.Output_Energies),18));
    MDP_Template = strrep(MDP_Template,'##CALCE##',pad(num2str(Settings.Calc_Energies),18));
    MDP_Template = strrep(MDP_Template,'##MET##',pad(Settings.Metal,3));
    MDP_Template = strrep(MDP_Template,'##HAL##',pad(Settings.Halide,3));
    
    % Ensure fast equilibration with Berendsen barostat + small time constant
    MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = strrep(MDP_Template,'##BAROSTAT##',pad('Berendsen',18));
    MDP_Template = strrep(MDP_Template,'##ISOTROPY##',pad(isotropy,18));
    MDP_Template = strrep(MDP_Template,'##PTIMECONST##',pad(num2str(tau_p),18));
    MDP_Template = strrep(MDP_Template,'##NSTPCOUPLE##',pad(num2str(nstpcouple),18));
    MDP_Template = strrep(MDP_Template,'##COMPRESS##',pad(Compresstxt,18));
    MDP_Template = strrep(MDP_Template,'##REFP##',pad(ref_p,18));
    
    % Pair it with velocity rescale thermostat + small time constant
    MDP_Template = strrep(MDP_Template,'##THERMOSTAT##',pad('v-rescale',18));
    MDP_Template = strrep(MDP_Template,'##TTIMECONST##',pad(num2str(tau_t),18));
    MDP_Template = strrep(MDP_Template,'##NSTTCOUPLE##',pad(num2str(nsttcouple),18));
    MDP_Template = strrep(MDP_Template,'##REFT##',pad(num2str(Settings.T0),18));
    
    % Remove annealing inputs from MDP
    MDP_Template = strrep(MDP_Template,'##ANNEALING##',pad('no',18));
    MDP_Template = regexprep(MDP_Template,'annealing-npoints +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'annealing-time +.+?\n','');
    MDP_Template = regexprep(MDP_Template,'annealing-temp +.+?\n','');
    
    MDP_Template = strrep(MDP_Template,'##TITLE##',[Settings.Salt ' BH_TestModel']);
    
    % Save MDP file
    MDP_Filename = fullfile(Settings.WorkDir,'Equil_Sol.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
    fclose(fidMDP);
    
    % Finish the topology file: add in title and the atom list
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##N##x##N##x##N##',num2str(nmol_solid));
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##GEOM##',['molecule ' Settings.Structure]);
    Atomlist = copy_atom_order(Settings.SuperCellFile);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
    Top_Filename = fullfile(Settings.WorkDir,'Equil_Sol.top');
    
    % Save topology file
    fidTOP = fopen(Top_Filename,'wt');
    fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
    fclose(fidTOP);
    
    TPR_File = fullfile(Settings.WorkDir,'Equil_Sol.tpr');
    MDPout_File = fullfile(Settings.WorkDir,'Equil_Sol_out.mdp');
    GrompLog_File = fullfile(Settings.WorkDir,'Equil_Sol_Grompplog.log');
    
    FEquil_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(Settings.SuperCellFile) ...
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
    Log_File = fullfile(Settings.WorkDir,'Equil_Sol.log');
    Energy_file = fullfile(Settings.WorkDir,'Equil_Sol.edr');
    TRR_File = fullfile(Settings.WorkDir,'Equil_Sol.trr');
    Final_Geom_File = fullfile(Settings.WorkDir,['Equil_Sol_out.' Settings.CoordType]);

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Final_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(Settings.WorkDir,'Equil_Sol')) ...
        Settings.mdrun_opts];
    
    if Table_Req
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Run solid Equilibration
    if Verbose
        disp(['Beginning Solid Equilibration for ' num2str(Settings.Solid_Test_Time) ' ps...'] )
    end
    mintimer = tic;
    [state,~] = system(mdrun_command);
    if state == 0
        if Verbose
            disp(['Solid Successfully Equilibrated! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
        end
    else
        try % Clean up
            [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.WorkDir) ' -iname "#*#" ^| xargs rm']);
            [~,~] = system([Settings.wsl 'find ' windows2unix(Settings.OuterDir) ' -iname "*core*" ^| xargs rm']);
        catch me
            disp(me.message)
        end
        SuperCellFile = Settings.SuperCellFile;
        WorkDir = Settings.WorkDir;
        Settings = Inp_Settings;
        Settings.SuperCellFile = SuperCellFile;
        Settings.WorkDir = WorkDir;
        Settings.Verbose = Verbose;
        if ~isfield(Settings,'QECompressibility_init')
            Settings.QECompressibility_init = Settings.QECompressibility;
        end
%         if ~isfield(Settings,'MinComplete')
%             Settings.MinComplete = false;
%         end
        
        if Settings.QECompressibility > 1e-8 % Retry until compressibility is very tight
            if Verbose
                disp('Solid Equilibration failed. Retrying with stiffer compressibility.')
            end
            Settings.QECompressibility = Settings.QECompressibility/2;
            Output = Calc_Solid_Properties_at_MP(Settings,'Verbose',Verbose,'Skip_Cell_Construction',true);
            return
%         elseif ~Settings.MinComplete
%             if Verbose
%                 disp('Solid Equilibration failed. Stiffer compressibility did not resolve.')
%                 disp('Running Pre-Minimization of Solid.')
%             end
%             Minimize_Solid(Settings);
%             Settings.QECompressibility = Settings.QECompressibility_init;
%             Settings.MinComplete = true;
%             Output = Calc_Solid_Properties_at_MP(Settings,'Verbose',Verbose,'Skip_Cell_Construction',true);
%             return
        elseif Settings.MDP.dt > 1e-4
            if Verbose
                disp('Solid Equilibration failed. Stiffer compressibility did not resolve.')
                disp('Reducing time step.')
            end
            Settings.QECompressibility = Settings.QECompressibility_init;
            Settings.MDP.dt = Settings.MDP.dt/2;
            Settings.Output_Coords = Settings.Output_Coords*2;
            Output = Calc_Solid_Properties_at_MP(Settings,'Verbose',Verbose,'Skip_Cell_Construction',true);
            return
        else
            if Verbose
                disp('Solid equilibration failed.')
                disp('Solid may be completely unstable!')
                disp(['WorkDir: ' WorkDir])
            end
            Output.Solid_V_MP = nan;
            Output.Solid_H_MP = nan;
            return
        end
    end
    
    % Check to ensure system remained in the correct solid structure
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(Settings.WorkDir, Settings.Salt, ...
        pyargs('SystemName','Equil_Sol',...
        'RefStructure',Settings.Structure,...
        'CheckFullTrajectory',true,...
        'FileType',Settings.CoordType,...
        'ML_TimeLength',10,...
        'ML_TimeStep',1,...
        'TimePerFrame',Settings.Output_Coords*Settings.MDP.dt,...
        'SaveTrajectory',true,...
        'SavePredictionsImage',true));
    Sol_Fraction = PyOut{4};
    
    if Sol_Fraction < 0.85
        if Verbose
            disp('Detected Solid Phase change at Experimental MP')
        end
        if Settings.Delete_Equil
            try
                rmdir(Settings.WorkDir,'s')
            catch
                disp(['Unable to remove directory: ' Settings.WorkDir])
            end
        end
        Output.Solid_V_MP = nan;
        Output.Solid_H_MP = nan;
        return
    end
    
    En_xvg_file = fullfile(Settings.WorkDir,'Equil_Sol_Energy.xvg');
    
    % Check energy options
    gmx_command = [strrep(Settings.gmx_loc,'gmx',['echo 0 ' Settings.pipe ' gmx']) ...
        ' energy -f ' windows2unix(Energy_file) ...
        ' -s ' windows2unix(TPR_File)];
    [~,outpt] = system(gmx_command);
    
    en_opts = regexp(outpt,'-+\n.+?-+\n','match','once');
    En_set = '';
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Volume','tokens','once'))];
    En_set = [En_set ' ' char(regexp(en_opts,'([0-9]{1,2})  Enthalpy','tokens','once'))];
    En_set = [En_set ' 0'];
    En_set = regexprep(En_set,' +',' ');
    
    % Grab second half of data from results
    startpoint = Settings.Solid_Test_Time*0.5; % ps
    gmx_command = [strrep(Settings.gmx_loc,'gmx',['echo' En_set ' ' Settings.pipe ' gmx']) ...
    ' energy -f ' windows2unix(Energy_file)...
    ' -o ' windows2unix(En_xvg_file) ' -s ' windows2unix(TPR_File) ...
    ' -b ' num2str(startpoint) ' -e ' num2str(Settings.Solid_Test_Time)];
    
    [err,~] = system(gmx_command);
    if err ~= 0
        warndlg('Failed to collect data.')
        return
    end
    
    Data = import_xvg(En_xvg_file); % Gather X,Y,Z lengths
    
    Output.Solid_V_MP = mean(Data(:,2))*(10^3)/nmol_solid; % A^3 / ion pair
    Output.Solid_H_MP = mean(Data(:,3))/nmol_solid; % kJ/mol
    
%     plot(Data(:,1),(10^3).*Data(:,2)./nmol_solid)
%     plot(Data(:,1),Data(:,3)./nmol_solid)
%     V = mean((10^3).*Data(timesteps/2:end,2)./nmol_solid) % A^3/molecule
%     stdevV = std((10^3).*Data(timesteps/2:end,2)./nmol_solid) % A^3/molecule
    if Verbose
        disp('*** Separate Equilibration of Solid Complete ***')
    end
    
    if Settings.Delete_Equil
        try
            rmdir(Settings.WorkDir,'s')
        catch
            disp(['Unable to remove directory: ' Settings.WorkDir])
        end
    end
end
function Output = Equilibrate_Solid(Settings,varargin)
    p = inputParser;
    p.FunctionName = 'Equilibrate_Solid';
    addOptional(p,'Skip_Cell_Construction',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
    
    parse(p,varargin{:});
    Settings.Skip_Cell_Construction =  p.Results.Skip_Cell_Construction;
    
    Output.StructureChange = false;
    Output.SolidMelted = false;
    Output.LiquidFroze = false;
    Output.LiquidAmorphous = false;
    Output.Aborted = false;
    
	if Settings.Verbose
    	disp('*** Separate Equilibration of Solid Selected ***')
	end
    
    Inp_Settings = Settings;
    WorkDir = fullfile(Settings.WorkDir,'Equil_Sol');
    if ~isfolder(WorkDir)
        mkdir(WorkDir)
    end
    
    % Generate a box with the structure of interest containing the smallest
    % possible number of atoms for the given cutoff, plus a buffer in case of contraction
    La = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_a; % nm, the minimum box dimension
    Lb = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_b; % nm, the minimum box dimension
    Lc = (2*Settings.Longest_Cutoff)*Settings.Cutoff_Buffer/Settings.Geometry.Skew_c; % nm, the minimum box dimension

    Na = ceil(La/(Settings.Geometry.a/10));
    Nb = ceil(Lb/(Settings.Geometry.b/10));
    Nc = ceil(Lc/(Settings.Geometry.c/10));
        
    SuperCell_File = fullfile(WorkDir,['Equil_Sol.' Settings.CoordType]);
    
    if ~Settings.Skip_Cell_Construction
        Supercell_command = [Settings.gmx_loc ' genconf -f ' windows2unix(Settings.UnitCellFile) ...
             ' -o ' windows2unix(SuperCell_File) ' -nbox ' num2str(Na) ' ' num2str(Nb) ' ' num2str(Nc)];
        [errcode,output] = system(Supercell_command);

        if errcode ~= 0
            if Settings.Verbose
                disp(output);
            end
            error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
        end
    end
    
    % Set the number of steps
    timesteps = Settings.Equilibrate_Solid/Settings.MDP.dt;
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
    
    % Ensure fast equilibration with Berendsen barostat + small time constant
    MDP_Template = regexprep(Settings.MDP_Template,'(nsteps += *)(.+?)( *);',['$1' num2str(timesteps) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = regexprep(MDP_Template,'(nstcalcenergy += *)(.+?)( *);','$11$3;');
    MDP_Template = regexprep(MDP_Template,'(nstxout += *)(.+?)( *);',['$1' num2str(round(0.1/Settings.MDP.dt)) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(pcoupl += *)(.+?)( *);','$1Berendsen$3;');
    MDP_Template = regexprep(MDP_Template,'(pcoupltype += *)(.+?)( *);',['$1' isotropy '$3;']);
    MDP_Template = regexprep(MDP_Template,'(tau-p += *)(.+?)( *);',['$1 ' num2str(tau_p) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nstpcouple += *)(.+?)( *);',['$1 ' num2str(nstpcouple) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(compressibility += *)(.+?)( *);',['$1 ' Compresstxt '$3;']);
    MDP_Template = regexprep(MDP_Template,'(ref-p += *)(.+?)( *);',['$1 ' ref_p '$3;']);
    MDP_Template = regexprep(MDP_Template,'(dt += *)(.+?)( *);',['$1' num2str(Settings.MDP.dt) '$3;']);    
    
    % Pair it with velocity rescale thermostat + small time constant
    MDP_Template = regexprep(MDP_Template,'(tcoupl += *)(.+?)( +);','$1v-rescale$3;');
    MDP_Template = regexprep(MDP_Template,'(tau-t += *)(.+?)( +);',['$1 ' num2str(tau_t) '$3;']);
    MDP_Template = regexprep(MDP_Template,'(nsttcouple += *)(.+?)( +);',['$1 ' num2str(nsttcouple) '$3;']);
    
    % Save MDP file
    MDP_Filename = fullfile(WorkDir,'Equil_Sol.mdp');
    fidMDP = fopen(MDP_Filename,'wt');
    fwrite(fidMDP,regexprep(MDP_Template,'\r',''));
    fclose(fidMDP);
    
    % Complete a topology file for the box to be minimized
    Atomlist = copy_atom_order(SuperCell_File);
    Settings.Topology_Text = strrep(Settings.Topology_Text,'##LATOMS##',Atomlist);
    Top_Filename = fullfile(WorkDir,'Equil_Sol.top');
    
    % Save topology file
    fidTOP = fopen(Top_Filename,'wt');
    fwrite(fidTOP,regexprep(Settings.Topology_Text,'\r',''));
    fclose(fidTOP);
    
    TPR_File = fullfile(WorkDir,'Equil_Sol.tpr');
    MDPout_File = fullfile(WorkDir,'Equil_Sol_out.mdp');
    GrompLog_File = fullfile(WorkDir,'Equil_Sol_Grompplog.log');
    
    FEquil_Grompp = [Settings.gmx_loc ' grompp -c ' windows2unix(SuperCell_File) ...
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
    Log_File = fullfile(WorkDir,'Equil_Sol.log');
    Energy_file = fullfile(WorkDir,'Equil_Sol.edr');
    TRR_File = fullfile(WorkDir,'Equil_Sol.trr');
    Final_Geom_File = fullfile(WorkDir,['Equil_Sol_out.' Settings.CoordType]);

    mdrun_command = [Settings.gmx ' mdrun -s ' windows2unix(TPR_File) ...
        ' -o ' windows2unix(TRR_File) ' -g ' windows2unix(Log_File) ...
        ' -e ' windows2unix(Energy_file) ' -c ' windows2unix(Final_Geom_File) ...
        ' -deffnm ' windows2unix(fullfile(WorkDir,'Equil_Sol')) ...
        Settings.mdrun_opts];

    if Settings.Table_Req
        mdrun_command = [mdrun_command ' -table ' windows2unix(Settings.TableFile_MX)];
    end

    % Final Equilibration
    if Settings.Verbose
        disp(['Beginning Solid Equilibration for ' num2str(Settings.Equilibrate_Solid) ' ps...'] )
    end
    mintimer = tic;
    [state,mdrun_output] = system(mdrun_command);
    if state == 0
        if Settings.Verbose
            disp(['Solid Successfully Equilibrated! Epalsed Time: ' datestr(seconds(toc(mintimer)),'HH:MM:SS')]);
        end
    else
        try % Clean up
            [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
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
%             Output = Equilibrate_Solid(Settings,'Skip_Cell_Construction',true);
%             return
        if Settings.MDP.dt > 1e-4
            if Settings.Verbose
                disp('Equilibration failed. Reducing time step.')
            end
            %Settings.QECompressibility = Settings.QECompressibility_init;
            Settings.MDP.dt = Settings.MDP.dt/2;
            Settings.Output_Coords = Settings.Output_Coords*2;
            Output = Equilibrate_Solid(Settings,'Skip_Cell_Construction',true);
            return
        else
            if Settings.Verbose
                disp('Equilibration failed. Reducing time step did not resolve.')
                disp(mdrun_output);
                disp(['Error running mdrun for solid equilibration. Problem command: ' newline mdrun_command]);
            end
            Output.Aborted = true;
            TDir = fullfile(strrep(Settings.WorkDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
            [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
            copyfile(WorkDir,TDir)
            try
                if Settings.Delete_Equil
                    rmdir(WorkDir,'s')
                end
            catch
                disp(['Unable to remove directory: ' WorkDir])
            end
            return
        end
    end
    
    % Check to ensure system remained in the correct solid structure
    PyOut = py.LiXStructureDetector.Calculate_Liquid_Fraction(WorkDir, Settings.Salt, ...
        pyargs('SystemName','Equil_Sol',...
        'RefStructure',Settings.Structure,...
        'CheckFullTrajectory',true,...
        'FileType',Settings.CoordType,...
        'ML_TimeLength',0,...
        'ML_TimeStep',0,...
        'SaveTrajectory',true,...
        'SavePredictionsImage',true,...
        'Verbose',Settings.Verbose));
    Sol_Fraction = PyOut{4};
    Liq_Fraction = PyOut{5};
    
    if Sol_Fraction < (1 - Settings.MeltFreezeThreshold)
        if Settings.Verbose
            disp('Detected Solid Phase change.')
        end
        
        if (1-Liq_Fraction-Sol_Fraction) >= Settings.MeltFreezeThreshold
            Output.StructureChange = true;
        else
            Output.SolidMelted = true;
        end
        
        Output.Aborted = true;
        TDir = fullfile(strrep(Settings.WorkDir,[filesep 'Minimization'],''),['T_' num2str(Settings.Target_T,'%.4f')]);
        [~,~] = system([Settings.wsl 'find ' windows2unix(WorkDir) ' -iname "#*#" ^| xargs rm']);
        copyfile(WorkDir,TDir)
        try
            if Settings.Delete_Equil
                rmdir(WorkDir,'s')
            end
        catch
            disp(['Unable to remove directory: ' WorkDir])
        end
        return
    end
    
    Box_xvg_file = fullfile(WorkDir,'Equil_Sol_Box.xvg');
    
    % Create gmx traj command
    startpoint = Settings.Equilibrate_Solid*0.50; % ps. Average over second half of equilibration period
    gmx_command = [strrep(Settings.gmx_loc,'gmx',['echo 0 ' Settings.pipe ' gmx']) ...
        ' traj -f ' windows2unix(TRR_File) ' -ob ' windows2unix(Box_xvg_file) ...
        ' -s ' windows2unix(TPR_File) ' -b ' num2str(startpoint) ...
        ' -e ' num2str(Settings.Equilibrate_Solid)];
    [err,~] = system(gmx_command);
    
    if err ~= 0
        error('Failed to collect box data.')
    end
    
    % Load data and calculate box vectors
    Data = import_xvg(Box_xvg_file); % Gather X,Y,Z lengths
    delete(Box_xvg_file) % remove temp output file
    
    a_vec = [Data(:,2) zeros(size(Data,1),1) zeros(size(Data,1),1)];
    b_vec = [Data(:,5) Data(:,3) zeros(size(Data,1),1)];
    c_vec = [Data(:,6) Data(:,7) Data(:,4)];

    a = vecnorm(a_vec,2,2);
    b = vecnorm(b_vec,2,2);
    c = vecnorm(c_vec,2,2);
    
%     hold on
%     plot(Data(:,1),10.*a./Na)
%     plot(Data(:,1),10.*b./Nb)
%     plot(Data(:,1),10.*c./Nc)
%     mean_a = mean(10.*a(ceil(length(a)/2):end)./Na)
%     std_a = std(10.*a(ceil(length(a)/2):end)./Na)
    
    switch isotropy
        case 'isotropic'
            Settings.Geometry.a = 10*mean(a)./Na; % Angstroms
            Settings.Geometry.b = Settings.Geometry.a; % Angstroms
            Settings.Geometry.c = Settings.Geometry.a; % Angstroms 
        case 'semiisotropic'
            Settings.Geometry.a = 10*mean(a)./Na; % Angstroms
            Settings.Geometry.b = Settings.Geometry.a; % Angstroms
            Settings.Geometry.c = 10*mean(c)./Nc; % Angstroms 
        case 'anisotropic'
            Settings.Geometry.a = 10*mean(a)./Na; % Angstroms
            Settings.Geometry.b = 10*mean(b)./Nb; % Angstroms
            Settings.Geometry.c = 10*mean(c)./Nc; % Angstroms 
    end
    
    % Load Coordinates text
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
    
    % Save the updated super cell file
    if sum([Settings.N_Supercell_a Settings.N_Supercell_b Settings.N_Supercell_c]) > 3
        
        fid = fopen(Settings.UnitCellFile,'wt');
        fwrite(fid,regexprep(Coordinate_Text,'\r',''));
        fclose(fid);
        
        Supercell_command = [Settings.gmx_loc ' genconf -f ' windows2unix(Settings.UnitCellFile) ...
             ' -o ' windows2unix(Settings.SuperCellFile) ' -nbox ' Settings.Na ' ' Settings.Nb ' ' Settings.Nc];
        
        % Remove old supercell file
        if isfile(Settings.SuperCellFile)
            delete(Settings.SuperCellFile)
        end
         
        [errcode,output] = system(Supercell_command);
        
        if errcode ~= 0
            if Settings.Verbose
                disp(output);
            end
            error(['Error creating supercell with genconf. Problem command: ' newline Supercell_command]);
        end
    else
        fid = fopen(Settings.SuperCellFile,'wt');
        fwrite(fid,regexprep(Coordinate_Text,'\r',''));
        fclose(fid);
    end
    
    try
        if Settings.Delete_Equil
            rmdir(WorkDir,'s')
        end
    catch
        if Settings.Verbose
            disp(['Unable to remove directory: ' WorkDir])
        end
    end
    
    if Settings.GenCluster
        Build_Cluster(Settings)
    end
    
    if Settings.Verbose
        disp('*** Separate Equilibration of Solid Complete ***')
    end
end